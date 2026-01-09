#!/usr/bin/env python3
"""
Merge two tabular files (A and B), duplicate rows from B to match A (A ids have suffixes),
compute rowwise product = A.value * B.value, print merged dataframe for debugging,
and write an output file with two columns:
  - sample (original id from file A)
  - product (numeric product)

Usage:
    python merge_and_multiply_and_write_fixed.py fileA.tsv fileB.tsv out.tsv
Options (examples):
    --sep '\t' --a-id sample --a-val value --b-id sample --b-val value
    --suffixes .genome1,.genome2 --float-format ".6f" --how left --print-all
"""
from typing import Sequence, Optional, Tuple, List
import re
import argparse
import sys
import pandas as pd


def _strip_colnames(df: pd.DataFrame) -> None:
    # Normalize column names: remove BOM and surrounding whitespace
    df.columns = df.columns.astype(str)
    df.columns = df.columns.str.replace(r"^\ufeff", "", regex=True).str.strip()


def _parse_suffixes(s: str) -> List[str]:
    return [tok for tok in (t.strip() for t in s.split(",")) if tok]


def _build_suffix_pattern(suffixes: Sequence[str]) -> str:
    escaped = "|".join(re.escape(s) for s in suffixes if s)
    return rf"(?:{escaped})$" if escaped else r"$^"  # match nothing if no suffixes provided


def _read_table(path: str, sep: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=[""])
    except Exception as e:
        raise RuntimeError(f"Failed to read '{path}': {e}")
    _strip_colnames(df)
    return df


def _infer_columns_if_needed(
    df: pd.DataFrame, id_col: str, val_col: str, role: str, suffixes: Sequence[str]
) -> Tuple[str, str]:
    """
    Ensure id_col and val_col exist in df. If not, attempt to infer:
      - if df has exactly 2 columns, use first -> id, second -> value
      - otherwise try to find a column whose values contain suffixes (useful for file A)
      - otherwise raise a helpful error
    Returns (id_col, val_col) actually used.
    """
    cols = list(df.columns)
    if id_col in df.columns and val_col in df.columns:
        return id_col, val_col

    # If file has exactly two columns, use them
    if len(cols) == 2:
        inferred_id, inferred_val = cols[0], cols[1]
        print(f"Notice: {role}: using inferred columns '{inferred_id}' (id) and '{inferred_val}' (value)")
        return inferred_id, inferred_val

    # Try to find id column by suffix presence (only sensible for file A)
    if role == "A" and suffixes:
        escaped = "|".join(re.escape(s) for s in suffixes if s)
        if escaped:
            pattern = re.compile(rf"(?:{escaped})$")
            for c in cols:
                try:
                    if df[c].astype(str).str.contains(pattern).any():
                        # pick another numeric-looking column as value if any
                        val_candidates = [x for x in cols if x != c]
                        if val_candidates:
                            # pick the first candidate that looks numeric
                            for vc in val_candidates:
                                try:
                                    sample_vals = pd.to_numeric(df[vc].astype(str).str.strip(), errors="coerce")
                                    if sample_vals.notna().any():
                                        return c, vc
                                except Exception:
                                    continue
                            # fallback: pick first other column
                            return c, val_candidates[0]
                except Exception:
                    continue

    # As last resort, if requested id exists but value doesn't, try to pick a numeric-like column
    if id_col in df.columns:
        num_like = None
        for c in cols:
            if c == id_col:
                continue
            try:
                if pd.to_numeric(df[c].astype(str).str.strip(), errors="coerce").notna().any():
                    num_like = c
                    break
            except Exception:
                continue
        if num_like:
            print(f"Notice: {role}: using '{id_col}' as id and inferred numeric column '{num_like}' as value")
            return id_col, num_like

    # If nothing worked, show available columns and fail
    raise KeyError(
        f"{role}: could not find columns '{id_col}' and '{val_col}' in file. Available columns: {cols}"
    )


def merge_and_write(
    file_a: str,
    file_b: str,
    output_file: str,
    a_id_col: str = "sample",
    a_val_col: str = "value",
    b_id_col: str = "sample",
    b_val_col: str = "value",
    suffixes: Sequence[str] = (".genome1", ".genome2"),
    sep: str = "\t",
    product_col: str = "product",
    how: str = "left",
    float_format: Optional[str] = None,
    print_all: bool = False,
) -> pd.DataFrame:
    # Read tables
    df_a = _read_table(file_a, sep)
    df_b = _read_table(file_b, sep)

    # Validate or infer columns
    a_id_col, a_val_col = _infer_columns_if_needed(df_a, a_id_col, a_val_col, role="A", suffixes=suffixes)
    b_id_col, b_val_col = _infer_columns_if_needed(df_b, b_id_col, b_val_col, role="B", suffixes=[])

    # Trim whitespace in ID columns (important for merges)
    df_a[a_id_col] = df_a[a_id_col].astype(str).str.strip()
    df_b[b_id_col] = df_b[b_id_col].astype(str).str.strip()

    # Convert value columns to numeric (coerce non-numeric to NaN)
    df_a[a_val_col] = pd.to_numeric(df_a[a_val_col].astype(str).str.strip(), errors="coerce")
    df_b[b_val_col] = pd.to_numeric(df_b[b_val_col].astype(str).str.strip(), errors="coerce")

    # Warn if conversion produced NaNs
    if df_a[a_val_col].isna().any():
        n = df_a[a_val_col].isna().sum()
        print(f"Warning: {n} non-numeric/missing values in file A column '{a_val_col}' -> coerced to NaN", file=sys.stderr)
    if df_b[b_val_col].isna().any():
        n = df_b[b_val_col].isna().sum()
        print(f"Warning: {n} non-numeric/missing values in file B column '{b_val_col}' -> coerced to NaN", file=sys.stderr)

    # Build base_id in A by stripping suffixes
    suffix_pattern = _build_suffix_pattern(suffixes)
    df_a = df_a.copy()
    df_a["base_id"] = df_a[a_id_col].astype(str).str.replace(suffix_pattern, "", regex=True)

    # Prepare B for merging
    df_b_merge = df_b.rename(columns={b_id_col: "base_id", b_val_col: "b_value"})

    # Merge (duplicates B rows to match A)
    merged = pd.merge(df_a, df_b_merge, on="base_id", how=how, sort=False)

    # Compute product
    merged[product_col] = merged[a_val_col] * merged["b_value"]

    # Print diagnostics
    print(f"File A rows: {len(df_a)}, File B rows: {len(df_b)}, Merged rows: {len(merged)}")
    if print_all:
        pd.set_option("display.max_rows", None)
        pd.set_option("display.max_columns", None)
        pd.set_option("display.width", None)
        print("Merged DataFrame (full):")
        print(merged)
    else:
        print("Merged DataFrame (first 20 rows):")
        print(merged.head(20))

    # Prepare output DataFrame with original A sample ids and product
    out_df = merged[[a_id_col, product_col]].copy()
    out_df.columns = ["sample", "product"]  # uniform output column names

    # Write output
    to_csv_kwargs = {"sep": sep, "index": False}
    if float_format:
        # user-supplied float format should be like ".6f" (no % in help to avoid argparse issues)
        try:
            # pandas expects a format string like "%.6f", but the user may pass ".6f" -> prepend % if needed
            fmt = float_format if float_format.startswith("%") else f"%{float_format}"
            to_csv_kwargs["float_format"] = fmt
        except Exception:
            print("Warning: invalid float_format ignored", file=sys.stderr)

    out_df.to_csv(output_file, **to_csv_kwargs)
    print(f"Wrote {len(out_df)} rows to '{output_file}'")

    return merged


def main():
    p = argparse.ArgumentParser(description="Merge fileA and fileB and write sample/product output.")
    p.add_argument("file_a", help="Path to file A (sample ids with suffixes).")
    p.add_argument("file_b", help="Path to file B (base sample ids).")
    p.add_argument("output_file", help="Path to write the two-column output (sample, product).")
    p.add_argument("--sep", default="\t", help="Separator for input/output files (default: tab).")
    p.add_argument("--a-id", default="sample", dest="a_id", help="ID column name in file A (default: sample).")
    p.add_argument("--a-val", default="value", dest="a_val", help="Numeric value column name in file A (default: value).")
    p.add_argument("--b-id", default="sample", dest="b_id", help="ID column name in file B (default: sample).")
    p.add_argument("--b-val", default="value", dest="b_val", help="Numeric value column name in file B (default: value).")
    p.add_argument("--suffixes", default=".genome1,.genome2", help="Comma-separated suffixes to strip from A ids (default: .genome1,.genome2).")
    p.add_argument("--float-format", default=None, help='Optional float format for output values, e.g. ".6f" (no percent sign).')
    p.add_argument("--how", default="left", choices=("left", "inner"), help="Merge type (left or inner).")
    p.add_argument("--print-all", action="store_true", help="Print the full merged DataFrame instead of only head().")
    args = p.parse_args()

    suffixes = _parse_suffixes(args.suffixes)
    try:
        merge_and_write(
            args.file_a,
            args.file_b,
            args.output_file,
            a_id_col=args.a_id,
            a_val_col=args.a_val,
            b_id_col=args.b_id,
            b_val_col=args.b_val,
            suffixes=suffixes,
            sep=args.sep,
            how=args.how,
            float_format=args.float_format,
            print_all=args.print_all,
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()