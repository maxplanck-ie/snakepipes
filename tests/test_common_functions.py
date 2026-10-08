import argparse
import pytest

from snakePipes.common_functions import (
    merge_dicts,
    sanity_dict_clean,
    namesOKinR,
    config_diff,
    get_sample_names,
    get_sample_names_bam,
    get_sample_names_suffix_bam,
    is_paired,
    check_replicates,
    isMultipleComparison,
    splitSampleSheet,
    returnComparisonGroups,
    sampleSheetGroups,
    check_sample_info_header,
    checkAlleleParams,
    tr,
)


class TestMergeDicts:
    def test_basic_merge(self):
        assert merge_dicts({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}

    def test_y_overrides_x(self):
        assert merge_dicts({"a": 1}, {"a": 2}) == {"a": 2}

    def test_y_none(self):
        assert merge_dicts({"a": 1}, None) == {"a": 1}

    def test_does_not_mutate_x(self):
        x = {"a": 1}
        merge_dicts(x, {"b": 2})
        assert x == {"a": 1}


class TestSanityDictClean:
    def test_removes_unwanted_keys(self):
        d = {"maindir": "/x", "workflow": "y", "keep": "z"}
        assert sanity_dict_clean(d) == {"keep": "z"}

    def test_none_input(self):
        assert sanity_dict_clean(None) is None

    def test_missing_keys_noop(self):
        assert sanity_dict_clean({"keep": "z"}) == {"keep": "z"}


class TestNamesOKinR:
    def test_valid_names_no_warning(self, capsys):
        namesOKinR(["sample1", ".sample2", "a_b.c"])
        assert capsys.readouterr().err == ""

    def test_starts_with_digit_warns(self, capsys):
        namesOKinR(["1sample"])
        assert "1sample" in capsys.readouterr().err

    def test_reserved_word_warns(self, capsys):
        namesOKinR(["NULL"])
        assert "reserved keyword" in capsys.readouterr().err

    def test_invalid_char_warns(self, capsys):
        namesOKinR(["sample-1"])
        assert "invalid" in capsys.readouterr().err


class TestConfigDiff:
    def test_changed_value(self):
        assert config_diff({"a": 1}, {"a": 2}) == {"a": 1}

    def test_missing_in_dict2(self):
        assert config_diff({"a": 1, "b": 2}, {"a": 1}) == {"b": 2}

    def test_identical_no_diff(self):
        assert config_diff({"a": 1}, {"a": 1}) == {}


class TestGetSampleNames:
    def test_paired_end(self):
        infiles = [
            "sample1_R1.fastq.gz",
            "sample1_R2.fastq.gz",
            "sample2_R1.fastq.gz",
            "sample2_R2.fastq.gz",
        ]
        assert get_sample_names(infiles, ".fastq.gz", ["_R1", "_R2"]) == [
            "sample1",
            "sample2",
        ]

    def test_bad_suffix_warns_and_skips(self, capsys):
        infiles = ["sample1_R1.fastq.gz", "sample1_bad.fastq.gz"]
        result = get_sample_names(infiles, ".fastq.gz", ["_R1", "_R2"])
        assert result == ["sample1"]
        assert "does not have" in capsys.readouterr().err

    def test_no_matching_files_exits(self):
        with pytest.raises(SystemExit):
            get_sample_names(["sample1_bad.fastq.gz"], ".fastq.gz", ["_R1", "_R2"])


class TestGetSampleNamesBam:
    def test_strips_ext(self):
        infiles = ["/path/sample1.bam", "/path/sample2.bam"]
        assert get_sample_names_bam(infiles, ".bam") == ["sample1", "sample2"]

    def test_dedups_and_sorts(self):
        infiles = ["b.bam", "a.bam", "a.bam"]
        assert get_sample_names_bam(infiles, ".bam") == ["a", "b"]


class TestGetSampleNamesSuffixBam:
    def test_strips_genome_suffix(self):
        infiles = ["/path/sample1.genome1.bam", "/path/sample1.genome2.bam"]
        assert get_sample_names_suffix_bam(infiles, ".bam") == ["sample1"]

    def test_no_match_returns_empty(self):
        assert get_sample_names_suffix_bam(["/path/sample1.bam"], ".bam") == []


class TestIsPaired:
    def test_paired(self):
        infiles = ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]
        assert is_paired(infiles, ".fastq.gz", ["_R1", "_R2"]) is True

    def test_single_end(self):
        infiles = ["sample1_R1.fastq.gz", "sample2_R1.fastq.gz"]
        assert is_paired(infiles, ".fastq.gz", ["_R1", "_R2"]) is False

    def test_mixed_exits(self):
        infiles = ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz", "sample2_R1.fastq.gz"]
        with pytest.raises(SystemExit):
            is_paired(infiles, ".fastq.gz", ["_R1", "_R2"])

    def test_no_files_exits(self):
        with pytest.raises(SystemExit):
            is_paired([], ".fastq.gz", ["_R1", "_R2"])


def write_tsv(path, lines):
    path.write_text("\n".join(lines) + "\n")


class TestCheckReplicates:
    def test_all_replicated(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA", "s2\tA", "s3\tB", "s4\tB"])
        assert check_replicates(str(f)) is True

    def test_missing_replicate_warns_false(self, tmp_path, capsys):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA", "s2\tB", "s3\tB"])
        assert check_replicates(str(f)) is False
        assert "no replicates" in capsys.readouterr().err

    def test_bad_header_exits(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tfoo", "s1\tA"])
        with pytest.raises(SystemExit):
            check_replicates(str(f))


class TestIsMultipleComparison:
    def test_no_group_col(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA"])
        assert isMultipleComparison(str(f)) is False

    def test_single_group(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition\tgroup", "s1\tA\tG1", "s2\tB\tG1"])
        assert isMultipleComparison(str(f)) is None

    def test_multiple_groups(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition\tgroup", "s1\tA\tG1", "s2\tB\tG2"])
        assert isMultipleComparison(str(f)) is True


class TestReturnComparisonGroups:
    def test_no_group_col(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA"])
        assert returnComparisonGroups(str(f)) is False

    def test_excludes_all(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition\tgroup", "s1\tA\tAll", "s2\tB\tG1"])
        assert set(returnComparisonGroups(str(f))) == {"G1"}


class TestSampleSheetGroups:
    def test_simple_grouping(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA", "s2\tB"])
        assert sampleSheetGroups(str(f), multipleComp=False) == {
            "A": ["s1"],
            "B": ["s2"],
        }

    def test_bad_header_exits(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tfoo", "s1\tA"])
        with pytest.raises(SystemExit):
            sampleSheetGroups(str(f), multipleComp=False)


class TestSplitSampleSheet:
    def test_splits_by_group(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "splitSampleSheets").mkdir()
        f = tmp_path / "sheet.tsv"
        write_tsv(
            f,
            [
                "name\tcondition\tgroup",
                "s1\tA\tG1",
                "s2\tB\tG1",
                "s3\tA\tG2",
                "s4\tB\tG2",
            ],
        )
        splitSampleSheet(str(f), "sheet")
        out1 = (tmp_path / "splitSampleSheets" / "sheet.G1.tsv").read_text()
        out2 = (tmp_path / "splitSampleSheets" / "sheet.G2.tsv").read_text()
        assert "s1\tA" in out1 and "s2\tB" in out1
        assert "s3\tA" in out2 and "s4\tB" in out2


class TestCheckSampleInfoHeader:
    def test_valid_header(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["name\tcondition", "s1\tA"])
        assert check_sample_info_header(str(f)) == str(f.resolve())

    def test_missing_file_exits(self, tmp_path):
        with pytest.raises(SystemExit):
            check_sample_info_header(str(tmp_path / "nope.tsv"))

    def test_bad_header_exits(self, tmp_path):
        f = tmp_path / "sheet.tsv"
        write_tsv(f, ["foo\tbar", "s1\tA"])
        with pytest.raises(SystemExit):
            check_sample_info_header(str(f))


class TestCheckAlleleParams:
    def _args(self, **kw):
        defaults = dict(
            mode="mapping",
            SNPfile="/nonexistent/snp",
            VCFfile="/nonexistent/vcf",
            strains="",
            NMaskedIndex="/nonexistent/idx/foo",
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_plain_mapping_returns_none(self):
        assert checkAlleleParams(self._args(mode="mapping")) is None

    def test_conflicting_modes_exits(self):
        with pytest.raises(SystemExit):
            checkAlleleParams(self._args(mode="allelic-mapping,mapping"))

    def test_allelic_no_snp_no_vcf_exits(self):
        with pytest.raises(SystemExit):
            checkAlleleParams(self._args(mode="allelic-mapping"))

    def test_allelic_vcf_without_strain_exits(self, tmp_path):
        vcf = tmp_path / "some.vcf"
        vcf.write_text("x")
        with pytest.raises(SystemExit):
            checkAlleleParams(
                self._args(mode="allelic-mapping", VCFfile=str(vcf), strains="")
            )

    def test_allelic_vcf_with_strain_ok(self, tmp_path):
        vcf = tmp_path / "some.vcf"
        vcf.write_text("x")
        result = checkAlleleParams(
            self._args(mode="allelic-mapping", VCFfile=str(vcf), strains="strainA")
        )
        assert result == "create_and_map"

    def test_allelic_snp_with_index_ok(self, tmp_path):
        snp = tmp_path / "some.snp"
        snp.write_text("x")
        idx_dir = tmp_path / "idxdir"
        idx_dir.mkdir()
        result = checkAlleleParams(
            self._args(
                mode="allelic-mapping",
                SNPfile=str(snp),
                NMaskedIndex=str(idx_dir / "foo"),
            )
        )
        assert result == "map_only"


class TestTr:
    def test_replaces_null_with_none(self):
        assert tr("value: null") == "value: None"

    def test_no_null_unchanged(self):
        assert tr("value: 1") == "value: 1"
