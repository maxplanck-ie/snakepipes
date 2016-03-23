import glob
import os

## from DNA-mapping workflow !!!
include: os.path.join(maindir, "workflows", "DNA-mapping", "internals.snakefile")


## ChIP-seq grouping of files from config.yaml: [[ChIP1, Input1],[ChIP2, Input2]]
groups = [[item["ChIP"], item["Input"]] for item in config["samples"]]
samples_from_groups = [item for sublist in groups for item in sublist]
## print(config["samples"])
## print("Groups:", groups)
## print("Samples_from_groups", samples_from_groups)
## print("Samples:", samples)

## Check if sample names form config.yaml match to files
unmatched_sample_names = list(set(samples_from_groups) - set(samples))
if unmatched_sample_names:
    print ("\nWARNING! Bad sample names:", ", ".join(unmatched_sample_names), "\n")
    with open(os.path.join(outdir, "bad_sample_names.txt"), "w") as f:
        for x in (unmatched_sample_names):
            f.write(x+"\n")
