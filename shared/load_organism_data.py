try:
    include: os.path.join(maindir, "shared", "organisms", genome+".py")
except:
    print("Error! Organism-specific configuration file NOT found for:", genome, "\n")
    exit(1)

if "verbose" in config and config["verbose"]==True:
    print("--- Genome ---------------------------------------------------------------------")
    print("Genome:", genome)
    print("Effective genome size:", genome_size)
    print("Bowtie2 index:", bowtie2_index)
    print("Gene annotation GTF:", genes_gtf)
    print("-" * 80)
