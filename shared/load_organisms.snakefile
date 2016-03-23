try:
    include: os.path.join(maindir, "shared", "organisms", genome+".py")
except:
    "Error! Organism specific configuration file NOT found for:", genome, "\n"

if "verbose" in config and config["verbose"]==True:
    print("-- Genome ----------------------------------------------------------------------")
    print("Gsize:", gsize)
    print("Bowtie2Index:", Bowtie2Index)
    print("-" * 80)
