### Read organism configuration file and set organism-specific variables #######
################################################################################
try:
    include: os.path.join(maindir, "shared", "organisms", genome+".py")
except:
    print("ERROR: Organism-specific configuration file NOT found for:", genome, "\n")
    exit(1)
