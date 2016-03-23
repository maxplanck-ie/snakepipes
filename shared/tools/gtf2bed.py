import sys

infile = sys.argv[1]

with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line.startswith("#"):
            columns = line.split("\t")
            if columns[2]=="gene":
                geneid = columns[8].split()[1].replace('"','').replace(';','')
                print "{}\t{}\t{}\t{}\t{}\t{}".format(columns[0], columns[3], columns[4], geneid, 0, columns[6])
