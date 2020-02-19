import os
import loompy
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-outf", "--output", dest = "outf", help = "output file", metavar = "FILE")
parser.add_argument('input', nargs='+', help = "input folders")

args = parser.parse_args()

input = args.input
outf = args.outf

filelist = []
for p in input:
    z = os.listdir(p)
    f = list(filter(lambda x: '.loom' in x, z))
    ifi = os.path.join(p,f[0])
    filelist.append(ifi)
print(filelist)
print(outf)
loompy.combine(files = filelist, output_file = outf, key = "Accession")
