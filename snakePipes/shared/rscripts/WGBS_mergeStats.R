.libPaths(R.home("library"))
sampleSheet = snakemake@params[['sampleSheet']]
groups = snakemake@params[['groups']]
blacklist = snakemake@params[['blacklist']]
minCoverage = snakemake@params[['minCoverage']]

require("bsseq")
require("rtracklayer")
require("BiocParallel")
BPPARAM = MulticoreParam(workers=snakemake@threads)

ss = read.delim(sampleSheet)
row.names(ss) = as.character(ss$name)
infiles = sprintf("MethylDackel/%s_CpG.bedGraph", ss$name)

bs = read.bismark(infiles, colData=ss, BPPARAM=BPPARAM)

if(length(blacklist)) {
    bl = import.bed(blacklist)
    bs = subsetByOverlaps(bs, bl, invert=TRUE)
}

# We require regions with a minimal coverage in all samples
RMV <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") < minCoverage) != 0)
if(length(RMV)) {
    bs = bs[-RMV,]
}

# The output format is: chrom position methylation values (between 0 and 100)
d = data.frame(chr=seqnames(bs),
               pos=start(bs))
m = 100 * getMeth(bs, type="raw")
d = cbind(d, m)
colnames(d)[3:ncol(d)] = sprintf("%s_%s", ss$condition, ss$name)

write.table(d, file=snakemake@output[["MetileneIN"]], sep="\t", row.names=FALSE, quote=FALSE)