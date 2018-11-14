library("sleuth")
library("dplyr")
#library("biomaRt")

args = commandArgs(trailingOnly=TRUE)

## Debug only !!!
## t2g_file = "snakemake_workflows/shared/organisms/dm6.t2g"
## sample_info_file = "sampleInfo.tsv"
## indir = "Salmon"
## outdir = "sleuth"

sample_info_file = args[[1]]
indir = args[[2]]
outdir = args[[3]]
fdr = args[[4]]
t2g_file = args[[5]]

setwd(outdir)
getwd()

sample_info = read.table(sample_info_file, header=T)[,1:2]
colnames(sample_info) = c("sample", "condition")
print(sample_info)
sample_info$sample

sample_id = list.dirs(file.path(indir), recursive=F, full.names=F)
sample_id = sort(sample_id[grep(".benchmark|SalmonIndex", sample_id, invert=T)])
sample_id = intersect(sample_info$sample, sample_id) # get only those sample that are defined in the sampleInfo!
sample_id

salmon_dirs = sapply(sample_id, function(id) file.path(indir, id))
salmon_dirs

s2c = mutate(sample_info, path=salmon_dirs)
## reorder conditions (for Wald test later on: order of comparison important for fold change)
if ( s2c$condition[[1]] != levels(s2c$condition)[[1]] ) {
  s2c$condition =  relevel(s2c$condition, as.character(s2c$condition[[1]]) )
}
print(s2c)


# ## get gene names from Biomart ###############################################
# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                          dataset = "dmelanogaster_gene_ensembl",
#                          host = 'ensembl.org')
#
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#                                      "external_gene_name"), mart = mart)
#
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#
# write.table(t2g, "dm6.t2g", col.names=T, row.names=F, sep="\t", quote=F)
#
################################################################################

## get gene names / symbol names
tryCatch( { t2g = read.table(t2g_file, header=F) },
          error = function(e) { print("No t2g file available!") },
          finally = {}
)

if (exists("t2g")) {
 	colnames(t2g) <- c("target_id","ens_gene","ext_gene")
	str(t2g)

  ## add gene names
  so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
} else {
  ## construct sleuth object
  so = sleuth_prep(s2c, ~ condition)
}

## model expression responding on condition
so = sleuth_fit(so)

# ## fit reduced model (not depending on any factor)
# so <- sleuth_fit(so, ~1, 'reduced')
#
# ## likelihood ratio test (LRT) between models to get transcripts affected by condition.
# ## Usually considered more stringent than Wald test. No fold change!
# so <- sleuth_lrt(so, 'reduced', 'full')
# results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
# head(results_table_lrt)

## Wald test (to get *fold change*)
wald_beta_name = paste("condition", tail(levels(s2c$condition), n=1), sep="")
so <- sleuth_wt(so, wald_beta_name, "full")
results_table_wt <- sleuth_results(so, wald_beta_name)
head(results_table_wt)
write.table(results_table_wt, "Wald-test.results.tsv", col.names=T, row.names=F, quote=F, sep="\t")


## view fitted models and tests
models(so)
tests(so)

saveRDS(so, file='so.rds')
so = readRDS("so.rds")

write(c("library(sleuth)",
        "so = readRDS('so.rds')",
        "sleuth_live(so)"),
      file="sleuth_live.R")

pdf("MA-plot.pdf")
plot_ma(so, wald_beta_name, sig_level = fdr)
dev.off()


# sleuth_live(so, "wt") # wt for Wald test
