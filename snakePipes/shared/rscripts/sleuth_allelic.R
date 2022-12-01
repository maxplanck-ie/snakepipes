.libPaths(R.home("library"))

library("sleuth")
library("dplyr")
#library("biomaRt")

args = commandArgs(trailingOnly=TRUE)

## Debug only !!!
## t2g_file = "snakemake_workflows/shared/organisms/dm6.t2g"
## sample_info_file = "sampleSheet.tsv"
## indir = "Salmon"
## outdir = "sleuth"

sample_info_file = args[[1]]
indir = args[[2]]
outdir = args[[3]]
fdr = args[[4]]
t2g_file = args[[5]]

setwd(outdir)

sample_info = read.table(sample_info_file, header=T)#[,1:2]
coldata_allelic <- data.frame(name = paste0(rep(sample_info$name,each=2),c(".genome1",".genome2")),
                   allele = rep(c("genome1", "genome2"), nrow(sample_info)),
                   condition = rep(sample_info$condition, each = 2) )
coldata_allelic$allele<-factor(coldata_allelic$allele,levels=c("genome1","genome2"))
coldata_allelic$condition<-factor(coldata_allelic$condition,levels=unique(coldata_allelic$condition))
sample_info<-coldata_allelic

if(length(unique(sample_info$condition))>1){
message("Fitting full model with allele, condition and allele:condition.")
    d<-formula(~allele + condition + allele:condition)
}else{
message("Fitting reduced model with allele only.")
    d<-formula(~allele)
}
colnames(sample_info)[colnames(sample_info) %in% "name"] ="sample"
print(sample_info)
sample_info$sample

sample_id = list.dirs(file.path(indir), recursive=F, full.names=F)
sample_id = sort(sample_id[grep('[^benchmark][^SalmonIndex]', sample_id, invert=F)])
print(sample_id)
#sample_id = intersect(sample_info$sample, sample_id) # get only those sample that are defined in the sampleInfo!
sample_id<-sample_id[match(sample_info$sample,sample_id)]
print(sample_id)

salmon_dirs = sapply(sample_id, function(id) file.path(indir, id))
print(salmon_dirs)

s2c = mutate(sample_info, path=salmon_dirs)
## reorder conditions (for Wald test later on: order of comparison important for fold change)
s2c$condition<-factor(s2c$condition)
if ( s2c$condition[[1]] != levels(s2c$condition)[[1]] ) {
  s2c$condition =  relevel(s2c$condition, as.character(s2c$condition[[1]]) )
}
print(s2c)

## get gene names / symbol names
tryCatch( { t2g = read.table(t2g_file, header=F) },error = function(e) { print('No t2g file available!') },finally = {})

if (exists('t2g')) {
  colnames(t2g) <- c("target_id","ens_gene","ext_gene")

  ## add gene names
  so <- sleuth_prep(s2c, full_model=d, target_mapping = t2g)
} else {
  ## construct sleuth object
  so = sleuth_prep(s2c, full_model=d)
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
if(length(unique(sample_info$condition))>1){
    wald_beta_name = paste0("allelegenome2.condition",unique(coldata$condition)[2])
}else{
    wald_beta_name = "allelegenome2"
}
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
