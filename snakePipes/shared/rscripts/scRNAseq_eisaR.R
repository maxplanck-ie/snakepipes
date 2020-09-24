#this is a modification of https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/generate_cdna_intron_fa_prepref.R , authored by C.Soneson
sink(snakemake@log[["out"]])
.libPaths(R.home("library"))

wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(GenomicFeatures)
  library(BiocGenerics)
  library(BSgenome)
  library(GenomicRanges)
})

gtf<-snakemake@input["gtf"]
genome<-snakemake@input["fasta"]
scriptdir<-snakemake@params["scriptdir"]
isoform_action<-snakemake@params["isoform_action"]
flanklength<-snakemake@params["flank_length"]
joint_fasta<-snakemake@output["joint_fasta"]
joint_t2g<-snakemake@output["joint_t2g"]

print(scriptdir)
print(gtf)
print(genome)
print(isoform_action)
print(flanklength)

source(file.path(scriptdir, "extractIntronSeqs.R"))
source(file.path(scriptdir, "extractTxSeqs.R"))

## Extract intronic sequences flanked by L-1 bases 
## of exonic sequences where L is the biological read length
genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
gtfdf <- as.data.frame(rtracklayer::import(gtf))

## Extract transcript and intron sequences
tx <- extractTxSeqs(gtf = gtf, genome = genome, type = "spliced")
intr <- extractIntronSeqs(gtf = gtf, genome = genome, type = isoform_action, 
                          flanklength = flanklength,
                          joinOverlappingIntrons = FALSE)

## Generate transcript/intron-to-gene mapping
t2gtx <- gtfdf %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()
if (isoform_action == "collapse") {
  ## Intron names already contain gene name
  t2gin <- data.frame(intr = names(intr),
                      gene = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE)
} else if (isoform_action == "separate") {
  ## Intron names contain transcript name
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
} else {
  stop("Unknown isoform_action")
}
colnames(t2gin) <- colnames(t2gtx)
t2g <- rbind(t2gtx, t2gin)

Biostrings::writeXStringSet(c(tx, intr), joint_fasta,  compress = FALSE)
write.table(names(tx), file = file.path(wdir, "cDNA_tx_to_capture.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(names(intr), file = file.path(wdir, "introns_tx_to_capture.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(t2g, file = joint_t2g, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




message('done all')
sink()

sink("sessionInfo.txt")
sessionInfo()
sink()


