#this is a modification of https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/generate_cdna_intron_fa_prepref.R , authored by C.Soneson
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

gtf<-snakemake@params[["gtf"]]
genome_fasta<-snakemake@input[["genome_fasta"]]
scriptdir<-snakemake@params[["scriptdir"]]
isoform_action<-snakemake@params[["isoform_action"]]
flanklength<-as.integer(snakemake@params[["flank_length"]])
joint_fasta<-snakemake@params[["joint_fasta"]]
joint_t2g<-snakemake@params[["joint_t2g"]]

print(scriptdir)
print(gtf)
print(genome_fasta)
print(isoform_action)
print(flanklength)

print("sourcing extractIntronSeqs.R ..")
source(file.path(scriptdir, "extractIntronSeqs.R"))
print("sourcing extractTxSeqs.R ..")
source(file.path(scriptdir, "extractTxSeqs.R"))
print("..done")

## Extract intronic sequences flanked by L-1 bases 
## of exonic sequences where L is the biological read length
print("loading genome ..")
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
print("..done")
print("loading gtf ..")
gtfdf <- as.data.frame(rtracklayer::import(gtf))
print("..done")

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
   t2gin$gene_id<-paste0(t2gin$gene_id,"-I")
} else if (isoform_action == "separate") {
  ## Intron names contain transcript name
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
  t2gin$gene_id<-paste0(t2gin$gene_id,"-I")
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

sink("sessionInfo.txt")
sessionInfo()
sink()


