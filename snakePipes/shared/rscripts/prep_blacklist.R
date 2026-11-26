#!/usr/bin/env Rscript
.libPaths(R.home("library"))

wdir<-snakemake@params[["wdir"]]
setwd(wdir)
input_file<-snakemake@params[["url"]]
output_file<-snakemake@output[["path"]]

gtf_file<-snakemake@input[["gtf"]]

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomeInfoDb)
})


# ── helpers ──────────────────────────────────────────────────────────────────

is_url <- function(x) grepl("^https?://|^ftp://", x)
is_rds <- function(x) grepl("\\.rds$", x, ignore.case = TRUE)


# ── load ─────────────────────────────────────────────────────────────────────

if (is_rds(input_file)) {

  if (is_url(input_file)) {
    tmp <- tempfile(fileext = ".rds")
    download.file(input_file, destfile = tmp, mode = "wb", quiet = TRUE)
    gr <- readRDS(tmp)
    unlink(tmp)
  } else {
    gr <- readRDS(input_file)
  }

} else {

  # rtracklayer::import handles URLs and local paths for BED, GTF, GFF, bigBed, etc.
  gr <- import(input_file)

}

stopifnot("Input must be a GRanges object" = is(gr, "GRanges"))


# ── deduplicate ───────────────────────────────────────────────────────────────

gr <- reduce(gr)

message(length(gr), " ranges after deduplication")

# check for chromosome style in gtf

gtf<-import(gtf_file)
lvls     <- seqlevels(gtf)
has_chr  <- grepl("^chr", lvls)

if (!any(has_chr)) {
  seqlevelsStyle(gr)<-"ucsc"}else{
  seqlevelsStyle(gr)<-"Ensembl"
}

# ── export ───────────────────────────────────────────────────────────────────

export(gr, output_file, format = "BED")
