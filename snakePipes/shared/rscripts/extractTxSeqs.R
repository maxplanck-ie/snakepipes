#this is a copy of https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/extractTxSeqs.R , authored by C. Soneson

#' Extract (spliced or unspliced) transcript sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'spliced' or 'unspliced'
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#'
extractTxSeqs <- function(gtf, genome, type = "spliced") {
  ## Construct TxDb from gtf file. 
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")

  ## Group exons by transcript. When using exonsBy with by = "tx", 
  ## the returned exons are ordered by ascending rank for each transcript, 
  ## that is, by their position in the transcript. 
  grl <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

  ## Extract transcript sequences.
  ## Here, it's important that for each transcript, the exons must be ordered 
  ## by ascending rank, that is, by ascending position in the transcript.
  if (type == "spliced") {
    txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = grl)
  } else if (type == "unspliced") {
    txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = range(grl))
  } else {
    stop("Unknown 'type' argument")
  }

  return(txout)
}
