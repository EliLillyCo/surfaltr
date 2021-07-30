#' Create csv and fasta files containing information about pairs of transcripts
#'
#' This function processes the input data to retrieve
#' information from ensembl and uniprot to generate a dataframe containing
#' the gene names, transcript IDs, APPRIS annotations, and protein sequences
#' for each pair of primary and alternative transcripts. Additionally, this
#' function creates a fasta file with the transcript ID followed by the
#' amino acid sequence for all inputted and associated primary transcripts. The
#' file is organized so that all transcripts from a gene are next to each other.
#' Finally, the function also produces a final table in csv form
#' containing the gene names, transcript IDs, APPRIS annotations, and amino acid
#' sequences for each transcript
#'
#' @usage get_pairs(data_file, if_aa = FALSE, organism = "human", temp = FALSE)
#' @param data_file Path to the input file
#' @param if_aa Boolean value indicating if the input file contains
#' amino acid sequences with TRUE indicating that sequences are present and
#' FALSE indicating that only IDs are present
#' @param organism String indicating if the transcripts are from a human or
#' a mouse
#' @param temp Boolean indicating if the fasta file should be deleted after the
#' function finishes running or not. Recommended to always be set to FALSE.
#' @return A data frame containing the gene names, transcript IDs, APPRIS
#' annotations,and protein sequences for each pair of primary and alternative
#' transcripts.
#' @note This function also creates a fasta file containing the transcript IDs
#' and associated amino acid sequences in the root directory. In addition to the
#' fasta file, a csv file containing the returned dataframe is saved to the
#' working directory.
#' @importFrom utils "write.csv"
#' @examples
#' \donttest{
#' AA_seq <- get_pairs(system.file("extdata", "CRB1.csv",
#' package = "surfaltr"), TRUE, "mouse", TRUE)
#' }
#' @export

get_pairs <- function(data_file, if_aa = FALSE, organism = "human", temp = FALSE){
  final_trans <- clean_data(data_file, if_aa, organism)
  princ <- ensembl_db_retrieval(organism)
  final_pairs <- merge_trans(princ, final_trans, if_aa)
  if (if_aa == FALSE){
    aa_trans <- format_ids(final_pairs)
    AA_seq <- get_prts(aa_trans, temp)
    expt <- data.frame(lapply(AA_seq, as.character), stringsAsFactors=FALSE)
    write.csv(expt, paste(getwd(), "/transcript_pairs.csv", sep = ""))
  }
  if (if_aa == TRUE){
    AA_seq <- get_aas(final_pairs, temp)
    expt <- data.frame(lapply(AA_seq, as.character), stringsAsFactors=FALSE)
    write.csv(expt, paste(getwd(), "/transcript_pairs.csv", sep = ""))
  }
  return(AA_seq)
}
