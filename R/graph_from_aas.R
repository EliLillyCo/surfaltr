#' Create a plot showing membrane locations of each protein
#' based on user provided amino acid sequences
#'
#' This function creates a ggplot figure showing the differences in membrane
#' location and length between primary and alternative transcripts from the
#' same gene. This process is performed based on input data
#' containing the gene names and amino acid sequences of the proteins
#' in question. Transcripts derived from the same gene are grouped together
#' to facilitate easy interpretation. The y axis lists the gene name and
#' transcript ID for each transcript and the x axis lists the length in
#' amino acids. Each fill color corresponds to a membrane location and either
#' principal or alternative isoform.
#'
#' @usage graph_from_aas(data_file, organism = "human", rank = "length",
#' n_prts = 20, mode = "phobius", size_txt = 2, space_left = -400, temp = FALSE)
#' @param data_file Path to the input file
#' @param organism String indicating if the transcripts are from a human or
#' a mouse
#' @param rank String indicating which method to use to rank proteins in graphical
#' output. Options include "Length", "TM", and "Combo".
#' @param n_prts Integer value indicating the number of genes that should be
#' displayed on the graphical output. Default value is 20.
#' @param mode String detailing whether TMHMM or Phobius should be used to
#' predict transmembrane regions. Input values include "phobius" or "tmhmm".
#' @param size_txt Integer value specifying the size of the row labels. Default
#' size is 2.
#' @param space_left Integer value specifying how far left the graph should
#' extend.
#' @param temp Boolean indicating if the fasta file should be deleted after the
#' function finishes running or not. Recommended to always be set to FALSE.
#' @return A ggplot figure showing the protein locations for each part of the
#' surface protein for each alternative and primary transcripts.
#' @examples graph_from_aas(system.file("extdata", "CRB1.csv", package = "surfaltr"),
#' "mouse", "combo", 10, "phobius", 2, -400, TRUE)
#' @export

graph_from_aas <- function(data_file, organism = "human", rank = "length", n_prts = 20, mode = "phobius", size_txt = 2, space_left = -400, temp = FALSE){
  old <- getwd()
  if_aa <- TRUE
  final_trans <- clean_data(data_file, if_aa, organism)
  princ <- ensembl_db_retrieval(organism)
  final_pairs <- merge_trans(princ, final_trans, if_aa)
  AA_seq <- get_aas(final_pairs, temp)
  if(mode == "phobius"){
    topo <- run_phobius(AA_seq, paste(getwd(), "/AA.fasta", sep = ""))
  }
  if(mode == "tmhmm"){
    topo <- get_tmhmm("AA.fasta")
  }
  counts <- process_tmhmm(topo, AA_seq)
  p <- graph_prots(counts, rank, n_prts, size_txt, space_left)
  setwd(old)
  p
}
