#' Create a plot showing where each amino acid is located within the cell for
#' each primary transcript compared to each alternative transcript
#'
#' This function creates a ggplot figure showing the differences in membrane
#' location and length between primary and alternative transcripts from the
#' same gene. Transcripts derived from the same gene are grouped together
#' to facilitate easy interpretation. The y axis lists the gene name and
#' transcript ID for each transcript and the x axis lists the length in
#' amino acids. Each fill color corresponds to a membrane location and either
#' principal or alternative isoform.
#'
#' @usage plot_isoforms(topo, AA_seq, rank = "length", n_prts = 20,
#' size_txt = 2, space_left = -400)
#' @param topo Outputted data frame from the run_phobius or get_tmhmm function
#' showing membrane locations of amino acids and transcript IDs
#' @param AA_seq A data frame outputted by the get_pairs function
#' containing the gene names, transcript IDs, APPRIS annotations, and protein
#' sequences for each transcript.
#' @param rank String indicating which method to use to rank proteins in graphicl
#' output. Options include "length", "TM", and "combo".
#' @param n_prts Integer value indicating the number of genes that should be
#' displayed on the graphical output. Default value is 20.
#' @param size_txt Integer value specifying the size of the row labels. Default
#' size is 2.
#' @return A ggplot figure showing the protein locations for each part of the
#' surface protein for each alternative and primary transcripts.
#' @param space_left Integer value specifying how far left the graph should
#' extend.
#' @examples
#' \donttest{
#' currwd <- getwd()
#' AA_seq <- get_pairs(system.file("extdata", "crb1_example.csv",
#' package = "surfaltr"), TRUE, "mouse", TRUE)
#' topo <- run_phobius(AA_seq, paste(getwd(), "/AA.fasta", sep = ""))
#' plot_isoforms(topo, AA_seq, "combo", 15, 3, -400)
#' setwd(currwd)
#' }
#' @export

plot_isoforms <- function(topo, AA_seq, rank = "length", n_prts = 20, size_txt = 2, space_left = -400){
  counts <- process_tmhmm(topo, AA_seq)
  p <- graph_prots(counts, rank, n_prts, size_txt, space_left)
  p
}
