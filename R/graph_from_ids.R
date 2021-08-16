#' Create a plot showing membrane locations of each protein based on
#' transcript IDs
#'
#' This function creates a ggplot figure showing the differences in membrane
#' location and length between primary and alternative transcripts from the
#' same gene. This process is performed based on input data
#' containing the gene names and transcript IDs of the proteins in question.
#' Transcripts derived from the same gene are grouped together
#' to facilitate easy interpretation. The y axis lists the gene name and
#' transcript ID for each transcript and the x axis lists the length in
#' amino acids. Each fill color corresponds to a membrane location and either
#' principal or alternative isoform.
#'
#' @usage graph_from_ids(data_file, organism = "human", rank = "length",
#' n_prts = 20, mode = "phobius", size_txt = 2, space_left = -400, temp = FALSE,
#' tmhmm_folder_name = NULL)
#' @param data_file Path to the input file
#' @param organism String indicating if the transcripts are from a human or
#' a mouse
#' @param rank String indicating which method to use to rank proteins in graphicl
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
#' @param tmhmm_folder_name Full path to folder containing installed TMHMM 2.0
#' software. This value should end in TMHMM2.0c and needs to be provided if
#' the mode used is TMHMM.
#' @return A ggplot figure showing the protein locations for each part of the
#' surface protein for each alternative and primary transcripts.
#' @examples
#' tmhmm_folder_name <- "~/TMHMM2.0c"
#' if (check_tmhmm_install(tmhmm_folder_name)) {
#'     graph_from_ids(
#'         system.file("extdata", "hpa_example.csv", package = "surfaltr"),
#'         "human", "length", 1, "tmhmm", 5, -300, TRUE
#'     )
#' }
#' @export

graph_from_ids <- function(data_file, organism = "human", rank = "length", 
    n_prts = 20, mode = "phobius", size_txt = 2, space_left = -400, temp = FALSE, 
    tmhmm_folder_name = NULL) {
    old <- getwd()
    if_aa <- FALSE
    final_trans <- clean_data(data_file, if_aa, organism)
    princ <- ensembl_db_retrieval(organism)
    final_pairs <- merge_trans(princ, final_trans, if_aa)
    aa_trans <- format_ids(final_pairs)
    AA_seq <- get_prts(aa_trans, temp)
    if (mode == "phobius") {
        topo <- run_phobius(AA_seq, paste(getwd(), "/AA.fasta", sep = ""))
    }
    if (mode == "tmhmm") {
        topo <- get_tmhmm("AA.fasta", tmhmm_folder_name)
    }
    counts <- process_tmhmm(topo, AA_seq)
    p <- graph_prots(counts, rank, n_prts, size_txt, space_left)
    setwd(old)
    p
}
