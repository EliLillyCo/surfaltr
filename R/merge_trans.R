#' Associate Inputted Transcripts with Corresponding Primary Transcripts
#'
#' This function matches each inputted transcript with its corresponding
#' primary transcripts and returns a data frame containing the gene name,
#' transcript ID and APPRIS annotation for each.
#'
#' @usage merge_trans(princ, final_trans, if_aa)
#' @param princ Data frame containing all primary transcripts and relevant
#' gene information for an organism
#' @param final_trans Data frame containing cleaned and formatted input data
#' @param if_aa Boolean value indicating if the input file contains
#' amino acid sequences with TRUE indicating that sequences are present and
#' FALSE indicating that only IDs are present
#' @return A data frame containing gene names, transcript IDs, and APPRIS
#' annotations for all inputted data and its corresponding primary transcripts.
#' If sequences were provided, the data frame will also contain
#' the amino acid sequences. If only IDs were provided, the data frame will
#' also contain the UniProt Swissprot ID, UniProt Swissprot isoform ID, and
#' UniProt TREMBL ID for both the inputted data and the primary transcripts.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter


merge_trans <- function(princ, final_trans, if_aa) {
    prim_alt_lst <- princ %>% dplyr::filter(princ$external_gene_name %in% 
    final_trans$external_gene_name)
    final_pairs <- merge(x = final_trans, y = prim_alt_lst, 
    by.x = "external_gene_name", by.y = "external_gene_name")
    if (if_aa == FALSE) {
        final_pairs <- dplyr::filter(final_pairs, 
            final_pairs$uniprotswissprot.x != "" |
            final_pairs$uniprot_isoform.x != "" |
            final_pairs$uniprotsptrembl.x != "")
        final_pairs <- dplyr::filter(final_pairs, 
            final_pairs$uniprotswissprot.y != "" |
            final_pairs$uniprot_isoform.y != "" |
            final_pairs$uniprotsptrembl.y != "")
        if (nrow(final_pairs) == 0) {
            stop("Please make sure your organism is correct and your transcript 
            IDs are valid and contain no decimal values before continuing.")
        }
    }
    if (if_aa == TRUE) {
        final_pairs <- dplyr::filter(final_pairs, final_pairs$uniprotswissprot != "" |
            final_pairs$uniprot_isoform != "" | final_pairs$uniprotsptrembl != "")
        if (nrow(final_pairs) == 0) {
            stop("Please make sure your organism is correct and your transcript 
            IDs are valid and contain no decimal values before continuing.")
        }
    }
    return(final_pairs)
}
