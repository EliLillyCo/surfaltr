#' Reformat transcripts to facilitate fasta file conversion
#'
#' Modify format of data to display all primary and alternative transcripts
#' from the same gene together and remove any duplicates.
#'
#' @usage format_ids(final_pairs)
#' @param final_pairs Data frame containing original row-wise pairings of
#' primary and alternative transcripts for inputted data without associated
#' sequences
#' @importFrom dplyr rename
#' @return A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, UniProt Swissprot IDs, UniProt Swissprot isoform IDs, and
#' UniProt TREMBL IDs for all the given and associated primary transcripts in
#' an alternating fashion


format_ids <- function(final_pairs) {
    prim_trans <- subset(final_pairs, select = c(
        "external_gene_name", "ensembl_transcript_id.y",
        "transcript_appris.y", "uniprotswissprot.y",
        "uniprot_isoform.y", "uniprotsptrembl.y"
    ))
    alt_trans <- subset(final_pairs, select = c(
        "external_gene_name", "ensembl_transcript_id.x",
        "transcript_appris.x", "uniprotswissprot.x",
        "uniprot_isoform.x", "uniprotsptrembl.x"
    ))
    alt_trans <- dplyr::rename(alt_trans, "ensembl_transcript_id" = 
                 "ensembl_transcript_id.x")
    alt_trans <- dplyr::rename(alt_trans, "transcript_appris" = 
                 "transcript_appris.x")
    alt_trans <- dplyr::rename(alt_trans, "uniprot_isoform" = 
                 "uniprot_isoform.x")
    alt_trans <- dplyr::rename(alt_trans, "uniprotsptrembl" = 
                 "uniprotsptrembl.x")
    alt_trans <- dplyr::rename(alt_trans, "uniprotswissprot" = 
                 "uniprotswissprot.x")
    prim_trans <- dplyr::rename(prim_trans, "ensembl_transcript_id" = 
                 "ensembl_transcript_id.y")
    prim_trans <- dplyr::rename(prim_trans, "transcript_appris" = 
                 "transcript_appris.y")
    prim_trans <- dplyr::rename(prim_trans, "uniprot_isoform" = 
                 "uniprot_isoform.y")
    prim_trans <- dplyr::rename(prim_trans, "uniprotsptrembl" = 
                 "uniprotsptrembl.y")
    prim_trans <- dplyr::rename(prim_trans, "uniprotswissprot" = 
                 "uniprotswissprot.y")
    prim_trans <- unique(prim_trans)
    alt_trans <- unique(alt_trans)
    aa_trans <- prim_trans[FALSE, ]
    gene_lst <- unique(subset(final_pairs, select = c("external_gene_name")))
    row.names(gene_lst) <- NULL
    for (row in seq_len(nrow(gene_lst))) {
    curr_gene <- gene_lst[row, "external_gene_name"]
    alt_data <- subset(alt_trans, alt_trans$external_gene_name == curr_gene)
    prim_data <- subset(prim_trans, prim_trans$external_gene_name == curr_gene)
    aa_trans <- rbind(aa_trans, alt_data)
    aa_trans <- rbind(aa_trans, prim_data)
    }
    aa_trans <- unique(aa_trans)
    aa_trans$id <- seq_len(nrow(aa_trans))
    return(aa_trans)
}
