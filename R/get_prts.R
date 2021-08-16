#' Create fasta file containing amino acid sequences based on IDs
#'
#' This function creates a fasta file with the transcript ID followed by the
#' amino acid sequence for all given alternative transcripts
#' and associated primary transcripts. The file is organized so that all
#' transcripts from a gene are next to each other. The function also returns
#' a final table containing the gene names, transcript IDs, APPRIS annotations,
#' and amino acid sequences for each transcript
#'
#' @usage get_prts(aa_trans, temp = FALSE)
#' @param aa_trans A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, UniProt Swissprot IDs, UniProt Swissprot isoform IDs, and
#' UniProt TREMBL IDs for all transcripts.
#' @param temp Boolean indicating if the fasta file should be deleted after the
#' function finishes running or not. Recommended to always be set to FALSE.
#' @importFrom seqinr write.fasta
#' @importFrom protr getUniProt
#' @return A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, UniProt IDs, and protein sequences for each transcript.
#' @note This function also creates a fasta file containing the transcript IDs
#' and associated amino acid sequences in the root directory.


get_prts <- function(aa_trans, temp = FALSE) {
    uniprot_aa <- aa_trans
    for (row in seq_len(nrow(aa_trans))) {
        if (aa_trans[row, "uniprot_isoform"] != "") {
            uniprot_aa[row, "uniprot_id"] <- uniprot_aa[row, "uniprot_isoform"]
        } else if (aa_trans[row, "uniprotswissprot"] != "") {
            uniprot_aa[row, "uniprot_id"] <- uniprot_aa[row, "uniprotswissprot"]
        } else {
            uniprot_aa[row, "uniprot_id"] <- uniprot_aa[row, "uniprotsptrembl"]
        }
    }
    uniprot_aa <- uniprot_aa[uniprot_aa$uniprot_id != "", ]
    seq <- protr::getUniProt(uniprot_aa$uniprot_id)
    AA_seq <- uniprot_aa
    AA_seq$protein_sequence <- seq
    if (temp == TRUE) {
        setwd(tempdir())
    }
    if (!(grepl("output", getwd()))) {
        dir.create(paste(getwd(), "/output", sep = ""), showWarnings = FALSE)
        setwd(paste(getwd(), "/output", sep = ""))
    }
    seqinr::write.fasta(
        sequences = as.list(AA_seq$protein_sequence),
        names = AA_seq$ensembl_transcript_id,
        file.out = paste(getwd(), "/AA.fasta", sep = ""), open = "w",
        nbchar = 80, as.string = FALSE
    )
    rownames(AA_seq) <- NULL
    return(AA_seq)
}
