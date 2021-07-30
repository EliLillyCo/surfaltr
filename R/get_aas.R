#' Create fasta file containing amino acid sequences based on user sequences
#'
#' This function creates a fasta file with the transcript ID followed by the
#' amino acid sequence for all inputted and associated primary transcripts. The
#' file is organized so that all transcripts from a gene are next to each other.
#' The function also returns a final table containing the gene names, transcript
#' IDs, APPRIS annotations, and amino acid sequences for each transcript
#'
#' @usage get_aas(final_pairs, temp = FALSE)
#' @param final_pairs A data frame containing gene names, transcript IDs, amino
#' acid sequences, and APPRIS annotations for all inputted data and its
#' corresponding primary transcripts.
#' @param temp Boolean indicating if the fasta file should be deleted after the
#' function finishes running or not. Recommended to always be set to FALSE.
#' @importFrom dplyr rename
#' @importFrom seqinr write.fasta
#' @importFrom protr getUniProt
#' @return A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, and protein sequences for each transcript.
#' @note This function also creates a fasta file containing the transcript IDs
#' and associated amino acid sequences in the root directory.


get_aas <- function(final_pairs, temp = FALSE){
  prim_trans <- subset(final_pairs, select = c("external_gene_name", "ensembl_transcript_id.y",
                                               "transcript_appris.y", "uniprotswissprot", "uniprot_isoform",
                                               "uniprotsptrembl"))
  alt_trans <- subset(final_pairs, select = c("external_gene_name", "ensembl_transcript_id.x",
                                              "transcript_appris.x", "protein_sequence"))
  alt_trans <- dplyr::rename(alt_trans, "ensembl_transcript_id" = "ensembl_transcript_id.x")
  alt_trans <- dplyr::rename(alt_trans, "transcript_appris" = "transcript_appris.x")
  prim_trans <- dplyr::rename(prim_trans, "ensembl_transcript_id" = "ensembl_transcript_id.y")
  prim_trans <- dplyr::rename(prim_trans, "transcript_appris" = "transcript_appris.y")
  uniprot_aa <- prim_trans
  for (row in 1:nrow(prim_trans)){
    if(prim_trans[row,"uniprot_isoform"] != ""){
      uniprot_aa[row, "uniprot_id"] = uniprot_aa[row,"uniprot_isoform"]
    }
    else if (prim_trans[row,"uniprotswissprot"] != ""){
      uniprot_aa[row, "uniprot_id"] = uniprot_aa[row,"uniprotswissprot"]
    }
    else {
      uniprot_aa[row, "uniprot_id"] = uniprot_aa[row,"uniprotsptrembl"]
    }
  }
  uniprot_aa <- uniprot_aa[uniprot_aa$uniprot_id != "", ]
  seq <- protr::getUniProt(uniprot_aa$uniprot_id)
  AA_prim_seq <- uniprot_aa
  AA_prim_seq$protein_sequence <- seq
  AA_prim_seq <- subset(AA_prim_seq, select = c("external_gene_name", "ensembl_transcript_id",
                                                "transcript_appris", "protein_sequence"))
  AA_prim_seq <- unique(AA_prim_seq)
  alt_trans <- unique(alt_trans)
  aa_trans <- prim_trans[FALSE,]
  gene_lst <- unique(subset(final_pairs, select = c("external_gene_name")))
  row.names(gene_lst) = NULL
  for(row in 1:(nrow(gene_lst))){
    curr_gene <- gene_lst[row, "external_gene_name"]
    alt_data <- subset(alt_trans, alt_trans$external_gene_name == curr_gene)
    prim_data <- subset(AA_prim_seq, AA_prim_seq$external_gene_name == curr_gene)
    aa_trans <- rbind(aa_trans, alt_data)
    aa_trans <- rbind(aa_trans, prim_data)
  }
  AA_seq <- aa_trans
  AA_seq$id  <- 1:nrow(aa_trans)
  if(temp == TRUE){
    setwd(tempdir())
  }
  if (!(grepl("output", getwd()))){
    dir.create(paste(getwd(), "/output", sep = ""), showWarnings = FALSE)
    setwd(paste(getwd(), "/output", sep = ""))
  }
  seqinr::write.fasta(sequences = as.list(AA_seq$protein_sequence),
                      names = AA_seq$ensembl_transcript_id,
                      file.out = paste(getwd(), "/AA.fasta", sep = ""), open = "w",
                      nbchar = 80, as.string = FALSE)
  rownames(AA_seq) <- NULL
  return(AA_seq)
}
