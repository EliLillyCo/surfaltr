#' Retrieve Transcript Information from Ensembl for all Primary Transcripts
#'
#' This function retrieves all the primary transcripts in the given organism
#' and their corresponding gene names, APPRIS annotations, and UniProt IDs.
#'
#' @usage ensembl_db_retrieval(organism)
#' @param organism String indicating if mouse or human transcripts should be
#' retrieved
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt getBM
#' @importFrom httr set_config
#' @importFrom httr config
#' @return A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, UniProt Swissprot IDs, UniProt Swissprot isoform IDs, and
#' UniProt TREMBL IDs for all the primary transcripts in an organism.

ensembl_db_retrieval <- function(organism){
  if (organism == "human"){
    httr::set_config(config(ssl_verifypeer = 0L))
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
    humandb <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                                             "ensembl_transcript_id", "transcript_appris",
                                             "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"),
                              mart = ensembl)
    princ <- humandb[ which(humandb$transcript_appris == "principal1" |
                              humandb$transcript_appris == "principal2" |
                              humandb$transcript_appris == "principal3" |
                              humandb$transcript_appris == "principal4" |
                              humandb$transcript_appris == "principal5"), ]
  }

  if (organism == "mouse"){
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
    mousedb <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                                             "ensembl_transcript_id", "transcript_appris",
                                             "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"),
                              mart = ensembl)
    princ <- mousedb[ which(mousedb$transcript_appris == "principal1" |
                              mousedb$transcript_appris == "principal2" |
                              mousedb$transcript_appris == "principal3" |
                              mousedb$transcript_appris == "principal4" |
                              mousedb$transcript_appris == "principal5"), ]
  }
  return(princ)
}
