#' Get aligned amino acid sequences for gene transcripts
#'
#' This function allows a user to specify genes of interest and subsequently
#' receive a pdf of all the corresponding aligned amino acid sequences in pdf
#' format.
#'
#' @usage align_prts(gene_names, data_file, if_aa = FALSE, organism = "human", temp = FALSE)
#' @param gene_names Vector containing names of genes of interest
#' (e.g. c(Crb1, Adgrl1))
#' @param data_file Path to the input file
#' @param if_aa Boolean value indicating if the input file contains
#' amino acid sequence. TRUE indicates that sequences are present and
#' FALSE indicates that only IDs are present.
#' @param organism String indicating if the transcripts are from a human or
#' a mouse
#' @param temp Boolean indicating if the fasta file should be saved to the
#' working directory or no
#' @return Nothing is returned.
#' @note Although the function returns nothing, it saves pdfs containing the
#' aligned sequences to the working directory under a file labeled with the
#' gene name.
#' @import msa
#' @importFrom seqinr write.fasta
#' @importFrom Biostrings readAAStringSet
#' @examples align_prts(c("Crb1"), system.file("extdata", "crb1_example.csv",
#' package = "surfaltr"), TRUE, "mouse", TRUE)
#' @export

align_prts <- function(gene_names, data_file, if_aa = FALSE, organism = "human", temp = FALSE){
  old <- getwd()
  if (if_aa == FALSE){
    final_trans <- clean_data(data_file, if_aa, organism)
    princ <- ensembl_db_retrieval(organism)
    final_pairs <- merge_trans(princ, final_trans, if_aa)
    aa_trans <- format_ids(final_pairs)
    AA_seq <- get_prts(aa_trans, temp)
  } else {
    final_trans <- clean_data(data_file, if_aa, organism)
    princ <- ensembl_db_retrieval(organism)
    final_pairs <- merge_trans(princ, final_trans, if_aa)
    AA_seq <- get_aas(final_pairs, temp)
  }
  gene_lst <- data.frame(gene_names)
  colnames(gene_lst) <- c("Gene_Name")
  rel_genes <- merge(x = gene_lst, y = AA_seq, by.x = "Gene_Name", by.y = "external_gene_name")
  gene_lst <- data.frame(unique(rel_genes$Gene_Name))
  for(row in 1:nrow(gene_lst)){
    curr_genes <- subset(rel_genes, rel_genes$Gene_Name == gene_lst[row,])
    temp_file_name <- paste(gene_lst[row,],".fasta", sep = "")
    temp_file_path <- paste(getwd(), "/", temp_file_name,  sep = "")
    seqinr::write.fasta(sequences = as.list(curr_genes$protein_sequence),
                        names = paste(curr_genes$Gene_Name, curr_genes$ensembl_transcript_id, sep = ", "),
                        file.out = temp_file_path, open = "w", nbchar = 80, as.string = FALSE)
    my_sequences <- Biostrings::readAAStringSet(temp_file_path)
    align_seq <- msa(my_sequences)
    pdf_name <- paste(gene_lst[row,],".pdf", sep = "")
    pdf_path <- paste(getwd(), "/", pdf_name,  sep = "")
    msa::msaPrettyPrint(align_seq, output="pdf", showNames="left",
                   showLogo="none", askForOverwrite=FALSE, verbose=FALSE, file = pdf_path)
  }
  setwd(old)
}
