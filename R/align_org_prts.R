#' Get aligned amino acid sequences for gene transcripts from multiple organisms
#'
#' This function allows a user to specify genes of interest and subsequently
#' receive a pdf of all the corresponding aligned amino acid sequences in pdf
#' format for both human and mice transcripts. In order for this to work,
#' transcripts for the same genes from both organisms need to be
#' provided in separate files.
#'
#' @usage align_org_prts(gene_names, hs_data_file, mm_data_file, if_aa = FALSE, temp = FALSE)
#' @param gene_names Vector containing names of genes of interest
#' (e.g. c("Crb1", "Adgrl1"))
#' @param hs_data_file Path to the input file containing the human transcripts
#' @param mm_data_file Path to the input file containing the mouse transcripts
#' @param if_aa Boolean value indicating if the input file contains
#' amino acid sequence. TRUE indicates that sequences are present and
#' FALSE indicates that only IDs are present.
#' @param temp Boolean indicating if the fasta file should be saved to the
#' working directory or no
#' @return Nothing is returned.
#' @note Although the function returns nothing, it saves pdfs containing the
#' aligned sequences to the working directory under a file labeled with the
#' gene name. It's also important to note that although the gene names will
#' be standardized to be fully capitalized, this may not match with the case
#' of the gene name for some organisms.
#' @import msa
#' @importFrom seqinr write.fasta
#' @importFrom Biostrings readAAStringSet
#' @export
#' @examples
#' align_org_prts(c("IGSF1", "TAPBP"),
#' system.file("extdata", "hpa_genes.csv", package = "surfaltr"),
#' system.file("extdata", "hpa_mouse_genes.csv", package = "surfaltr"),
#' FALSE, TRUE)

align_org_prts <- function(gene_names, hs_data_file, mm_data_file, if_aa = FALSE, temp = FALSE){
  old <- getwd()
  if (if_aa == FALSE){
    hs_final_trans <- clean_data(hs_data_file, if_aa, "human")
    hs_princ <- ensembl_db_retrieval("human")
    hs_final_pairs <- merge_trans(hs_princ, hs_final_trans, if_aa)
    hs_aa_trans <- format_ids(hs_final_pairs)
    hs_AA_seq <- get_prts(hs_aa_trans, temp)
    mm_final_trans <- clean_data(mm_data_file, if_aa, "mouse")
    mm_princ <- ensembl_db_retrieval("mouse")
    mm_final_pairs <- merge_trans(mm_princ, mm_final_trans, if_aa)
    mm_aa_trans <- format_ids(mm_final_pairs)
    mm_AA_seq <- get_prts(mm_aa_trans, temp)
    AA_seq <- rbind(hs_AA_seq, mm_AA_seq)
    AA_seq$external_gene_name <- toupper(AA_seq$external_gene_name)
  } else {
    hs_final_trans <- clean_data(hs_data_file, if_aa, "human")
    hs_princ <- ensembl_db_retrieval("human")
    hs_final_pairs <- merge_trans(hs_princ, hs_final_trans, if_aa)
    hs_AA_seq <- get_aas(hs_final_pairs, temp)
    mm_final_trans <- clean_data(mm_data_file, if_aa, "mouse")
    mm_princ <- ensembl_db_retrieval("mouse")
    mm_final_pairs <- merge_trans(mm_princ, mm_final_trans, if_aa)
    mm_AA_seq <- get_aas(mm_final_pairs, temp)
    AA_seq <- rbind(hs_AA_seq, mm_AA_seq)
    AA_seq$external_gene_name <- toupper(AA_seq$external_gene_name)
  }
  gene_lst <- data.frame(gene_names)
  colnames(gene_lst) <- c("Gene_Name")
  rel_genes <- merge(x = gene_lst, y = AA_seq, by.x = "Gene_Name", by.y = "external_gene_name")
  gene_lst <- data.frame(unique(rel_genes$Gene_Name))
  for(row in 1:nrow(gene_lst)){
    curr_genes <- subset(rel_genes, rel_genes$Gene_Name == gene_lst[row,])
    curr_genes <- curr_genes[order(curr_genes$id),]
    temp_file_name <- paste(gene_lst[row,],"_multi_org.fasta", sep = "")
    seqinr::write.fasta(sequences = as.list(curr_genes$protein_sequence),
                        names = paste(curr_genes$Gene_Name, curr_genes$ensembl_transcript_id, sep = ", "),
                        file.out = temp_file_name, open = "w", nbchar = 80, as.string = FALSE)
    my_sequences <- Biostrings::readAAStringSet(temp_file_name)
    align_seq <- msa(my_sequences)
    pdf_name <- paste(gene_lst[row,],"_multi_org.pdf", sep = "")
    msa::msaPrettyPrint(align_seq, output="pdf", showNames="left",
                        showLogo="none", askForOverwrite=FALSE, verbose=FALSE, file = pdf_name)
  }
  setwd(old)
}
