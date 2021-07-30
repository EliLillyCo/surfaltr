#' Rank the surface proteins by differences in principal and alternative isoforms
#'
#' This function creates a data frame containing the primary and alternative
#' transcripts of each gene ranked by how different the resultant surface proteins
#' are. Transcripts can be ranked by length, number of transmembrane domains,
#' or a combo metric that multiplied the difference in length by the number of
#' transmembrane domains and ranks accordingly. This function can also be set
#' to restrict the number of genes that are returned to the user to show only
#' the most significant gene transcripts.
#'
#' @usage rank_prts(counts, rank, n_prts)
#' @param counts A data frame containing the overall length and individual lengths of
#' each section of the surface protein corresponding to a certain transcript.
#' @param rank String indicating which method to use to rank proteins in graphicl
#' output. Options include "Length", "TM", and "Combo".
#' @param n_prts Integer value indicating the number of genes that should be
#' displayed on the graphical output. Default value is 20.
#' @return A data frame containing the overall length and individual lengths of
#' each section of the surface protein corresponding to a certain transcript
#' ranked by how different the primary and alternative transcripts are
#' functionally.
#' @importFrom dplyr distinct
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr rename



rank_prts <- function(counts, rank, n_prts){
  ranks <- counts
  ranks$Position <- 1:nrow(ranks)
  genes <- as.data.frame(dplyr::distinct(ranks, ranks$Gene_Name))
  if(nrow(genes) == 1){
    return(ranks)
  }
  if (rank == "length"){
    for(row in 1:nrow(ranks)){
      if(substr(ranks[row, "Membrane_Location"], 1, 1) == "E"){
        curr_length = as.numeric(ranks[row, "End"])
      }
      ranks[row, "Length"] = curr_length
    }
    genes <- as.data.frame(unique(ranks$Gene_Name))
    len_genes <- genes
    len_genes$Position <- 1:nrow(genes)
    len_genes <- dplyr::rename(len_genes, "Gene_Name" = "unique(ranks$Gene_Name)")
    for(row in 1:nrow(genes)){
      curr_gene <- genes[row,]
      curr_gene_trans <- subset(ranks, ranks$Gene_Name == curr_gene)
      princ_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "p" )
      alt_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "a")
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "N"))
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == " "))
      princ_length <- mean(princ_vals$Length)
      alt_length_low <- min(alt_vals$Length)
      alt_length_high <- max(alt_vals$Length)
      len_diff_low <- abs(princ_length-alt_length_low)
      len_diff_high <- abs(princ_length-alt_length_high)
      len_diff <- max(len_diff_high,len_diff_low)
      for (count in 1:nrow(ranks)){
        if(ranks[count, "Gene_Name"] == curr_gene){
          ranks[count,"Length_Diff"] <-len_diff
        }
      }
      len_genes[row, "Length_Diff"] <- len_diff
    }
    ranks <- ranks[with(ranks, order(-Length_Diff, Position)),]
    len_genes <- len_genes[with(len_genes, order(-Length_Diff, Position)),]
    if(nrow(len_genes) < n_prts){
      n_prts <- nrow(len_genes)
    }
    len_genes <- len_genes[1:n_prts,]
    len_genes <- subset(len_genes, select = c("Gene_Name"))
    final_ranks <- merge( x= ranks, y = len_genes, by = "Gene_Name")
    final_ranks <- final_ranks[with(final_ranks, order(Length_Diff, Position)),]
    ord_num <- 0
    for(row in 1:nrow(final_ranks)){
      if(substr(final_ranks[row, "Membrane_Location"], 1, 1) == "E"){
        ord_num <- ord_num + 1
      }
      final_ranks[row,"Order"] <- ord_num
    }
  }
  if(rank == "TM"){
    genes <- as.data.frame(unique(ranks$Gene_Name))
    len_genes <- genes
    len_genes$Position <- 1:nrow(genes)
    len_genes <- dplyr::rename(len_genes, "Gene_Name" = "unique(ranks$Gene_Name)")
    pos <- 0
    for(row in 1:nrow(genes)){
      pos <- pos + 1
      curr_gene <- genes[row,]
      curr_gene_trans <- subset(ranks, ranks$Gene_Name == curr_gene)
      p_tm_count <- 0
      a_tm_count <- 0
      princ_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "p" )
      alt_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "a")
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "N"))
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == " "))
      alt_vals <- alt_vals %>% dplyr::filter(alt_vals$Membrane_Location == "M")
      princ_vals <- princ_vals %>% dplyr::filter(princ_vals$Membrane_Location == "M")
      p_trans <- as.data.frame(unique(princ_vals$Transcript_ID))
      a_trans <- as.data.frame(unique(alt_vals$Transcript_ID))
      if (nrow(princ_vals) == 0 | nrow(p_trans) == 0){
        p_tm_count <- 0
      } else {
        p_tm_count <- nrow(princ_vals)/nrow(p_trans)
      }
      a_counts <- a_trans
      for(row in 1:nrow(a_trans)){
        tx_curr <- subset(alt_vals, alt_vals$Transcript_ID == a_trans[row,])
        a_tm_count <- nrow(tx_curr)
        a_counts[row, "TM_Domains"] <- a_tm_count
      }
      max_a <- max(a_counts$TM_Domains)
      min_a <- min(a_counts$TM_Domains)
      top_diff <- max_a - p_tm_count
      bottom_diff <- min_a - p_tm_count
      if (top_diff <= 0 && bottom_diff <= 0){
        tm_diff <- bottom_diff
      } else {
        tm_diff <- top_diff
      }
      for (count in 1:nrow(ranks)){
        if(ranks[count, "Gene_Name"] == curr_gene){
          ranks[count,"TM_Domains_Diff"] <-tm_diff
        }
      }
      len_genes[pos, "TM_Domains_Diff"] <- tm_diff
    }
    if(n_prts > nrow(len_genes)){
      n_prts <- nrow(len_genes)
    }
    ranks <- ranks[with(ranks, order(-TM_Domains_Diff, Position)),]
    len_genes <- len_genes[with(len_genes, order(-TM_Domains_Diff, Position)),]
    len_genes_great <- subset(len_genes, len_genes$TM_Domains_Diff > 0)
    len_genes_less <- subset(len_genes, len_genes$TM_Domains_Diff <= 0)
    len_genes_less <- len_genes_less[with(len_genes_less,
                                          order(-abs(len_genes_less$TM_Domains_Diff), Position)),]
    len_genes <- rbind(len_genes_great, len_genes_less)
    len_genes <- len_genes[1:n_prts,]
    len_genes <- subset(len_genes, select = c("Gene_Name"))
    final_ranks <- merge( x= ranks, y = len_genes, by = "Gene_Name")
    great_zero <- subset(final_ranks, final_ranks$TM_Domains_Diff > 0)
    less_zero <- subset(final_ranks, final_ranks$TM_Domains_Diff <= 0)
    great_zero <- great_zero[with(great_zero, order(great_zero$TM_Domains_Diff, Position)),]
    less_zero <- less_zero[with(less_zero, order(abs(less_zero$TM_Domains_Diff), Position)),]
    final_ranks <- rbind(less_zero, great_zero)
    ord_num <- 0
    for(nums in 1:nrow(final_ranks)){
      if(substr(final_ranks[nums, "Membrane_Location"], 1, 1) == "E"){
        ord_num <- ord_num + 1
      }
      final_ranks[nums,"Order"] <- ord_num
    }
  }
  if(rank == "combo"){
    pos <- 0
    for(row in 1:nrow(ranks)){
      if(substr(ranks[row, "Membrane_Location"], 1, 1) == "E"){
        curr_length = as.numeric(ranks[row, "End"])
      }
      ranks[row, "Length"] = curr_length
    }
    genes <- as.data.frame(unique(ranks$Gene_Name))
    len_genes <- genes
    len_genes$Position <- 1:nrow(genes)
    len_genes <- dplyr::rename(len_genes, "Gene_Name" = "unique(ranks$Gene_Name)")
    for(row in 1:nrow(genes)){
      curr_gene <- genes[row,]
      curr_gene_trans <- subset(ranks, ranks$Gene_Name == curr_gene)
      princ_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "p" )
      alt_vals <- subset(curr_gene_trans, substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "a")
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == "N"))
      alt_vals <- rbind(alt_vals, subset(curr_gene_trans,
                                         substr(curr_gene_trans$APPRIS_Annotation, 1, 1) == " "))
      princ_length <- mean(princ_vals$Length)
      alt_length_low <- min(alt_vals$Length)
      alt_length_high <- max(alt_vals$Length)
      len_diff_low <- abs(princ_length-alt_length_low)
      len_diff_high <- abs(princ_length-alt_length_high)
      len_diff <- max(len_diff_high,len_diff_low)
      pos <- pos + 1
      p_tm_count <- 0
      a_tm_count <- 0
      alt_vals <- alt_vals %>% dplyr::filter(alt_vals$Membrane_Location == "M")
      princ_vals <- princ_vals %>% dplyr::filter(princ_vals$Membrane_Location == "M")
      p_trans <- as.data.frame(unique(princ_vals$Transcript_ID))
      a_trans <- as.data.frame(unique(alt_vals$Transcript_ID))
      if (nrow(princ_vals) == 0 | nrow(p_trans) == 0){
        p_tm_count <- 0
      } else {
        p_tm_count <- nrow(princ_vals)/nrow(p_trans)
      }
      a_counts <- a_trans
      for(row in 1:nrow(a_trans)){
        tx_curr <- subset(alt_vals, alt_vals$Transcript_ID == a_trans[row,])
        a_tm_count <- nrow(tx_curr)
        a_counts[row, "TM_Domains"] <- a_tm_count
      }
      max_a <- max(a_counts$TM_Domains)
      min_a <- min(a_counts$TM_Domains)
      top_diff <- max_a - p_tm_count
      bottom_diff <- min_a - p_tm_count
      if (top_diff <= 0 && bottom_diff <= 0){
        tm_diff <- bottom_diff
      } else {
        tm_diff <- top_diff
      }
      combo_diff <- len_diff * tm_diff
      if(combo_diff == 0){
        combo_diff <- -len_diff/100
      }
      for (count in 1:nrow(ranks)){
        if(ranks[count, "Gene_Name"] == curr_gene){
          ranks[count,"Combo_Value"] <-combo_diff
        }
      }
      len_genes[pos, "Combo_Value"] <- combo_diff
    }
    if(n_prts > nrow(len_genes)){
      n_prts <- nrow(len_genes)
    }
    ranks <- ranks[with(ranks, order(-Combo_Value, Position)),]
    len_genes <- len_genes[with(len_genes, order(-Combo_Value, Position)),]
    len_genes_great <- subset(len_genes, len_genes$Combo_Value > 0)
    len_genes_less <- subset(len_genes, len_genes$Combo_Value <= 0)
    len_genes_less <- len_genes_less[with(len_genes_less,
                                          order(-abs(len_genes_less$Combo_Value), Position)),]
    len_genes <- rbind(len_genes_great, len_genes_less)
    len_genes <- len_genes[1:n_prts,]
    len_genes <- subset(len_genes, select = c("Gene_Name"))
    final_ranks <- merge( x= ranks, y = len_genes, by = "Gene_Name")
    great_zero <- subset(final_ranks, final_ranks$Combo_Value > 0)
    less_zero <- subset(final_ranks, final_ranks$Combo_Value <= 0)
    great_zero <- great_zero[with(great_zero, order(great_zero$Combo_Value, Position)),]
    less_zero <- less_zero[with(less_zero, order(abs(less_zero$Combo_Value), Position)),]
    final_ranks <- rbind(less_zero, great_zero)
    ord_num <- 0
    for(nums in 1:nrow(final_ranks)){
      if(substr(final_ranks[nums, "Membrane_Location"], 1, 1) == "E"){
        ord_num <- ord_num + 1
      }
      final_ranks[nums,"Order"] <- ord_num
    }
  }
  return(final_ranks)
}
