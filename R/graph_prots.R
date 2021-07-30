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
#' @usage graph_prots(counts, rank = "length", n_prts = 20, size_txt = 2,
#' space_left = -400)
#' @param counts A data frame containing the overall length and individual lengths of
#' each section of the surface protein corresponding to a certain transcript.
#' @param rank String indicating which method to use to rank proteins in graphicl
#' output. Options include "Length", "TM", and "Combo".
#' @param n_prts Integer value indicating the number of genes that should be
#' displayed on the graphical output. Default value is 20.
#' @param size_txt Integer value specifying the size of the row labels. Default
#' size is 2.
#' @param space_left Integer value specifying how far left the graph should
#' extend.
#' @import ggplot2
#' @return A ggplot figure showing the protein locations for each part of the
#' surface protein for each alternative and primary transcripts.
#'


graph_prots <- function(counts, rank = "length", n_prts = 20, size_txt = 2, space_left = -400){
  counts <- rank_prts(counts, rank, n_prts)
  counts$Order<- as.numeric(counts$Order)
  counts$Begin <- as.numeric(counts$Begin)
  counts$End <- as.numeric(counts$End)
  Begin <- End <- Order <- NULL
  cols <- c()
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0.5, max(counts$Order)+0.5)
  p <- p + ggplot2::xlim(space_left,
                max(counts$End) + max(counts$End)*0.1)
  p <- p + ggplot2::labs(x = "Amino Acid Number")
  p <- p + ggplot2::labs(y = "Transcript ID")
  p <- p + ggplot2::annotate("text", x = -10,
                    y = counts[counts$Membrane_Location != "M" | counts$Membrane_Location != "O" |
                                 counts$Membrane_Location != "o" | counts$Membrane_Location != "i",]$Order,
                    label = paste(counts$Gene_Name, counts$Transcript_ID, sep = "; "),
                    hjust = 1,
                    size = size_txt)
  if(nrow(counts[counts$Membrane_Location == "M" & substr(counts$APPRIS_Annotation, 1, 2) == "pr",]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location == "M" &
                                                substr(counts$APPRIS_Annotation, 1, 2) == "pr",],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "M, Pr",),
                     size = 0.1,)
    cols <- append(cols, "#004754")
  }
  if(nrow(counts[counts$Membrane_Location == "M" & (substr(counts$APPRIS_Annotation, 1, 2) == " _"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "al"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location ==
                                                "M" & (substr(counts$APPRIS_Annotation, 1, 2) == " _"
                                                       | substr(counts$APPRIS_Annotation, 1, 2) == "al"
                                                       | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "M, Alt"),
                     size = 0.1)
    cols <- append(cols, "#540000")
  }
  if(nrow(counts[(counts$Membrane_Location == "O" | counts$Membrane_Location == "o") &
                 substr(counts$APPRIS_Annotation, 1, 2) == "pr",]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[(counts$Membrane_Location == "O" |
                                                 counts$Membrane_Location == "o") &
                                                substr(counts$APPRIS_Annotation, 1, 2) == "pr",],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "O, Pr"),
                     size = 0.1)
    cols <- append(cols, "#ABC3C9")
  }
  if(nrow(counts[(counts$Membrane_Location == "O" | counts$Membrane_Location == "o")
                 & (substr(counts$APPRIS_Annotation, 1, 2) == " _" | substr(counts$APPRIS_Annotation, 1, 2) == "al" | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[(counts$Membrane_Location == "O" |
                                                 counts$Membrane_Location == "o") &
                                                (substr(counts$APPRIS_Annotation, 1, 2) == " _" |
                                                   substr(counts$APPRIS_Annotation, 1, 2) == "al" |
                                                   substr(counts$APPRIS_Annotation, 1, 2) == "NA"),],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "O, Alt"),
                     size = 0.1)
    cols <- append(cols,"#63ACBE")
  }
  if(nrow(counts[counts$Membrane_Location == "i" & substr(counts$APPRIS_Annotation, 1, 2) == "pr",]) != 0) {
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location == "i" &
                                                substr(counts$APPRIS_Annotation, 1, 2) == "pr",],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "i, Pr"),
                     size = 0.1)
    cols <- append(cols, "#B5D99C")
  }
  if(nrow(counts[counts$Membrane_Location == "i" & (substr(counts$APPRIS_Annotation, 1, 2) == " _"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "al"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location == "i" &
                                                (substr(counts$APPRIS_Annotation, 1, 2) == " _" |
                                                   substr(counts$APPRIS_Annotation, 1, 2) == "al" |
                                                   substr(counts$APPRIS_Annotation, 1, 2) == "NA"),],
                     ggplot2::aes(xmin=as.numeric(Begin),
                         xmax=as.numeric(End),
                         ymin=as.numeric(Order)-0.2,
                         ymax=as.numeric(Order)+0.2,
                         fill = "i, Alt"),
                     size = 0.1)
    cols <- append(cols, "#EE442F")
  }
  if(nrow(counts[counts$Membrane_Location == "S" & substr(counts$APPRIS_Annotation, 1, 2) == "pr",]) != 0){
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location == "S" &
                                                substr(counts$APPRIS_Annotation, 1, 2) == "pr",],
                       ggplot2::aes(xmin=as.numeric(Begin),
                           xmax=as.numeric(End),
                           ymin=as.numeric(Order)-0.2,
                           ymax=as.numeric(Order)+0.2,
                           fill = "S, Pr",),
                       size = 0.1,)
    cols <- append(cols, "#FDB261")
  }
  if(nrow(counts[counts$Membrane_Location == "S" & (substr(counts$APPRIS_Annotation, 1, 2) == " _"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "al"
                                                    | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),]) != 0) {
    p <- p + ggplot2::geom_rect(data = counts[counts$Membrane_Location == "S" &
                                                (substr(counts$APPRIS_Annotation, 1, 2) == " _"
                                                 | substr(counts$APPRIS_Annotation, 1, 2) == "al"
                                                 | substr(counts$APPRIS_Annotation, 1, 2) == "NA"),],
                       ggplot2::aes(xmin=as.numeric(Begin),
                           xmax=as.numeric(End),
                           ymin=as.numeric(Order)-0.2,
                           ymax=as.numeric(Order)+0.2,
                           fill = "S, Alt"),
                       size = 0.1)
    cols <- append(cols, "#B3A78C")
  }
  p <- p + ggplot2::theme_light()
  p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p <- p + ggplot2::theme(legend.text=ggplot2::element_text(size=10))
  p <- p + ggplot2::theme(axis.text=ggplot2::element_text(size=15),
                 axis.title=ggplot2::element_text(size=15))
  cols <- as.data.frame(cols)
  colnames(cols) <- "colors"
  order_col = as.data.frame(c("#EE442F","#B5D99C","#540000", "#004754","#63ACBE" ,
                              "#ABC3C9", "#FDB261", "#B3A78C"))
  colnames(order_col) <- "colors"
  order_col$ord <- 1:nrow(order_col)
  col_corr <- merge(x = cols, y = order_col, by = "colors")
  col_corr <- col_corr[order(col_corr$ord),]
  fin_col <- as.vector(col_corr$colors)
  p <- p + ggplot2::scale_fill_manual(name = "", values = fin_col)
  p
  return(p)
}
