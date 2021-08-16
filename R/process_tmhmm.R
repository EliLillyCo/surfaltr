#' Create a data frame with the membrane locations of each amino acid in a
#' sequence
#'
#' This function creates a data frame with columns containing transcript
#' IDs and corresponding output from tmhmm. The tmhmm output includes a location
#' for each amino acid, with O and o representing extracellular, M representing
#' transmembrane, and i representing intracellular. The data frame includes
#' columns with the transcript ID, membrane location, gene name, starting amino
#' acid, and ending amino acid for a certain transcript. The first row for each
#' transcript contains the overall length of the amino acid sequence.
#'
#' @usage process_tmhmm(topo, AA_seq)
#' @param topo A data frame containing each transcript ID and the corresponding
#' membrane location for each amino acid in its sequence formatted as a string.
#' @param AA_seq A data frame containing the gene names, transcript IDs, APPRIS
#' annotations, and protein sequences for each transcript.
#' @importFrom dplyr rename
#' @importFrom stringr str_split
#' @return A data frame containing the overall length and individual lengths of
#' each section of the surface protein corresponding to a certain transcript.


process_tmhmm <- function(topo, AA_seq) {
    counts <- data.frame("Membrane_Location" = character(0), "Count" = integer(0), "Transcript_ID" = character(0))
    comma_data <- topo
    for (row in seq_len(nrow(comma_data))) {
        comma_data[row, "Output"] <- gsub("Si", "S,i", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("So", "S,o", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("SO", "S,O", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("Mo", "M,o", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("oM", "o,M", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("MO", "M,O", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("OM", "O,M", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("iM", "i,M", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("Mi", "M,i", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("Oi", "O,i", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("iO", "i,O", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("oi", "o,i", comma_data[row, "Output"])
        comma_data[row, "Output"] <- gsub("io", "i,o", comma_data[row, "Output"])
        split_let <- unlist(stringr::str_split(comma_data[row, "Output"], ","))
        let_counts <- data.frame(split_let)
        let_counts <- dplyr::rename(let_counts, "Membrane_Location" = "split_let")
        let_counts$counts <- c(nchar(split_let))
        let_counts <- dplyr::rename(let_counts, "Count" = "counts")
        for (val in seq_len(nrow(let_counts))) {
            let_counts[val, "Membrane_Location"] <- substr(let_counts[val, 
            "Membrane_Location"], 1, 1)
            let_counts[val, "Transcript_ID"] <- topo[row, "Transcript_ID"]
            let_counts[val, "Gene_Name"] <- AA_seq[row, "external_gene_name"]
            let_counts[val, "APPRIS_Annotation"] <- paste(AA_seq[row, 
            "transcript_appris"], "_")
            let_counts[val, "Order"] <- row
            if (val == 1) {
                let_counts[val, "Begin"] <- 1
                let_counts[val, "End"] <- let_counts[val, "Count"]
            } else {
                let_counts[val, "Begin"] <- let_counts[val - 1, "End"]
                let_counts[val, "End"] <- let_counts[val - 1, "End"] + 
                let_counts[val, "Count"]
            }
        }
        let_counts <- rbind(
            c(
                topo[row, "Transcript_ID"], 100, topo[row, "Transcript_ID"],
                AA_seq[row, "external_gene_name"], "", row, 1, 
                let_counts[nrow(let_counts), "End"]
            ),
            let_counts
        )
        counts <- rbind(counts, let_counts)
    }
    counts <- subset(counts, select = c(
        "Gene_Name", "Transcript_ID", "Membrane_Location",
        "Begin", "End", "Order", "APPRIS_Annotation"
    ))
    return(counts)
}
