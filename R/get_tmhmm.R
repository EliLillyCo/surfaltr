#' Create a data frame with the membrane locations of each amino acid in a
#' protein using TMHMM
#'
#' This function creates a data frame with columns containing transcript
#' IDs and corresponding output from TMHMM. The TMHMM output includes a location
#' for each amino acid, with O and o representing extracellular, M representing
#' transmembrane, and i representing intracellular.
#'
#' @usage get_tmhmm(fasta_file_name, tmhmm_folder_name)
#' @param fasta_file_name Name of .fasta file containing amino acid
#' sequences
#' @param tmhmm_folder_name Full path to folder containing installed TMHMM 2.0
#' software. This path should end in TMHMM2.0c
#' @return A data frame containing each transcript ID and the corresponding
#' membrane location for each amino acid in its sequence formatted as a string
#' @note In order for this function to work, there needs to be a .fasta file
#' containing the amino acid sequences for each transcript called "AA.fasta"
#' saved to a folder called output within the working directory.
#' Additionally, the file saves a copy of the
#' returned data frame in csv format to the output folder in the working
#' directory.
#' @importFrom utils write.csv
#' @import stringr
#' @examples
#' tmhmm_folder_name <- "~/TMHMM2.0c"
#' if (check_tmhmm_install(tmhmm_folder_name)) {
#'     AA_seq <- get_pairs(system.file("extdata", "crb1_example.csv",
#'         package = "surfaltr"
#'     ), TRUE, "mouse", TRUE)
#'     topo <- get_tmhmm("AA.fasta", tmhmm_folder_name)
#' }
#' @export


get_tmhmm <- function(fasta_file_name, tmhmm_folder_name) {
    if (!(check_tmhmm_install(tmhmm_folder_name))) {
        stop("Please check your TMHMM installation and the provided path.")
    }
    topo_fasta <- tmhmm_fix_path(fasta_file_name, tmhmm_folder_name)
    topo_fasta <- data.frame(topo_fasta)
    topo_fasta_clean <- data.frame("data" = character(0))
    for (row in seq_len(nrow(topo_fasta))) {
        curr_sub <- substr(topo_fasta[row, ], 1, 1)
        if (curr_sub == ">" | curr_sub == "O" | curr_sub == "i" | curr_sub == "M" ) {
            topo_fasta_clean[nrow(topo_fasta_clean) + 1, "data"] <- substr(topo_fasta[row, ], 
                                                            1, nchar(topo_fasta[row, ]))
        }
    }
    topo <- data.frame("Transcript_ID" = character(0), "Output" = character(0))
    id_val <- 0
    for (row in seq_len(nrow(topo_fasta_clean))) {
        curr_sub <- substr(topo_fasta_clean[row, ], 1, 1)
        if (curr_sub == ">") {
            topo[nrow(topo) + 1, "Transcript_ID"] <- substr(topo_fasta_clean[row, ], 
                                                            2, nchar(topo_fasta_clean[row, ]))
            id_val <- id_val + 1
        } else {
            if (is.na(nchar(topo[id_val, "Output"], keepNA = TRUE))) {
                topo[id_val, "Output"] <- topo_fasta_clean[row, ]
            } else {
                topo[id_val, "Output"] <- paste(topo[id_val, "Output"], 
                                                topo_fasta_clean[row, ], sep = "")
            }
        }
    }
    write.csv(topo, paste(getwd(), "/mem_topo_tmhmm.csv", sep = ""))
    return(topo)
}
