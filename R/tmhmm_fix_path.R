#' Retrieve Data from TMHMM and Fix Functionality of TMHMM R Package
#'
#' This function retrieves the raw data from tmhmm containing information
#' about the membrane location of each amino acid in a transcript.
#' In order to set a standard path that allows tmhmm to run, the path
#' is set to match that of the fasta file contining the amino acids.
#'
#' @usage tmhmm_fix_path(fasta_filename, folder_name)
#' @param fasta_filename Parameter containing input fasta file to be run on
#' tmhmm
#' @param folder_name Path to folder containing installed tmhmm software
#' @return Raw results from tmhmm containing membrane locations for
#' each transcript
#' @note In order for this function to work, there needs to be a .fasta file
#' containing the amino acid sequences for each transcript called "AA.fasta"
#' saved to a folder called output within the working directory.
#' @importFrom testthat expect_equal
#' @importFrom stringr str_remove

tmhmm_fix_path <- function(fasta_filename, folder_name) {
    bin_path <- normalizePath(file.path(folder_name, "tmhmm-2.0c", "bin", 
                                        "decodeanhmm.Linux_x86_64"))
    options_path <- file.path(
        folder_name, "tmhmm-2.0c", "lib",
        "TMHMM2.0.options"
    )
    model_path <- file.path(
        folder_name, "tmhmm-2.0c", "lib",
        "TMHMM2.0.model"
    )
    cmds <- c("-f", options_path, "-modelfile", model_path)
    text <- NA
    suppressWarnings(text <- system2(
        command = bin_path, args = cmds,
        stdout = TRUE, stderr = NULL, stdin = fasta_filename
    ))
    if (length(text) == 0) {
        if (is.null(attr(text, "status"))) {
            stop("Protein sequence must have at least one character")
        }
        testthat::expect_equal(attr(text, "status"), 100)
        suppressWarnings(text <- system2(
            command = bin_path,
            args = cmds, stdout = NULL, stderr = TRUE, stdin = fasta_filename
        ))
        stop(text[6])
    }
    text <- text[text != ""]
    stringr::str_remove(string = text, pattern = "^\\?0 ")
}
