#' Check to make sure TMHMM 2.0 is installed in the file path specified
#'
#' This function checks to make sure that TMHMM is installed correctly
#' at the file path specified by the user. If TMHMM is not installed correctly,
#' then the function will output an error message telling the user to check 
#' their installation.
#'
#' @usage check_tmhmm_install(tmhmm_folder_name)
#' @param tmhmm_folder_name Full path to folder containing installed TMHMM 2.0
#' software. This value should end in TMHMM2.0c
#' @return A Boolean stating if TMHMM is installed correctly, will be TRUE
#' if TMHMM 2.0 is located at the path specified and FALSE if it is not.
#' @note This function also prints a helpful method providing tips on how
#' to fix the installation if TMHMM is not found at the folder path specified.
#' @importFrom utils write.csv
#' @importFrom utils assignInNamespace
#' @examples
#' tmhmm_folder_name <- "~/TMHMM2.0c"
#' install_correct <- check_tmhmm_install(tmhmm_folder_name)
#' @export

check_tmhmm_install <- function(tmhmm_folder_name) {
    final_path <- paste(tmhmm_folder_name, "/bin/tmhmm", sep = "")
    if (!(file.exists(final_path))) {
        print("Please check and make sure that TMHMM 2.0 is installed at the 
              file path specified. To install TMHMM, please go to 
              https://services.healthtech.dtu.dk/service.php?TMHMM-2.0 and 
              follow the instructions.")
        return(FALSE)
    } else {
        return(TRUE)
    }
}
