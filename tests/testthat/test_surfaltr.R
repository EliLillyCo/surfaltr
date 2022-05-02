#' Test the functionality of surfaltr
#'
#' This function runs all of surfaltr's other functions on the CRB1 data set
#' to ensure that the function output matches the expected output. An incorrect
#' output or error indicates that something went wrong in installation.
#'
#' @usage test_surfaltr()
#' @return Nothing is returned.
#' @note If the results from the test match the expected
#' results, a message stating that the test worked will be printed. If not, the
#' user will be prompted to check the installation
#' @examples
#' \donttest{
#' test_surfaltr()
#' }
#' @importFrom dplyr all_equal
#' @export

test_surfaltr <- function() {
    old <- getwd()
    AA_seq <- get_pairs(system.file("extdata", "CRB1.csv", package = "surfaltr"), TRUE, "mouse", TRUE)
    topo <- run_phobius(AA_seq, paste(getwd(), "/AA.fasta", sep = ""))
    counts <- process_tmhmm(topo, AA_seq)
    final_test <- rank_prts(counts, "combo", 50)
    test_equal <- dplyr::all_equal(final_ranks, final_test)
    if (test_equal == TRUE) {
        print("The test was successful. SurfaltR works as expected!")
    } else {
        print("The test was unsuccessful. Please check your installation.")
    }
    setwd(old)
}
