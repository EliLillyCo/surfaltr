#' Create a data frame with the membrane locations of each amino acid in a
#' protein using Phobius
#'
#' This function creates a data frame with columns containing transcript
#' IDs and corresponding output from Phobius. The Phobius output includes a
#' location for each amino acid, with O representing extracellular, M
#' representing transmembrane, S representing signal, and i representing
#' intracellular.
#'
#' @usage run_phobius(AA_seq, fasta_file_name)
#' @param AA_seq A data frame outputted by the get_pairs function
#' containing the gene names, transcript IDs, APPRIS annotations, and
#' protein sequences for each transcript.
#' @param fasta_file_name Path to fasta file containing amino acid sequences
#' @return A data frame containing each transcript ID and the corresponding
#' membrane location for each amino acid in its sequence formatted as a string
#' @note In order for this function to work, there needs to be a .fasta file
#' containing the amino acid sequences for each transcript called "AA.fasta"
#' saved to the working directory. Additionally, the file saves a copy of the
#' returned data frame in csv format to the working directory.
#' @importFrom utils write.csv
#' @importFrom stringr str_split
#' @importFrom ragp get_phobius
#' @importFrom readr parse_number
#' @examples
#' \donttest{
#' AA_seq <- get_pairs(system.file("extdata", "CRB1.csv",
#' package = "surfaltr"), TRUE, "mouse", TRUE)
#' topo <- run_phobius(AA_seq, paste(getwd(), "/AA.fasta", sep = ""))
#' }
#' @export


run_phobius <- function(AA_seq, fasta_file_name){
  phobius <- ragp::get_phobius(fasta_file_name)
  topo <- data.frame("Transcript_ID" = character(nrow(phobius)),
                     "Output" = character(nrow(phobius)))
  topo$Transcript_ID <- phobius$Name
  for(row in 1:nrow(phobius)){
    regs <- as.data.frame(stringr::str_split(phobius[row, "prediction"], "/"))
    phobius[row,"Signal"] <- regs[1,]
    phobius[row, "Surf_Protein"] <- regs[2,]
    if (is.na(phobius[row, "Surf_Protein"])){
      phobius[row, "Surf_Protein"] <- regs[1,]
      phobius[row, "Signal"] <- NA
    }
    if (is.na(phobius[row, "Signal"]) == FALSE){
      num_sig <- substr(phobius[row, "Signal"], nchar(phobius[row, "Signal"])-1,
                        nchar(phobius[row, "Signal"]))
      phobius[row, "Mem_String"] <- paste(replicate(num_sig, "S"), collapse = "")
    }
    surf_prot <- as.data.frame(stringr::str_split(phobius[row, "Surf_Protein"], "-"))
    colnames(surf_prot) <- c("Split")
    for(val in 1:nrow(surf_prot)){
      if((substr(surf_prot[val, "Split"], 1, 1) == "i" | substr(surf_prot[val, "Split"], 1, 1) == "o")
         && nchar(surf_prot[val, "Split"]) == 1){
        curr_let <- substr(surf_prot[val, "Split"], 1, 1)
        num_lets <- nchar(AA_seq[row,"protein_sequence"])
        phobius[row, "Mem_String"] <- paste(replicate(num_lets, curr_let), collapse = "")
      }
      if((substr(surf_prot[val, "Split"], 1, 1) == "i" | substr(surf_prot[val, "Split"], 1, 1) == "o")
         && nchar(surf_prot[val, "Split"]) != 1){
        curr_let <- substr(surf_prot[val, "Split"], 1, 1)
        num_lets <- readr::parse_number(surf_prot[val, "Split"])
        phobius[row, "Mem_String"] <- paste(replicate(num_lets, curr_let), collapse = "")
        nxt_len <- stringr::str_extract_all(surf_prot[val+1, "Split"],"\\(?[0-9,.]+\\)?")[[1]]
        mem_len <- as.numeric(nxt_len[1])-num_lets
        phobius[row, "Mem_String"] <- paste(phobius[row, "Mem_String"], paste(replicate(mem_len, "M"),
                                                                              collapse = ""), sep = "")
      }
      if(suppressWarnings(!is.na(as.numeric(substr(surf_prot[val, "Split"], 1, 1))))){
        len_nums <- stringr::str_extract_all(surf_prot[val, "Split"],"\\(?[0-9,.]+\\)?")[[1]]
        if (length(len_nums) != 1){
          seg_len <- as.numeric(len_nums[2])-as.numeric(len_nums[1])
          curr_let <- gsub("[^a-zA-Z]", "", surf_prot[val, "Split"])
          phobius[row, "Mem_String"] <- paste(phobius[row, "Mem_String"],
                                              paste(replicate(seg_len, curr_let), collapse = ""), sep = "")
          nxt_len <- stringr::str_extract_all(surf_prot[val+1, "Split"],"\\(?[0-9,.]+\\)?")[[1]]
          mem_len <- as.numeric(nxt_len[1])-as.numeric(len_nums[2])
          phobius[row, "Mem_String"] <- paste(phobius[row, "Mem_String"],
                                              paste(replicate(mem_len, "M"), collapse = ""), sep = "")
        } else {
          seg_len <- nchar(AA_seq[row,"protein_sequence"])-as.numeric(len_nums)
          curr_let <- gsub("[^a-zA-Z]", "", surf_prot[val, "Split"])
          phobius[row, "Mem_String"] <- paste(phobius[row, "Mem_String"],
                                              paste(replicate(seg_len, curr_let), collapse = ""), sep = "")
        }
      }
    }
  }
  topo$Output <- phobius$Mem_String
  write.csv(topo, paste(getwd(), "/mem_topo_phobius.csv", sep = ""))
  return(topo)
}
