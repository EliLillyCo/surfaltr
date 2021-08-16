#' Retrieve, Clean, and Format Input Data
#'
#' This function cleans and formats input data. The cleaning and formatting
#' portion involves removing any non-protein coding transcripts, removing any
#' principal transcripts, and standardizing all column names.
#' If the sequence is provided directly, the function also extracts the APPRIS
#' annotation and UniProt IDs of each transcript from Ensembl. Provided data can
#' follow 2 formats — the first option only contain transcript IDs and gene
#' names and the second option contains a unique transcript identifier, gene
#' names, and amino acid sequences. The function will return a data frame
#' containing the transcript IDs, gene names, and APPRIS Annotation for each
#' inputted transcript. If the amino acid sequence is included in the input 
#' data, this will also be included in the data frame. If only gene names and
#' transcript IDS are provided, UniProt IDs will be included in the data frame.
#'
#' @usage clean_data(data_file, if_aa, organism)
#' @param data_file Path to the input file
#' @param if_aa Boolean value indicating if the input file contains
#' amino acid sequences with TRUE indicating that sequences are present and
#' FALSE indicating that only IDs are present
#' @param organism String indicating if the transcripts are from a human or
#' a mouse
#' @importFrom utils read.csv
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
#' @importFrom httr set_config
#' @importFrom httr config
#' @return A data frame containing gene names, transcript IDs, and APPRIS
#' annotations for the given data. If sequences were provided, the data frame
#' will also contain amino acid sequences. If only IDs were provided, the data
#' frame will also contain the UniProt Swissprot ID, UniProt Swissprot
#' isoform ID, and UniProt TREMBL ID.



clean_data <- function(data_file, if_aa, organism) {
    all_genome_data <- read.csv(file = data_file, header = TRUE)
    if (if_aa == TRUE) {
        if (!("external_gene_name" %in% colnames(all_genome_data)) |
            !("transcript_id" %in% colnames(all_genome_data)) |
            !("protein_sequence" %in% colnames(all_genome_data))) {
            stop("Your input data has incorrect column names. Please rename 
            columns according to this function's help page before proceeding.")
        }
        final_trans <- all_genome_data
        final_trans$transcript_appris <- ""
        final_trans <- dplyr::rename(final_trans, "ensembl_transcript_id" = 
                                    "transcript_id")
        return(final_trans)
    }
    if (if_aa == FALSE) {
        if (!("gene_name" %in% colnames(all_genome_data)) |
            !("transcript" %in% colnames(all_genome_data))) {
            stop("Your input data has incorrect column names. Please rename 
            columns according to this function's help page before proceeding.")
        }
        all_genome_data <- dplyr::distinct(all_genome_data, 
        all_genome_data$transcript, .keep_all = TRUE)
        if (organism == "human") {
            httr::set_config(httr::config(ssl_verifypeer = 0L))
            ensembl <- biomaRt::useEnsembl(
                biomart = "genes", dataset = "hsapiens_gene_ensembl",
                mirror = "useast"
            )
            compare <- biomaRt::getBM(
                attributes = c(
                    "external_gene_name", "ensembl_transcript_id",
                    "transcript_appris", "transcript_biotype",
                    "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"
                ),
                mart = ensembl
            )
            merge_trans <- merge(
                x = compare, y = all_genome_data, by.x = "ensembl_transcript_id",
                by.y = "transcript"
            )
            if (nrow(merge_trans) == 0) {
        stop("Please make sure your organism is correct and your transcript IDs 
        are valid and contain no decimal values before continuing. ")
            }
            merge_trans <- merge_trans[which(
                merge_trans$transcript_appris != "principal1" &
                merge_trans$transcript_appris != "principal2" &
                merge_trans$transcript_appris != "principal3" &
                merge_trans$transcript_appris != "principal4" &
                merge_trans$transcript_appris != "principal5"), ]
            for (row in seq_len(nrow(merge_trans))) {
                type <- merge_trans[row, "transcript_biotype"]
                trans <- merge_trans[row, "ensembl_transcript_id"]
                if (type != "protein_coding") {
                    print("The transcript")
                    print(trans)
                    print("is non-protein coding. It has been removed 
                          from further analysis.")
                }
            }
            final_trans <- merge_trans[grep("protein_coding", 
                           merge_trans$transcript_biotype), ]
            final_trans <- subset(final_trans, select = c(
                "ensembl_transcript_id",
                "external_gene_name", "transcript_appris",
                "uniprotswissprot", "uniprotsptrembl",
                "uniprot_isoform"
            ))
        }
        if (organism == "mouse") {
            ensembl <- biomaRt::useEnsembl(biomart = "genes", 
                                           dataset = "mmusculus_gene_ensembl")
            compare <- biomaRt::getBM(
                attributes = c(
                    "external_gene_name",
                    "ensembl_transcript_id",
                    "transcript_appris", "transcript_biotype",
                    "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"
                ),
                mart = ensembl
            )
            merge_trans <- merge(x = compare, y = all_genome_data, 
            by.x = "ensembl_transcript_id", by.y = "transcript")
            if (nrow(merge_trans) == 0) {
                stop("Please make sure your organism is correct and your 
            transcript IDs are validand contain no decimal values 
            before continuing. ")
            }
            merge_trans <- merge_trans[which(
                merge_trans$transcript_appris != "principal1" &
                merge_trans$transcript_appris != "principal2" &
                merge_trans$transcript_appris != "principal3" &
                merge_trans$transcript_appris != "principal4" &
                merge_trans$transcript_appris != "principal5"), ]
            for (row in seq_len(nrow(merge_trans))) {
                type <- merge_trans[row, "transcript_biotype"]
                trans <- merge_trans[row, "ensembl_transcript_id"]
                if (type != "protein_coding") {
                print("The transcript")
                print(trans)
                print("is non-protein coding. It has been removed from further 
                analysis.")}
            }
            final_trans <- merge_trans[grep("protein_coding", 
                                            merge_trans$transcript_biotype), ]
            final_trans <- subset(final_trans, select = c(
                "ensembl_transcript_id",
                "external_gene_name", "transcript_appris",
                "uniprotswissprot", "uniprotsptrembl",
                "uniprot_isoform"
            ))
        }
        return(final_trans)
    }
}
