## Pattern search in parallel
## This should be switched to mclapply and foreach\
## This function takes in fastq files and searches for telomere
## repeats in the sequences. It returns a list of sequences which
## contain telomere repeats and a list of statistics for each file.
#' @title Search for Telomeres
#' @description \code{telo_search} takes as input a \code{list} of fastq files and 
#' searches for telomere repeats in the sequences. It returns a 
#' list of sequences which contain telomere repeats (if \code{return_telomeres}: TRUE),
#' a \code{data.table} of statistics for each file and a \code{list} of files.
#' \code{return_telomeres} should only be TRUE when you have an exssesive amount of
#' memory (RAM) to store the telomere sequences in the R session, roughly > 2 TB.
#' @details The 'telo_search' function searches for telomere repeats in FASTQ files and 
#' generates related statistics. It begins by loading necessary libraries and configuring 
#' the TensorFlow environment, followed by cleaning up memory. The function reads the FASTQ 
#' files, searches for motifs from the 'grep_list', and, based on the 'return_telomeres' 
#' parameter, either returns the sequences containing telomere repeats or just the 
#' statistics. It supports parallel processing with a specified number of threads and 
#' can process multiple files at once. Additionally, 'verbose' and 'progress' options provide 
#' user feedback during execution. The function returns a list with the telomere sequences 
#' (if requested), a 'data.table' of statistics, and a list of the processed files.
#'
#' @param fastq_files A character vector of fastq files to search.
#' @param grep_list A character vector of telomere motifs to search for.
#' @param names_prefix A character string to prefix the names of the telomere sequences.
#' @param threads An integer of the number of threads to use.
#' @param out_dir A character string of the output directory.
#' @param process_file_number An integer of the number of files to process at once.
#' @param return_telomeres A logical of whether to return the telomere sequences.
#' @param verbose A logical of whether to print messages.
#' @param progress A logical of whether to show a progress bar.
#' @return A \code{list} of telomere sequences (if \code{return_telomere}: is TRUE), 
#' \code{data.table} of statistics, and a file \code{list}.
telo_search <- function(fastq_files = NULL,
                        grep_list = NULL,
                        names_prefix = "telomere_",
                        threads = 1,
                        out_dir = getwd(),
                        process_file_number = 1,
                        return_telomeres = FALSE,
                        verbose = FALSE,
                        progress = FALSE) {
    # Load libraries
    library(keras)
    library(reticulate)
    use_condaenv(condaenv = "r-tensorflow")

    # Clean up
    rm(results)
    gc()
}
