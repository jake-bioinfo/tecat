## Isolate reads which contain telomere repeats
#' @title Search FASTQ Files for Telomere Repeats and Collect Statistics
#'
#' @description telomere_matching takes in a fastq file and searches for telomere
#' repeats in the sequences. It returns a list of sequences which
#' contain telomere repeats and a list of statistics for each file.
#' 
#' @details The telomere_matching function extracts sequences containing telomere repeats 
#' from a FASTQ file by searching for specified motifs and computes relevant statistics. 
#' It reads the FASTQ file, identifies sequences with the motifs, and saves these 
#' sequences to a file while also providing a summary of telomere-related statistics. 
#' The function can return either just the statistics or both the statistics and the 
#' extracted sequences, depending on the return_telomeres parameter.
#' 
#' @param fastq_file A character string of the fastq file to search.
#' @param names_prefix A character string to prefix the names of the telomere sequences.
#' @param out_dir A character string of the output directory.
#' @param grep_list A character vector of telomere motifs to search for.
#' @param return_telomeres A logical of whether to return the telomere sequences.
#' @param verbose A logical of whether to print messages.
#' @return A list of telomere sequences (if return_telomere: TRUE) and a list of statistics.
#' @importFrom Biostrings readDNAStringSet writeXStringSet vcountPattern
#' @importFrom IRanges width
#' @export
telomere_matching <- function(fastq_file,
                              names_prefix = "telomere_",
                              out_dir = getwd(),
                              grep_list = c("TTAGGG", "CCCTAA"),
                              return_telomeres = FALSE,
                              verbose = FALSE) {

    if (verbose) {
        message("Processing: ", fastq_file, "\n")
    }

    # File type check
    file_type <- get_mime_type(fastq_file)
    if (file_type != "application/fasta" && file_type != "application/fastq") {
        stop("File type is not a fastq or fasta file. Please provide a fastq or fasta file.")
    }

    # Read in sequences
    if(file_type == "application/fasta") {
        sequences <- Biostrings::readDNAStringSet(
            filepath = fastq_file,
            format = "fasta"
        )
    } else {
        sequences <- Biostrings::readDNAStringSet(
            filepath = fastq_file,
            format = "fastq"
        )
    }

    # Try Biostings::matchPattern, fixed = FALSE
    g_ls <- Biostrings::DNAStringSet(as.character(grep_list))

    indices <- unlist(lapply(g_ls, function(pat) {
        rngs <- Biostrings::vcountPattern(
            pattern = pat,
            subject = sequences,
            algorithm = "auto",
            fixed = FALSE
        )
        which(rngs > 0)
    }))
    telomere_sequences <- sequences[unique(indices)]

    # File Names
    nm <- gsub(".fastq.gz|.fastq|.fq|.fa|.fasta.gz|.fasta",
         "", basename(fastq_file))

    # Collect telomere information
    telomere_info <- list(
        mean_telomere_read_length = sum(IRanges::width(telomere_sequences)) / length(telomere_sequences),
        total_telomere_length = sum(IRanges::width(telomere_sequences)),
        telomere_count = length(telomere_sequences),
        telomere_percent_of_reads = (length(telomere_sequences) /
            length(sequences)) * 100,
        telomere_percent_of_bases = (sum(IRanges::width(telomere_sequences)) /
            sum(IRanges::width(sequences))) * 100,
        file_name = nm
    )

    # Clean up
    rm(sequences)
    gc()

    # Add names
    if (length(telomere_sequences) > 0) {
        names(telomere_sequences) <- paste0(
            nm, "_",
            names_prefix,
            seq(001, length(telomere_sequences), 1)
        )
    } else {
        telomere_sequences <- NULL
    }

    # Save fasta
    Biostrings::writeXStringSet(
        x = telomere_sequences,
        file = file.path(
            out_dir,
            paste0(
                "telomeres_",
                nm, ".fasta.gz")), 
                format = "fasta", 
                compress = TRUE)

    # Return
    if (return_telomeres) {
        return(list(
            telomeres = telomere_sequences,
            telomere_stats = telomere_info
        ))
    } else {
        return(list(telomere_stats = telomere_info))
    }

    # Clean up
    rm(telomere_sequences)
    gc()
}

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
#' 
#' @details The telo_search function efficiently searches for telomere repeats across 
#' a list of FASTQ files using parallel processing. It leverages the telomere_matching 
#' function to analyze each file for telomere motifs, collecting sequences and statistics 
#' such as average read length and percentage of telomere bases. The function runs these 
#' searches concurrently with pbmcapply::pbmclapply, handles multiple threads, and can 
#' optionally show progress updates. After processing, it compiles a data.table of statistics 
#' and, based on the return_telomeres parameter, returns either a list of telomere sequences, 
#' statistics, and file paths, or just the statistics and file paths. The function also manages 
#' verbosity and output directory options to control the search process and results.
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
#' @export
telo_search <- function(fastq_files = NULL,
                        grep_list = NULL,
                        names_prefix = "telomere_",
                        threads = 1,
                        out_dir = getwd(),
                        process_file_number = 1,
                        return_telomeres = FALSE,
                        verbose = FALSE,
                        progress = FALSE) {

    # Status message
    message(
        "Searching for telomeres in ", length(fastq_files), " files\n",
        "Using ", threads, " threads\n"
    )

    # Create output directory if not all ready created
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    # Run parallel search
    results <- pbmcapply::pbmclapply(
        fastq_files,
        telomere_matching,
        grep_list = grep_list,
        names_prefix = names_prefix,
        out_dir = out_dir,
        verbose = verbose,
        return_telomeres = return_telomeres,
        mc.cores = threads, 
        mc.style = "ETA", 
        ignore.interactive = TRUE
    )

    # Add file list to output
    file_list <- list.files(
        path = out_dir,
        pattern = "telomeres_",
        full.names = TRUE
    )

    # Combine fasta files
    if (length(file_list) > 1) {
        combined_fasta <- combine_fasta(
            list_of_fasta = file_list,
            write_files = TRUE,
            out_dir = out_dir,
            return_reads = FALSE,
            verbose = verbose
        )
    }

    # Split results
    telomeres <- lapply(results, function(x) {
        x$telomeres
    })
    telomeres <- bio_ul(telomeres)

    # Formatting results
    stats_df <- data.table::rbindlist(
        lapply(results, function(x) {
            x$telomere_stats
        })
    )
    colnames(stats_df) <- c(
        "mean_telomere_read_length",
        "total_telomere_length",
        "telomere_count",
        "telomere_percent_of_reads",
        "telomere_percent_of_bases",
        "file_name"
    )

    # Return telomeres
    if (return_telomeres) {
        return(list(
            telomeres = telomeres,
            telomere_stats = stats_df,
            combined_fasta = combined_fasta
        ))
    } else {
        return(list(
            telomere_stats = stats_df,
            combined_fasta = combined_fasta
        ))
    }

    # Clean up
    rm(results)
    gc()
}