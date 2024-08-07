## Create sliding window of reads
# Sliding window function
# !-- Depends on foreach, IRanges, doSNOW, parallel --!
# !-- fixed Frequencies, this needs the same treatment, as well as pattern searching -
#' @title Sliding Window Telomere Sequence Processing
#' @description The 'sliding_window' function generates overlapping windows
#' of specified length from telomere sequences, facilitating the analysis
#' of smaller, manageable segments of the sequences.
#' @details The 'sliding_window' function processes a set of telomere
#' sequences by creating overlapping windows of a specified length.
#' For each sequence in the input list, it generates a series of windows
#' starting from the beginning of the sequence, moving in steps defined
#' by the 'step' parameter. Each window has a length specified by the
#' 'window_length' parameter. The windows are stored in an
#' 'IRanges::Views' object, with each window named according to the
#' original sequence. After processing all sequences, the function
#' cleans up memory and returns the list of windowed views, allowing
#' for detailed analysis of these smaller sequence segments.
#'
#' @param telomere_sequences A DNAStringSet object containing the telomere
#' sequences to be processed.
#' @param window_length An integer specifying the length of each sliding
#' window. Default is 200.
#' @param step An integer specifying the step size for moving the window.
#' Default is 100.
#' @returns A list of IRanges::Views objects, each containing overlapping
#' windows of the original telomere sequences.
#' @importFrom IRanges Views
#' @export
sliding_window <- function(telomere_sequence,
                           window_length = 200,
                           step = 100) {

  # Create sliding window of telomere sequences
  names <- rep(
    names(telomere_sequence),
    ceiling(nchar(telomere_sequence) / window_length)
  )
  result <- IRanges::Views(telomere_sequence,
      start = seq(
        from = 1,
        to = length(telomere_sequence),
        by = step
      ),
      width = window_length,
      names = names
    )

  # Clean up
  rm(telomere_sequence)
  gc()

  return(result)
}

## Sliding windows parallel function
# Sliding window parallel function
# !-- Depends on sliding_window --!
#' @title Sliding Window Telomere Sequence Processing in Parallel
#' @description The 'sliding_window_parallel' function generates overlapping
#' windows of specified length from telomere sequences in parallel, facilitating
#' the analysis of smaller, manageable segments of the sequences.
#' @details The 'sliding_window_parallel' function processes a set of telomere
#' sequences by creating overlapping windows of a specified length in parallel.
#' For each sequence in the input list, it generates a series of windows starting
#' from the beginning of the sequence, moving in steps defined by the 'step'
#' parameter. Each window has a length specified by the 'window_length' parameter.
#' The windows are stored in an 'IRanges::Views' object, with each window named
#' according to the original sequence. After processing all sequences, the function
#' cleans up memory and returns the list of windowed views, allowing for detailed
#' analysis of these smaller sequence segments.
#' @param telomere_file A DNAStringSet object containing the telomere sequences
#' to be processed.
#' @param window_length An integer specifying the length of each sliding window.
#' Default is 200.
#' @param step An integer specifying the step size for moving the window. Default
#' is 100.
#' @param environment A character string specifying the environment in which to
#' run the function. Default is 'linux'.
#' @param threads An integer specifying the number of threads to use for parallel
#' processing. Default is 1.
#' @param verbose A logical specifying whether to print messages. Default is FALSE.
#' @returns A list of IRanges::Views objects, each containing overlapping windows
#' of the original telomere sequences.
#' @importFrom Biostrings readDNAStringSet
#' @import parallel pbmcapply doSNOW foreach
#' @export
sliding_window_parallel <- function(telomere_file,
                                    window_length = 200,
                                    step = 100,
                                    environment = "linux",
                                    threads = 1,
                                    verbose = FALSE) {
  # Status message
  if (verbose) {
    cat(
      "\nCreating sliding windows from reads: \n\t",
      paste(telomere_file, collapse = "\n\t"), "\n\t",
      "using ", threads, " threads.\n"
    )
  }

  # Read in files
  X <- Biostrings::readDNAStringSet(telomere_file)

  # Run in parallel
  if (environment == "linux") {
    if (verbose) {
      cat(
        "\n\tRunning in parallel with ", threads,
        " threads, in ", environment, " environment."
      )
    }
    # Run in parallel in linux
  ret <- pbmcapply::pbmclapply(X, sliding_window, 
    window_length = window_length, 
    step = step, mc.cores = threads, 
    ignore.interactive = TRUE)

  } else {
    if (verbose) {
      message("\n\tRunning in parallel with ", threads, " threads,
              in ", environment, ".")
    }
    # Running in parallel in windows
    library(doSNOW)
    library(foreach)
    cl <- makeCluster(threads)
    registerDoSNOW(cl)

    # Run parallel
    ret <- foreach(
      x = X,
      .combine = "c",
      .export = "sliding_window"
    ) %dopar% {
      sliding_window(x,
        window_length = window_length, 
        step = step
      )
    }
    stopCluster(cl)
  }
  
  # Clean up
  rm(X)
  gc()

  return(ret)
}