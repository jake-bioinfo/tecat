## Get telomere frequencies from telomeres cut into sliding windows
# Create telomere frequencies based on each sliding window in an
# individual telomere.
#' @title Generate telomere frequencies
#' @description Generate telomere frequencies based on each sliding window in an
#' individual telomere.
#' @param telomere_windows An XStringViews object containing telomere windows.
#' @param motifs A character vector of telomere motifs to search for.
#' @return A data frame containing the start, end, window length, telomere
#' count, motif frequency, and hit motif.
#' @importFrom Biostrings DNAStringSetList DNAStringSet countPattern
#' @export
generate_telomere_frequencies <- function(telomere_windows, 
                                          motifs) {

    # Check variables
    stopifnot(
        class(telomere_windows) == "XStringViews",
        is.character(motifs),
        !is.null(telomere_windows),
        !is.null(motifs)
    )
    
    # Create grep list for counting telomere repeats
    motifs <- Biostrings::DNAStringSetList(motifs)
    motifs <- unlist(motifs)

    # Create result matrix to store all the results
    result_df <- data.frame(matrix(0, nrow = length(telomere_windows), ncol = 7))
    
    for (bpwin in seq_len(length(telomere_windows))) {
        start <- Biostrings::start(telomere_windows[bpwin])
        seq_length <- Biostrings::nchar(telomere_windows[bpwin])
        end <- start + seq_length - 1
        window_length <- end - start + 1
        subject <- Biostrings:::fromXStringViewsToStringSet(
            telomere_windows[bpwin], 
            out.of.limits = "ok")
                
        out <- try({
            sapply(motifs, function(mot) {
            Biostrings::countPattern(mot, 
                subject,
                fixed = FALSE)
        })}, silent = TRUE)

        if (inherits(out, "try-error")) {
            subject <- Biostrings:::fromXStringViewsToStringSet(
                telomere_windows[bpwin], 
                out.of.limits = "ok")[[1]]
                
            out <- sapply(motifs, function(mot) {
                Biostrings::countPattern(mot, 
                    subject,
                    fixed = FALSE)
            })
        }

        # Assign hits
        chars <- out * nchar(motifs)
        mx_char <- max(chars)
        mx_char <- ifelse(mx_char > window_length, window_length, mx_char)
        hit_motif <- as.character(motifs[which.max(chars)][[1]])
        hit_motif <- ifelse(mx_char == 0, "None", hit_motif)
        motif_frequency <- mx_char / window_length * 100
        telomere_count <- out[which.max(chars)]
        
        # Assign values to result matrix
        result_df[bpwin, ] <- data.frame(as.numeric(start), 
            as.numeric(end), 
            as.numeric(window_length), 
            as.numeric(telomere_count), 
            as.numeric(motif_frequency), 
            as.character(hit_motif))
    }
    
    # Return result
    colnames(result_df) <- c("start", "end", "window_length", 
        "telomere_count", "motif_frequency", "hit_motif")
    return(result_df)
}

# Tie together all the functions into an easy to use function
# that will process the telomere files and return the results
# in a list.
#' @title Generate telomere frequencies
#' @description Generate telomere frequencies based on each sliding window in an
#' individual telomere.
#' @param telomere_windows An XStringViews object containing telomere windows.
#' @param motifs A character vector of telomere motifs to search for.
#' @return A data frame containing the start, end, window length, telomere
#' count, motif frequency, and hit motif.
#' @importFrom Biostrings DNAStringSetList DNAStringSet countPattern
#' @export
frequencies <- function(windows,
                        motifs,
                        environment = "linux",
                        parallel = FALSE,
                        threads = 0,
                        save_files = TRUE,
                        out_dir = file.path(getwd(), "telomere_frequencies"),
                        progress = FALSE,
                        verbose = FALSE) {

    # Check variables
    stopifnot(
        is.list(windows),
        is.character(motifs),
        environment %in% c("linux", "windows"),
        is.logical(parallel),
        is.numeric(threads)
    )

    nms <- names(windows)
    freqs <- pbmcapply::pbmclapply(windows, function(w) {
        generate_telomere_frequencies(w, motifs)
    }, mc.cores = threads)

    # Assign names
    names(freqs) <- nms
    for (i in 1:length(freqs)) {
        freqs[[i]]$sequence_name <- nms[i]
    }

    # Save files
    if (save_files) {
        if (verbose) {
            message("\tSaving files to ", out_dir)
        }
        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }
        file_name <- file.path(
            out_dir,
            paste0("telomere_frequencies", ".rds")
        )
        saveRDS(freqs, file_name)
    }

    # Clean up
    rm(list = setdiff(ls(), "freqs"))
    gc()

    # Reformat frequencies list and return
    return(freqs)
}