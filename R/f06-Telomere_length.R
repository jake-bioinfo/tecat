## This script is used to calculate the telomere length of each sample
## in a given dataset)
# Draw to Smash Logic Puzzle

# This function is used to calculate the lenght of the telomere
# repeat in a given read.
#' @title findTelLength
#' @description This function is used to calculate the length of the
#' telomere repeat in a given read.
#' @param df A data frame containing the telomere repeat information.
#' @param st.thresh A numeric value specifying the start threshold for the
#' telomere length calculation. Default is 50.
#' @param en.thresh A numeric value specifying the end threshold for the
#' telomere length calculation. Default is 50.
#' @param windows_size A numeric value specifying the size of the sliding
#' window. Default is 5.
#' @return A data frame containing the telomere length, read length, start
#' and end truncation points, and the telomere end.
#' @export
findTelLength <- function(df, 
    st.thresh = 50, 
    en.thresh = 50, 
    windows_size = 5) {

    # Assign frequency vector
    freq_vec <- df$motif_frequency

    # Create initial values
    start.tel <- mean(head(freq_vec, windows_size)) > st.thresh
    end.tel <- mean(tail(freq_vec, windows_size)) > st.thresh
    read.length <- max(df$end)

    # Initialized tel length to 0
    tel.length <- 0
    trunc.s <- tel.length + 1
    trunc.e <- read.length
    tel.end <- NA

    # First condidtion, telomere starts at the beginning of the read
    if (start.tel & !end.tel) {
        tel.end <- "L"

        # Reassign first row to the first of the 3 rows that are above the start threshold
        start.ind <- 1
        
        # Calculate on a 5 window size sliding window with a step of 1
        ## At some point I would like to flag the very low values within telomeres iot 
        ## retain model based on motifs found here
        end.ind <- which(zoo::rollapply(freq_vec, windows_size, mean, by = 1) < en.thresh)[1] + 1

        # Determine telomere length      
        ml <- nchar(df[end.ind, "hit_motif"])  
        tc <- df[end.ind, "telomere_count"]
        tel.length <- df[end.ind, "start"] + tc * ml

        # Return telomere value
        trunc.s <- tel.length + 1
        trunc.e <- read.length
    }

    # Second condition, telomere starts at the end of the read
    if (end.tel & !start.tel) {
        tel.end <- "R"

        # Reassign first row to the first of the 3 rows that are above the start threshold
        start.ind <- length(freq_vec)
        
        # Calculate on a 5 window size sliding window with a step of 1
        ## At some point I would like to flag the very low values within telomeres iot 
        ## retain model based on motifs found here
        end.ind <- length(freq_vec) - which(zoo::rollapply(rev(freq_vec), windows_size, mean, by = 1) < en.thresh)[1]

        # Determine telomere length      
        ml <- nchar(df[end.ind, "hit_motif"])  
        tc <- df[end.ind, "telomere_count"]
        tel.length <- read.length - df[end.ind, "end"] + tc * ml

        # Return telomere value
        trunc.s <- 1
        trunc.e <- read.length - tel.length
    }

    # Assign value of 0, if entire read is telomere repeat
    if (start.tel & end.tel) {
        tel.length <- 0
        tel.end <- NA
    }

    # If everything else fails, assign 0
    if (is.na(tel.length) || tel.length == 0) {
        tel.length <- 0
        tel.end <- NA
    }

    ret_df <- data.frame(
        telomere_length = as.numeric(round(tel.length)),
        read_length = as.numeric(read.length),
        trunc_start = as.numeric(trunc.s),
        trunc_end = as.numeric(trunc.e),
        telomere_end = as.character(tel.end)
    )
    return(ret_df)
}

# This function is used to determine the start and stop thresholds
# for the telomere length calculation.
#' @title determine_threshold
#' @description Determine the start and stop thresholds for the telomere
#' length calculation.
#' @param telomere_list A list of data frames containing the telomere
#' sequences to be processed.
#' @param platform A character string specifying the sequencing platform
#' used to generate the data. Default is "pb" for PacBio.
#' @param sample_ratio A numeric value specifying the proportion of the
#' data to sample for threshold determination. Default is 0.15.
#' @param threads An integer specifying the number of parallel processors
#' to use. Default is 0.
#' @return A list of data frames containing the start and end thresholds
#' for the telomere length calculation.
#' @export
determine_threshold <- function(telomere_list = NULL,
                                platform = "pb",
                                sample_ratio = 0.15,
                                threads = 0) {

    # Check that telomere_list is not NULL
    if (is.null(telomere_list)) {
        stop("telomere_list is NULL")
    }

    # Determining number of threads
    if (threads == 0) {
        cat("\nNumber of parallel processors automatically set to max-2.\n")
        threads <- parallel::detectCores() - 2
    }

    # Take sample of the data
    sample_list <- sample(telomere_list,
        size = ceiling(length(telomere_list) *
            sample_ratio)
    )

    # Clean up some
    rm(telomere_list)
    gc()

    # Calculate sample lengths based on platform
    if (platform == "pb") {
        cat("\nDetermining thresholds for full \n
           data set based off of downsampled data, \n
           and start threshold of 50.\n")
        # Set start and end thresholds
        start <- 1
        end <- 100

        # Generate temporary result data
        #end_result_list <- parallel::parLapply(cl, start:end, function(en) {
        end_result_list <- parallel::mclapply(start:end, function(en) {
            start <- 50
            result <- data.frame(matrix(nrow = 0, ncol = 8))
            for (telomere in sample_list) {
                row <- findTelLength(
                    df = telomere,
                    st.thresh = start,
                    en.thresh = en
                )
                row$start <- start
                row$end <- en
                row$threshold <- paste(c(start, "-", en), collapse = "")
                result <- rbind(result, row)
            }

            # Clean up
            rm(list = setdiff(ls(), "result"))
            gc()

            # Return result
            return(result)
        }, mc.cores = threads)

        cat("\nDetermining thresholds for full \n
            data set based off of downsampled data, \n
            and start threshold of 50.\n")

        # Generate temporary result data
         start_result_list <- parallel::mclapply(start:end, function(st) {
            end <- 50
            result <- data.frame(matrix(nrow = 0, ncol = 8))
            for (telomere in sample_list) {
                row <- findTelLength(
                    df = telomere,
                    st.thresh = st,
                    en.thresh = end
                )
                row$start <- st
                row$end <- end
                row$threshold <- paste(c(st, "-", end), collapse = "")
                result <- rbind(result, row)
            }

            # Clean up
            rm(list = setdiff(ls(), "result"))
            gc()

            # Return result
            return(result)
        }, mc.cores = threads)
    }

    # Clean up
    rm(list = setdiff(ls(), c(
        "start_result_list",
        "end_result_list"
    )))
    gc()

    # Combine results
    return(list(
        start = start_result_list,
        end = end_result_list
    ))
}

# This function is used to determine the optimal thresholds for the
# telomere length calculation.
#' @title determine_optimal_thresholds
#' @description Determine the optimal thresholds for the telomere length
#' calculation.
#' @param threshold_dataframe A data frame containing the start and end
#' thresholds for the telomere length calculation.
#' @return A list containing the optimal thresholds for the telomere length
#' calculation.
#' @export
optimal_thresholds <- function(threshold_dataframe) {
    # Start sensitivity
    thresh_df <- threshold_dataframe
    start_number_of_telos <- unlist(lapply(thresh_df$start, function(x) {
        length(which(!is.na(x$telomere_end)))
    }))
    sensitivity_df <- data.frame(
        start = seq(1, 100, 1),
        end = rep(5, 100),
        number_of_telos = start_number_of_telos
    )

    # Sens threshold 
    sens_thresh <- sensitivity_df[which.max(sensitivity_df$number_of_telos), "start"]

    # Telomere length
    # End telomere length
    telomere_lengths_end <- lapply(thresh_df$end, function(x) {
        vals <- x$telomere_length[!is.na(x$telomere_end)]
        x <- mean(vals)
        if (x != 0 & length(vals) > 1) {
            std_dev <- sd(vals)
        } else {
            std_dev <- 0
        }
        ret_frame <- data.frame(
            mean = x,
            std_dev = std_dev
        )
        return(ret_frame)
    })
    end_df <- do.call(rbind, telomere_lengths_end)
    telomere_length_df <- cbind(
        start = rep(20, 100),
        end = seq(1, 100, 1),
        end_df
    )

    # Telomere threshold, find inflection
    tel_thresh <- round(inflection::ede(telomere_length_df$end, 
                                  telomere_length_df$mean, 
                                  index = 0)[,3] * 0.5)
    # We multiply by 0.5 to get the very start of the flattening out of the 
    # curve, which is the portion we are insterested in.

    # Return values
    return(list(
        sensitivity = sens_thresh,
        telomere_length = tel_thresh
    ))
}

# Plot thresholds data with thresholds
#' @title plot_thresholds
#' @description Plot the thresholds data with the thresholds.
#' @param threshold_dataframe A data frame containing the start and end
#' thresholds for the telomere length calculation.
#' @param optimal_thresholds A list containing the optimal thresholds for
#' the telomere length calculation.
#' @return A ggplot object containing the plot of the thresholds data.
#' @export
#' @import ggplot2 cowplot
plot_thresholds <- function(threshold_dataframe = NULL,
                            optimal_thresholds = NULL) {
    # Check that threshold_dataframe is not NULL
    if (is.null(threshold_dataframe)) {
        stop("threshold_dataframe is NULL")
    }
    # Check that optimal_thresholds is not NULL
    if (is.null(optimal_thresholds)) {
        stop("optimal_thresholds is NULL")
    }

    # Reformat data
    thresh_df <- threshold_dataframe
    # Start sensitivity
    start_number_of_telos <- unlist(lapply(thresh_df$start, function(x) {
        length(which(!is.na(x$telomere_end)))
    }))
    sensitivity_df <- data.frame(
        start = seq(1, 100, 1),
        end = rep(5, 100),
        number_of_telos = start_number_of_telos
    )

    # Telomere length
    # End telomere length
    telomere_lengths_end <- lapply(thresh_df$end, function(x) {
        vals <- x$telomere_length[!is.na(x$telomere_end)]
        x <- mean(vals)
        if (x != 0 & length(vals) > 1) {
            std_dev <- sd(vals)
        } else {
            std_dev <- 0
        }
        ret_frame <- data.frame(
            mean = x,
            std_dev = std_dev
        )
        return(ret_frame)
    })
    end_df <- do.call(rbind, telomere_lengths_end)
    telomere_length_df <- cbind(
        start = rep(20, 100),
        end = seq(1, 100, 1),
        end_df
    )

    # Plot the thresholds data
    library(cowplot)

    # Plot sensitivity and telomere line graphs next to each other
    sens_plot <- ggplot(sensitivity_df, aes(x = start, y = number_of_telos)) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept = optimal_thresholds$sensitivity, color = "red") +
        labs(
            x = "Start Threshold",
            y = "Number of Telomeres",
            title = "Telomere Sensitivity"
        ) +
        theme_bw()

    telomere_plot <- ggplot(telomere_length_df, aes(x = end, y = mean)) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept = optimal_thresholds$telomere_length, color = "red") +
        labs(
            x = "End Threshold",
            y = "Telomere Length",
            title = "Telomere Length"
        ) +
        theme_bw()

    return(plot_grid(sens_plot, telomere_plot, ncol = 2))
}
