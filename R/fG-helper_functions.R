#' TECAT: Telomere End Chromosome Assaying Tool
#' 
#' @description 
#' A package for determining chromosome end specific telomere lengths using pattern 
#' matching with 3rd generation sequencing data.
#'
#' @details
#' The main functions of this package are:
#' \itemize{
#'   \item Telomere motif identification
#'   \item Pattern matching and analysis
#'   \item Length calculation
#'   \item Visualization
#' }
#'
#' @docType package
#' @name TECAT
NULL

## This creates a mime map for other functions etc to use
# Define a mapping of file extensions to corresponding MIME types
library(mime)
library(dplyr)
mimemap[["rds"]] <- "application/x-rds"
mimemap[["fastq"]] <- "application/fastq"
mimemap[["fq"]] <- "application/fastq"
mimemap[["fasta"]] <- "application/fasta"
mimemap[["fa"]] <- "application/fasta"
mimemap[["fastq.gz"]] <- "application/fastq-gz"
mimemap[["fq.gz"]] <- "application/fastq-gz"
mimemap[["fasta.gz"]] <- "application/fasta-gz"
mimemap[["fa.gz"]] <- "application/fasta-gz"

## Get the type of mime based on a separate function
#' @title Get MIME type
#' @description Get the MIME type of a file based on its extension.
#' @param file_path A character string specifying the path to the file.
#' @return A character string specifying the MIME type of the file.
#' @export
get_mime_type <- function(file_path) {
    # Correction for biologic data types
    if((grepl("\\.f", file_path) & grepl("gz", file_path))) {
        file_path <- tools::file_path_sans_ext(file_path)
    }
    tools::file_ext(file_path) %>% 
        tolower() %>%
        mimemap[.]
}

mime_read_map <- list(
    "application/x-rds" = function(file_path) {
        readRDS(file_path)
    },
    "text/csv" = function(file_path) {
        read.csv(file_path)
    },
    # Biology relevant files
    "application/fastq" = function(file_path) {
        readDNAStringSet(file_path, format = "fastq")
    },
    "application/fasta" = function(file_path) {
        readDNAStringSet(file_path)
    },
    "application/fastq-gz" = function(file_path) {
        readDNAStringSet(file_path, format = "fastq")
    },
    "application/fasta-gz" = function(file_path) {
        readDNAStringSet(file_path)
    }
)

#' @title Get system RAM
#' @description Get the total amount of RAM on the system and suggest an amount
#' of RAM to use for processing.
#' @param environment A character string specifying the operating system. Default
#' is "linux".
#' @param verbose A logical value specifying whether to display verbose messages.
#' Default is FALSE.
#' @return A list containing the total amount of RAM on the system and the suggested
#' amount of RAM to use for processing.
#' @export
get_system_ram <- function(environment, verbose = FALSE) {
    # Determine amount of ram on system
    if (environment == "linux") {
        lin_ram <- system("free -h | awk 'NR==2{print $2}'", 
            intern = TRUE)
        lin_ram <- gsub("Gi", "", lin_ram)
    } else {
        win_ram <- system("(Get-CimInstance Win32_PhysicalMemory | Measure-Object -Property capacity -Sum).sum /1gb", 
            intern = TRUE)
    }
    # Max ram usage
    ram <- ifelse(environment == "linux", lin_ram, win_ram)
    suggested_ram <- floor(as.numeric(gsub("G", "", ram)) * 0.7)
    
    if(verbose) {
        message("Total RAM: ", ram, " GB\n\tSuggested RAM: ", suggested_ram, 
            " GB.")
    }

    return(list(total_ram = ram, suggested_ram = suggested_ram))
}

## This function takes in a list of files and estimates the RAM usage of reading in each file
#' @title Estimate RAM usage
#' @description Estimate the RAM usage of reading in a list of files.
#' @param files A character vector specifying the paths to the files.
#' @param verbose A logical value specifying whether to display verbose messages.
#' Default is FALSE.
#' @return A list containing the mean RAM usage per file, a data frame of RAM usage
#' per file, and the suggested number of files to process.
#' @export
estimate_ram_usage <- function(files, verbose = FALSE) {

    sample_num <- max(floor(length(files) * 0.1), 2)
    sl <- sample(files, sample_num)

    # Generate ram estimates
    ram_ls <- lapply(sl, function(f) {
        nm <- basename(f)
        st_r <- as.numeric(pryr::mem_used()) / 10E8
        read <- mime_read_map[[get_mime_type(f)]]
        f <- read(f)
        end_r <- as.numeric(pryr::mem_used()) / 10E8
        dif <- end_r - st_r
        rm(f)
        gc()
        data.frame(file = nm, ram = dif, unit = "GB")
    })
    file_ram_df <- do.call(rbind, ram_ls)
    file_ram_mu <- mean(file_ram_df$ram)

    # Get system info
    system_ram <- get_system_ram("linux", verbose = verbose)
    
    # Calculate suggest file number
    suggested_num_files <- floor(system_ram$suggested_ram / file_ram_mu * 0.75)

    if(verbose) {
        message("\tMean RAM per file: ", file_ram_mu, " GB\n\t",
                "Suggested number of files: ", suggested_num_files)
    }

    return(list(mean_ram_per_file = file_ram_mu, 
                file_ram_df = file_ram_df, 
                suggested_num_files = suggested_num_files))
}

## This function returns size of objects in environment in Mb in descending order
#' @title Find object sizes
#' @description Find the sizes of objects in the environment in descending order.
#' @param clean A logical value specifying whether to clean the environment. Default
#' is FALSE.
#' @param which A character string specifying which objects to remove. Default is "all".
#' @param environ An environment object specifying the environment to search. Default
#' is .GlobalEnv.
#' importFrom Biostrings DNAStringSet nchar substring vcountPattern
#' @export
find_sizes <- function(clean = FALSE,
    which = "all",
    environ = .GlobalEnv) {
    # List all loaded objects
    loaded_objects <- ls(envir = environ)

    # Iterate through each object and get its size
    sizes <- lapply(loaded_objects, function(x) {
        object.size(get(x, envir = environ)) / 1024^2
    })

    # Combine object names and their sizes into a data frame
    object_sizes <- data.frame(Object = loaded_objects, Size_MB = as.numeric(unlist(sizes)))
    obj_sizes <- object_sizes[order(object_sizes$Size_MB, decreasing = TRUE), ]

    # Print the object sizes
    print(obj_sizes)
    cat("\n")

# Clean up the environment
    if (clean) {
        if (which == "all") {
            rm(list = ls(envir = environ), envir = environ)
        } else {
            print(paste("Removing objects: ", obj_sizes$Object[1:which]))
            rm(list = c(as.character(obj_sizes$Object[1:which])),
                envir = environ)
        }
        gc()
    }
    cat("\n")
}


## This function takes a list of `Biostrings::DNAStringSet` and creates a single list from 
## all the sequences
bio_ul <- function(rds) {
    ret <- Biostrings::DNAStringSet()
    for (i in seq_along(rds)) {
        ret <- Biostrings::DNAStringSet(c(ret, rds[[i]]))
    }
    return(ret)
}


# Encode read into a vector of integers
encode_match <- function(read, motifs) {
    if (Biostrings::nchar(read) < 400) {
        return(NULL)
    }

    # Define the window size and step size
    window_size <- max(nchar(motifs))
    step_size <- ceiling(0.5 * window_size)

    # Function to extract sliding windows from a single sequence
    extract_windows <- function(seq) {
        start_positions <- seq(1, Biostrings::nchar(seq) - window_size + 1, by = step_size)
        windows <- lapply(start_positions, function(start) {
            end <- start + window_size - 1
            if (end > Biostrings::nchar(seq)) {
                end <- Biostrings::nchar(seq)
            }
            Biostrings::substring(seq, start, end)
        })
        return(Biostrings::DNAStringSet(windows))
    }
    windows <- extract_windows(read)

    # Match each window to motifs
    out <- lapply(motifs, function(mot) {
        Biostrings::vcountPattern(mot,
            windows,
            fixed = FALSE
        )
    })
    out <- do.call(rbind, out)
    out <- colSums(out)

    # Assign vector
    out <- ifelse(out > 0, 1, 0)
    return(out)
}

#' @title Combine FASTA files
#' @description Combine multiple FASTA files into a single FASTA file.
#' @param list_of_fasta A character vector specifying the paths to the FASTA files.
#' @param write_files A logical value specifying whether to write the combined FASTA
#' file to disk. Default is TRUE.
#' @param out_dir A character string specifying the output directory. Default is
#' "combined_reads".
#' @param prefix A character string specifying the prefix for the output file. Default
#' is "combined_telo".
#' @param return_reads A logical value specifying whether to return the combined reads.
#' Default is FALSE.
#' @param verbose A logical value specifying whether to display verbose messages. Default
#' is TRUE.
#' @return A list containing the combined reads.
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @export
combine_fasta <- function(list_of_fasta = NULL, 
                          write_files = TRUE,
                          out_dir = file.path(getwd(), "combined_reads"),
                          prefix = "combined_telo",
                          return_reads = FALSE, 
                          verbose = TRUE) {

  # Check inputs
  stopifnot(is.character(list_of_fasta),
            is.logical(write_files),
            is.character(out_dir),
            is.character(prefix))

  # Create output directory if it doesn't exist
  dir.create(out_dir, showWarnings = FALSE)

  # Check if file already exists
  check <- list.files(out_dir, pattern = prefix)
  if(!is.null(file.exists(check)) & length(check) > 0){
    cat("File already exists. Please choose a different prefix.")
    return(invisible())
  }

  # Load files
  fasta <- lapply(list_of_fasta, Biostrings::readDNAStringSet)
  fasta <- bio_ul(fasta)

  # Combine
  if(write_files) {
    if(verbose) {
      message("\nWriting files to: ", out_dir)
    }
    if(!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    Biostrings::writeXStringSet(fasta, 
                                file.path(out_dir, paste0(prefix, ".fa.gz")), 
                                format = "fasta", 
                                compress = TRUE)
  }
  fp <- file.path(out_dir, paste0(prefix, ".fa.gz"))
  
  # Return
  return(file.path = fp)
  if(return_reads) return(list(combined_reads = fasta), file.path = fp)
}