## This R script contains helper functions for processing
## fastq files
## First function will be a wrapper for the seqkit split2
## Call fastqsplitter with options with system call
#' @title Divide a FASTQ File into Smaller Files
#'
#' @description Split fastq file into smaller files by number of reads
#' in order to allow R instances to process the data without
#' segmentation faults 
#'
#' @details The 'split_fastq' function is designed to divide a large 
#' FASTQ file into smaller files based on a specified size and includes 
#' options for multithreading. It constructs a 
#' command-line call to the 'seqkit split2' tool, which is a part of the 
#' SeqKit toolkit for processing sequencing data. The function takes as 
#' input the path to the FASTQ file to be split, the desired size 
#' for each split file, a prefix for naming the split files, the 
#' number of threads to use for parallel processing, and the directory 
#' where the split files will be saved. An optional file extension can 
#' also be specified for the split files. The function constructs the 
#' appropriate command with these parameters, prints the command for user 
#' verification, and then executes it using the system2 function to run the 
#' command on the system’s command line. The resulting output consists of 
#' multiple smaller FASTQ files, each with the specified size, which are 
#' saved in the designated output directory.
#'
#' @param in_fastq character path to input fastq file.
#' @param size numeric number of reads to split larger 
#' fastq file into smaller files.
#' @param prefix character prefix for output files.
#' @param threads numeric number of threads to use.
#' @param out_dir character path to output directory.
#' @param extension character extension for output files.
#' @return NULL, does not have a return value
#' @export
split_fastq <- function(in_fastq, size = 284289, prefix = "split_",
                        threads = 1, out_dir = paste0(dirname(in_fastq), "/split"),
                        extension = NULL) {
  # Create output directory if it doesn't exist
  dir.create(out_dir, showWarnings = FALSE)

  # Create command
  cmd <- paste0("seqkit")
  args <- c("split2", in_fastq, "-s", size, 
            "--threads", threads, "--by-size-prefix", prefix,
            "--out-dir", out_dir)
  if (!is.null(extension)) args <- c(args, "--extension", extension)
  # Print command
  cat("\nThis is the command for system call: \n\t", cmd, args, "\n\n")
  # Run command
  system2(command = cmd, args = args, stdout = TRUE)
}

## Need to determine basic stats for submitted fastq file
# Use system2 call to seqkit stats
#' @title Generate FASTQ File Statistics and Splitting Plan
#' @description Generate stats about intial fastq file and determine
#' how to split the file into smaller files based on RAM constraints
#' @details The 'fastq_stats' function analyzes a FASTQ file 
#' to generate basic statistics and estimates a file-splitting 
#' strategy based on target RAM size. It begins by constructing 
#' and executing a command-line call to the seqkit stats tool, 
#' which collects statistics such as the total number of sequences 
#' (reads) and the total length of sequences (bases). The function 
#' parses this information to extract the number of reads and 
#' the sum of base pairs, then calculates the total gigabases of 
#' data. Using the target RAM size and a predefined ratio for RAM 
#' usage efficiency, the function estimates the optimal number of 
#' smaller FASTQ files required to stay within the specified RAM 
#' limit and determines how many reads each of these files should 
#' contain. Finally, the function returns a list with statistics on 
#' the total number of bases, reads, the estimated number of files, 
#' the number of reads per file, and the total gigabases of data, 
#' which helps in planning the processing of large FASTQ files.
#' 
#' @param in_fastq Path to the input FASTQ file.
#' @param threads Number of threads to use for processing (default is 1).
#' @param target_ram_size_gb Desired RAM size in GB for processing (default is 10 GB).
#' @return Returns a list of base stats, estimated files, and read counts.
#' @export
fastq_stats <- function(in_fastq, threads = 1, 
                        target_ram_size_gb = 10) {
  # Call seqkit stats to get basic statistics
  cmd <- paste0("seqkit")
  args <- c("stats", in_fastq, "--threads ", threads, "--tabular")
  stats <- system2(command = cmd, args = args, stdout = TRUE)
  list <- strsplit(stats, "\t")
  stats_df <- as.data.frame(list[[2]], row.names = list[[1]], 
                            stringsAsFactors = FALSE, 
                            col.names = c("value"))
  bases <- as.numeric(stats_df["sum_len", ])
  reads <- as.numeric(stats_df["num_seqs", ])
  # Bases to Ram size, 1.5 gigabases per gigabyte of ram
  ratio <- 1.07
  gig_bases <- bases / 1e9
  number_files <- floor((gig_bases * ratio) / target_ram_size_gb)
  number_of_reads_per_file <- ceiling(reads / number_files)
  # Return list
  return(list(bases = bases, reads = reads, 
              number_files = number_files, 
              number_of_reads_per_file = number_of_reads_per_file, 
              gigabases = gig_bases))
}

## Create an autosplit function that will automatically determine
## values and cut
# Function to automatically split fastq files
# Returns list of fastq files for further processing
#' @title Divide FASTQ Files According to RAM Constraints
#' @description Automatically split fastq file into smaller files by number of reads
#' in order to allow R instances to process the data without
#' segmentation faults
#' @details The 'auto_split' function automatically splits a FASTQ file into smaller 
#' files based on the available RAM size. It first calculates the number of files needed 
#' and the number of reads per file to ensure that each file's size fits within the 
#' specified 'target_ram_size_gb'. It uses the 'fastq_stats' function to obtain statistics 
#' on the total number of bases and reads in the input FASTQ file and estimates the 
#' number of files required for the target RAM size. If more than one file is needed, 
#' it calls the 'split_fastq' function to perform the actual splitting, applying the 
#' specified 'size', 'prefix', and 'threads' for the operation. The function 
#' finally returns a list containing the file statistics and the paths to the created 
#' split files.
#'
#' @param in_fastq character path to input fastq file.
#' @param size numeric number of reads to split larger file into
#' smaller files.
#' @param prefix character prefix for output files 
#' (file label, default "split_").
#' @param threads numeric number of threads to use.
#' @param out_dir character path to output directory.
#' @param extension character extension for output files (fastq | fastq.gz).
#' @param target_ram_size_gb numeric target size of ram in gigabytes
#' which each loaded file will occupy (threads * target_ram_size_gb) < ram.
#' @return full path list of files generated by auto_split
#' @export
auto_split <- function(in_fastq, size = 284289, prefix = "split_",
                       threads = 1,
                       out_dir = paste0(dirname(in_fastq), "/split"),
                       extension = NULL, target_ram_size_gb = 10) {
  # Print progress
  cat("\nDetermining number of files to split into\n")
  # Get stats
  stats <- fastq_stats(in_fastq, threads = threads, 
                       target_ram_size_gb = target_ram_size_gb)

  # Check if number of files is greater than 1
  if (stats$number_files > 1) {
    # Print message
    cat("\nSplitting fastq file into", stats$number_files, "files\n")
    # Split fastq
    split_fastq(in_fastq, size = stats$number_of_reads_per_file,
                prefix = prefix, 
                threads = threads,
                out_dir = out_dir,
                extension = ".gz")
  } else {
    # Print message
    cat("\nOnly one file would be created, no need to split\n")
  }
  # Return list of files
  file_list <- list.files(path = out_dir, full.names = TRUE,
                          pattern = "split_")
  return(list(stats = stats, file_list = file_list))
}
