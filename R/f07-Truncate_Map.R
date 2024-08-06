## Truncate reads at the read level
#' @title Truncate reads
#' @description Truncate reads based on the telomere length and hit motif.
#' @param reads A DNAStringSet object containing reads to be truncated
#' @param result_df A data frame containing the telomere length and hit motif
#' information.
#' @return A DNAStringSet object containing truncated reads.
#' @param reads A DNAStringSet object containing reads to be truncated
#' @importFrom Biostrings subseq 
#' @export
truncation <- function(reads, result_df) {

  # Get initial values
  nms <- names(reads)
  rows <- result_df[rownames(result_df) %in% nms,]
  st <- rows$trunc_start
  en <- rows$trunc_end

  truncated_sequences <- Biostrings::subseq(reads, start = as.numeric(st), 
                           end = as.numeric(en))
  names(truncated_sequences) <- nms
 
  return(truncated_sequences)
}


## Truncate reads at the file level
#' @title Truncate files
#' @description Truncate reads based on the telomere length and hit motif.
#' @param combined_telomere_file A character string specifying the combined
#' telomere file.
#' @param results_data_frame A data frame containing the telomere length and hit
#' motif information.
#' @param out_dir A character string specifying the output directory.
#' @param prefix A character string specifying the prefix for the output files.
#' @param write_file A logical value specifying whether to write the output
#' files to disk. Default is TRUE.
#' @param progress A logical value specifying whether to display progress
#' messages. Default is TRUE.
#' @param return A logical value specifying whether to return the truncated
#' reads. Default is FALSE.
#' @param verbose A logical value specifying whether to display verbose
#' messages. Default is TRUE.
#' @return A list containing the truncated reads and the truncated file list.
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @export
truncate_files <- function(combined_telomere_file = NULL,
                           results_data_frame = NULL,
                           out_dir = file.path(getwd(), "truncated_reads"),
                           prefix = "truncated",
                           write_file = TRUE,
                           progress = TRUE,
                           return = FALSE,
                           verbose = TRUE) {

  # Check inputs
  stopifnot(
    file.exists(combined_telomere_file),
    is.data.frame(results_data_frame),
    is.character(out_dir),
    is.character(prefix),
    is.logical(progress),
    is.logical(return),
    is.logical(verbose)
  )

  # Read in the combined telomere file
  telomeres <- Biostrings::readDNAStringSet(combined_telomere_file)
  # telomeres <- bio_ul(telomeres)

  # Cross telos to throw out NAs
  ind <- which(!is.na(results_data_frame$telomere_end))
  names <- rownames(results_data_frame)[ind]
  rdf <- results_data_frame[ind, ]
  select <- names(telomeres) %in% names
  telomeres <- telomeres[select]

  ret_ls <- list()
  for (iter in seq_along(telomeres)) {
    ret_ls[[iter]] <- truncation(telomeres[iter], results_data_frame)
  }

  ret_ls <- bio_ul(ret_ls)

  if (write_file) {
    # Save files
    if (verbose) {
      message("\nWriting files to: ", out_dir)
    }
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    Biostrings::writeXStringSet(ret_ls,
      file.path(out_dir, paste0(
        prefix, "_", "reads",
        ".fa"
      )),
      format = "fasta",
      compress = TRUE
    )
  }

  # Create files_list
  files_list <- list.files(out_dir, full.names = TRUE)

  # Return
  if (return) {
    return(list(
      truncated_reads = ret,
      truncated_file_list = files_list
    ))
  }
  return(list(truncated_files_list = files_list))
}

## Map reads to reference
#' @title Map reads to reference
#' @description Map reads to a reference genome using minimap2.
#' @param fasta A character string specifying the path to the fasta file.
#' @param results_data_frame A data frame containing the telomere length and hit
#' motif information.
#' @param reference_file A character string specifying the path to the reference
#' genome.
#' @param preset_string A character string specifying the preset to use for
#' minimap2. Default is "map-hifi".
#' @param out_dir A character string specifying the output directory.
#' @param prefix A character string specifying the prefix for the output files.
#' @param threads An integer specifying the number of parallel processors to use.
#' Default is 1.
#' @param return_mapped A logical value specifying whether to return the mapped
#' reads. Default is TRUE.
#' @param underscore A logical value specifying whether to use underscores in
#' the chromosome names. Default is TRUE.
#' @param verbose A logical value specifying whether to display verbose messages.
#' Default is TRUE.
#' @return A list containing the mapped reads and the results data frame.
#' @importFrom dplyr distinct left_join
#' @importFrom Biostrings readDNAStringSet width
#' @import stringr
#' @export
map <- function(fasta = NULL, 
                results_data_frame = NULL, 
                reference_file = NULL,
                preset_string = "map-hifi",
                out_dir = file.path(getwd(), "mapped_reads"),
                prefix = "mapped",
                threads = 1, 
                return_mapped = TRUE,
                underscore = TRUE,
                verbose = TRUE) {

  # Check inputs
  stopifnot(is.character(fasta), 
            is.data.frame(results_data_frame),
            is.character(reference_file),
            is.character(preset_string),
            is.character(out_dir),
            is.character(prefix),
            is.numeric(threads),
            is.logical(return_mapped),
            is.logical(verbose))

  # Create output directory if it doesn't exist
  dir.create(out_dir, showWarnings = FALSE)
  
  mapped_reads <- minimapR::minimap2(reference = reference_file, 
                                     query = fasta, 
                                     preset = preset_string, 
                                     output_file_prefix = prefix,
                                     return = return_mapped,
                                     threads = threads, 
                                     verbose = verbose, 
                                     ... = paste0("--secondary=no"))

  # Mapped reads, dplyr deduplication based on qname
  mapped_reads <- mapped_reads %>% 
    dplyr::distinct(qname, .keep_all = TRUE)

  # Add telo_name to mapped reads
  results_data_frame$telo_name <- rownames(results_data_frame)

  # Create a bam data frame with only pertinent columns
  mapped_reads <- as.data.frame(mapped_reads)
  bam_df <- data.frame(telo_name = mapped_reads$qname,
                       ref_name = mapped_reads$rname,
                       strand = mapped_reads$strand,
                       start = mapped_reads$pos,
                       end = mapped_reads$pos + mapped_reads$qwidth - 1,
                       qwidth = mapped_reads$qwidth,
                       score = mapped_reads$mapq)
                       
  # Merge with results data frame
  results <- dplyr::left_join(bam_df, results_data_frame, by = "telo_name")

  # Add chromosomal end
  ## Load reference to get lengths
  ref <- Biostrings::readDNAStringSet(reference_file)
  # Only take the character string before the first space
  if (underscore) {
    chromosomes <- gsub("_.*", "", names(ref))
    results$ref_name <- gsub("_.*", "", results$ref_name)
    lengths <- data.frame(chromosome = chromosomes,
      length = Biostrings::width(ref))
    lengths <- distinct(lengths, chromosome, .keep_all = TRUE)
    ord <- order(as.numeric(stringr::str_extract(lengths$chromosome, "\\d+")))
    lengths <- lengths[ord,]
  } else {
    chromosomes <- gsub(" .*", "", names(ref))
    lengths <- data.frame(chromosome = chromosomes,
      length = Biostrings::width(ref))
  }


  # Clean up some
  rm(bam_df, ref)
  gc()

  # Assign chromosomal end
  # Remove NAs
  results <- results[!is.na(results$ref_name),]
  for(i in seq_len(nrow(results))) {
    chrom <- results$ref_name[i]
    chrom_length <- lengths[lengths$chromosome == chrom, "length"]
    pos <- results$start[i]
    end <- NULL
    l <- abs(pos - 0)
    r <- abs(pos - chrom_length)
    if(l < r) {
      end <- "5'"
    } else {
      end <- "3'"
    }
    results$chromEnd[i] <- end
  }

  # Return
  if(return_mapped) {
    return(list(mapped_reads = mapped_reads, 
                results = results))
  } else {
    return(list(results = results))
  }

}