## This R script contains processing fasta reference
## in order to determine telomere motif for pattern
## matching in fastq files

## Determine telomere motifs
## This function will determine telomere motifs from a reference
## fasta file. The function will use MEME/STREME to determine
## the motifs.
#' @title
#' Determine telomere motifs from reference fasta file
#' 
#' @description 
#' Identifies telomere motifs from a reference FASTA file and generates patterns for matching in FASTQ files.
#'
#' @details
#' The 'telomere_motif' function identifies telomere motifs from a reference DNA sequence dataset and 
#' generates search patterns for FASTQ files. It first validates the input data and extracts the telomere 
#' regions from the start and end of each sequence. The function then uses the MEME tool to discover telomere 
#' motifs and processes these motifs into various forms (complement, reverse, etc.), adding known ONT 
#' basecalling errors if applicable. It repeats and extends the motifs to create search patterns for 
#' FASTQ files. The function returns a list with the search patterns ('grep_list'), identified motifs 
#' ('motifs'), and MEME analysis results ('results').
#' @param data_ref DNAStringSet reference data.
#' @param reference_telo_length numeric predicted length of telomeres
#' present in the reference file.
#' @param telo_motif_length numeric length of telomere motif.
#' @param number_of_motifs numeric number of motifs to identify.
#' @param number_of_repeats numeric number of repeats to identify.
#' @param threads numeric number of threads MEME will use.
#' @param platform character sequencing platform used (PB | ONT).
#' @param seed numeric seed for random number generation.
#' @return list object which contains the list elements grep_list and
#' motifs
#' @importFrom Biostrings DNAString DNAStringSet DNAStringSetList complement reverse reverseComplement readDNAStringSet width
#' @import stringi
#' @import memes
#' @import universalmotif
#' @export
telomere_motif <- function(reference_file = NULL,
                           reference_telo_length = 2000,
                           telo_motif_length = 6,
                           number_of_motifs = 5,
                           number_of_repeats = 3, 
                           threads = 1,
                           platform = "PB",
                           seed = 9273) {
  # Check inputs
  stopifnot(is.character(reference_file),
            is.numeric(reference_telo_length),
            is.numeric(telo_motif_length),
            is.numeric(number_of_motifs),
            is.numeric(number_of_repeats),
            is.numeric(threads),
            is.character(platform),
            is.numeric(seed))

  # Assign reference sequence after check
  data_ref <- Biostrings::readDNAStringSet(reference_file)
  ref.seq <- data_ref

  # Print data is a DNAStringSet
  message("\nSuccessfully loaded DNAStringSet reference.\n")

  # Check if reference data was loaded
  if (length(DNAStringSet(ref.seq)) == 0) {
    stop("No reference data was loaded")
  }
  # Print reference was not empty
  message("\nReference data is not empty.\n")

  # Create only 1000bp ends of each sequence in reference file
  ref.seq.l <- DNAStringSet(ref.seq,
    start = 1,
    end = reference_telo_length,
    use.names = TRUE
  )
  ref.seq.r <- DNAStringSet(ref.seq,
    start = (width(ref.seq) -
      (reference_telo_length - 1)),
    end = width(ref.seq),
    use.names = TRUE
  )

  ref.seq.ends <- xscat(ref.seq.l, ref.seq.r)
  names(ref.seq.ends) <- names(data_ref)

  # MEME/STREME method
  # conda install conda-forge::ghostscript also dependency
  ## Maybe try installing meme from source! also uninstall from conda
  # Run motif analysis
  seed <- 3873
  ctrl <- universalmotif::create_sequences(
    alphabet = "DNA",
    seqnum = 1000, 
    seqlen = reference_telo_length,
    rng.seed = seed, 
    nthreads = threads
  )

  message("\nRunning MEME analysis with ", threads, 
          " threads.\n")
  meme_out <- memes::runMeme(
    input = ref.seq.ends,
    control = ctrl,
    seed = seed,
    objfun = "de",
    parse_genomic_coord = FALSE,
    dna = TRUE,
    nmotifs = 100,
    evt = 0.1,
    minw = telo_motif_length,
    maxw = telo_motif_length * 2,
    minsites = 5,
    maxsites = 50000,
    norand = TRUE,
    mod = "anr",
    revcomp = TRUE,
    p = paste(threads, "--oversubscribe"), 
    silent = FALSE
  )

  # Motifs
  ord <- meme_out[order(meme_out$eval), ]
  motifs <- ord$consensus[1:number_of_motifs]
  motifs <- motifs[!is.na(motifs)]
  if (platform == "ONT") {
      motifs <- c(motifs, "CCCTGG")
  }

  message("Motifs identified: ", 
         "\n\t", paste(motifs, collapse = "\n\t"))

  # Garbage Collection
  gc()

  ## !-- Probably don't need stringi package --! ##
  comp <- lapply(motifs, function(mot) {
    dna <- Biostrings::DNAString(mot)
    comp <- Biostrings::complement(dna)
    return(as.character(comp))
  })

  revComp <- lapply(motifs, function(mot) {
    dna <- Biostrings::DNAString(mot)
    rev <- Biostrings::reverseComplement(dna)
    return(as.character(rev))
  })

  rev <- lapply(motifs, function(mot) {
    dna <- Biostrings::DNAString(mot)
    rev <- Biostrings::reverse(dna)
    return(as.character(rev))
  })

  revRevComp <- lapply(revComp, function(mot) {
    dna <- Biostrings::DNAString(mot)
    rev <- Biostrings::reverseComplement(Biostrings::reverse(dna))
    return(as.character(rev))
  })

  # Correct for L telomeres, I don't believe Biostrings::countPattern is
  # actually searching for all motifs
  all_motifs <- c(motifs, comp, rev, revComp, revRevComp)

  mult <- function(in_mot) {
    if (nchar(in_mot) <= telo_motif_length) {
      paste0(rep(in_mot, number_of_repeats), collapse = "")
    } else if (nchar(in_mot) > telo_motif_length) {
      chars <- number_of_repeats * telo_motif_length
      m <- ceiling(chars / nchar(in_mot))
      paste0(rep(in_mot, m), collapse = "")
    }
  }

  orig <- lapply(motifs, mult)
  comp <- lapply(comp, mult) 
  rev_comp <- lapply(revComp, mult)
  rev <- lapply(rev, mult)
  rev_rev_comp <- lapply(revRevComp, mult)
  
  list <- c(original = orig,
            complement = comp, 
            reverse = rev, 
            reverse_complement = rev_comp, 
            rev_reverse_complement = rev_rev_comp)
  grep_list <- list

  # Return grep pattern
  return(list(grep_list = grep_list,
              motifs = as.character(all_motifs),
              results = meme_out))
}
