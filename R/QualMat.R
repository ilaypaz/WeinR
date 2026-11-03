#' Compute and store per-file quality-control metrics
#'
#' Internal helper that computes read-level and file-level summary metrics from a
#' `DNAStringSet`/`QualityScaledDNAStringSet`-like input and stores them into a
#' QC object. Poly-A and base-modification results that were already computed
#' elsewhere (the chunked BAM pass) can be passed in and are stored alongside.
#'
#' @param qc_obj An S4 QC object with at least `@metrics` and `@summary_metrics`
#'   slots. This function appends metrics for `filename` into these slots.
#' @param stringset A `Biostrings` string set containing the sequences (and
#'   optionally qualities, depending on object type).
#' @param filename Character scalar naming the file being summarized (used as a
#'   key in `qc_obj@metrics` and in the summary table).
#' @param polyA_lengths Optional numeric vector of valid poly-A tail lengths
#'   (output of the chunking code).
#' @param modification_density Optional named list (one element per modification
#'   type) of per-read modification density, calls per 100 bp.
#' @param read_means Optional named list of per-read mean ML probabilities (0-1)
#'   per modification type.
#' @param type_values Optional named list of raw ML values per modification type.
#' @param average_ML_file Optional named list of file-level mean ML probability
#'   (0-1) per modification type.
#'
#' @return The updated `qc_obj`.
#'
#' @keywords internal
QualMat <- function(qc_obj, stringset, filename,
                    polyA_lengths = NULL,
                    modification_density = NULL,
                    read_means = NULL,
                    type_values = NULL,
                    average_ML_file = NULL){

  seqNames <- names(stringset)

  #calculate metrics

  #LENGTH
  lengths <- Biostrings::width(Biostrings::DNAStringSet(stringset))
  yield <- sum(lengths)

  #calculate base frequency
  af <- Biostrings::alphabetFrequency(Biostrings::DNAStringSet(stringset))

  #N content per read
  N <- af[,'N']

  #GC CONTENT
  #per read
  perReadGC <- (af[,'G'] + af[,'C'])/lengths

  #per file
  perFileGC <- sum(af[,'G']+ af[,'C'])/yield

  rm(af)

  #per position gc content
  # dataframe, 4 columns for each base: each has proportion of given base
  # row for each read

  #Quality scores
  #per file
  qual_list <- lapply(as.character(Biostrings::quality(stringset)), function(q) as.integer(utf8ToInt(q) - 33))

  avgQscore <- round(mean(unlist(qual_list)),2)

  #per read
  perReadQscore <- lapply(qual_list, function(q) mean(q))

  #per position
  q_stats <- chunked_quality_per_position(qual_list,lengths, chunk_size = 1000)

  rm(qual_list)

  #calculate summary metrics
  #N50
  decreasinglengths <- sort(lengths, decreasing = TRUE)
  tmp <- cumsum(decreasinglengths)
  N50 <- decreasinglengths[which(tmp >= yield/2)[1]]



  #fill object

  # add metrics
  qc_obj@metrics[[filename]] <- list(
    readLengths = lengths,
    meanprQscore = perReadQscore,
    perPosQuality = q_stats,
    Ncount = as.list(N),
    prGCcontent = as.list(perReadGC),
    # poly-A and base-modification results from the chunking code (NULL if not supplied)
    polyALengths = polyA_lengths,
    modDensity = modification_density,
    modReadMeans = read_means,
    modValues = type_values,
    avgML = average_ML_file
  )

  # Add summary metrics
  summary_df <- data.frame(
    File = filename,
    `Read Count` = length(lengths),
    `Yield (bp)` = yield,
    N50 = N50,
    `Average Per Read Q Score` = avgQscore,
    `N Count` = sum(N),
    stringsAsFactors = FALSE
  )
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)

  #return the object
  base::return(qc_obj)
}
