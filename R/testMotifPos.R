#' Test for a Uniform Distribution across a set of best matches
#'
#' @details
#' This function tests for an even positional spread of motif matches across a
#' set of sequences, using the assumption (i.e. H~0~) that if there is no
#' positional bias, matches will be evenly distributed across all positions
#' within a set of sequences.
#' Conversely, if there is positional bias, typically but not necessarily near
#' the centre of a range, this function intends to detect this signal, as a
#' rejection of the null hypothesis.
#'
#' Input can be provided as the output from \link{getPwmMatches} setting
#' `best_only = TRUE` if these matches have already been identified.
#' If choosing to provide this object to the argument `bm`, nothing is required
#' for the arguments `pwm`, `stringset`, `rc`, `min_score` or `break_ties`
#' Otherwise, a Position Weight Matrix (PWM) and an `XStringSet` are required,
#' along with the relevant arguments, with best matches identified within the
#' function.
#'
#' The set of best matches are then grouped into bins along the range, with the
#' central bin containing zero, and tallied.
#' Setting `abs` to `TRUE` will set all positions from the centre as
#' *absolute values*, returning counts purely as bins with distances from zero,
#' marking this as an inclusive lower bound.
#'
#' The \link[stats]{binom.test} is performed on each bin using the alternative
#' hypothesis, with the returned p-values across all bins combined using the
#' Harmonic Mean p-value (HMP) (See \link[harmonicmeanp]{p.hmp}).
#' All bins with raw p-values below the HMP are identified and the returned
#' values for start, end, centre, width, matches in region, expected and
#' enrichment are across this set of bins.
#' The expectation is that where a positional bias is evident, this will be a
#' narrow range containing a non-trivial proportion of the total matches.
#'
#'
#' @return
#' A data.frame with columns `start`, `end`, `centre`, `width`, `total_matches`,
#' `matches_in_region`, `expected`, `enrichment`, `prop_total`, `p`
#' and `consensus_motif`
#' The total matches represent the total number of matches within the set of
#' sequences, whilst the number observed in the final region are also given,
#' along with the proportion of the total this represents.
#' Enrichment is simply the ratio of observed to expected based on the
#' expectation of the null hypothesis
#'
#' The consensus motif across all matches is returned as a Position Frequency
#' Matrix (PFM) using \link[Biostrings]{consensusMatrix}.
#'
#'
#' @param pwm A Position Weight Matrix. Not required if `bm` is supplied.
#' @param stringset An XStringSet. Not required if `bm` is supplied
#' @param bm An optional set of 'best matches' as returned by
#' \link{getPwmMatches} setting `best_only = TRUE`. Any sequences with multiple
#' 'best matches' will have each match weighted. If provided, will override
#' anything passed via `pwm` or `stringset.`
#' @param binwidth Width of bins across the range to group data into
#' @param abs Use absolute positions around zero to find symmetrical enrichment
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param break_ties Choose how to resolve matches with tied scores
#' @param alt Alternative hypothesis for the binomial test
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#'
#' @importFrom harmonicmeanp p.hmp
#' @importFrom stats binom.test
#' @export
testMotifPos <- function(
    pwm, stringset, bm, binwidth = 10, abs = FALSE, rc = TRUE,
    min_score = "80%", break_ties = "random",
    alt = c("greater", "less", "two.sided"), ...
) {

  if (missing(bm)) {
    ## This will perform all checks as well
    ## pwm & stringset are not required beyond this initial call
    bm <- getPwmMatches(
      pwm, stringset, rc, min_score, best_only = TRUE, break_ties, ...
    )
  }
  stopifnot(is(bm, "DataFrame"))

  ## Setup a well formed output object for easy combining with other tests
  out_cols <- c(
    "start", "end", "centre", "width", "total_matches", "matches_in_region",
    "expected", "enrichment", "prop_total", "p", "consensus_motif"
  )
  out <- lapply(out_cols, \(x) integer())
  out$consensus_motif <- list()
  names(out) <- out_cols
  n_matches <- nrow(bm)
  if (n_matches == 0) return(data.frame(out))

  ## The required columns for all downstream analysis are
  reqd_cols <- c("seq", "from_centre", "seq_width")
  stopifnot(all(reqd_cols %in% colnames(bm)))
  bm <- .makeBmBins(bm, binwidth, abs)
  bins <- levels(bm$bin)

  ## Get the counts & form a df with every bin as a row
  counts <- vapply(splitAsList(bm, bm$bin), \(x) sum(x$weight), numeric(1))
  df <- data.frame(bin = names(counts), matches_in_bin = as.integer(counts))
  bin_coords <- do.call("rbind", strsplit(df$bin, split = ","))
  bin_coords <- apply(bin_coords, 2, \(x) as.integer(gsub("\\[|\\(|\\]", "", x)))
  df$start <- bin_coords[,1]
  df$end <- bin_coords[,2]
  df$bin <- factor(df$bin, levels = bins)
  df$width <- df$end - df$start + 1
  df$bin_prob <- df$width / sum(df$width)
  df$expected <- n_matches * df$bin_prob
  alt <- match.arg(alt)
  df$p <- vapply(
    seq_along(df$bin),
    \(i) {
      x <- as.list(df[i,])
      binom.test(x$matches_in_bin, n_matches, x$bin_prob, alt)$p.value
    }, numeric(1)
  )

  ## Summarise the output selecting rows with p < hmp
  hmp <- p.hmp(df$p, L = length(bins))
  df <- subset(df, df$p < hmp | df$p == min(df$p))
  out <- data.frame(start = min(df$start), end = max(df$end))
  out$centre <- (out$start + out$end) / 2
  out$width <- out$end - out$start
  out$total_matches <- n_matches
  out$matches_in_region <- sum(df$matches_in_bin)
  out$expected <- sum(df$expected)
  out$enrichment <- out$matches_in_region / out$expected
  out$prop_total <- out$matches_in_region / out$total_matches
  out$p <- as.numeric(hmp)

  ## Setup the consensus motif to return
  motif_cols <- strsplit(consensusString(bm$match), "")[[1]]
  consensus <- consensusMatrix(bm$match)[c("A", "C", "G", "T"),]
  colnames(consensus) <- motif_cols
  out$consensus_motif <- list(consensus)
  out[, out_cols]

}

#' @keywords internal
.makeBmBins <- function(bm, binwidth, abs){

  ## Set sequence weights for multiple matches
  reps <- table(as.character(bm$seq)) # Counts reps
  bm$weight <- as.numeric(1 / reps[as.character(bm$seq)]) # 1 / reps

  ## Make the bins
  longest_seq <- max(bm$seq_width)
  stopifnot(binwidth < longest_seq)
  if (abs) {
    ## Set the dist from centre to be +ve only
    bm$from_centre <- abs(bm$from_centre)
    all_breaks <- c(seq(0, longest_seq / 2, by = binwidth), longest_seq / 2)
    seq_per_bin <- vapply(all_breaks, \(x) sum(bm$seq_width >= x), integer(1))
  } else {
    ## To define bins, make sure the central bin is around zero and that
    ## bins extends symetrically from this point to the widest sequence
    pos_breaks <- c(
      ## Ensure the upper limit is included, even if it leaves a smaller bin
      seq(binwidth / 2, longest_seq / 2, by = binwidth), longest_seq / 2
    )
    all_breaks <- c(-1 * pos_breaks, pos_breaks)
    pos_bins <- vapply(pos_breaks, \(x) sum(bm$seq_width >= x), integer(1))
    seq_per_bin <- c(pos_bins, pos_bins)
  }

  ## Cut into bins
  all_breaks <- unique(sort(all_breaks))
  bm$bin <- cut(bm$from_centre, breaks = all_breaks, include.lowest = TRUE)
  bm
}
