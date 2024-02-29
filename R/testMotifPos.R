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
#' Motif alignments are assigned into bins based on the central position of the
#' match, as provided in the column `from_centre` when calling
#' \link{getPwmMatches}.
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
#' @param pwm A Position Weight Matrix, or list of PWMs. Not required if `bm`
#' is supplied.
#' @param stringset An XStringSet. Not required if `bm` is supplied
#' @param bm An optional set of 'best matches' as returned by
#' \link{getPwmMatches} setting `best_only = TRUE`. Alternatively, a list of
#' returned 'best matches' can be used. Any individual sequence with multiple
#' 'best matches' will have each match weighted. If provided, will override
#' anything passed via `pwm` or `stringset.`
#' @param binwidth Width of bins across the range to group data into
#' @param abs Use absolute positions around zero to find symmetrical enrichment
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param break_ties Choose how to resolve matches with tied scores
#' @param alt Alternative hypothesis for the binomial test
#' @param mc.cores Passed to \link[parallel]{mclapply}
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#'
#' ## Load the example Peaks
#' data("ar_er_peaks")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' seq <- getSeq(genome, ar_er_peaks)
#'
#' ## Get the best match and use this data
#' bm <- getPwmMatches(esr1, seq, best_only = TRUE)
#' ## Test for enrichment in any position
#' testMotifPos(bm = bm)
#'
#' ## Provide a list of PWMs, testing for distance from zero
#' testMotifPos(ex_pwm, seq, abs = TRUE, binwidth = 10)
#'
#'
#' @importFrom parallel mclapply
#' @importFrom stats p.adjust
#' @export
testMotifPos <- function(
    pwm, stringset, bm, binwidth = 10, abs = FALSE, rc = TRUE,
    min_score = "80%", break_ties = "random",
    alt = c("greater", "less", "two.sided"), mc.cores = 1, ...
) {

  alt <- match.arg(alt)
  args <- c(as.list(environment()), list(...))
  args <- args[!names(args) %in% c("pwm", "stringset", "mc.cores")]
  ## Handle any situation without a bm df or list
  if (missing(bm)) {
    ## This will perform all checks as well
    ## pwm & stringset are not required beyond this initial call
    bm <- getPwmMatches(
      pwm, stringset, rc, min_score, best_only = TRUE, break_ties, mc.cores, ...
    )
  }
  .checkInputBM(bm)
  cols <- c(
    "start", "end", "centre", "width", "total_matches", "matches_in_region",
    "expected", "enrichment", "prop_total", "p", "consensus_motif"
  )
  if (is(bm, "DataFrame")) {
    args$bm <- bm
    out <- do.call(.testSingleMotifPos, args)
  } else {
    out <- mclapply(
      bm, .testSingleMotifPos, binwidth = binwidth, abs = abs, rc = rc,
      min_score = min_score, break_ties = break_ties, alt = alt, ...,
      mc.cores = mc.cores
    )
    cols <- c(cols[1:10], "fdr", "consensus_motif")
    out <- do.call("rbind", out)
    out$fdr <- p.adjust(out$p, "fdr")
  }

  out[, cols]

}

#' @importFrom harmonicmeanp p.hmp
#' @importFrom stats binom.test
#' @keywords internal
.testSingleMotifPos <- function(
    bm, binwidth , abs, rc, min_score, break_ties, alt, ...
) {

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

