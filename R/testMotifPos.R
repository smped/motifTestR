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
#' If choosing to provide this object to the argument `matches`, nothing is required
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
#' @param x A Position Weight Matrix, universalmotif object or list thereof.
#' Alternatively can be a single DataFrame or list of DataFrames as returned
#' by \link{getPwmMatches} with `best_only = TRUE`
#' @param stringset An XStringSet. Not required if matches are supplied as x
#' @param binwidth Width of bins across the range to group data into
#' @param abs Use absolute positions around zero to find symmetrical enrichment
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param break_ties Choose how to resolve matches with tied scores
#' @param alt Alternative hypothesis for the binomial test
#' @param sort_by Column to sort results by
#' @param mc.cores Passed to \link[parallel]{mclapply}
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#'
#' ## Load the example sequences
#' data("ar_er_seq")
#'
#' ## Get the best match and use this data
#' matches <- getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
#' ## Test for enrichment in any position
#' testMotifPos(matches)
#'
#' ## Provide a list of PWMs, testing for distance from zero
#' testMotifPos(ex_pwm, ar_er_seq, abs = TRUE, binwidth = 10)
#'
#'
#' @importFrom parallel mclapply
#' @importFrom stats p.adjust
#' @export
testMotifPos <- function(
    x, stringset, binwidth = 10, abs = FALSE, rc = TRUE, min_score = "80%",
    break_ties = "all", alt = c("greater", "less", "two.sided"),
    sort_by = c("p", "none"), mc.cores = 1,
    ...
) {

  ## Handle single objects by coercing to a list where appropriate
  valid_cl <- c("matrix", "universalmotif", "DataFrame")
  if (any(vapply(valid_cl, \(cl) is(x, cl), logical(1)))) x <- list(x)
  if (!is.list(x))
    stop("Input should be provided as a list, PWM/matrix or DF of matches")

  ## Handle a list of PWMs
  isPWM <- vapply(
    x, \(el) is(el, "matrix") | is(el, "universalmotif"), logical(1)
  )
  if (all(isPWM)) {
    matches <- getPwmMatches(
      x, stringset, rc, min_score, best_only = TRUE, break_ties, mc.cores, ...
    )
  } else {
    ## Now for a list of matches
   .checkMatches(x) ## Will fail if not valid
    matches <- x
  }
  if (missing(matches)) stop("Provided input is not in a recognised format")
  alt <- match.arg(alt)

  cols <- c(
    "start", "end", "centre", "width", "total_matches", "matches_in_region",
    "expected", "enrichment", "prop_total", "p", "fdr", "consensus_motif"
  )

  out <- mclapply(
    matches, .testSingleMotifPos, binwidth = binwidth, abs = abs, rc = rc,
    min_score = min_score, break_ties = break_ties, alt = alt, ...,
    mc.cores = mc.cores
  )
  out <- do.call("rbind", out)
  out$fdr <- p.adjust(out$p, "fdr")
  o <- seq_len(nrow(out))
  sort_by <- match.arg(sort_by)
  if (sort_by != "none") o <- order(out[[sort_by]])
  out[o, cols]

}

#' @importFrom harmonicmeanp p.hmp
#' @importFrom stats binom.test
#' @keywords internal
.testSingleMotifPos <- function(
    matches, binwidth, abs, rc, min_score, break_ties, alt, ...
) {

  ## Setup a well formed output object for easy combining with other tests
  out_cols <- c(
    "start", "end", "centre", "width", "total_matches", "matches_in_region",
    "expected", "enrichment", "prop_total", "p", "consensus_motif"
  )
  out <- lapply(out_cols, \(x) integer())
  out$consensus_motif <- list()
  names(out) <- out_cols
  n_matches <- nrow(matches)
  if (n_matches == 0) return(data.frame(out))
  matches <- .makeBmBins(matches, binwidth, abs)
  bins <- levels(matches$bin)

  ## Get the counts & form a df with every bin as a row
  counts <- vapply(splitAsList(matches, matches$bin), \(x) sum(x$weight), numeric(1))
  df <- data.frame(bin = names(counts), matches_in_bin = as.integer(counts))
  bin_coords <- do.call("rbind", strsplit(df$bin, split = ","))
  bin_coords <- apply(bin_coords, 2, \(x) as.integer(gsub("\\[|\\(|\\]", "", x)))
  df$start <- bin_coords[, 1]
  df$end <- bin_coords[, 2]
  df$bin <- factor(df$bin, levels = bins)
  ## Given matches are measured from the centre, only bases > d_from_edge
  ## away from the ends of the range are viable
  motif_width <- max(width(matches$match))
  d_from_edge <- (motif_width - 1) / 2
  valid_pos <- seq(min(df$start) + d_from_edge, max(df$end) - d_from_edge)
  if (abs) valid_pos <- seq(d_from_edge %% 1, max(df$end) - d_from_edge)
  df$valid_width <- vapply(
    seq_len(nrow(df)),
    \(i) sum(valid_pos > df$start[i] & valid_pos <= df$end[i]), numeric(1)
  )
  df$bin_prob <- df$valid_width / sum(df$valid_width)
  df$expected <- n_matches * df$bin_prob
  df$p <- vapply(
    seq_along(df$bin),
    \(i) {
      x <- df[i,]
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
  consensus <- consensusMatrix(matches$match)[c("A", "C", "G", "T"),]
  out$consensus_motif <- list(consensus)
  out[, out_cols]
}

