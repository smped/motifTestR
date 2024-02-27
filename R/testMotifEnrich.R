#' @title Test motif enrichment using a background set of sequences
#'
#' @description
#' Test for motif enrichment within a set of sequences using a background set to
#' derive a NULL hypothesis
#'
#' @details
#' This function offers two alternatives for assessing the enrichment of a motif
#' within a set of sequences, in comparison to a background set of sequences.
#' Enrichment using an iterative model assumes that the background sequences
#' can be subset using a column in their `mcols()` element and these are
#' iterated through to derive a mean and sd operating under the null hypothesis.
#' These values are then used in a Z test to provide an estimate of motif
#' enrichment.
#' Whilst relying on few distributional assumptions and relying on the Central
#' Limit Theorem for the Z score, this can be a time consuming and
#' computationally demanding task.
#'
#' Testing on real datasets has revealed the the background set of counts
#' commonly follows a Poisson distribution across multiple iterations.
#' As a less time-consuming and computationally-demanding approach, is also
#' offered as a testing strategy.
#' No iterations are required, but a background set of sequences which is far
#' greater than the sequences of interest is advised, with background sets
#' >1000 times larger being preferred.
#'
#'
#' @return
#' A data.frame with columns: `sequences`, `matches`, `expected`, `enrichment`,
#' `Z` and `p`, with additional columns `est_bg_rate` (Poisson) or `sd_bg`,
#' `n_iter` and `perm_p` (Iterations).
#' The numbers of sequences and matches refer to the test set of sequences,
#' whilst expected is the expected number of matches under the Poisson or
#' iterative null distribution. The ratio of matches to expected is given as
#' `enrichment`, along with the Z score an p-value. The Z score is only
#' representative under the Poisson model, but is used to estimate p-values
#' under the iterative approach.
#' Under this latter approach, the sd of the matches found in the background
#' sequences is also given, along with the number of iterations and the p-values
#' from permutations testing the one-sided hypothesis hypothesis for enrichment.
#'
#' It may also be worth noting that if producing background sequences using
#' \link{makeRMRanges} with `replace = TRUE` and `force_ol = TRUE`, the
#' iterative model corresponds to a bootstrap, given that the test sequences
#' will overlap the background sequences and background ranges are able to be
#' sampled with replacement.
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet with equal sequence widths
#' @param bg An XStringSet with the same sequence widths as the test XStringset
#' @param model The model used for analysis
#' @param var A column in the mcols element of bg, usually denoting an iteration
#' number
#' @param mc.cores Passed to \link[parallel]{mclapply}
#' @param ... Passed internally to \link[Biostrings]{countPWM}
#'
#' @examples
#' ## Load the example peaks
#' data("ar_er_peaks")
#' sq <- seqinfo(ar_er_peaks)
#' ## Now sample size-matched ranges 10 times larger. In real-world analyses,
#' ## this set should be sampled as at least 1000x larger, ensuring features
#' ## are matched to your requirements
#' set.seed(305)
#' bg_ranges <- makeRMRanges(ar_er_peaks, GRanges(sq)[1], n_iter = 10)
#'
#' ## Convert ranges to DNAStringSets
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' test_set <- getSeq(genome, ar_er_peaks)
#' bg_set <- getSeq(genome, bg_ranges)
#' ## Remove sequences with Ns to avoid annoying messages
#' bg_set <- bg_set[vapply(bg_set, hasOnlyBaseLetters, logical(1))]
#'
#' ## Test for enrichment of the ESR1 motif
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#' testMotifEnrich(esr1, test_set, bg_set, model = "poisson")
#'
#'
#' @export
testMotifEnrich <- function(
    pwm, stringset, bg, model = c("poisson", "iteration"), var = "iteration",
    mc.cores = 1, ...
) {

  ## Checks
  stopifnot(is(bg, "XStringSet"))
  model <- match.arg(model)
  args <- c(as.list(environment()), list(...))
  out <- NULL
  if (is.matrix(pwm)) return(do.call(".testSingleMotifEnrich", args))
  if (is.list(pwm)) {
    pwm <- .cleanMotifList(pwm)
    out <- lapply(
      pwm, .testSingleMotifEnrich, stringset = stringset, bg = bg,
      model = model, var = var, mc.cores = mc.cores
    )
    out <- do.call("rbind", out)
  }

  out

}

#' @keywords internal
.testSingleMotifEnrich <- function(
    pwm, stringset, bg, model, var, mc.cores = 1, ...
){

  ## Find the matches in the test set. This will also perform tests on inputs
  matches <- countPwmMatches(pwm, stringset, ...)
  n <- length(stringset)
  w <- unique(width(stringset))
  if (length(w) != 1) stop("All sequences are required to be the same length")
  bg_w <- width(bg)
  if (!all(bg_w == w))
    stop("All background sequences must be the same length as the test sequences")

  if (model == "poisson") {
    out <- .testPois(pwm, bg, matches, n, ...)
  }
  if (model == "iteration") {
    out <- .testIter(pwm, bg, var, mc.cores, matches, n, ...)
  }

  out

}

#' @importFrom stats poisson.test
#' @keywords internal
.testPois <- function(pwm, bg, matches, n, ...){

  bg_rate <- countPwmMatches(pwm, bg, ...) / length(bg)
  p <- poisson.test(matches, n, bg_rate)$p.value
  expected <- bg_rate * n
  Z  <- (matches - expected) / sqrt(expected) ## Assumes a strict poisson
  data.frame(
    sequences = n, matches, expected, enrichment = matches / expected,
    Z, p, est_bg_rate = bg_rate
  )

}

#' @importFrom parallel mclapply
#' @importFrom stats pchisq sd
#' @keywords internal
.testIter <- function(pwm, bg, var, mc.cores, matches, n, ...) {

  stopifnot(var %in% colnames(mcols(bg)))
  bglist <- split(bg, mcols(bg)[[var]])
  if (!all(vapply(bglist, length, integer(1)) == n))
    stop("All iterations should be the same size as the test sequences")
  bg_matches <- mclapply(
    bglist, \(x) countPwmMatches(pwm, x, ...), mc.cores = mc.cores
  )
  bg_matches <- as.integer(unlist(bg_matches))
  mean_bg <- mean(bg_matches)
  sd_bg <- sd(bg_matches)
  n_iter <- length(bg_matches)
  Z <- (matches - mean_bg) / sd_bg
  perm_p <- (sum(bg_matches > matches) + 1) / n_iter
  p <- 1 - pchisq(Z^2, 1)
  data.frame(
    sequences = n, matches, expected = mean_bg, enrichment = matches / mean_bg,
    Z, p, sd_bg, n_iter, perm_p
  )

}
