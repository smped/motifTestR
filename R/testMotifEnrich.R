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
#' can often follow a Poisson distribution across multiple iterations.
#' As a less time-consuming and computationally-demanding approach, is also
#' offered as a testing strategy.
#' No iterations are required, but a background set of sequences which is far
#' greater than the sequences of interest is advised, with background sets
#' >1000 times larger being preferred.
#'
#' Alternatively, a quasipoisson model is also available to account for
#' overdispersed counts, and this may be more widely applicable than using a
#' Poisson test. Under this approach, background sequences are divided into
#' sets of exactly the same size as the test set (i.e. iterations) and these
#' are used to estimate the parameters of the underlying model.
#'
#' When comparing a test set of sequences against a directly and specifically
#' matched background set of sequences, a one-sided hypergeometric test is also
#' implemented.
#'
#'
#' @return
#' A data.frame with columns: `sequences`, `matches`, `expected`, `enrichment`,
#' and `p`, with additional columns `Z`, `est_bg_rate` (Poisson), `odds_ratio`
#' (Hypergeometric) or `Z`, `sd_bg`, `n_iter` and `iter_p` (Iterations).
#' The numbers of sequences and matches refer to the test set of sequences,
#' whilst expected is the expected number of matches under the Poisson or
#' iterative null distribution. The ratio of matches to expected is given as
#' `enrichment`, along with the Z score and p-value. Whilst the Z score is only
#' representative under the Poisson model, it is used to directly estimate
#' p-values under the iterative approach.
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
#' @param pwm A Position Weight Matrix or list of PWMs
#' @param stringset An XStringSet with equal sequence widths
#' @param bg An XStringSet with the same sequence widths as the test XStringset
#' @param model The model used for analysis
#' @param var A column in the mcols element of bg, usually denoting an iteration
#' number
#' @param sort_by Column to sort results by
#' @param mc.cores Passed to \link[parallel]{mclapply}
#' @param ... Passed internally to \link[Biostrings]{countPWM}
#'
#' @examples
#' ## Load the example peaks & the sequences
#' data("ar_er_peaks")
#' data("ar_er_seq")
#' sq <- seqinfo(ar_er_peaks)
#' ## Now sample size-matched ranges 10 times larger. In real-world analyses,
#' ## this set should be sampled as at least 1000x larger, ensuring features
#' ## are matched to your requirements. This example masks regions with known N
#' ## content, including centromeres & telomeres
#' data("hg19_mask")
#' set.seed(305)
#' bg_ranges <- makeRMRanges(
#'   ar_er_peaks, GRanges(sq)[1], exclude = hg19_mask, n_iter = 10
#' )
#'
#' ## Convert ranges to DNAStringSets
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' bg_seq <- getSeq(genome, bg_ranges)
#'
#' ## Test for enrichment of the ESR1 motif
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#' testMotifEnrich(esr1, ar_er_seq, bg_seq, model = "poisson")
#'
#' ## Test all motifs
#' testMotifEnrich(ex_pwm, ar_er_seq, bg_seq, model = "poisson")
#'
#'
#' @importFrom stats p.adjust
#' @export
testMotifEnrich <- function(
        pwm, stringset, bg, var = "iteration",
        model = c("quasipoisson", "hypergeometric", "poisson", "iteration"),
        sort_by = c("p", "none"), mc.cores = 1, ...
) {

    ## Checks
    stopifnot(is(bg, "XStringSet"))
    model <- match.arg(model)
    args <- c(as.list(environment()), list(...))
    ## Prepare the output
    cols <- c("sequences", "matches", "expected", "enrichment", "Z", "p", "fdr")
    mod_cols <- list(
        poisson = "est_bg_rate", iteration = c("iter_p", "n_iter", "sd_bg"),
        hypergeometric = "odds_ratio", quasipoisson = c("n_iter", "sd_bg")
    )
    cols <- c(cols, mod_cols[[model]])
    if (model == "hypergeometric") cols <- setdiff(cols, "Z")
    ## Run the analysis
    if (is.matrix(pwm)) pwm <- list(pwm)
    pwm <- .cleanMotifList(pwm)
    if (model == "poisson") out <- .testPois(pwm, stringset, bg, mc.cores, ...)
    if (model == "iteration")
        out <- .testIter(pwm, stringset, bg, var, mc.cores, ...)
    if (model == "quasipoisson")
        out <- .testQuasi(pwm, stringset, bg, var, mc.cores, ...)
    if (model == "hypergeometric")
        out <- .testHyper(pwm, stringset, bg, mc.cores, ...)

    out$fdr <- p.adjust(out$p, "fdr")
    o <- seq_len(nrow(out))
    sort_by <- match.arg(sort_by)
    if (sort_by != "none") o <- order(out[[sort_by]])
    out[o,cols]

}

#' @importFrom parallel mclapply
#' @importFrom stats phyper
.testHyper <- function(pwm, stringset, bg, mc.cores, ...){

    ## Check there's no overlap between the two by removing shared sequences
    ss <- setdiff(stringset, bg)
    bg <- setdiff(bg, stringset)
    if (length(ss) != length(stringset)) {
        msg <- paste(
            "Shared sequences found in the test and background set.",
            "These have been excluded from both"
        )
        message(msg)
        stopifnot(length(ss) > 0)
    }
    n_ss <- length(ss)
    n_bg <- length(bg)

    ## Find & count the sequences with a match in both sets
    matches_ss <- getPwmMatches(pwm, stringset, mc.cores = mc.cores, ...)
    n_matches_ss <- mclapply(
        matches_ss, \(x) length(unique(x$seq)), mc.cores = mc.cores
    )
    n_matches_ss <- unlist(n_matches_ss)
    matches_bg <- getPwmMatches(pwm, bg, mc.cores = mc.cores, ...)
    n_matches_bg <- mclapply(
        matches_bg, \(x) length(unique(x$seq)), mc.cores = mc.cores
    )
    n_matches_bg <- unlist(n_matches_bg)
    or_denom <- (n_ss - n_matches_ss) / (n_bg - n_matches_bg)
    or <- (n_matches_ss / n_matches_bg) / or_denom

    ## Perform a 1-sided hypergeometric test
    p_vals <- mclapply(
        seq_along(pwm),
        \(i) {
            q <- n_matches_ss[[i]]
            k <- q + n_matches_bg[[i]]
            phyper(q - 1, m = n_ss, n = n_bg, k = k, lower.tail = FALSE)
        },
        mc.cores = mc.cores
    )

    ## Return the output
    data.frame(
        sequences = n_ss, matches = n_matches_ss,
        expected = n_ss * n_matches_bg / n_bg,
        enrichment = (n_matches_ss / n_ss) / (n_matches_bg / n_bg),
        p = unlist(p_vals),
        odds_ratio = or
    )

}

#' @importFrom parallel mclapply
#' @importFrom stats glm quasipoisson
#' @importFrom matrixStats colSds
#' @keywords internal
.testQuasi <- function(pwm, stringset, bg, var, mc.cores, ...) {

    n <- length(stringset)
    stopifnot(var %in% colnames(mcols(bg)))
    splitbg <- split(bg, mcols(bg)[[var]])
    if (!all(vapply(splitbg, length, integer(1)) == n))
        stop("All iterations must be the same size as the test sequences")

    matches <- countPwmMatches(pwm, stringset, mc.cores = mc.cores, ...)
    bg_matches <- mclapply(
        splitbg, \(x) countPwmMatches(pwm, x, mc.cores = 1, ...),
        mc.cores = mc.cores
    )
    bg_mat <- do.call("rbind", bg_matches)
    n_iter <- nrow(bg_mat)
    stopifnot(n_iter > 1)
    mean_bg <- colMeans(bg_mat)
    sd_bg <- colSds(bg_mat)
    Z <- (matches - mean_bg) / sd_bg

    p <- vapply(
        seq_along(pwm),
        \(i) {
            df <- data.frame(
                x = c(matches[[i]], bg_mat[,i]),
                type = c("test", rep_len("control", n_iter))
            )
            fit <- glm(x~type, family = quasipoisson(), data = df)
            summary(fit)$coef[2, 4]
        }, numeric(1)
    )

    data.frame(
        sequences = n, matches, expected = mean_bg,
        enrichment = matches / mean_bg, Z, p, n_iter, sd_bg
    )
}

#' @importFrom parallel mclapply
#' @importFrom stats pchisq
#' @importFrom matrixStats colSds
#' @keywords internal
.testIter <- function(pwm, stringset, bg, var, mc.cores, ...) {

    stopifnot(var %in% colnames(mcols(bg)))
    n <- length(stringset)
    matches <- countPwmMatches(pwm, stringset, mc.cores = mc.cores, ...)
    splitbg <- split(bg, mcols(bg)[[var]])
    if (!all(vapply(splitbg, length, integer(1)) == n))
        stop("All iterations must be the same size as the test sequences")

    bg_matches <- mclapply(
        splitbg, \(x) countPwmMatches(pwm, x, mc.cores = 1, ...),
        mc.cores = mc.cores
    )
    bg_mat <- do.call("rbind", bg_matches)
    mean_bg <- colMeans(bg_mat)
    sd_bg <- colSds(bg_mat)
    n_iter <- nrow(bg_mat)
    stopifnot(n_iter > 1)
    diff <- bg_mat - matrix(
        matches, nrow = n_iter, ncol = length(matches), byrow = TRUE
    )
    iter_p <- (colSums(diff > 0) + 1) / n_iter
    Z <- (matches - mean_bg) / sd_bg
    p <- 1 - pchisq(Z^2, 1)

    data.frame(
        sequences = length(stringset), matches, expected = mean_bg,
        enrichment = matches / mean_bg, Z, p, iter_p, n_iter, sd_bg
    )

}


#' @importFrom stats poisson.test
#' @keywords internal
.testPois <- function(pwm, test_seq, bg_seq, mc.cores, ...){

    n_seq <- length(test_seq)
    matches <- countPwmMatches(pwm, test_seq, mc.cores = mc.cores)
    n_bg <- countPwmMatches(pwm, bg_seq, mc.cores = mc.cores)
    est_bg_rate <- n_bg / length(bg_seq)
    expected <- est_bg_rate * n_seq
    ## Running vapply seems faster than mclappy here
    p <- vapply(
        seq_along(pwm),
        \(i) poisson.test(matches[i], n_seq, est_bg_rate[i])$p.value,
        numeric(1)
    )
    data.frame(
        sequences = n_seq, matches, expected, enrichment = matches / expected,
        Z  = (matches - expected) / sqrt(expected), p, est_bg_rate
    )

}

# @importFrom parallel mclapply
# @importFrom MASS glm.nb
# @importFrom matrixStats colSds
# @importFrom stats glm.control
# @keywords internal
#.testNB <- function(pwm, stringset, bg, var, mc.cores, ...) {
#
#     ## Currently disabled while I check in more detail
#     ## Convergence issues may be tricky to deal with...
#
#     stopifnot(var %in% colnames(mcols(bg)))
#     n <- length(stringset)
#     matches <- countPwmMatches(pwm, stringset, mc.cores = mc.cores, ...)
#     splitbg <- split(bg, mcols(bg)[[var]])
#     if (!all(vapply(splitbg, length, integer(1)) == n))
#         stop("All iterations must be the same size as the test sequences")
#
#     bg_matches <- mclapply(
#         splitbg, \(x) countPwmMatches(pwm, x, mc.cores = 1, ...),
#         mc.cores = mc.cores
#     )
#     bg_mat <- do.call("rbind", bg_matches)
#     n_iter <- nrow(bg_mat)
#     stopifnot(n_iter > 1)
#     mean_bg <- colMeans(bg_mat)
#     sd_bg <- colSds(bg_mat)
#     Z <- (matches - mean_bg) / sd_bg
#
#     p <- vapply(
#         seq_along(pwm),
#         \(i) {
#             df <- data.frame(
#                 x = c(matches[[i]], bg_mat[,i]),
#                 type = c("test", rep_len("control", n_iter))
#             )
#             ## Return NA if there is an error. This mostly derives from:
#             ## Error in while ((it <- it + 1) < limit && abs(del) > eps) { :
#             ## missing value where TRUE/FALSE needed
#             ## Currently ignoring warnings. May need to add a handler later
#             tryCatch(
#                 summary(
#                     ## Increasing maxit doesn't add much time overhead but
#                     ## helps with quite a few outliers
#                     glm.nb(x~type, data = df, control = glm.control(maxit = 100))
#                 )$coef[2,4],
#                 error = function(e){NA_real_}
#             )
#         }, numeric(1)
#     )
#
#     data.frame(
#         sequences = n, matches, expected = mean_bg,
#         enrichment = matches / mean_bg, Z, p, n_iter, sd_bg
#     )
# }
