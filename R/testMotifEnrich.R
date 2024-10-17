#' @title Test motif enrichment using a background set of sequences
#'
#' @description
#' Test for motif enrichment within a set of sequences using a background set to
#' derive a NULL hypothesis
#'
#' @details
#' This function offers four alternative models for assessing the enrichment of
#' a motif within a set of sequences, in comparison to a background set of
#' sequences.
#' Selection of the BG sequences plays an important role and, in conjunction
#' with the question being addressed, determines the most appropriate model to
#' use for testing, as described below.
#'
#' It should also be noted that the larger the BG set of sequences, the larger
#' the computational burden, and results can take far longer to return. For
#' many millions of background sequences, this may run beyond an hour
#'
#' ## Descriptions of Models and Use Cases
#'
#' ### Hypergeometric Tests
#'
#' Hypergeometric tests are best suited to the use case where the test set of
#' sequences represents a subset of a larger set, with a specific feature or
#' behaviour, whilst the BG set may be the remainder of the set without that
#' feature. For example, the test set may represent ChIP-Seq binding sites
#' where signal changes in response to treatment, whilst the BG set represents
#' the sites where no changed signal was observed. Testing is one-sided, for
#' enrichment of motifs within the test set.
#'
#' Due to these relatively smaller sized datasets, setting
#' model = "hypergeometric", will generally return results quickly
#'
#' ### Poisson Tests
#'
#' This approach requires a set of background sequences which should be much
#' larger than the test set of sequences.
#' The parameters for a Poisson model are estimated in a per-sequence manner on
#' the set of BG sequences, and the observed rate of motif-matches within the
#' test set is then tested using \link[stats]{poisson.test}.
#'
#' This approach assumes that all matches follow a Poisson distribution, which
#' is often true, but data can also be over-dispersed. Given that this model can
#' also return results relatively quickly, is it primarily suitable for data
#' exploration, such as quickly checking for expected behaviours, but not for
#' final results.
#'
#' ### Quasi-Poisson Test
#'
#' The quasipoisson model allows for over-dispersion and will return more
#' conservative results than using the standard Poisson model.
#' Under the method currently implemented here, BG sequences should be divided
#' into blocks (i.e. iterations), identical in size to the test set of sequences.
#' Model parameters are estimated per iteration across the BG set of sequences,
#' with the rate of matches in the test set being compared against these blocks.
#' This ensures more conservative results that if analysing test and bg
#' sequences as collections of individual sequences.
#'
#' It is expected that the BG set will matched for the features of interest and
#' chosen using \link{makeRMRanges} with a large number of iterations, e.g.
#' n_iter = 1000.
#' Due to this parameterisation, quasipoisson approaches can be computationally
#' time-consuming, as this is effectively an iterative approach.
#'
#' ### Iteration
#'
#' Setting the model as "iteration" performs a non-parametric analysis, with
#' the exception of returning Z-scores under the Central Limit Theorem.
#' Mean and SD of matches is found for each iteration, and used to return Z
#' scores, with p-values returned from both a Z-test and from comparing
#' observed values directly to sampled values obtained from the BG sequences.
#' Sampled values are calculated directly and as such, are limited in precision.
#'
#' As for the QuasiPoisson model, a very large number of iterations is expected
#' to be used, to ensure the CLT holds, again making this a computationally
#' demanding test.
#' Each iteration/block is expected to be identically-sized to the test set,
#' and matched for any features as appropriate using [makeRMRanges()].
#'
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
#' @param ... Passed to \link{getPwmMatches} or \link{countPwmMatches}
#'
#' @seealso [makeRMRanges()], [getPwmMatches()], [countPwmMatches()]
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
#' data("ex_pfm")
#' esr1 <- ex_pfm$ESR1
#' testMotifEnrich(esr1, ar_er_seq, bg_seq, model = "poisson")
#'
#' ## Test all motifs
#' testMotifEnrich(ex_pfm, ar_er_seq, bg_seq, model = "poisson")
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
.testHyper <- function(
        x, stringset, bg, mc.cores, type = c("pwm", "cluster"), ...
){

    ## Check there's no overlap between the two by removing shared sequences
    cl <- class(stringset)
    ss <- as(setdiff(stringset, bg), cl)
    bg <- as(setdiff(bg, stringset), cl)
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

    type <- match.arg(type)

    if (type == "pwm") {
        matches_ss <- mclapply(
            x, .hasPwmMatch, stringset = ss, ..., mc.cores = mc.cores
        )
        matches_bg <- mclapply(
            x, .hasPwmMatch, stringset = bg, ..., mc.cores = mc.cores
        )
    }
    if (type == "cluster") {
        matches_ss <- mclapply(
            x, .hasClusterMatch, stringset = ss, ..., mc.cores = mc.cores
        )
        matches_bg <- mclapply(
            x, .hasClusterMatch, stringset = bg, ..., mc.cores = mc.cores
        )
    }
    n_matches_ss <- vapply(matches_ss, sum, integer(1))
    n_matches_bg <- vapply(matches_bg, sum, integer(1))

    or_denom <- (n_ss - n_matches_ss) / (n_bg - n_matches_bg)
    or <- (n_matches_ss / n_matches_bg) / or_denom

    ## Perform a 1-sided hypergeometric test
    p_vals <- mclapply(
        seq_along(x),
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
.testQuasi <- function(
        x, stringset, bg, var, mc.cores, type = c("pwm", "cluster"), ...
) {

    n <- length(stringset)
    stopifnot(var %in% colnames(mcols(bg)))
    splitbg <- split(bg, mcols(bg)[[var]])
    if (!all(vapply(splitbg, length, integer(1)) == n))
        stop("All iterations must be the same size as the test sequences")

    type <- match.arg(type)
    if (type == "pwm") {
        matches <- countPwmMatches(x, stringset, mc.cores = mc.cores, ...)
        bg_matches <- mclapply(
            splitbg, \(i) countPwmMatches(x, i, mc.cores = 1, ...),
            mc.cores = mc.cores
        )
    }
    if (type == "cluster") {
        matches <- countClusterMatches(x, stringset, mc.cores = mc.cores, ...)
        bg_matches <- mclapply(
            splitbg, \(i) countClusterMatches(x, i, mc.cores = 1, ...),
            mc.cores = mc.cores
        )
    }
    bg_mat <- do.call("rbind", bg_matches)
    n_iter <- nrow(bg_mat)
    stopifnot(n_iter > 1)
    mean_bg <- colMeans(bg_mat)
    sd_bg <- colSds(bg_mat)
    Z <- (matches - mean_bg) / sd_bg

    p <- vapply(
        seq_along(x),
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
.testIter <- function(
        x, stringset, bg, var, mc.cores, type = c("pwm", "cluster"), ...
) {

    stopifnot(var %in% colnames(mcols(bg)))
    n <- length(stringset)
    splitbg <- split(bg, mcols(bg)[[var]])
    if (!all(vapply(splitbg, length, integer(1)) == n))
        stop("All iterations must be the same size as the test sequences")
    type <- match.arg(type)
    if (type == "pwm") {
        matches <- countPwmMatches(x, stringset, mc.cores = mc.cores, ...)
        bg_matches <- mclapply(
            splitbg, \(i) countPwmMatches(x, i, mc.cores = 1, ...),
            mc.cores = mc.cores
        )
    }
    if (type == "cluster") {
        matches <- countClusterMatches(x, stringset, mc.cores = mc.cores, ...)
        bg_matches <- mclapply(
            splitbg, \(i) countClusterMatches(x, i, mc.cores = 1, ...),
            mc.cores = mc.cores
        )
    }
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
.testPois <- function(
        x, test_seq, bg_seq, mc.cores, type = c("pwm", "cluster"), ...
){

    n_seq <- length(test_seq)
    type <- match.arg(type)
    if (type == "pwm") {
        matches <- countPwmMatches(x, test_seq, mc.cores = mc.cores, ...)
        n_bg <- countPwmMatches(x, bg_seq, mc.cores = mc.cores, ...)
    }
    if (type == "cluster") {
        matches <- countClusterMatches(x, test_seq, mc.cores = mc.cores, ...)
        n_bg <- countClusterMatches(x, bg_seq, mc.cores = mc.cores, ...)
    }
    est_bg_rate <- n_bg / length(bg_seq)
    expected <- est_bg_rate * n_seq
    ## Running vapply seems faster than mclappy here
    p <- vapply(
        seq_along(x),
        \(i) poisson.test(matches[i], n_seq, est_bg_rate[i])$p.value,
        numeric(1)
    )
    data.frame(
        sequences = n_seq, matches, expected, enrichment = matches / expected,
        Z  = (matches - expected) / sqrt(expected), p, est_bg_rate
    )

}

