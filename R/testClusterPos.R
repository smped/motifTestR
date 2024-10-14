#' Test positional bias motifs within a cluster
#'
#' @description
#' Test positional bias for all motifs within a given cluster
#'
#' @details
#' This is a reimplmentation of \link{testMotifPos} for sets of motifs which
#' have been clustered for similarity.
#' The positions test the bias of any motifs within the cluster given that
#' overlapping matches are only counted once, and with the match retained being
#' the one with the highest relative score.
#'
#' It should also be noted that some motif clusters will contain PWMs of
#' varying length. When finding positional bias, the widest motif is taken as
#' the width for all, and any matches from narrower motifs outside of the range
#' allowed by wider motifs are discarded. This reduction in signal will make a
#' small difference in the outer bins, but is not considered to be problematic
#' for the larger analysis.
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
#' @param x A Position Weight Matrix, universalmotif object or list thereof.
#' Alternatively can be a single DataFrame or list of DataFrames as returned
#' by \link{getClusterMatches} with `best_only = TRUE`
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
#' @examples
#' ## Load the example PWM
#' data("ex_pwm")
#' ## Load the example sequences
#' data("ar_er_seq")
#'
#' ## Cluster the motifs
#' cl <- list(A = ex_pwm[1], B = ex_pwm[2:3])
#'
#' ## Get the best match and use this data
#' matches <- getClusterMatches(cl, ar_er_seq, best_only = TRUE)
#' ## Test for enrichment in any position
#' testClusterPos(matches)
#'
#' ## Or just pass the clustered matrices
#' ## Here we've set abs = TRUE to test absolute distance from the centre
#' testClusterPos(cl, ar_er_seq, abs = TRUE, binwidth = 10)
#'
#' @importFrom parallel mclapply
#' @importFrom stats p.adjust
#' @export
testClusterPos <- function(
        x, stringset, binwidth = 10, abs = FALSE, rc = TRUE, min_score = "80%",
        break_ties = "all", alt = c("greater", "less", "two.sided"),
        sort_by = c("p", "none"), mc.cores = 1,
        ...
) {

    cols <- c(
        "start", "end", "centre", "width", "total_matches", "matches_in_region",
        "expected", "enrichment", "prop_total", "p", "fdr", "consensus_motif"
    )
    alt <- match.arg(alt)
    sort_by <- match.arg(sort_by)

    if (is(x, "DataFrame")) x <- list(x)
    is_DF <- vapply(x, is, logical(1), "DataFrame")
    matches <- NULL
    if (all(is_DF)) {
        ## If we have a list of best matches
        check <- .checkMatches(x)
        matches <- x
    } else {
        ## Otherwise, obtain the list of best match DFs
        is_list <- vapply(x, is, logical(1), "list")
        if (!all(is_list)) stop("All clusters should contain lists of motifs")
        clust_lens <- vapply(x, length, integer(1))
        if (all(clust_lens == 1)) stop("All clusters cannot contain a single motif")
        cl <- lapply(x, .cleanMotifList)
        matches <- getClusterMatches(
            cl, stringset, rc, min_score, TRUE, break_ties, mc.cores, ...
        )
    }
    if (is.null(matches)) stop("Couldn't determine structure of input")

    out <- mclapply(
        matches, .testSingleClusterPos, binwidth = binwidth, abs = abs, rc = rc,
        min_score = min_score, break_ties = break_ties, alt = alt, ...,
        mc.cores = mc.cores
    )
    out <- do.call("rbind", out)
    out$fdr <- p.adjust(out$p, "fdr")
    o <- seq_len(nrow(out))
    if (sort_by != "none") o <- order(out[[sort_by]])
    out[o, cols]

}

#' @keywords internal
.testSingleClusterPos <- function(
        matches, binwidth, abs, rc, min_score, break_ties, alt, ...
) {

    ## This differs from the single motif case in that different width motifs
    ## may exist within a cluster. As a result, the boundary points need to be
    ## recalculated using the widest motif and all matches outside these
    ## boundary points will be excluded. This small loss of information
    ## shouldn't be profoundly deleterious to the process

    ## Setup a well formed output object for easy combining with other tests
    out_cols <- c(
        "start", "end", "centre", "width", "total_matches", "matches_in_region",
        "expected", "enrichment", "prop_total", "p", "consensus_motif"
    )
    out <- lapply(out_cols, \(x) integer())
    out$consensus_motif <- list()
    names(out) <- out_cols
    if (!nrow(matches)) return(data.frame(out))

    ## Remove those where the centre is outside of the range
    max_width <- max(width(matches$match))
    seq_width <- min(matches$seq_width)
    lim <- (seq_width - max_width) / 2
    matches <- matches[abs(matches$from_centre) <= lim, ]

    .testSingleMotifPos(
        matches, binwidth, abs, rc, min_score, break_ties, alt, ...
    )
}
