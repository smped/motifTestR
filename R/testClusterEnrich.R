#' @title Test enrichment across a cluster of motifs using a background set of sequences
#'
#' @description
#' Test for enrichment of any motif within a cluster across a set of sequences
#' using a background set to
#' derive a NULL hypothesis
#'
#' @details
#' This extends the analytic methods offered by \link{testMotifEnrich} using
#' PWMs grouped into a set of clusters.
#' As with all cluster-level approaches, hits from multiple PWMs which overlap
#' are counted as a single hit ensuring that duplicated matches are not
#' double-counted, and that only individual positions within the sequences are.
#'
#'
#' @return
#' See \link{testMotifEnrich}
#'
#' @param cl  A list of Position Weight Matrices, universalmotifs, with each
#' element representing clusters of related matrices
#' @param stringset An XStringSet with equal sequence widths
#' @param bg An XStringSet with the same sequence widths as the test XStringset
#' @param model The model used for analysis
#' @param var A column in the mcols element of bg, usually denoting an iteration
#' number
#' @param sort_by Column to sort results by
#' @param mc.cores Passed to \link[parallel]{mclapply}
#' @param ... Passed to \link{getPwmMatches} or \link{countPwmMatches}
#'
#' @seealso [makeRMRanges()], [getClusterMatches()], [countClusterMatches()], [testMotifEnrich()]
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
#' ## Test for enrichment of clustered motifs
#' data("ex_pfm")
#' cl <- list(A = ex_pfm[1], B = ex_pfm[2:3])
#' testClusterEnrich(cl, ar_er_seq, bg_seq, model = "poisson")
#'
#'
#' @importFrom stats p.adjust
#' @export
testClusterEnrich <- function(
        cl, stringset, bg, var = "iteration",
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
    cl <- lapply(cl, .cleanMotifList)
    if (model == "poisson")
        out <- .testPois(cl, stringset, bg, mc.cores, type = "cluster", ...)
    if (model == "iteration")
        out <- .testIter(cl, stringset, bg, var, mc.cores, type = "cluster", ...)
    if (model == "quasipoisson")
        out <- .testQuasi(cl, stringset, bg, var, mc.cores, type = "cluster", ...)
    if (model == "hypergeometric")
        out <- .testHyper(cl, stringset, bg, mc.cores, type = "cluster", ...)

    out$fdr <- p.adjust(out$p, "fdr")
    o <- seq_len(nrow(out))
    sort_by <- match.arg(sort_by)
    if (sort_by != "none") o <- order(out[[sort_by]])
    out[o,cols]

}



