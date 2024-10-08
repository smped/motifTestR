#' @title Test Enrichment of motifs within a cluster
#'
#' @description Test enrichment of related motifs assigned to a cluster
#'
#' @details
#' This function uses a similar strategy to \link{testMotifEnrich}, however
#' motif detection is applied to all motifs within a cluster instead of a
#' single motif
#'

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
    ## Run the analysis
    is_mat <- vapply(cl, is.matrix, logical(1))
    cl[is_mat] <- lapply(cl[is_mat], \(x) list(x))
    cl <- lapply(cl, .cleanMotifList)
    ## Set these up one at a time ASAP
    # if (model == "poisson") out <- .testPois(pwm, stringset, bg, mc.cores, ...)
    # if (model == "iteration")
    #     out <- .testIter(pwm, stringset, bg, var, mc.cores, ...)
    # if (model == "quasipoisson")
    #     out <- .testQuasi(pwm, stringset, bg, var, mc.cores, ...)
    # if (model == "hypergeometric") {
    #     cols <- setdiff(cols, "Z")
    #     out <- .testHyper(pwm, stringset, bg, mc.cores, ...)
    # }
    #
    # out$fdr <- p.adjust(out$p, "fdr")
    # o <- seq_len(nrow(out))
    # sort_by <- match.arg(sort_by)
    # if (sort_by != "none") o <- order(out[[sort_by]])
    # out[o,cols]

}


#' @import Biostrings
#' @importFrom IRanges Views IRangesList DataFrameList findOverlaps
#' @importFrom S4Vectors mcols<- queryHits subjectHits
#' @keywords internal
.getClusterPwmMatches <- function(
        cl, stringset, rc, min_score, best_only = FALSE, break_ties, nm_type,
        ...
){

    ## Checks
    cl <- lapply(cl, .checkPWM)
    stopifnot(is(stringset, "XStringSet"))
    ## Handle empty stringsets
    empty_df <- .emptyPwmDF(nm_type)
    empty_df$motif <- character()
    if (!length(stringset)) return(empty_df)

    # Form the entire XStringSetList into a Views object & find all hits
    map <- .viewMapFromXStringset(stringset)
    views <- Views(
        unlist(stringset), map$start, width = map$width, names = map$names
    )
    hits <- lapply(
        cl, matchPWM, subject = views, min.score = min_score, with.score = TRUE,
        ...
    )
    .add_dir <- function(x, dir) {
        mcols(x)$direction <- rep_len(dir, length(x))
        x
    }
    hits <- lapply(hits, .add_dir, dir = "F")
    hits_rev <- list()
    if (rc) {
        rev_cl <- lapply(cl, reverseComplement)
        hits_rev <- lapply(
            rev_cl, matchPWM, subject = views, min.score = min_score,
            with.score = TRUE, ...
        )
        hits_rev <- lapply(hits_rev, .add_dir, dir = "R")
    }

    ## Combine Hits
    ir <- lapply(c(hits, hits_rev), slot, "ranges")
    all_hits <- unlist(IRangesList(ir))
    if (length(all_hits) == 0) return(empty_df)

    ## Add mcols if we have hits
    mcols(all_hits) <- unlist(DataFrameList(lapply(c(hits, hits_rev), mcols)))
    mcols(all_hits)$motif <- names(all_hits)
    all_hits <- unname(all_hits)
    ## Sort by decreasing relative score to enable selecting the highest scoring
    max_scores <- vapply(cl, maxScore, numeric(1))
    prop <- mcols(all_hits)$score / max_scores[mcols(all_hits)$motif]
    mcols(all_hits)$prop <- prop
    all_hits <- all_hits[order(prop, decreasing = TRUE)]
    merged_hits <- IRanges::reduce(all_hits) # Will be in positional order
    ## For overlapping hits, choose the one with highest relative score (first)
    ol <- findOverlaps(all_hits, merged_hits)
    best <- queryHits(ol)[!duplicated(subjectHits(ol))]
    final_hits <- sort(all_hits[best])

    ## Map back to the original Views
    hits_to_map <- findInterval(start(final_hits), map$start)
    w <- width(final_hits)

    ## Form the output
    cols <- c(
        "seq", "score", "direction", "start", "end", "from_centre",
        "seq_width", "motif", "match"
    )
    out <- mcols(final_hits)
    out$direction <- factor(out$direction, levels = c("F", "R"))
    out$seq <- hits_to_map
    out$start <- as.integer(start(final_hits) - c(0, map$end)[hits_to_map])
    out$end <- as.integer(out$start + w - 1)
    out$seq_width <- width(stringset[out$seq])
    out$from_centre <- (out$start + out$end - out$seq_width) / 2

    ## The match itself
    to_rev <- out$direction == "R"
    out$match <- as(
        Views(unlist(stringset), start = out$start, end = out$end), "XStringSet"
    )
    out$match[to_rev] <- reverseComplement(out$match[to_rev])

    ## Setup any named strings to appear in the same order
    if (nm_type == "character") {
        out$seq <- factor(map$names[hits_to_map], map$names)
        out$seq <- droplevels(out$seq)
    }

    ## The final object
    if (best_only) out <- .getBestMatch(out, break_ties, score_col = "prop")
    o <- order(out$seq, out$start)
    if (nm_type == "character") out$seq <- as.character(out$seq)
    out[o, cols]

}

