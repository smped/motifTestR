#' @title Find matches from a PWM cluster within an XStringSet
#'
#' @description
#' Find matches from a PWM cluster within a set of sequences
#'
#' @details
#' This function extends \link{getPwmMatches} by returning a single set of
#' results for set of clustered motifs.
#' This can help remove some of the redundancy in results returned for highly
#' similar PWMs, such as those in the GATA3 family.
#'
#' Taking a set of sequences as an XStringSet, find all matches above the
#' supplied score (i.e. threshold) for a list of Position Weight Matrices
#' (PWMs), which have been clustered together as highly-related motifs.
#' By default, matches are performed using the PWMs as provided and the reverse
#' complement, however this can easily be disabled by setting `rc = FALSE`.
#'
#' The function relies heavily on \link[Biostrings]{matchPWM} and
#' \link[IRanges]{Views} for speed.
#'
#' Where overlapping matches are found for the PWMs within a cluster, only a
#' single match is returned.
#' The motif with the highest relative score (score / maxScore(PWM)) is selected.
#'
#' When choosing to return the best match (`best_only = TRUE`), only the match
#' with the highest relative score is returned for each sequence.
#' Should there be tied scores, the best match can be chosen as either the first,
#' last, most central, all tied matches, or choosing one at random (the default).
#'
#' @param cl A list of Position Weight Matrices, universalmotifs, with each
#' element representing clusters of related matrices
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param best_only logical(1) Only return the best match
#' @param break_ties Method for breaking ties when only returning the best match
#' Ignored when all matches are returned (the default)
#' @param mc.cores Passed to \link[parallel]{mclapply} if passing multiple PWMs
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#' @return A list of DataFrames with columns: `seq`, `score`, `direction`,
#' `start`, `end`, `from_centre`, `seq_width`, `motif` and `match`
#'
#' The first three columns describe the sequence with matches, the score of
#' the match and whether the match was found using the forward or reverse PWM.
#' The columns `start`, `end` and `width` describe the where the match was found
#' in the sequence, whilst `from_centre` defines the distance between the centre
#' of the match and the centre of the sequence being queried.
#' The motif column denotes which individual motif was found to match in this
#' position, again noting that when matches overlap, only the one with the
#' highest relative score is returned.
#' The final column contains the matching fragment of the sequence as an
#' `XStringSet`.
#'
#'
#' @importFrom parallel mclapply
#' @importFrom methods slot
#' @importClassesFrom universalmotif universalmotif
#' @export
getClusterMatches <- function(
        cl, stringset, rc = TRUE, min_score = "80%", best_only = FALSE,
        break_ties = c("all", "random", "first", "last", "central"),
        mc.cores = 1, ...
) {

    stopifnot(is.list(cl))
    break_ties <- match.arg(break_ties)
    args <- c(as.list(environment()), list(...))
    args <- args[names(args) != "mc.cores"]
    nm_type <- "integer"
    if (!is.null(names(stringset))) nm_type <- "character"
    out <- NULL
    if (
        all(vapply(cl, is, logical(1), "list")) |
        all(vapply(cl, is, logical(1), "universalmotif"))
    ) {
        ret_cols <- c(
            "seq", "score", "direction", "start", "end", "from_centre",
            "seq_width", "motif", "match"
        )
        cl <- lapply(cl, .cleanMotifList)
        nm <- names(cl)
        if (is.null(nm)) nm <- seq_along(cl)
        single_pwm <- vapply(cl, length, integer(1)) == 1
        out <- vector("list", length(cl))
        out[single_pwm] <- mclapply(
            cl[single_pwm],
            \(x) {
                DF <- .getSinglePwmMatches(
                    x[[1]], stringset = stringset, rc = rc,
                    min_score = min_score, best_only = best_only,
                    break_ties = break_ties, nm_type = nm_type,
                )
                DF$motif <- names(x)
                DF[ret_cols]
            },  mc.cores = mc.cores
        )
        out[!single_pwm] <- mclapply(
            cl[!single_pwm], .getClusterPwmMatches, stringset = stringset,
            rc = rc, min_score = min_score, best_only = best_only,
            break_ties = break_ties, nm_type = nm_type, mc.cores = mc.cores
        )
    }
    if (all(vapply(cl, is, logical(1), "matrix"))) {
        args$nm_type <- nm_type
        out <- do.call(".getClusterPwmMatches", args)
    }

    if (is.null(out)) message("Could not determine clusters")
    out

}



#' @import Biostrings
#' @importFrom IRanges Views IRangesList DataFrameList findOverlaps
#' @importFrom S4Vectors mcols<- queryHits subjectHits
#' @keywords internal
.getClusterPwmMatches <- function(
        cl, stringset, rc, min_score, best_only = FALSE, break_ties, nm_type,
        counts_only = FALSE, ...
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
    if (counts_only) return(length(merged_hits))

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

