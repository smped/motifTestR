#' @title Find all PWM matches within an XStringSet
#'
#' @description
#' Find all PWM matches within a set of sequences
#'
#' @details
#' Taking a set of sequences as an XStringSet, find all matches above the
#' supplied score (i.e. threshold) for a single Position Weight Matrix (PWM),
#' generally representing a transcription factor binding motif.
#' By default, matches are performed using the PWM as provided and the reverse
#' complement, however this can easily be disabled by setting `rc = FALSE`.
#'
#' The function relies heavily on \link[Biostrings]{matchPWM} and
#' \link[IRanges]{Views} for speed.
#'
#' When choosing to return the best match (`best_only = TRUE`), only the match
#' with the highest score is returned for each sequence.
#' Should there be tied scores, the best match can be chosen as either the first,
#' last, most central, all tied matches, or choosing one at random (the default).
#'
#' @param pwm A Position Weight Matrix, list of PWMs or universalmotif list
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param best_only logical(1) Only return the best match
#' @param break_ties Method for breaking ties when only returning the best match
#' Ignored when all matches are returned (the default)
#' @param mc.cores Passed to \link[parallel]{mclapply} if passing multiple PWMs
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#' @return A DataFrame with columns: `seq`, `score`, `direction`, `start`,
#' `end`, `from_centre`, `seq_width`, and `match`
#'
#' The first three columns describe the sequence with matches, the score of
#' the match and whether the match was found using the forward or reverse PWM.
#' The columns `start`, `end` and `width` describe the where the match was found
#' in the sequence, whilst `from_centre` defines the distance between the centre
#' of the match and the centre of the sequence being queried.
#' The final column contains the matching fragment of the sequence as an
#' `XStringSet`.
#'
#' When passing a list of PWMs, a list of the above DataFrames will be returned.
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pfm")
#' esr1 <- ex_pfm$ESR1
#'
#' ## Load the example Peaks
#' data("ar_er_seq")
#'
#' ## Return all matches
#' getPwmMatches(esr1, ar_er_seq)
#'
#' ## Just the best match
#' getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
#'
#' ## Apply multiple PWMs as a list
#' getPwmMatches(ex_pfm, ar_er_seq, best_only = TRUE)
#'
#' @importFrom universalmotif scan_sequences convert_motifs convert_type to_list
#' @importFrom universalmotif to_list to_df
#' @importFrom Biostrings width
#' @importClassesFrom universalmotif universalmotif
#' @export
getPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "50%", pvalue = 1e-4,
        qvalue = 0.05, threshold_type = "pvalue", best_only = FALSE,
        break_ties = c("all", "random", "first", "last", "central"),
        mc.cores = 1, ...
){

    break_ties <- match.arg(break_ties)
    threshold_type <- match.arg(
        threshold_type, eval(formals(scan_sequences)$threshold.type)
    )
    pwm_was_matrix <- is.matrix(pwm)

    ## Check the min_score
    if (!is.null(min_score)){
        if (is.character(min_score)) {
            stopifnot(grepl("%$", min_score))
            min_score <- as.numeric(gsub("%$", "", min_score)) / 100
        }
        if (is.numeric(min_score)) stopifnot(min_score > 0 & min_score < 1)
    }
    stopifnot(pvalue > 0 & pvalue < 1) # Cannot be NULL
    thresh <- pvalue
    if (grepl("logodds", threshold_type)) thresh <- min_score
    if (threshold_type == "qvalue") {
        stopifnot(qvalue > 0 & qvalue < 1)
        thresh <- qvalue
    }
    stopifnot(thresh > 0) # Would lead to silly results

    ## Now run the checks on the PWM
    if (is.list(pwm)) {
        if (all(vapply(pwm, is.matrix, logical(1)))) {
            um <- lapply(pwm, convert_motifs)
            um_df <- to_df(um)
            stopifnot(length(unique(names(pwm))) == length(pwm)) # Uniq names
            um_df$name <- names(pwm)
            um_df$nsites <- 100 # Avoids an annoying message in the next step
            pwm <- to_list(um_df)
        }
        stopifnot(all(vapply(pwm, is, logical(1), class2 = "universalmotif")))
        pwm <- convert_type(pwm, type = "PWM")
    }

    matches <- scan_sequences(
        pwm, stringset, threshold = thresh, threshold.type = threshold_type,
        RC = rc, nthreads = mc.cores, ...
    )
    if (!is.null(min_score)) matches <- matches[matches$score.pct > min_score*100,]
    names(matches)[c(3, 6)] <- c("seq", "end")
    matches <- matches[order(matches$motif, matches$sequence.i, matches$start),]
    matches$direction <- c("+" = "F", "-" = "R")[matches$strand]
    matches$direction <- factor(matches$direction, levels = c("F", "R"))
    matches$seq_width <- width(stringset)[matches$sequence.i]
    matches$from_centre <- (matches$start + matches$end - matches$seq_width) / 2
    matches$match <- as(matches$match, is(stringset)[[1]])
    keep_cols <- c(
        "seq", "score", "direction", "start", "end", "from_centre", "seq_width",
        "match", "pvalue", "qvalue"
    )
    out <- matches[, keep_cols]
    matches$motif <- factor(matches$motif, levels = to_df(pwm)$name)
    out <- split(out, matches$motif)

    if (best_only) out <- lapply(out, .getBestMatch, ties = break_ties)
    if (pwm_was_matrix) out <- DataFrame(out[[1]])
    out
}

#' @keywords internal
.getBestMatch <- function(df, ties, score_col = "score"){

    if (nrow(df) == 0) return(df)

    ## Decide how to break ties
    pos <- sample.int(nrow(df), nrow(df)) # Default to random
    if (ties == "first") pos <- df$start
    if (ties == "last") pos <- (-1) * df$end
    if (ties == "central") pos <- abs(df$from_centre)

    ## Pick the best match using the values we have
    if (ties != "all") {
        o <- order(df$seq, -df[[score_col]], pos)
        df <- df[o,]
        df <- df[!duplicated(df$seq), ]
    } else {
        split_scores <- split(df[[score_col]], cumsum(!duplicated(df$seq)))
        keep <- as.logical(unlist(lapply(split_scores, \(x) x == max(x))))
        df <- df[keep,]
    }
    df

}

