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
#' When choosing to return the best match (`best_only = TRUE`), only the match
#' with the highest score is returned for each sequence.
#' Should there be tied scores, the best match can be chosen as either the first,
#' last, most central, all tied matches, or choosing one at random (the default).
#'
#' Importantly, prior to Bioconductor Release 3.25, this relied on
#' \link[Biostrings]{matchPWM} but now utilises
#' \link[universalmotif]{scan_sequences} for greater flexibility.
#' There may be a slight decrease in speed as a result, especially when
#' incorporating pvalues or qvalues for thresholding matches.
#'
#' @param pwm A Position Weight Matrix, list of PWMs or universalmotif list
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match if matching by logodds
#' or logadds.abs. Can be passed as a percentage or as a value between 0 & 1
#' @param pvalue,qvalue thresholds for returning amatch when using either as the
#' threshold type
#' @param threshold_type See \link[universalmotif]{scan_sequences}. Defaults to
#' "logodds" for approximate compatiblity with earlier versions
#' @param best_only logical(1) Only return the best match
#' @param break_ties Method for breaking ties when only returning the best match
#' Ignored when all matches are returned (the default)
#' @param mc.cores Passed as `nthreads` to  \link[universalmotif]{scan_sequences}
#' @param calc.pvals See \link[universalmotif]{scan_sequences}
#' @param ... Passed to See \link[universalmotif]{scan_sequences}
#'
#' @return A DataFrame with columns: `seq`, `score`, `direction`, `start`,
#' `end`, `from_centre`, `seq_width`, and `match`.
#' Additional columns contain `pvalue` and `qvalue` will also be included if
#' requested via `calc.pvals`, or if using either of these methods.
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
#' @importFrom methods as
#' @importClassesFrom universalmotif universalmotif
#' @export
getPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "50%", pvalue = 1e-4,
        qvalue = 0.05, threshold_type = "logodds", best_only = FALSE,
        break_ties = c("all", "random", "first", "last", "central"),
        mc.cores = 1, calc.pvals = NULL, ...
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
    }
    thresh <- pvalue
    if (grepl("logodds", threshold_type)) {
        thresh <- min_score
        if (is.null(calc.pvals)) calc.pvals <- FALSE
    }
    if (is.null(calc.pvals)) calc.pvals <- TRUE
    if (threshold_type == "qvalue") thresh <- qvalue
    stopifnot(thresh > 0 & thresh < 1) # Would lead to silly results

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
        RC = rc, nthreads = mc.cores, calc.pvals = calc.pvals, ...
    )
    matches$motif <- factor(matches$motif, levels = to_df(pwm)$name)
    if (!is.null(min_score)) matches <- matches[matches$score.pct > min_score*100,]
    if (is.null(names(stringset))) matches$sequence <- as.integer(matches$sequence)
    out <- list(
        seq = matches$sequence, score = matches$score,
        direction = factor(
            c("+" = "F", "-" = "R")[matches$strand], levels = c("F", "R")
        ),
        start = pmin(matches$start, matches$stop),
        end = pmax(matches$start, matches$stop),
        from_centre = numeric(),
        seq_width = width(stringset)[matches$sequence.i],
        match =  as(matches$match, is(stringset)[[1]])
    )
    out$from_centre <- (out$start + out$end - out$seq_width) / 2
    if (calc.pvals) {
        out$pvalue <- matches$pvalue
        out$qvalue <- matches$qvalue
    }
    o <- order(out$seq)
    out <- split(as(out, "DataFrame")[o,], matches$motif[o])

    if (best_only) out <- lapply(out, .getBestMatch, ties = break_ties)
    if (pwm_was_matrix) out <- out[[1]]
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

