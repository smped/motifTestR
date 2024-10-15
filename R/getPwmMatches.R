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
#' @importFrom parallel mclapply
#' @importFrom methods slot
#' @importClassesFrom universalmotif universalmotif
#' @export
getPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "80%", best_only = FALSE,
        break_ties = c("all", "random", "first", "last", "central"),
        mc.cores = 1, ...
) {

    break_ties <- match.arg(break_ties)
    args <- c(as.list(environment()), list(...))
    args <- args[names(args) != "mc.cores"]
    nm_type <- "integer"
    if (!is.null(names(stringset))) nm_type <- "character"
    if (is.list(pwm)) {
        pwm <- .cleanMotifList(pwm)
        out <- mclapply(
            pwm, .getSinglePwmMatches, stringset = stringset, rc = rc,
            min_score = min_score, best_only = best_only,
            break_ties = break_ties, nm_type = nm_type, mc.cores = mc.cores
        )
    } else {
        args$nm_type <- nm_type
        out <- do.call(".getSinglePwmMatches", args)
    }
    out

}

#' @import Biostrings
#' @importFrom IRanges Views IRanges
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom methods as is
#' @importFrom stats setNames
#' @keywords internal
.getSinglePwmMatches <- function(
        pwm, stringset, rc, min_score, best_only = FALSE, break_ties, nm_type,
        ...
){

    ## Checks
    pwm <- .checkPWM(pwm)
    stopifnot(is(stringset, "XStringSet"))
    ## Handle empty stringsets
    if (!length(stringset)) return(.emptyPwmDF(nm_type))

    # Form the entire XStringSetList into a Views object
    map <- .viewMapFromXStringset(stringset)
    views <- Views(
        unlist(stringset), start = map$start, width = map$width,
        names = map$names
    )
    hits <- matchPWM(pwm, views, min.score = min_score, with.score = TRUE)
    mcols(hits)$direction <- rep_len("F", length(hits))
    if (rc) {
        rev_pwm <- reverseComplement(pwm)
        hits_rev <- matchPWM(
            rev_pwm, views, min.score = min_score, with.score = TRUE
        )
        mcols(hits_rev)$direction <- rep_len("R", length(hits_rev))
        hits <- c(hits, hits_rev)
        hits <- hits[!duplicated(IRanges(hits))] # Remove strict palindromic hits
    }

    ## Setup the return object as empty for any cases without matches
    if (length(hits) == 0) return(.emptyPwmDF(nm_type))

    ## Form it as a list, the eventually wrap into a DF
    cols <- c(
        "seq", "score", "direction", "start", "end", "from_centre",
        "seq_width", "match"
    )
    out <- mcols(hits)
    out$direction <- factor(out$direction, levels = c("F", "R"))

    ## Map back to the original Views
    hits_to_map <- findInterval(start(hits), map$start)
    w <- width(hits)

    ## Form the output
    out$seq <- hits_to_map
    out$start <- as.integer(start(hits) - c(0, map$end)[hits_to_map])
    out$end <- as.integer(out$start + w - 1)
    out$seq_width <- width(stringset[out$seq])
    out$from_centre <- (out$start + out$end - out$seq_width) / 2

    ## The match itself
    to_rev <- out$direction == "R"
    out$match <- as(hits, "XStringSet")
    out$match[to_rev] <- reverseComplement(out$match[to_rev])

    ## Setup any named strings to appear in the same order
    if (nm_type == "character") {
        out$seq <- factor(map$names[hits_to_map], map$names)
        out$seq <- droplevels(out$seq)
    }

    ## The final object
    if (best_only) out <- .getBestMatch(out, break_ties)
    o <- order(out$seq, out$start)
    if (nm_type == "character") out$seq <- as.character(out$seq)
    out[o, cols]

}

#' @importFrom S4Vectors splitAsList
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

#' @importFrom S4Vectors DataFrame
#' @keywords internal
.emptyPwmDF <- function(nm_type){
    DataFrame(
        seq = vector(nm_type, 0),
        score = numeric(),
        direction = factor(NULL, levels = c("F", "R")),
        start = integer(),
        end = integer(),
        from_centre = numeric(),
        seq_width = integer(),
        match = DNAStringSet()
    )
}
