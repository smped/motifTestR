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
#' complement, however this can easily be disable by setting `rc = FALSE`.
#'
#' The function relies heavily on \link[Biostrings]{matchPWM} and
#' \link[IRanges]{Views} for speed.
#'
#' When choosing to return the best match (`best_only = TRUE`), only the match
#' with the highest score is returned for each sequence.
#' Should there be tied scores, the best match can be chosen as either the first,
#' last, most central, all tied matches, or choosing one at random (the default).
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param best_only logical(1) Only return the best match
#' @param break_ties Method for breaking ties when only returning the best match
#' Ignored when all matches are returned (the default)
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#' @return A DataFrame with columns: `seq`, `score`, `direction`, `start`,
#' `end`, `fromCentre`, `seq_width`, and `match`
#'
#' The first three columns describe the sequence with matches, the score of
#' the match and whether the match was found using the forward or reverse PWM.
#' The columns `start`, `end` and `width` describe the where the match was found
#' in the sequence, whilst `from_centre` defines the distance between the centre
#' of the match and the centre of the sequence being queried.
#' The final column contains the matching fragment of the sequence as an
#' `XStringSet`.
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#'
#' ## Load the example Peaks
#' data("ar_er_peaks")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' seq <- getSeq(genome, ar_er_peaks)
#'
#' ## Return all matches
#' getPwmMatches(esr1, seq)
#'
#' ## Just the best match
#' getPwmMatches(esr1, seq, best_only = TRUE)
#'
#'
#' @import Biostrings
#' @importFrom IRanges Views
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom methods as is
#' @importFrom stats setNames
#' @export
getPwmMatches <- function(
    pwm, stringset, rc = TRUE, min_score = "80%", best_only = FALSE,
    break_ties = c("random", "first", "last", "central", "all"), ...
){

  ## Checks & the map
  .checkPWM(pwm)
  stopifnot(is(stringset, "XStringSet"))
  map <- .viewMapFromXStringset(stringset)

  # Form the entire XStringSetList into a Views object
  views <- Views(
    unlist(stringset), start = map$start, width = map$width, names = map$names
  )
  hits <- matchPWM(pwm, views, min.score = min_score, with.score = TRUE, ...)
  mcols(hits)$direction <- rep_len("F", length(hits))
  if (rc) {
    rev_pwm <- reverseComplement(pwm)
    hits_rev <- matchPWM(
      rev_pwm, views, min.score = min_score, with.score = TRUE, ...
    )
    mcols(hits_rev)$direction <- rep_len("R", length(hits_rev))
    hits <- c(hits, hits_rev)
  }

  ## Setup the return object as empty for any cases without matches
  cols <- c(
    "seq", "score", "direction", "start", "end", "from_centre", "seq_width",
    "match"
  )
  ## Form it as a list, the eventually wrap into a DF
  out <- as.list(mcols(hits))
  int_cols <- c("start", "end", "seq_width")
  out[int_cols] <- lapply(int_cols, \(x) integer())
  out$from_centre <- numeric()
  if (!"score" %in% names(out)) out$score <- numeric()
  out$seq <- character()
  out$match <- DNAStringSet()
  if (length(hits) == 0) return(DataFrame(out)[,cols])

  ## Map back to the original Views
  hits_to_map <- findInterval(start(hits), map$start)
  w <- width(hits)

  ## Form the output
  out$seq <- hits_to_map
  out$start <- as.integer(start(hits) - c(0, map$end)[hits_to_map])
  out$end <- as.integer(out$start + w - 1)
  out$seq_width <- width(stringset[out$seq])
  out$from_centre <- 0.5 * (out$start + out$end) - out$seq_width / 2

  ## The match itself
  out$direction <- factor(mcols(hits)$direction, levels = c("F", "R"))
  to_rev <- out$direction == "R"
  out$match <- as(hits, "XStringSet")
  out$match[to_rev] <- reverseComplement(out$match[to_rev])

  ## Setup any named strings to appear in the same order
  names_are_char <- is.character(names(stringset))
  if (names_are_char) {
    out$seq <- factor(map$names[hits_to_map], map$names)
    out$seq <- droplevels(out$seq)
  }

  ## The final object
  out <- DataFrame(out)[, cols]
  if (best_only) {
    break_ties <- match.arg(break_ties)
    out <- .getBestMatch(out, break_ties)
  }
  o <- order(out$seq, out$start)
  if (names_are_char) out$seq <- as.character(out$seq)
  out[o, ]

}

#' @importFrom S4Vectors splitAsList
#' @keywords internal
.getBestMatch <- function(df, ties){

  if (nrow(df) == 0) return(df)

  ## Decide how to break ties
  pos <- sample.int(nrow(df), nrow(df)) # Default to random
  if (ties == "first") pos <- df$start
  if (ties == "last") pos <- (-1) * df$end
  if (ties == "central") pos <- abs(df$from_centre)
  if (ties == "all") pos <- rep_len(1, nrow(df))

  ## Pick the best match using the values we have
  o <- order(df$seq, 1 / df$score, pos)
  df <- df[o,]
  if (ties != "all") {
    df <- df[!duplicated(df$seq), ]
  } else {
   df_list <- splitAsList(df, df$seq)
   df_list <- lapply(df_list, \(x) subset(x, score == max(score)))
   df <- do.call("rbind", df_list)
  }
  df

}

