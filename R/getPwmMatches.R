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
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param min_score The minimum score to return a match
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#' @return A DataFrame with columns: `seq`, `start`, `end`, `width`, `match`,
#' `score` and `direction`
#' The first four columns indicate the position of the match within the query
#' stringset, whilst `match` returns the matching sequence as an `XStringSet`.
#' The `direction` column denotes `F` and `R` as a forward or reverse match
#' respectively
#'
#' @import Biostrings
#' @importFrom IRanges Views
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom methods as is
#' @importFrom stats setNames
#' @export
getPwmMatches <- function(pwm, stringset, rc = TRUE, min_score = "80%", ...){

  ## Checks
  .checkPWM(pwm)
  stopifnot(is(stringset, "XStringSet"))

  # Form a map from the original sequences to a views object
  map <- data.frame(end = cumsum(width(stringset)), width = width(stringset))
  map$start <- map$end - map$width + 1
  map$names <- names(stringset)

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
    "seq", "score", "direction", "start", "end", "width", "from_centre", "match"
  )
  ## Form it as a list, the eventually wrap into a DF
  out <- as.list(mcols(hits))
  int_cols <- c("start", "end", "width")
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
  out$width <- w
  seq_widths <- width(stringset[out$seq])
  out$from_centre <- 0.5 * (out$start + out$end) - seq_widths / 2

  ## The match itself
  out$match <- as(hits, "XStringSet")
  to_rev <- mcols(hits)$direction == "R"
  out$match[to_rev] <- reverseComplement(out$match[to_rev])
  out$direction <- factor(mcols(hits)$direction, levels = c("F", "R"))

  ## Setup any named strings to appear in the same order
  names_are_char <- is.character(names(stringset))
  if (names_are_char) {
    out$seq <- factor(map$names[hits_to_map], map$names)
    out$seq <- droplevels(out$seq)
  }
  o <- order(out$seq, out$start)
  if (names_are_char) out$seq <- as.character(out$seq)

  ## The final object
  DataFrame(out)[o, cols]
}
