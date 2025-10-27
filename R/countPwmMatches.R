#' @title Count the matches to a PWM within an XStringSet
#'
#' @description Count the matches to a PWM within an XStringSet
#'
#' @details
#' Will simply count the matches within an XStringSet and return an integer.
#' All matches are included.
#'
#' @return An integer vector
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param mc.cores Passed to \link[parallel]{mclapply} when analysing a list of
#' PWMs
#' @param ... Passed to \link[Biostrings]{countPWM}
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pfm")
#' esr1 <- ex_pfm$ESR1
#'
#' ## Load the example Peaks
#' data("ar_er_seq")
#' countPwmMatches(esr1, ar_er_seq)
#'
#' ## Count all PWMs
#' countPwmMatches(ex_pfm, ar_er_seq)
#'
#' @importFrom parallel mclapply
#' @importFrom matrixStats rowSums2
#' @import Biostrings
#' @export
countPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "80%", mc.cores = 1, ...
) {

    if (!is.list(pwm)) pwm <- list(pwm)
    mc.cores <- min(mc.cores, length(pwm))
    pwm <- mclapply(pwm, .checkPWM, mc.cores = mc.cores)
    pwm <- .cleanMotifList(pwm)
    nm <- names(pwm)
    ## Append any reverse matrices to the existing list
    if (rc) pwm <- c(pwm, lapply(pwm, reverseComplement))

    # Form the entire XStringSetList into a Views object
    map <- .viewMapFromXStringset(stringset)
    views <- Views(
        unlist(stringset), start = map$start, width = map$width,
        names = map$names
    )

    ## Now the counts
    counts <- mclapply(
        pwm,
        \(x) countPWM(x, views, min.score = min_score, ...),
        mc.cores = mc.cores
    )
    out <- unlist(counts)
    if (rc) out <- rowSums2(matrix(out, ncol = 2, dimnames = list(nm, NULL)))
    names(out) <- nm
    out

}

