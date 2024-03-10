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
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#'
#' ## Load the example Peaks
#' data("ar_er_seq")
#' countPwmMatches(esr1, ar_er_seq)
#'
#' ## Count all PWMs
#' countPwmMatches(ex_pwm, ar_er_seq)
#'
#' @importFrom parallel mclapply
#' @export
countPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "80%", mc.cores = 1, ...
) {

    args <- c(as.list(environment()), list(...))
    args <- args[names(args) != "mc.cores"]
    if (is.list(pwm)) {
        pwm <- .cleanMotifList(pwm)
        out <- mclapply(
            pwm, .countSinglePwmMatches, stringset = stringset, rc = rc,
            min_score = min_score, mc.cores = mc.cores
        )
        out <- unlist(out)
    }
    if (is.matrix(pwm)) out <- do.call(".countSinglePwmMatches", args)
    out

}

#' @import Biostrings
#' @keywords internal
.countSinglePwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "80%", ...
){
    ## Checks & the map
    pwm <- .checkPWM(pwm)
    map <- .viewMapFromXStringset(stringset)

    # Form the entire XStringSetList into a Views object
    views <- Views(
        unlist(stringset), start = map$start, width = map$width,
        names = map$names
    )
    n_matches <- countPWM(pwm, views, ...)
    if (rc)
        n_matches <- c(n_matches, countPWM(reverseComplement(pwm), views, ...))

    as.integer(sum(n_matches))
}
