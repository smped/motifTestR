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
#' @param ... Passed to \link{getPwmMatches}
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
#' @export
countPwmMatches <- function(
        pwm, stringset, rc = TRUE, min_score = "50%", mc.cores = 1, ...
) {

    matches <- getPwmMatches(
        pwm = pwm, stringset = stringset, rc = rc, min_score = min_score,
        mc.cores = mc.cores, ...
    )
    if (is(matches, "SplitDataFrameList")){
        n <- vapply(matches, nrow, integer(1))
    } else {
        n <- nrow(matches)
    }
    n

}

