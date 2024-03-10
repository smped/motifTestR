#' motifTestR: Perform Key Analyses on Transcription Factor Binding Motifs
#'
#' The package `motifTestR` has been designed for two primary analyses of TFBMs,
#' testing for positional bias and overall enrichment.
#'
#' The package `motifTestR` provides two primary functions for testing TFBMs
#' within a set of sequences
#'
#' * [testMotifPos()] for detecting positional bias within a set of test
#' sequences
#' * [testMotifEnrich()] for testing overall enrichment of a TFBM within a set
#' of test sequences
#'
#' The main functions rely on lower-level functions such as:
#'
#' * [countPwmMatches()] simply counts the number of matches within an
#' `XStringSet`
#' * [getPwmMatches()] returns the position of matches within an `XStringSet`
#' * [makeRMRanges()] which produces a set of random, matching ranges based on
#' key characteristics of the set of test sequences/ranges
#'
#' A simple utility function is provided yo enable visualisation of results
#'
#' * [plotMatchPos()] enables visualisation of the matches within a set of
#' sequences using multiple strategies
#'
#' @author
#' Stevie Pederson
#'
#' @keywords internal
"_PACKAGE"
