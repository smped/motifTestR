#' @title Simulate sequences using optional TFBMs
#'
#' @description
#' Simulate a set of fixed-width sequences using optional TFBMs
#'
#' @details
#' Using the nucleotide and probabilities provided as set of sequences can be
#' simulated. By default, this will effectively be a set of 'background'
#' sequences, with letters effectively chosen at random.
#'
#' If a PWM/PFM is supplied, the shape parameters are first passed to
#' \link[VGAM]{rbetabinom.ab} to determine the random positions the motif will
#' be placed, with the default parameters representing a discrete uniform
#' distribution.
#' Once positions for the TFBM have been selected, nucleotides will be randomly
#' sampled using the probabilities provided in the PWM and these motifs will be
#' placed at the randomly sample positions
#'
#' @return
#' By default a DNAStringSet will be returned.
#' If possible, the position of any randomly sampled motifs will be included
#' in the mcols element of the returned object.
#'
#' @param n The number of sequences to simulate
#' @param width Width of sequences to simulate
#' @param pfm Probability Weight/Frequency Matrix
#' @param nt Nucleotides to include
#' @param prob Sampling probablities for each nucleotide
#' @param shape1,shape2 Passed to \link[VGAM]{rbetabinom.ab}
#' @param as ObjectClass to return objects as. Defaults to DNAStringSet, but
#' other viable options may include 'character', 'CharacterList' or any
#' other class from which a character vector may be coerced.
#' @param ... Not used
#'
#' @examples
#' ## Randomly generate 10x50nt sequences without any TFBMs present
#' simSeq(10, 50)
#'
#' ## Now place a motif at random positions
#' data('ex_pfm')
#' sim_seq <- simSeq(10, width = 20, pfm = ex_pfm$ESR1)
#' sim_seq
#' ## The position of the motif within each sequence is included in the mcols
#' mcols(sim_seq)
#'
#' ## Use this to extract the random motifs from the random sequences
#' library(IRanges)
#' i <- mcols(sim_seq)$pos + cumsum(width(sim_seq)) - width(sim_seq)
#' Views(unlist(sim_seq), start = i, width = 10)
#'
#'
#' @importFrom S4Vectors mcols<-
#' @importFrom methods slot is
#' @export
simSeq <- function(
        n, width, pfm = NULL, nt = c("A", "C", "G", "T"), prob = rep(0.25, 4),
        shape1 = 1, shape2 = 1, as = "DNAStringSet", ...
){

    ## Assuming an even frequency, create a vector randomly
    prob <- rep_len(prob, length(nt))
    bg <- sample(nt, n * width, replace = TRUE, prob = prob)
    pos <- NULL

    ## If a PWM is provided, now sample using the motifs
    if (!is.null(pfm)) {

        if (!requireNamespace('VGAM', quietly = TRUE))
            stop("Please install 'VGAM' to inject TFBMs into the sequences.")

        ## Check we have PFMs, not PWMs or any other format
        if (is(pfm, "universalmotif")) pfm <- slot(pfm, "motif")
        col_sums <- colSums(pfm)
        is_pfm <- isTRUE(
            all.equal(unname(col_sums), rep(1, length(col_sums)))
        )
        if (!is_pfm) {
            ## Convert to a PWM, then to PFM using existing checks
            pwm <- .checkPWM(pfm)
            pfm <- slot(create_motif(pwm, type = "PPM"), "motif")
        }
        stopifnot(all(rownames(pfm) %in% nt)) # Same alphabet

        ## Sample the random motifs as a matrix, then coerce to a vector
        pfm_width <- ncol(pfm)
        stopifnot(pfm_width <= width)
        rnd_mot <- replicate(
            n, apply(pfm, MARGIN = 2, FUN = \(p) sample(nt, 1, prob = p))
        )

        ## Determine the valid positions given the width of the PFM.
        max_start <- width - pfm_width
        pos <- VGAM::rbetabinom.ab(n, max_start, shape1, shape2) + 1

        ## Inject into the existing sequences. This is faster treating
        ## bg as a vector, not a matrix to be iterated through.
        i <- seq(0, n - 1) * width + pos
        vec_pos <- vapply(
            i, \(i) seq(i, length.out = pfm_width), numeric(pfm_width)
        )
        vec_pos <- as.integer(vec_pos)
        bg[vec_pos] <- as.character(rnd_mot)

    }
    seq <- apply(matrix(bg, ncol = n), MARGIN = 2, paste, collapse = "")
    seq <- as(seq, as)
    if (is(seq, "Vector")) mcols(seq)$pos <- pos ## The base class with mcols
    seq

}
