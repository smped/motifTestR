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
#'
#' The sequences to have a motif inserted will be selected, along with the
#' number of motifs, using the rate and theta parameters.
#' If both are NULL, every sequence will have a single motif inserted.
#' If the rate is > 0 and theta is NULL, sequences will be selected to have
#' motifs inserted using a poisson distribution.
#' If theta is also provided, sequences will be selected to contain motifs
#' using a negative binomial distribution
#'
#' Once positions and sequences for the TFBM have been selected, nucleotides
#' will be randomly sampled using the probabilities provided in the PWM and
#' these motifs will be placed at the randomly sampled positions
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
#' @param prob Sampling probabilities for each nucleotide
#' @param shape1,shape2 Passed to \link[VGAM]{rbetabinom.ab}
#' @param rate The expected rate of motifs per sequence. Is equivalent to
#' \eqn{ \lambda } in \link[stats]{rpois}. If set to NULL, all sequences will
#' be simulated with a single motif, otherwise a Poisson distribution will be used
#' @param theta Overdispersion parameter passed to \link[MASS]{rnegbin}.
#' If set to NULL the rate parameter will be passed to \link[stats]{rpois}.
#' However if this value is set, the rate and theta parameters are passed to
#' \link[MASS]{rnegbin} to simulate overdispersed counts
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
        shape1 = 1, shape2 = 1, rate = NULL, theta = NULL, as = "DNAStringSet",
        ...
){

    ## Assuming an even frequency, create a vector randomly
    prob <- rep_len(prob, length(nt))
    bg <- sample(nt, n * width, replace = TRUE, prob = prob)
    seq_starts <- seq(1, n*width, by = width) # Where each sequence starts
    pos_vec <- NULL

    ## If a PWM is provided, now sample using the motifs
    if (!is.null(pfm)) {

        if (!requireNamespace('VGAM', quietly = TRUE))
            stop("Please install 'VGAM' to insert TFBMs into the sequences.")

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

        pfm_width <- ncol(pfm)
        stopifnot(pfm_width <= width)
        ## Determine the valid positions given the width of the PFM.
        max_start <- width - pfm_width

        ## We really need to just choose the positions using the
        ## different distributions. Everything else can follow.
        ## Placing them in the mcols at the end will take some thought though

        if (is.null(rate)) {
            pos <- VGAM::rbetabinom.ab(n, max_start, shape1, shape2) + seq_starts
            ## Inject into the existing sequences. This is faster treating
            ## bg as a vector, not a matrix to be iterated through.
            vec_pos <- vapply(
                pos, \(i) seq(i, length.out = pfm_width), numeric(pfm_width)
            )
            vec_pos <- as.integer(vec_pos)
        } else {
            stopifnot(rate > 0)
            if (is.null(theta)) {
                pos <- .simPoisSeq(n, rate, shape1, shape2, max_start, seq_starts)
            } else {
                stopifnot(theta > 0)
                pos <- .simNBSeq(n, rate, theta, shape1, shape2, max_start, seq_starts)
            }
            vec_pos <- vapply(
                pos, \(i) seq(i, length.out = pfm_width), numeric(pfm_width)
            )
            vec_pos <- as.integer(vec_pos)
        }

        ## Sample the random motifs as a matrix, then coerce to a vector
        rnd_mot <- replicate(
            length(pos),
            apply(pfm, MARGIN = 2, FUN = \(p) sample(nt, 1, prob = p))
        )

        ## Handle any overlapping positions by treating them as duplicates
        ## This will result in some partial motifs being placed
        dups <- duplicated(vec_pos)
        bg[vec_pos[!dups]] <- as.character(rnd_mot)[!dups]

        ## Now setup for inclusion in the mcols
        pos_vec <- rep_len(NA, n)
        temp_pos <- pos[apply(matrix(!dups, ncol = length(pos)), MARGIN = 2, all)]
        i <- ceiling(temp_pos / width)
        pos_vec[i] <- temp_pos %% width

    }
    seq <- apply(matrix(bg, ncol = n), MARGIN = 2, paste, collapse = "")
    seq <- as(seq, as)
    if (is(seq, "Vector")) mcols(seq)$pos <- pos_vec ## The base class with mcols
    seq

}


#' @keywords internal
.simPoisSeq <- function(n, rate, shape1, shape2, max_start, seq_starts) {

    ## Unlike simulating for all sequences, we need a list here.
    ## First get the number of motifs per sequences
    n_pois <- stats::rpois(n, rate)
    which_seq <- which(n_pois > 0)
    pos <- lapply(
        n_pois[which_seq],
        \(x) VGAM::rbetabinom.ab(x, max_start, shape1, shape2)
    )

    ## Expand to the global position with the random sequence
    seq_starts <- seq_starts[which_seq]
    motif_starts <- lapply(
        seq_along(seq_starts), \(i) seq_starts[i] + pos[[i]] ## may be variable length
    )
    unlist(motif_starts)

}

#' @keywords internal
.simNBSeq <- function(n, rate, theta, shape1, shape2, max_start, seq_starts) {

    if (!requireNamespace('MASS', quietly = TRUE))
        stop("Please install 'MASS' to inject over-dispersed TFBMs into the sequences.")

    ## Unlike simulating for all sequences, we need a list here.
    ## First get the number of motifs per sequences
    n_negbin <- MASS::rnegbin(n, rate, theta)
    which_seq <- which(n_negbin > 0)
    pos <- lapply(
        n_negbin[which_seq],
        \(x) VGAM::rbetabinom.ab(x, max_start, shape1, shape2)
    )

    ## Expand to the global position with the random sequence
    seq_starts <- seq_starts[which_seq]
    motif_starts <- lapply(
        seq_along(seq_starts), \(i) seq_starts[i] + pos[[i]] ## may be variable length
    )
    unlist(motif_starts)

}
