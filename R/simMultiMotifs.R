#' @title Simulate sequences with multiple motifs
#'
#' @description
#' Simulate a set of sequences incorporating multiple motifs
#'
#' @details
#' Simulate a set of sequences with multiple motifs inserted using different
#' rates and distributions, as specified.
#' All shape, rate and theta parameters are recycled to match the length of the
#' supplied motif list, and can be supplied as vectors to tailor these
#' parameters to each provided element of the list of matrices
#'
#' @return
#' A DNAStringSet with mcols denoting the positions of all inserted motifs
#'
#' @param n The number of sequences to simulate
#' @param width Width of sequences to simulate
#' @param pfm List of Probability Weight/Frequency Matrices
#' @param bg Optional, pre-defined set of background sequences. Can be passed as
#' an XStringSet or character vector. All sequences must be the same width
#' @param nt Nucleotides to include
#' @param prob Sampling probabilities for each nucleotide
#' @param shape1,shape2 Passed to \link[VGAM]{rbetabinom.ab}
#' @param rate The expected rate of motifs per sequence. Is equivalent to
#' \eqn{ \lambda } in \link[stats]{rpois}. If set to NULL or NA, all sequences will
#' be simulated with a single motif, otherwise a Poisson distribution will be used
#' @param theta Overdispersion parameter passed to \link[MASS]{rnegbin}.
#' If set to NULL or NA the rate parameter will be passed to \link[stats]{rpois}.
#' However if this value is set, the rate and theta parameters are passed to
#' \link[MASS]{rnegbin} to simulate overdispersed counts
#' @param as ObjectClass to return objects as. Defaults to DNAStringSet, but
#' other viable options may include 'character', 'CharacterList' or any
#' other class from which a character vector may be coerced.
#' @param ol When randomly simulated positions overlap, choose one either at
#' random, by the first occurring PFM in the list of PFMs, or by the last.
#' @param ... Not used
#'
#' @examples
#' data("ex_pfm")
#' ## Simulate sequences including both ESR1 and ANDR, but with
#' ## ESR1 being included at a higher rate
#' seq <- simMultiMotifs(10, 100, ex_pfm[1:2], rate = c(2, 1))
#' seq
#' ## The positions of the motifs are included in the mcols
#' mcols(seq)
#'
#'
#' @export
simMultiMotifs <- function(
    n, width, pfm = NULL, bg = NULL, nt = c("A", "C", "G", "T"),
    prob = rep(0.25, 4), shape1 = 1, shape2 = shape1, rate = NA, theta = NA,
    as = "DNAStringSet", ol = c("random", "first", "last"), ...
){

  args <- c(as.list(environment()), list(...))
  args$bg <- NULL
  ## Pass to simSeq if there are no motifs
  if (is.null(pfm)) return(do.call(motifTestR::simSeq, args))
  prob <- rep_len(prob, length(nt))

  ## Setup the original (BG) set of sequences to be character(1)
  bg <- .defineBG(bg, n, width, nt, prob)
  seq_starts <- 1 + c(0, cumsum(nchar(bg)))[seq_along(bg)]
  n <- length(seq_starts)
  bg <- strsplit(paste(bg, collapse = ""), "")[[1]]
  width <- as.integer(length(bg) / n)

  ## Check the PFMs. Should be a list here! And all need to be PFMs not PWMs
  pfm <- .checkPfmList(pfm)

  ## Now expand the relevant params to match the nbr of motifs
  n_motifs <- length(pfm)
  shape1 <- rep_len(shape1, n_motifs)
  shape2 <- rep_len(shape2, n_motifs)
  if (!is.null(rate)) rate <- rep_len(rate, n_motifs)
  if (!is.null(theta)) theta <- rep_len(theta, n_motifs)
  pfm_width <- vapply(pfm, ncol, integer(1))
  max_start <- width - pfm_width
  stopifnot(all(max_start > 0))

  ## From here we step through each motif, or do we lapply the positions
  ## and sampled motifs and inject all at once? Not sure how to handle
  ## overlaps either
  pos_list <- mapply(
    .samplePos, rate = rate, theta = theta,
    shape1 = shape1, shape2 = shape2, max_start = max_start,
    MoreArgs = list(n = n, seq_starts = seq_starts), SIMPLIFY = FALSE
  )
  ## We need to ensure pfm is a named list!!!
  vec_pos <- mapply(
    \(x, y) {
      vapply(x, \(i) seq(i, length.out = y), numeric(y)) |> as.integer()
    }, x = pos_list, y = pfm_width, SIMPLIFY = FALSE
  )
  ## Now sample the new nucleotides
  new_nuc <- mapply(
    \(x, y) {
      replicate(
        length(x), apply(y, MARGIN = 2, FUN = \(p) sample(nt, 1, prob = p))
      ) |> as.character()
    }, x = pos_list, y = pfm, SIMPLIFY = FALSE
  )
  index_df <- data.frame(i = unlist(vec_pos), nt = unlist(new_nuc))

  ## Handle overlapping positions
  ol <- match.arg(ol)
  if (ol == "random") # Randomise the order
    index_df <- index_df[sample.int(nrow(index_df), nrow(index_df)),]
  if (ol == "last") # Reverse the order
    index_df <- index_df[seq(nrow(index_df), 1),]
  ## Keep the first of any duplicate rows using the above re-ordering
  index_df <- index_df[!duplicated(index_df[["i"]]),]

  ## Create the final set of sequences
  bg[index_df$i] <- index_df$nt
  seq <- apply(matrix(bg, ncol = n), MARGIN = 2, paste, collapse = "")
  seq <- as(seq, as)
  ## Now update the mcols if required
  if (is(seq, "Vector")) seq <- .addPosList(seq, pos_list, width, pfm)
  seq

}

#' @keywords internal
#' @importFrom IRanges IntegerList
#' @importFrom S4Vectors DataFrame endoapply
.addPosList <- function(seq, pos_list, width, pfm) {
  df_list <- lapply(
    pos_list,
    \(pos) {
      which_seq <- ceiling(pos / width)
      int_list <- IntegerList(vector("list", length(seq)))
      int_list[unique(which_seq)] <- split(pos %% width, f = which_seq)
      out <- endoapply(endoapply(int_list, unique), sort)
      n_motifs <- vapply(out, length, integer(1))
      if (all(n_motifs == 1)) out <- unlist(out)
      out
    }
  )
  nm <- names(pfm)
  if (is.null(nm)) nm <- paste0("X", seq_along(pfm))
  names(df_list) <- nm
  counts <- lapply(df_list, \(x) vapply(x, length, integer(1)))
  df_list$n_motifs <- rowSums(do.call("cbind", counts))
  DF <- DataFrame(df_list)
  mcols(seq) <- DF
  seq
}

#' @importFrom methods is
#' @importFrom IRanges width
.defineBG <- function(bg, n, wd, nt, prob) {
  ## Return a character vector which will be collapsed
  if (is.null(bg))
    bg <- simSeq(n, wd, nt = nt, prob = prob, as = "character")

  if (is(bg, "XStringSet")) {
    stopifnot(length(unique(width(bg))) == 1) # Fixed width check
    bg <- as.character(bg)
  }

  ## Must be a character beyond this point
  stopifnot(is.character(bg))

  ## Ensure that only the passed nucleotides are present
  pat <- paste0("[^", paste(nt, collapse = ""), "].+")
  has_non_nuc <- vapply(bg, \(x) grepl(pat, x), logical(1))
  stopifnot(all(!has_non_nuc))

  ## These need to be the same width for sampling motif positions
  w <- unique(nchar(bg))
  stopifnot(length(w) == 1)
  bg

}

#' @importFrom universalmotif create_motif
.checkPfmList <- function(pfm) {
  stopifnot(is.list(pfm))
  if (length(pfm) == 1) stop("Please use simSeq for a single motif")
  is_uni <- vapply(pfm, is, logical(1), "universalmotif")
  if (all(is_uni)) {
    nm <- vapply(pfm, slot, character(1), "name")
    pfm <- convert_type(pfm, "PPM")
    pfm <- lapply(pfm, slot, "motif")
    names(pfm) <- nm
  } else {
    col_sums <- lapply(pfm, \(x) unname(colSums(x)))
    is_pfm <- lapply(col_sums, \(x) all.equal(sum(x), length(x), tolerance = 1e-6))
    is_pfm <- unlist(is_pfm)
    if (is.character(is_pfm)) {
      ## Perform the conversion using universal motif
      nm <- names(pfm)
      pfm <- lapply(pfm, create_motif, type = "PPM")
      pfm <- lapply(pfm, slot, "motif")
      names(pfm) <- nm
    }
  }
  stopifnot(all(vapply(pfm, is.matrix, logical(1))))
  pfm
}
