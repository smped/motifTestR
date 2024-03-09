#' @title Form a set of random, matching ranges for bootstrapping or permuting
#'
#' @description
#' Form a set of ranges from y which (near) exactly match those in x for use as
#' a background set requiring matching
#'
#' @details
#' This function uses the width distribution of the 'test' ranges (i.e. `x`) to
#' randomly sample a set of ranges with matching width from the ranges provided
#' in `y`. The width distribution will clearly be exact when a set of
#' fixed-width ranges is passed to `x`, whilst random sampling may yield some
#' variability when matching ranges of variable width.
#'
#' When both x and y are GRanges objects, they are implcitly assumed to both
#' represent similar ranges, such as those overlapping a promoter or enhancer.
#' When passing two GRangesList objects, both objects are expected to contain
#' ranges annotated as belonging to key features, such that the list elements in
#' y must encompass all elements in x.
#' For example if `x` contains two elements named 'promoter' and 'intron', `y`
#' should also contain elements named 'promoter' and 'intron' and these will
#' be sampled as matching ranges for the same element in `x`. If elements of
#' `x` and `y` are not named, they are assumed to be in matching order.
#'
#' The default behaviour is to assume that randomly-generated ranges are for
#' iteration, and as such, ranges are randomly formed in multiples of the number
#' of 'test' ranges provided in `x`. The column `iteration` will be added to the
#' returned ranges.
#' Placing any number into the `n_total` argument will instead select a total
#' number of ranges as specified here. In this case, no `iteration` column will
#' be included in the returned ranges.
#'
#' Sampling is assumed to be with replacement as this is most suitable for
#' bootstrapping and related procedures, although this can be disabled by
#' setting `replace = FALSE`
#'
#' @return A GRanges or GRangesList object
#'
#'
#' @param x GRanges with ranges to be matched
#' @param y GRanges with ranges to select random matching ranges from
#' @param exclude GRanges of ranges to omit from testing
#' @param n_iter The number of times to repeat the random selection process
#' @param n_total Setting this value will over-ride anything set using n_iter
#' @param replace logical(1) Sample with our without replacement when creating
#' the set of random ranges.
#' @param mc.cores Passsed to \link[parallel]{mclapply}
#' @param ... Not used
#' @param force_ol logical(1) Enforce an overlap between every site in x and y
#' @param unlist logical(1) Return as a sorted GRanges object, or leave as a
#' GRangesList
#'
#' @examples
#' ## Load the example peaks
#' data("ar_er_peaks")
#' sq <- seqinfo(ar_er_peaks)
#' ## Now sample size-matched ranges for two iterations from chr1
#' makeRMRanges(ar_er_peaks, GRanges(sq)[1], n_iter = 2)
#'
#' ## Or simply sample 100 ranges if not planning any iterative analyses
#' makeRMRanges(ar_er_peaks, GRanges(sq)[1], n_total = 100)
#'
#'
#' @rdname makeRMRanges-methods
#' @aliases makeRMRanges
#' @export
setGeneric(
  "makeRMRanges",
  function(x, y, ...) standardGeneric("makeRMRanges")
)
#' @import GenomicRanges
#' @importFrom IRanges overlapsAny width
#' @importFrom GenomeInfoDb seqinfo
#' @rdname makeRMRanges-methods
#' @aliases makeRMRanges
#' @export
setMethod(
  "makeRMRanges", signature = signature(x = "GRanges", y = "GRanges"),
  function(
    x, y, exclude = GRanges(), n_iter = 1, n_total = NULL, replace = TRUE, ...,
    force_ol = TRUE
  ) {

    ## Check some basic features, e.g. y must be bigger than x
    stopifnot(sum(width(y)) > sum(width(x)))
    ## There must be overlap if bootstrapping.
    if (force_ol) stopifnot(all(overlapsAny(x, y)))

    # Make sure any new ranges don't fall off the end of chromosomes
    stopifnot(is(exclude, "GRanges"))
    sq_gr <- GRanges(seqinfo(x))
    w <- width(x)
    max_w <- max(w)
    chr_ends <- c(
      resize(sq_gr, width = max_w / 2, fix = "start"),
      resize(sq_gr, width = max_w / 2, fix = "end")
    )
    # Chop the ends off so resizing won't spit an error
    exclude <- GenomicRanges::setdiff(exclude, chr_ends)
    exclude <- resize(exclude, width = width(exclude) + max_w, fix = "center")
    exclude <- c(granges(exclude), chr_ends)
    y <- GenomicRanges::setdiff(y, exclude)

    ## Set the iterations
    has_iter <- is.null(n_total)
    if (has_iter) {
      iter_size <- length(x)
      n_total <- as.integer(n_iter * iter_size)
    }

    ## Now the random sampling
    rand_w <- sample(w, n_total, replace = TRUE)
    too_big <- sum(width(y)) > .Machine$integer.max
    if (too_big) {
      ## The human genome is ~3.1e9 & the usual max integer is 2.1e9
      ## Splitting in half should solve problems for all genomes < 4e9
      ## Sorry to wheat people...
      i <- floor(length(y) / 2)
      n1 <- floor(n_total / 2)
      y1 <- y[seq_len(i)]
      i1 <- sample(sum(width(y1)), n1, replace = replace)
      gr1 <- GPos(y1)[i1]

      y2 <- GenomicRanges::setdiff(y, y1)
      n2 <- n_total - n1
      i2 <- sample(sum(width(y2)), n2, replace = replace)
      gr2 <- GPos(y2)[i2]

      ## Randomise to make up for chunking
      i <- sample(n_total, n_total, replace = FALSE)
      gr <- GRanges(c(gr1, gr2))[i]

    } else {
      gpos <- GPos(y)
      i <- sample(seq_along(gpos), n_total, replace = replace)
      gr <- GRanges(gpos[i])
    }
    gr <- resize(gr, width = rand_w, fix = "center")
    ## Ranges will be randomised so adding iterations in order is fine
    if (has_iter) gr$iteration <- rep(seq_len(n_iter), times = iter_size)
    sort(gr)

  }
)
#' @import GenomicRanges
#' @importFrom parallel mclapply
#' @rdname makeRMRanges-methods
#' @aliases makeRMRanges
#' @export
setMethod(
  "makeRMRanges",
  signature = signature(x = "GRangesList", y = "GRangesList"),
  function(
    x, y, exclude = GRanges(), n_iter = 1, n_total = NULL, replace = TRUE,
    mc.cores = 1, ..., force_ol = TRUE, unlist = TRUE
  ){

    stopifnot(length(x) <= length(y))
    nm_x <- names(x)
    nm_y <- names(y)
    if (!is.null(nm_x) & !is.null(nm_y)) {
      ## Set GRLists in matching order if named
      stopifnot(all(nm_x %in% nm_y))
      y <- y[nm_x]
    }
    i <- seq_along(x)
    n_iter <- rep_len(n_iter, length(x))
    n_total <- lapply(i, \(j) n_total[[j]])
    out <- mclapply(
      i,
      \(i) makeRMRanges(
        x[[i]], y[[i]], exclude = exclude,  n_iter = n_iter[[i]],
        n_total = n_total[[i]], replace = replace, ..., force_ol = force_ol
      ), mc.cores = mc.cores
    )
    if (!is.null(nm_x)) names(out) <- nm_x
    out <- GRangesList(out)
    if (unlist) out <- sort(unlist(out))
    out
  }
)
