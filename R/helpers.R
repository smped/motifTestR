#' @importFrom methods slot is
#' @importClassesFrom universalmotif universalmotif
#' @keywords internal
.checkPWM <- function(pwm){
  if (is(pwm, "universalmotif")) pwm <- slot(pwm, "motif")
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
  pwm
}

#' @keywords internal
.viewMapFromXStringset <- function(stringset){
  stopifnot(is(stringset, "XStringSet"))
  # Form a map from the original sequences to a views object
  map <- data.frame(end = cumsum(width(stringset)), width = width(stringset))
  map$start <- map$end - map$width + 1
  map$names <- names(stringset)
  map
}

#' @keywords internal
#' @importClassesFrom universalmotif universalmotif
.cleanMotifList <- function(x) {
  if (all(vapply(x, is, logical(1), "universalmotif"))) {
    nm <- vapply(x, slot, character(1), "name")
    x <- lapply(x, slot, "motif")
    names(x) <- nm
  }
  ## Now check everything is a matrix
  all_mat <- vapply(
    x,
    \(x) all(is.matrix(x), rownames(x) == c("A", "C", "G", "T")),
    logical(1)
  )
  stopifnot(all_mat)
  x
}

#' @keywords internal
.checkMatches <- function(matches) {
  colTypes <- c(
    score = "numeric", direction = "factor", start = "integer", end = "integer",
    from_centre = "numeric", seq_width = "integer", match = "DNAStringSet"
  )
  reqdCols <- c("seq", names(colTypes)) # seq can be any type

  if (is(matches, "DataFrame")) {
    stopifnot(all(reqdCols %in% colnames(matches)))
    types <- vapply(matches, \(x) is(x)[[1]], character(1))
    stopifnot(
      identical(types[names(colTypes)], colTypes)
    )
  } else {
    stopifnot(is(matches, "list") & !is(matches, "data.frame"))
    hasCols <- vapply(matches, \(x) all(reqdCols %in% colnames(x)), logical(1))
    if (any(!hasCols))
      stop("All objects must contain the columns: \n", paste(reqdCols, "\n"))
    correctTypes <- vapply(
      matches,
      \(x) {
        types <- vapply(x, \(cols) is(cols)[[1]], character(1))
        identical(types[names(colTypes)], colTypes)
      }, logical(1)
    )
    if (any(!correctTypes)) stop("All columns are not of the correct type")
  }

  invisible(TRUE)
}

#' @keywords internal
.makeBmBins <- function(matches, binwidth, abs){

  ## Set sequence weights for multiple matches
  reps <- table(as.character(matches$seq)) # Counts reps
  matches$weight <- as.numeric(1 / reps[as.character(matches$seq)]) # 1 / reps

  ## Make the bins
  longest_seq <- max(matches$seq_width)
  stopifnot(binwidth < longest_seq)
  if (abs) {
    ## Set the dist from centre to be +ve only
    matches$from_centre <- abs(matches$from_centre)
    all_breaks <- c(seq(0, longest_seq / 2, by = binwidth), longest_seq / 2)
    seq_per_bin <- vapply(all_breaks, \(x) sum(matches$seq_width >= x), integer(1))
  } else {
    ## To define bins, make sure the central bin is around zero and that
    ## bins extends symetrically from this point to the widest sequence
    pos_breaks <- c(
      ## Ensure the upper limit is included, even if it leaves a smaller bin
      seq(binwidth / 2, longest_seq / 2, by = binwidth), longest_seq / 2
    )
    all_breaks <- c(-1 * pos_breaks, pos_breaks)
    pos_bins <- vapply(pos_breaks, \(x) sum(matches$seq_width >= x), integer(1))
    seq_per_bin <- c(pos_bins, pos_bins)
  }

  ## Cut into bins
  all_breaks <- unique(sort(all_breaks))
  matches$bin <- cut(matches$from_centre, breaks = all_breaks, include.lowest = TRUE)
  matches
}
