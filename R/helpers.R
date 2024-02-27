#' @keywords internal
.checkPWM <- function(pwm){
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
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
