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
