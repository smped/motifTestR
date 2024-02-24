#' @keywords internal
.checkPWM <- function(pwm){
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
}
