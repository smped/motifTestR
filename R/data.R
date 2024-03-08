#' @title Example Position Weight Matrices
#'
#' @description Example Position Weight Matrices
#'
#' @details
#' This object contains 5 PWMs taken from HOCOMOCOv11-coreA for examples and testing
#'
#' @usage data("ex_pwm")
#'
#' @examples
#' data("ex_pwm")
#' ex_pwm$ESR1
#'
#' @name ex_pwm
#' @rdname ex_pwm
"ex_pwm"

#' @title A set of peaks with AR and ER detected
#'
#' @description
#' A set of ChIP-Seq peaks where AR and ER were both detected
#'
#' @details
#' The subset of peaks found on chr1, and which contained signal from AR and ER,
#' along with H3K27ac signal were taken from GSE123767.
#' Peaks were resized to a uniform width of 400bp after downloading
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123767}
#'
#' @usage data("ar_er_peaks")
#'
#' @examples
#' data("ar_er_peaks")
#' ar_er_peaks
#' @name ar_er_peaks
#' @rdname ar_er_peaks
"ar_er_peaks"

#' @title Regions from hg19 with high N content
#'
#' @description
#' A GRanges object with regions annotated as telomeres or centromeres
#'
#' @details
#' The regions defined as centromeres or telomeres in hg19, taken from
#' AnnotationHub objects "AH107360" and "AH107361". These were combined with
#' regions containing Ns from the UCSC 2bit file, and regions with Ns in the
#' BSgenome.Hsapiens.UCSC.hg19 were retained.
#'
#' @source The package AnnotationHub and
#' \url{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.masked.gz}
#'
#' @usage data("hg19_mask")
#' @examples
#' data("hg19_mask")
#' hg19_mask
#' @name hg19_mask
#' @rdname hg19_mask
"hg19_mask"
