#' @title Example Position Weight Matrices
#'
#' @description Example Position Weight Matrices
#'
#' @details
#' This object contains 5 PWMs taken from HOCOMOCOv11-coreA for examples and testing
#'
#' Generation of this motif list is documented in
#' `system.file("scripts/ex_pwm.R", package = "motifTestR")`
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
#' Generation of these ranges is documented in
#' `system.file("scripts/ar_er_peaks.R", package = "motifTestR")`
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

#' @title Sequences from peaks with AR and ER detected
#'
#' @description
#' The genomic sequences obtained from the ar_er_peaks
#'
#' @details
#' These sequences represent the sequences obtained from
#' BSgenome.Hsapiens.UCSC.hg19 for thw peaks supplied as `ar_er_peaks`
#'
#' Generation of these sequences is documented in
#' `system.file("scripts/ar_er_peaks.R", package = "motifTestR")`
#'
#' @usage data("ar_er_seq")
#'
#' @examples
#' data("ar_er_seq")
#' ar_er_seq
#' @name ar_er_seq
#' @rdname ar_er_seq
"ar_er_seq"

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
#' Generation of these ranges is documented in
#' `system.file("scripts/hg19_mask.R", package = "motifTestR")`
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

#' @title Candidate Enhancer Regions from ZR-75-1 Cells
#'
#' @description
#' The chr1 subset of candidate enhancers for ZR-75-1 cells
#'
#' @details
#' These enhancers are the chr1 subset of enhancer regions for ZR-75-1 cells as
#' identified by EnhancerAtlas 2.0
#'
#' #' Generation of these ranges is documented in
#' `system.file("scripts/zr75_enh.R", package = "motifTestR")`
#'
#'
#' @source \url{http://www.enhanceratlas.org/index.php}
#'
#' @usage data("zr75_enh")
#' @examples
#' data("zr75_enh")
#' zr75_enh
#' @name zr75_enh
#' @rdname zr75_enh
"zr75_enh"
