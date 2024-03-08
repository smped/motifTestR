#' Form a set of masked regions for hg19 using telomeres & centromeres
library(AnnotationHub)
ah <- AnnotationHub()
centro <- granges(ah[["AH107360"]])
centro$type <- "centromere"
telo <- ah[["AH107361"]]
telo <- trim(telo)
telo <- granges(telo)
telo$type <- "telomere"
hg19_mask <- sort(c(centro, telo))
names(hg19_mask) <- c()

#' Now load the N.bed regions from the 2bit file.
#' These were produced on an HPC by taking the hard-masked fasta file from UCSC
#' converting to 2bit and extracting the regions with Ns
#'
#' faToTwoBit hg19.fa.masked hg19.fa.masked.2bit
#' twoBitInfo -nBed hg19.fa.masked.2bit N.bed
#' gzip N.bed
library(rtracklayer)
masks <- import.bed("inst/scripts/N.bed.gz") # Hidden in .gitignore & .Rbuildignore
masks <- subset(masks, seqnames %in% seqlevels(hg19_mask))
masks <- sortSeqlevels(masks)
masks <- sort(masks)
seqlevels(masks) <- seqnames(seqinfo(hg19_mask))
seqinfo(masks) <- seqinfo(hg19_mask)
masks$type <- 'N-content'

#' Get all the sequences corresponding to these masks from the actual BSgenome
library(BSgenome.Hsapiens.UCSC.hg19)
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, masks)
n_probs <- letterFrequency(seq, letters = "N", as.prob = TRUE)
#" Only retain those with any Ns at all
hg19_mask <- c(hg19_mask, masks[n_probs > 0]) |>
  sort() |>
  reduce()
save(hg19_mask, file = "data/hg19_mask.RData")
