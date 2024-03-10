#' Script for creating example peaks.
#'
#' These represent peaks where ER, AR are found along with H3K27ac under DHT treatment


library(rtracklayer)
library(extraChIPs)

sq <- defineSeqinfo(build = "GRCh37")
genome(sq) <- "hg19"

## ER
er_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3511nnn/GSM3511085/suppl/GSM3511085%5FER%5Fpeaks%5FED.bed.gz"
er <- import.bed(er_url)
seqlevels(er) <- seqnames(sq)
seqinfo(er) <- sq
er <- sort(er)
er

## AR
ar_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3511nnn/GSM3511083/suppl/GSM3511083%5FAR%5Fpeaks%5FED.bed.gz"
ar <- import.bed(ar_url)
seqlevels(ar) <- seqnames(sq)
seqinfo(ar) <- sq
ar <- sort(ar)
ar

## H3K27ac
h3_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3511nnn/GSM3511087/suppl/GSM3511087%5FH3K27ac%5Fpeaks%5FED.bed.gz"
h3 <- import.bed(h3_url)
seqlevels(h3) <- seqnames(sq)
seqinfo(h3) <- sq
h3 <- sort(h3)
h3

ar_er_peaks <- GRangesList(
  AR = ar, ER = er, H3K27ac = h3
) |>
  makeConsensus(p = 2/3, method = "coverage", min_width = 200) |>
  subset(seqnames == "chr1") |>
  subset(n == 3) |>
  resize(width = 400, fix = 'center') |>
  granges()

save(ar_er_peaks, file = "data/ar_er_peaks.RData")

library(BSgenome.Hsapiens.UCSC.hg19)
ar_er_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ar_er_peaks)
save(ar_er_seq, file = "data/ar_er_seq.RData")
