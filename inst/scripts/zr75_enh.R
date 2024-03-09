library(rtracklayer)
library(extraChIPs)
sq <- defineSeqinfo("GRCh37", chr = TRUE)
genome(sq) <- "hg19"
enh <- import.bed("http://www.enhanceratlas.org/data/download/enhancer/hs/ZR75-1.bed")
enh <- granges(enh)
seqlevels(enh) <- seqnames(sq)
seqinfo(enh) <- sq
zr75_enh <- subset(enh, seqnames == "chr1")

save(zr75_enh, file = "data/zr75_enh.RData")
