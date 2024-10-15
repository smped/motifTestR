#' Script used to create example motifs
#' All taken from the coreA set of HOCOMOCOv11
library(MotifDb)
ex_pfm <- MotifDb |>
  subset(organism == "Hsapiens") |>
  query(andStrings = "HOCOMOCOv11-core-A", orStrings = c("ESR1", "ANDR", "FOXA1", "ZNF")) |>
  as.list()
names(ex_pfm) <- gsub(".+-core-A-([A-Z0-9]+)_HUMAN.+", "\\1", names(ex_pfm))

save(ex_pfm, file = "data/ex_pfm.RData")
