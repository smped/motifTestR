#' Script used to create example motifs
#' All taken from the coreA set of HOCOMOCOv11
library(MotifDb)
ex_pwm <- MotifDb |>
  subset(organism == "Hsapiens") |>
  query(andStrings = "HOCOMOCOv11-core-A", orStrings = c("ESR1", "ANDR", "FOXA1", "ZNF")) |>
  as.list()
names(ex_pwm) <- gsub(".+-core-A-([A-Z0-9]+)_HUMAN.+", "\\1", names(ex_pwm))
