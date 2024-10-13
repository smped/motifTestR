## Try clustering motifs
## The idea is to test for cluster enrichment/bias instead of individual
## motifs.
## A function could be written to cluster motifs and then all motifs in the
## cluster are tested as a single unit, allowing for matches to any cluster
## member.
## Obviously, existing methods in TFBSTools and PWMEnrich need checking
library(universalmotif)
library(tidyverse)

motif_db <- read_meme(
    "~/github/motifTestR_paper/data/motif_db/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
) %>%
    to_df() %>%
    as_tibble()

comp <- motif_db %>%
    to_list() %>%
    compare_motifs(method = "EUCL")

d <- motif_db %>%
    dplyr::filter(str_detect(name, "GATA|FOX|E[SR]R|ANDR|PRGR|ZNF|AP2|STAT")) %>%
    to_list() %>%
    compare_motifs(method = "WPCC")
d[d < 0.8] <- 0
(1 -d) %>%
    as.dist() %>%
    hclust() %>%
    plot(cex = 0.7)
abline(a = 0.2, b = 0)
