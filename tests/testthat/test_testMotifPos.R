esr1 <- structure(
  ## Taken from HOCOMOCO v12
  c(
    0.226744, 0.26938, 0.275194, 0.228682, 0.275194, 0.275194, 0.313953,
    0.135659, 0.023256, 0.005814, 0.005814, 0.965116, 0.011628, 0.003876,
    0.965116, 0.01938, 0.903101, 0.087209, 0.003876, 0.005814, 0.027132,
    0.94186, 0.013566, 0.017442, 0.236434, 0.738372, 0.001938, 0.023256,
    0.079457, 0.273256, 0.017442, 0.629845, 0.290698, 0.114341, 0.337209,
    0.257752
    ),
  dim = c(4L, 9L),
  dimnames = list(
    c("A", "C", "G", "T"), c("N", "V", "T", "G", "A", "C", "C", "Y", "D")
  )
)
stringset <- DNAStringSet(
  c(seq1 = "ATCCGTCTCCGGTGACTCACAACACTCCACTGACTCATGAAATCACCCCAGGTGTCAGGACACCGCCAGGAAGAATCGCTCCTCCTCCCACACCTCCCAACTTCAGTTATTCACATTCATGATTTCTGCCTCACCCAACAGCAACTGTACTAATATTTACAGAGCGTTTATCTTTAAATCTGTGACCAGCCCAAACCGAGCCTTGGTGTCTCTTAAAAGAACAAAGGGAGCACAGTGTTCCTGAATGGCACTGTCACTCACCCACTGCCTGTCCACCAACCTCACCCCTGGCCTCTCCATCTCCCCAAGGAGCCCCAGGGTGCAGCCTCTGAGACTCCGGAGCACCTCAGCATGACCATAAAGGCCTGGCTTCCTTCACATCATTTCCCATGTTCCTGGG" )
)

exp_cols <- c(
  "start", "end", "centre", "width", "total_matches", "matches_in_region",
  "expected", "enrichment", "prop_total", "p", "consensus_motif"
)

## Bodgy up a set of matches
bm <- getPwmMatches(esr1, stringset)
bm$seq <- c(1, seq_len(nrow(bm) - 1))

test_that("testMotifPos returns symmetrical peaks as expected",{
  res <- testMotifPos(bm = bm)
  expect_true(nrow(res) == 1)
  expect_equal(colnames(res), exp_cols)
  expect_equal(
    res[,1:6],
    data.frame(
      start = -155, end = 155, centre = 0, width = 310,
      total_matches = 7, matches_in_region = 5
    )
  )
  expect_true(is(res$consensus_motif[[1]], "matrix"))
  expect_equal(dim(res$consensus_motif[[1]]), c(4, 9))
})

test_that("setting abs = TRUE behaves as expected", {
  bins <- motifTestR:::.makeBmBins(bm, binwidth = 10, abs = TRUE)
  expect_equal(bins$weight, c(0.5, 0.5, rep_len(1, 5)))
  exp_bins <- c(
    "[0,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", "(50,60]",
    "(60,70]", "(70,80]", "(80,90]", "(90,100]", "(100,110]", "(110,120]",
    "(120,130]", "(130,140]", "(140,150]", "(150,160]", "(160,170]",
    "(170,180]", "(180,190]", "(190,200]"
  )
  expect_equal(levels(bins$bin), exp_bins)

})

test_that("Passing to getPwmMatches works", {
  res <- testMotifPos(esr1, stringset, abs = TRUE, binwidth = 100)
  expect_true(nrow(res) == 1)
  expect_equal(
    res[,1:6],
    data.frame(
      start = 0, end = 100, centre = 50, width = 100,
      total_matches = 1, matches_in_region = 1
    )
  )
})

test_that("testMotifPos errors correctly", {
  empty <- testMotifPos(bm = bm[0,])
  expect_true(nrow(empty) == 0)
  expect_true(is(empty, "data.frame"))
  expect_equal(colnames(empty), exp_cols)
  expect_error(testMotifPos(esr1, stringset, binwidth = 500))
})
