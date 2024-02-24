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
  "seq", "score", "direction", "start", "end", "width", "from_centre", "match"
)

test_that("getPwmMatches works as expected",{
  all <- getPwmMatches(esr1, stringset)
  fwd_only <- getPwmMatches(esr1, stringset, rc = FALSE)
  expect_true(is(all, "DataFrame"))
  expect_true(nrow(all) == 7)
  expect_true(nrow(fwd_only) == 5)
  expect_equal(colnames(all), colnames(fwd_only))
  expect_equal(colnames(all), exp_cols)
})

test_that("No matches are returned correctly", {
  no_hit <- DNAStringSet(c(no_hit = stringset[[1]][41:100]))
  empty <- getPwmMatches(esr1, no_hit, rc = FALSE)
  expect_true(nrow(empty) == 0)
  expect_equal(colnames(empty), exp_cols)
})
