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
  "seq", "score", "direction", "start", "end", "from_centre", "seq_width", "match"
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

test_that("Choosing only the best matches works", {
  ## Test the default first
  best <- getPwmMatches(esr1, stringset, best_only = TRUE)
  expect_true(nrow(best) == 1)
  expect_true(best$start == 181)
  ## Now modify & check other options in the private function
  all <- getPwmMatches(esr1, stringset)
  all$score <- 5
  expect_true(.getBestMatch(all, "first")$start == 11)
  expect_true(.getBestMatch(all, "central")$start == 181)
  expect_true(.getBestMatch(all, "last")$start == 351)
  expect_equal(nrow(.getBestMatch(all, "all")), nrow(all))

})

test_that("Fails are correct",{
  expect_error(getPwmMatches(letters, stringset))
  expect_error(getPwmMatches(esr1 * 2, stringset))
  expect_error(getPwmMatches(esr1, stringset[[1]]))
  expect_true(nrow(.getBestMatch(data.frame(start = integer()))) == 0)
})
