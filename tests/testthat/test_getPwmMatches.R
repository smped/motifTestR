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
