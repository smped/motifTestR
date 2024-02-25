test_that("countPwmMatches behaves as expected",{
  expect_true(countPwmMatches(esr1, stringset) == 7)
  expect_true(countPwmMatches(esr1, stringset, rc = FALSE) == 5)
})
