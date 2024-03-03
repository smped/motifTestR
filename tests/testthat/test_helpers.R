test_that(".cleanMotifList works as expected", {
  ml <- .cleanMotifList(um_db)
  expect_equal(names(ml), vapply(um_db, slot, character(1), "name"))
  expect_true(all(vapply(ml, is, logical(1), "matrix")))
})

test_that("inputBM checks work", {
  expect_error(.checkMatches(data.frame()))
  expect_error(.checkMatches(DataFrame()))
  matches <- getPwmMatches(esr1, stringset)
  expect_true(.checkMatches(matches))
  expect_true(.checkMatches(list(matches)))
  expect_error(.checkMatches(matches[1:2]))
})
