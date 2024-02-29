test_that(".cleanMotifList works as expected", {
  ml <- .cleanMotifList(um_db)
  expect_equal(names(ml), vapply(um_db, slot, character(1), "name"))
  expect_true(all(vapply(ml, is, logical(1), "matrix")))
})

test_that("inputBM checks work", {
  expect_error(.checkInputBM(data.frame()))
  expect_error(.checkInputBM(DataFrame()))
  bm <- getPwmMatches(esr1, stringset)
  expect_true(.checkInputBM(bm))
  expect_true(.checkInputBM(list(bm)))
  expect_error(.checkInputBM(bm[1:2]))
})
