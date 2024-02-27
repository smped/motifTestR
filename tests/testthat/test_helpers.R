test_that(".cleanMotifList works as expected", {
  ml <- .cleanMotifList(um_db)
  expect_equal(names(ml), vapply(um_db, slot, character(1), "name"))
  expect_true(all(vapply(ml, is, logical(1), "matrix")))
})
