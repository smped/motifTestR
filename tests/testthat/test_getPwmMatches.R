exp_cols <- c(
    "seq", "score", "direction", "start", "end", "from_centre", "seq_width", "match"
)

test_that("getPwmMatches works as expected",{
    all <- getPwmMatches(esr1, stringset)
    fwd_only <- getPwmMatches(esr1, stringset, rc = FALSE)
    expect_true(is(all, "DataFrame"))
    expect_true(nrow(all) == 1)
    expect_true(nrow(fwd_only) == 0)
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
    stringset <- ar_er_seq[17]
    best <- getPwmMatches(esr1, stringset, best_only = TRUE)
    expect_true(nrow(best) == 1)
    expect_true(best$start == 148)
    ## Now modify & check other options in the private function
    all <- getPwmMatches(esr1, stringset)
    all$score <- 5
    expect_true(.getBestMatch(all, "first")$start == 148)
    expect_true(.getBestMatch(all, "central")$start == 221)
    expect_true(.getBestMatch(all, "last")$start == 371)
    expect_equal(nrow(.getBestMatch(all, "all")), nrow(all))

})

test_that("Passing a list also produces a list",{
    multi <- getPwmMatches(ex_pwm[1:2], stringset)
    expect_true(is(multi, "list"))
    expect_true(length(multi) == 2)
    expect_true(all(vapply(multi, is, logical(1), "DataFrame")))
})

test_that("Fails are correct",{
    expect_error(getPwmMatches(letters, stringset))
    expect_error(getPwmMatches(esr1, stringset[[1]]))
    expect_true(nrow(.getBestMatch(data.frame(start = integer()))) == 0)
})
