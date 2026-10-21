exp_cols <- c(
    "seq", "score", "direction", "start", "end", "from_centre", "seq_width",
    "match", "pvalue", "qvalue"
)

test_that("getPwmMatches works as expected",{
    all <- getPwmMatches(esr1, stringset, min_score = "80%", threshold_type = "logodds")
    fwd_only <- getPwmMatches(
        esr1, stringset, rc = FALSE, min_score = "80%",
        threshold_type = "logodds"
    ) |> suppressMessages()
    expect_true(is(all, "DataFrame"))
    expect_true(nrow(all) == 1)
    expect_true(nrow(fwd_only) == 0)
    expect_equal(colnames(all), colnames(fwd_only))
    expect_equal(colnames(all), exp_cols)
})

test_that("No matches are returned correctly", {
    empty <- getPwmMatches(
        esr1, stringset, rc = FALSE, min_score = "80%",
        threshold_type = "logodds"
    ) |> suppressMessages()
    expect_true(nrow(empty) == 0)
    expect_equal(colnames(empty), exp_cols)
})

test_that("Choosing only the best matches works", {
    ## Test the default first
    stringset <- ar_er_seq[2]
    best <- getPwmMatches(esr1, stringset, best_only = TRUE, threshold_type = "logodds")
    expect_true(nrow(best) == 1)
    expect_true(best$end == 99)
    ## Now modify & check other options in the private function
    all <- getPwmMatches(
        esr1, stringset, min_score = "80%", threshold_type = "logodds"
    )
    all$score <- 5 # Bodge up a tie
    expect_true(.getBestMatch(all, "first")$start == 107)
    expect_true(.getBestMatch(all, "central")$start == 194)
    expect_true(.getBestMatch(all, "last")$start == 194)
    expect_equal(nrow(.getBestMatch(all, "all")), nrow(all))

})

test_that("Passing a list also produces a list",{
    multi <- getPwmMatches(ex_pfm[1:2], stringset)
    expect_true(is(multi, "SplitDataFrameList"))
    expect_true(length(multi) == 2)
})

test_that("Fails are correct",{
    expect_error(getPwmMatches(letters, stringset))
    expect_error(getPwmMatches(esr1, stringset[[1]]))
    expect_true(nrow(.getBestMatch(data.frame(start = integer()))) == 0)
})
