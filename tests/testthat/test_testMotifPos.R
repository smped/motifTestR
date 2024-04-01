exp_cols <- c(
    "start", "end", "centre", "width", "total_matches", "matches_in_region",
    "expected", "enrichment", "prop_total", "p", "fdr", "consensus_motif"
)

## Bodgy up a set of matches
matches <- getPwmMatches(esr1, stringset)
matches$seq <- c(1, seq_len(nrow(matches) - 1))

test_that("testMotifPos returns symmetrical peaks as expected",{
    res <- testMotifPos(matches)
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
    bins <- .makeBmBins(matches, binwidth = 10, abs = TRUE)
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

test_that("List input works", {
    res <- testMotifPos(ex_pwm, seq, sort_by = "none")
    expect_true(is(res, "data.frame"))
    expect_equal(rownames(res), names(ex_pwm))
    expect_true(all(exp_cols %in% colnames(res)))
    expect_true("fdr" %in% colnames(res))
})

test_that("testMotifPos errors correctly", {
    empty <- testMotifPos(matches[0,])
    expect_true(nrow(empty) == 0)
    expect_true(is(empty, "data.frame"))
    expect_equal(colnames(empty), exp_cols)
    expect_error(testMotifPos(esr1, stringset, binwidth = 500))
    expect_error(testMotifPos(""))
    expect_error(testMotifPos(c(list(""), ex_pwm)))
})

test_that("Same width sequences are enforced", {
    matches$seq_width[1] <- 100
    expect_error(testMotifPos(matches))
    mixed_width <- c(stringset, as(stringset[[1]][1:200], "DNAStringSet"))
    expect_error(testMotifPos(ex_pwm, mixed_width))
})
