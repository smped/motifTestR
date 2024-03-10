set.seed(305)
data("hg19_mask")
bg_ranges <- makeRMRanges(ar_er_peaks, GRanges(sq)[1], exclude = hg19_mask, n_iter = 10)

## Convert ranges to DNAStringSets
test_set <- ar_er_seq
bg_set <- getSeq(genome, bg_ranges)
mcols(bg_set) <- mcols(bg_ranges)


test_that("Basic Poisson analysis works", {
    res <- testMotifEnrich(esr1, test_set, bg_set, model = "poisson")
    expect_true(is(res, "data.frame"))
    expect_true(nrow(res) == 1)
    expect_true("est_bg_rate" %in% colnames(res))
})

test_that("Lists are tested", {
    res <- testMotifEnrich(um_db, test_set, bg_set, model = "poisson")
    expect_true(is(res, "data.frame"))
    expect_equal(rownames(res), vapply(um_db, slot, character(1), "name"))
})

test_that("Iteration works", {
    mcols(bg_set) <- mcols(bg_ranges)
    iter <- testMotifEnrich(esr1, test_set, bg_set, model = "iter")
    expect_true(is(iter, "data.frame"))
    expect_true(nrow(iter) == 1)
    expect_true("iter_p" %in% colnames(iter))
    expect_true(iter$n_iter == 10)

})

test_that("Quasipoisson works", {
    mcols(bg_set) <- mcols(bg_ranges)
    iter <- testMotifEnrich(esr1, test_set, bg_set, model = "quasi")
    expect_true(is(iter, "data.frame"))
    expect_true(nrow(iter) == 1)
    expect_true(iter$n_iter == 10)

})

test_that("NegBinom works", {
    mcols(bg_set) <- mcols(bg_ranges)
    iter <- suppressWarnings(
        # 10 is too few iterations in reality
        testMotifEnrich(esr1, test_set, bg_set, model = "neg")
    )
    expect_true(is(iter, "data.frame"))
    expect_true(nrow(iter) == 1)
    expect_true(iter$n_iter == 10)

})

test_that("Errors are caught",{
    expect_error(.testIter(list(esr1), test_set, bg_set[-1]))
})
