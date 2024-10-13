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
    res <- testMotifEnrich(um_db, test_set, bg_set, model = "poisson", sort_by = "none")
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

test_that("Hypergeometric works", {
    mcols(bg_set) <- mcols(bg_ranges)
    hg <- testMotifEnrich(esr1, test_set, bg_set, model = "hyper")
    expect_true(is(hg, "data.frame"))
    expect_true(nrow(hg) == 1)
    expect_true("odds_ratio" %in% colnames(hg))
    expect_message(
        testMotifEnrich(
            esr1, test_set, c(bg_set, test_set[1]), model = "hyper"
        ),
        "Shared.+"
    )

})

test_that("Errors are caught",{
    expect_error(.testIter(list(esr1), test_set, bg_set[-1]))
})
