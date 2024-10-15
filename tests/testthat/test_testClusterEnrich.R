## Convert ranges to DNAStringSets
test_set <- ar_er_seq
bg_set <- getSeq(genome, bg_ranges)
mcols(bg_set) <- mcols(bg_ranges)

cl <- list(A = ex_pfm[1], B = ex_pfm[2:3], C = ex_pfm[4:5])

test_that("Basic Poisson analysis works", {
    res <- testClusterEnrich(cl, test_set, bg_set, model = "poisson")
    expect_true(is(res, "data.frame"))
    expect_true(nrow(res) == 3)
    expect_true("est_bg_rate" %in% colnames(res))
})


test_that("Iteration works", {
    mcols(bg_set) <- mcols(bg_ranges)
    iter <- testClusterEnrich(cl, test_set, bg_set, model = "iter")
    expect_true(is(iter, "data.frame"))
    expect_true(nrow(iter) == 3)
    expect_true("iter_p" %in% colnames(iter))
    expect_true(all(iter$n_iter == 10))

})

test_that("Quasipoisson works", {
    mcols(bg_set) <- mcols(bg_ranges)
    iter <- testClusterEnrich(cl, test_set, bg_set, model = "quasi")
    expect_true(is(iter, "data.frame"))
    expect_true(nrow(iter) == 3)
    expect_true(all(iter$n_iter == 10))

})

test_that("Hypergeometric works", {
    mcols(bg_set) <- mcols(bg_ranges)
    hg <- testClusterEnrich(cl, test_set, bg_set, model = "hyper")
    expect_true(is(hg, "data.frame"))
    expect_true(nrow(hg) == 3)

})

test_that("clusterMotifPos also works as expected", {
    res <- testClusterPos(cl, test_set)
    expect_equal(nrow(res), length(cl))
    expect_equal(nrow(testClusterPos(cl, test_set[0])), 0)
    matches <- getClusterMatches(cl, ar_er_seq, best_only = TRUE)
    res <- testClusterPos(matches, abs = TRUE)
    expect_equal(nrow(res), length(cl))
    expect_true(all(res$start >= 0))

    expect_equal(nrow(testClusterPos(matches[[1]])), 1L)
    expect_error(testClusterPos(ex_pfm, ar_er_seq), "All clusters should.+")
    expect_error(testClusterPos(cl[1], ar_er_seq), "All clusters cannot.+")
})
