test_that("clusterMotifs works as expected", {
    cl <- clusterMotifs(ex_pfm)
    expect_true(is(cl, "integer"))
    expect_equal(names(cl), names(ex_pfm))
    expect_error(
        clusterMotifs(ex_pfm, type = "ICM", method = "ALLR"),
        "Cannot use ICM with ALLR or ALLR_LL"
    )
    expect_error(clusterMotifs(c("cat", "dog")))
})
