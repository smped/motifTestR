test_that("clusterMotifs works as expected", {
    cl <- clusterMotifs(ex_pfm, power = 2, plot = TRUE)
    expect_true(is(cl, "integer"))
    expect_equal(names(cl), names(ex_pfm))
    expect_error(
        clusterMotifs(ex_pfm, type = "ICM", method = "ALLR"),
        "Cannot use ICM with ALLR or ALLR_LL"
    )
    expect_error(clusterMotifs(c("cat", "dog")))
})

test_that("returning the distance matrix works", {
    cl <- clusterMotifs(ex_pfm, return_d = TRUE)
    expect_equal(c("cl", "d"), names(cl))
    expect_equal(max(cl$cl), length(cl$d))
})
