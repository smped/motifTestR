test_that("clusterMotifs works as expected", {
    cl <- clusterMotifs(ex_pwm)
    expect_true(is(cl, "integer"))
    expect_equal(names(cl), names(ex_pwm))
    expect_error(
        clusterMotifs(ex_pwm, type = "ICM", method = "ALLR"),
        "Cannot use ICM with ALLR or ALLR_LL"
    )
    expect_error(clusterMotifs(c("cat", "dog")))
})
