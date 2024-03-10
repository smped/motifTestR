bm <- getPwmMatches(ex_pwm[1:2], seq, best_only = TRUE, break_ties = "all")

test_that("Defaults work", {
    p <- plotMatchPos(bm, se = FALSE)
    expect_true(is(p$layers[[1]]$geom, "GeomSmooth"))
    expect_equal(levels(p$data$name), names(bm))
    expect_equal(p$labels, list(x = "bin_centre", y = "p", colour = "name"))
    expect_true(!p$layers[[1]]$geom_params$se)
})

test_that("cdf plots work", {
    p <- plotMatchPos(bm, abs = TRUE, type = "cdf", geom = "line")
    expect_true(is(p$layers[[1]]$geom, "GeomLine"))
    expect_equal(p$labels, list(x = "bin_start", y = "cdf", colour = "name"))
})

test_that("col plots work", {
    p <- plotMatchPos(bm, abs = TRUE, geom = "col", position = "dodge")
    expect_true(is(p$layers[[1]]$geom, "GeomCol"))
    expect_true(is(p$layers[[1]]$position, "PositionDodge"))
    expect_equal(p$labels, list(x = "bin_start", y = "p", fill = "name"))
})

test_that("heatmaps work", {
    p <- plotMatchPos(bm, type = "heat")
    expect_true(is(p$layers[[1]]$geom, "GeomTile"))
    expect_equal(p$labels, list(x = "bin", y = "name", fill = "p"))
})

test_that("single input works", {
    p <- plotMatchPos(bm[[1]])
    expect_true(is(p$layers[[1]]$geom, "GeomSmooth"))
    expect_equal(p$labels, list(x = "bin_centre", y = "p", colour = "colour"))
})

test_that("errors are as expected",{
    expect_error(plotMatchPos(bm[[1]], type = "heat"))
    expect_error(plotMatchPos(bm[c(1, 1)]))
})
