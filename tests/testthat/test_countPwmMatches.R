test_that("countPwmMatches behaves as expected",{
    expect_true(countPwmMatches(esr1, stringset) == 7)
    expect_true(countPwmMatches(esr1, stringset, rc = FALSE) == 5)
})

test_that("countPwmMatches accepts a list", {
    n <- countPwmMatches(ex_pwm, stringset)
    expect_true(length(n) == 5)
    expect_equal(names(n), names(ex_pwm))
})
