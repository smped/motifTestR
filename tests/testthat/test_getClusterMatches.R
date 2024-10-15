cl <- list(A = ex_pfm[1], B = ex_pfm[2:3], C = ex_pfm[4:5])

test_that("countClusterMatches returnes expexted values", {
    counts <- countClusterMatches(cl, ar_er_seq)
    expect_true(is(counts, "integer"))
    expect_equal(names(counts), c("A", "B", "C"))
    expect_null(countClusterMatches(cl[0], ar_er_seq))

    counts <- countClusterMatches(ex_pfm, ar_er_seq)
    expect_equal(length(counts), 1L)

    expect_message(
        countClusterMatches(list("cat"), ar_er_seq),
        "Could not determine clusters"
    )

})

test_that("getClusterMatches returns expected values", {
    res <- getClusterMatches(cl, ar_er_seq)
    expect_true(is(res, "list"))
    expect_true(all(vapply(res, is, integer(1), "DataFrame")))
    expect_equal(names(res), c("A", "B", "C"))
    expect_equal(getClusterMatches(cl[0], ar_er_seq), list())

    res <- getClusterMatches(ex_pfm, ar_er_seq)
    expect_true(is(res, "DataFrame"))
    nm <- c(
        "seq", "score", "direction", "start", "end", "from_centre",
        "seq_width", "motif", "match"
    )
    expect_equal(colnames(res), nm)

    names(ar_er_seq) <- paste("Seq", seq_along(ar_er_seq))
    best <- getClusterMatches(cl, ar_er_seq, best_only = TRUE)
    expect_true(is(best, "list"))
    expect_true(all(vapply(best, is, integer(1), "DataFrame")))
    expect_true(is(best$A$seq, "character"))

    expect_message(
        getClusterMatches(list("cat"), ar_er_seq),
        "Could not determine clusters"
    )

    empty <- getClusterMatches(ex_pfm, ar_er_seq[0])
    expect_equal(colnames(empty), nm)
    expect_equal(nrow(empty), 0)

})
