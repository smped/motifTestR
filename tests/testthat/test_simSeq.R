pfm <- matrix(0, ncol = 4, nrow = 4)
diag(pfm) <- 1
rownames(pfm) <- c("A", "C", "G", "T")
um <- universalmotif::create_motif(pfm, type = "PWM", pseudocount = 0.01)

test_that("Basic sequence simulation works", {
    sim_seq <- simSeq(10, 4, pfm)
    expect_true(is(sim_seq, "DNAStringSet"))
    expect_equal(mcols(sim_seq)$pos, rep_len(1, 10))
    expect_equal(unique(as.character(sim_seq)), paste(rownames(pfm), collapse = ""))

    sim_seq <- simSeq(10, 4, um, as = "character")
    expect_true(mean(sim_seq ==  paste(rownames(pfm), collapse = "")) > 0.5)
    expect_true(is(sim_seq, "character"))
})

test_that("Poisson simulation works", {
    set.seed(1000)
    sim_seq <- simSeq(10, 5, pfm, rate = 0.5)
    expect_true(is(mcols(sim_seq)$pos, "IntegerList"))
    expect_true(is(mcols(sim_seq)$n_motifs, "integer"))
    ## For a width of 5, can only be pos 1 or 2
    expect_true(max(unlist(mcols(sim_seq)$pos)) <= 2)

})

test_that("NegBin simulation works", {
    set.seed(1000)
    sim_seq <- simSeq(10, 5, pfm, rate = 1, theta = 1)
    expect_true(is(mcols(sim_seq)$pos, "IntegerList"))
})

test_that("simSeq errors where expected", {
    expect_error(simSeq(1, 3, pfm))
    expect_error(simSeq(1, 4, pfm, rate = -1))
    expect_error(simSeq(1, 4, pfm, rate = 1, theta = -1))
})
