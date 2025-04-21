test_that(".checkPfmList functions as expected", {

    mat_list <- .checkPfmList(um_db)
    expect_true(all(vapply(mat_list, is.matrix, logical(1))))
    nm <- vapply(um_db, slot, character(1), "name")
    expect_equal(names(mat_list), nm)

    mat_list <- .checkPfmList(ex_pfm)
    expect_equal(names(mat_list), names(ex_pfm))

    ## Errors
    expect_error(.checkPfmList(um_db[[1]]))
    expect_error(.checkPfmList(um_db[1]), "Please use simSeq.+")
})

test_that(".defineBG functions as expected", {

    bg <- .defineBG(stringset, nt = c("A", "C", "G", "T"))
    expect_true(is.character(bg))

    bg <- .defineBG(NULL, 2, 5, c("A", "C"), c(0.5, 0.5))
    expect_true(is.character(bg))
    expect_true(length(bg) == 2)
    expect_true(all(vapply(bg, nchar, integer(1)) == 5))

    err_bg <- as(list(stringset[[1]], stringset[[1]][1:5]), "DNAStringSet")
    expect_error(.defineBG(err_bg, nt = c("A", "C", "G", "T")))

    expect_error(.defineBG(1))
    expect_error(defineBG("cat", 1, 5, nt = c("A", "C"), c(0.5, 0.5)))
})

test_that("Main function simulates data as expected", {

    pfm_list <- ex_pfm[1:2]
    ## Check character output
    seq <- simMultiMotifs(2, 100, pfm_list, as = "character")
    expect_true(is.character(seq))
    expect_true(sum(nchar(seq)) == 2 * 100)

    ## Check DNAStringSet output
    seq <- simMultiMotifs(2, 100, pfm_list, ol = "last")
    expect_true(is(seq, "DNAStringSet"))
    expect_true(sum(width(seq)) == 2 * 100)
    expect_equal(colnames(mcols(seq)), c(names(pfm_list), "n_motifs"))
    ## These should be simulated with 1 motif/sequence
    expect_true(all(vapply(mcols(seq), is.numeric, logical(1))))
    ## Re-simulate using the rate & more sequences
    set.seed(10)
    seq <- simMultiMotifs(10, 100, pfm_list, rate = c(2, 0.5))
    mc <- mcols(seq)[names(pfm_list)]
    expect_true(all(vapply(mc, \(x) is(x, "IntegerList"), logical(1))))
    counts <- lapply(mc, \(x) vapply(x, length, integer(1)))
    expect_true(mean(counts[[1]]) > mean(counts[[2]]))


})
