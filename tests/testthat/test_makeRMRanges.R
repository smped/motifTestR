sq <- Seqinfo(seqnames = "chr1", seqlengths = 500, FALSE, "test")
rng <- GRanges(c("chr1:101-150", "chr1:201-250"), seqinfo = sq)
chr <- GRanges("chr1:51-500", seqinfo = sq)
grl <- splitAsList(rng, c("promoter", "other"))
y <- GRangesList(
  promoter = GRanges("chr1:1-200"), other = GRanges("chr1:201-500")
)
seqinfo(y) <- sq

test_that("Basic functions of makeRMRanges work", {
  bg <- makeRMRanges(rng, chr, n_iter = 2)
  expect_true(length(bg) == length(rng) * 2)
  expect_true("iteration" %in% colnames(mcols(bg)))
  expect_true(all(width(bg) == 50))

  bg <- makeRMRanges(rng, chr, n_total = 3)
  expect_true(length(bg) == 3)
  expect_true(length(colnames(mcols(bg))) == 0)

})

test_that("GRangesList behaves as expected", {

  bg <- makeRMRanges(grl, y, n_iter = 2)
  expect_true("iteration" %in% colnames(mcols(bg)))
  nm <- structure(
    c(other = 2L, promoter = 2L), dim = 2L,
    dimnames = structure(list(c("other", "promoter")), names = ""),
    class = "table"
  )
  expect_equal(table(names(bg)), nm)

  bg <- makeRMRanges(grl, y, n_total = c(2, 1))
  nm <- structure(
    c(other = 2L, promoter = 1L), dim = 2L,
    dimnames = structure(list(c("other", "promoter")), names = ""),
    class = "table"
  )
  expect_equal(table(names(bg)), nm)

  bg <- makeRMRanges(grl, y, unlist = FALSE)
  expect_equal(names(bg), c("other", "promoter"))
  expect_true(all(vapply(bg, length, integer(1)) == 1))

})

test_that("makeRMRanges errors as expected", {
  expect_error(makeRMRanges(rng, chr, "chr1:1"))
  expect_error(makeRMRanges(rng, GRanges("chr1:150-120"))) # BG too short
  expect_error(makeRMRanges(rng, GRanges("chr1:1-200"))) # Range 2 excluded
  expect_error(makeRMRanges(grl, y[1]))
})

test_that("splitting entire genomes works", {
  hg19 <- extraChIPs::defineSeqinfo("GRCh37")
  genome(hg19) <- "hg19"
  rng <- GRanges(c("chr1:101-150", "chr1:201-250"), seqinfo = hg19)
  rm_ranges <- makeRMRanges(rng, GRanges(hg19), n_total = 2)
  expect_true(length(rm_ranges) == 2)
  expect_equal(dim(mcols(rm_ranges)), c(2, 0))
})
