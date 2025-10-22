# Candidate Enhancer Regions from ZR-75-1 Cells

The chr1 subset of candidate enhancers for ZR-75-1 cells

## Usage

``` r
data("zr75_enh")
```

## Format

An object of class `GRanges` of length 5237.

## Source

<http://www.enhanceratlas.org/index.php>

## Details

These enhancers are the chr1 subset of enhancer regions for ZR-75-1
cells as identified by EnhancerAtlas 2.0

\#' Generation of these ranges is documented in
`system.file("scripts/zr75_enh.R", package = "motifTestR")`

## Examples

``` r
data("zr75_enh")
zr75_enh
#> GRanges object with 5237 ranges and 0 metadata columns:
#>          seqnames              ranges strand
#>             <Rle>           <IRanges>  <Rle>
#>      [1]     chr1         28481-29320      *
#>      [2]     chr1       234931-236780      *
#>      [3]     chr1       440761-443340      *
#>      [4]     chr1       459151-460990      *
#>      [5]     chr1       462821-464300      *
#>      ...      ...                 ...    ...
#>   [5233]     chr1 248800421-248800510      *
#>   [5234]     chr1 249132751-249134450      *
#>   [5235]     chr1 249151701-249152020      *
#>   [5236]     chr1 249166051-249169210      *
#>   [5237]     chr1 249239441-249240740      *
#>   -------
#>   seqinfo: 24 sequences from hg19 genome
```
