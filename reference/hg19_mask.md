# Regions from hg19 with high N content

A GRanges object with regions annotated as telomeres or centromeres

## Usage

``` r
data("hg19_mask")
```

## Format

An object of class `GRanges` of length 345.

## Source

The package AnnotationHub and
<https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.masked.gz>

## Details

The regions defined as centromeres or telomeres in hg19, taken from
AnnotationHub objects "AH107360" and "AH107361". These were combined
with regions containing Ns from the UCSC 2bit file, and regions with Ns
in the BSgenome.Hsapiens.UCSC.hg19 were retained.

Generation of these ranges is documented in
`system.file("scripts/hg19_mask.R", package = "motifTestR")`

## Examples

``` r
data("hg19_mask")
hg19_mask
#> GRanges object with 345 ranges and 0 metadata columns:
#>         seqnames            ranges strand
#>            <Rle>         <IRanges>  <Rle>
#>     [1]     chr1           1-11447      *
#>     [2]     chr1     177418-227417      *
#>     [3]     chr1     267720-317719      *
#>     [4]     chr1     471369-521368      *
#>     [5]     chr1   2634221-2684220      *
#>     ...      ...               ...    ...
#>   [341]     chrY 22369680-22419679      *
#>   [342]     chrY 23901429-23951428      *
#>   [343]     chrY 28819362-58819361      *
#>   [344]     chrY 58917657-58967656      *
#>   [345]     chrY 59363566-59373566      *
#>   -------
#>   seqinfo: 24 sequences from hg19 genome
```
