# A set of peaks with AR and ER detected

A set of ChIP-Seq peaks where AR and ER were both detected

## Usage

``` r
data("ar_er_peaks")
```

## Format

An object of class `GRanges` of length 849.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123767>

## Details

The subset of peaks found on chr1 which contained signal from at least
two of AR, ER and H3K27ac, taken from GSE123767. Peaks were resized to a
uniform width of 400bp after downloading

Generation of these ranges is documented in
`system.file("scripts/ar_er_peaks.R", package = "motifTestR")`

## Examples

``` r
data("ar_er_peaks")
ar_er_peaks
#> GRanges object with 849 ranges and 0 metadata columns:
#>         seqnames              ranges strand
#>            <Rle>           <IRanges>  <Rle>
#>     [1]     chr1     1008982-1009381      *
#>     [2]     chr1     1014775-1015174      *
#>     [3]     chr1     1051296-1051695      *
#>     [4]     chr1     1299561-1299960      *
#>     [5]     chr1     2179886-2180285      *
#>     ...      ...                 ...    ...
#>   [845]     chr1 246771887-246772286      *
#>   [846]     chr1 246868678-246869077      *
#>   [847]     chr1 246873126-246873525      *
#>   [848]     chr1 247095351-247095750      *
#>   [849]     chr1 247267507-247267906      *
#>   -------
#>   seqinfo: 24 sequences from hg19 genome
```
