# Example Position Frequency Matrices

Example Position Frequency Matrices

## Usage

``` r
data("ex_pfm")
```

## Format

An object of class `list` of length 5.

## Details

This object contains 5 PFMs taken from HOCOMOCOv11-coreA for examples
and testing

Generation of this motif list is documented in
`system.file("scripts/ex_pfm.R", package = "motifTestR")`

## Examples

``` r
data("ex_pfm")
ex_pfm$ESR1
#>       1     2     3     4     5     6     7     8     9    10    11    12    13
#> A 0.638 0.074 0.046 0.094 0.002 0.856 0.108 0.396 0.182 0.104 0.054 0.618 0.040
#> C 0.048 0.006 0.018 0.072 0.888 0.006 0.442 0.604 0.376 0.078 0.034 0.198 0.884
#> G 0.260 0.808 0.908 0.178 0.048 0.112 0.312 0.000 0.286 0.044 0.908 0.070 0.014
#> T 0.054 0.112 0.028 0.656 0.062 0.026 0.138 0.000 0.156 0.774 0.004 0.114 0.062
#>      14    15
#> A 0.090 0.058
#> C 0.822 0.330
#> G 0.008 0.066
#> T 0.080 0.546
```
