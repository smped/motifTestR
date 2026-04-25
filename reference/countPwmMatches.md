# Count the matches to a PWM within an XStringSet

Count the matches to a PWM within an XStringSet

## Usage

``` r
countPwmMatches(
  pwm,
  stringset,
  rc = TRUE,
  min_score = "50%",
  mc.cores = 1,
  ...
)
```

## Arguments

- pwm:

  A Position Weight Matrix

- stringset:

  An XStringSet

- rc:

  logical(1) Also find matches using the reverse complement of pwm

- min_score:

  The minimum score to return a match

- mc.cores:

  Passed to [mclapply](https://rdrr.io/r/parallel/mclapply.html) when
  analysing a list of PWMs

- ...:

  Passed to [countPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html)

## Value

An integer vector

## Details

Will simply count the matches within an XStringSet and return an
integer. All matches are included.

## Examples

``` r
## Load the example PWM
data("ex_pfm")
esr1 <- ex_pfm$ESR1

## Load the example Peaks
data("ar_er_seq")
countPwmMatches(esr1, ar_er_seq)
#> [1] 199

## Count all PWMs
countPwmMatches(ex_pfm, ar_er_seq)
#>  ESR1  ANDR FOXA1 ZN143 ZN281 
#>   199   290  1041    76   213 
```
