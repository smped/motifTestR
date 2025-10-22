# Simulate sequences with multiple motifs

Simulate a set of sequences incorporating multiple motifs

## Usage

``` r
simMultiMotifs(
  n,
  width,
  pfm = NULL,
  bg = NULL,
  nt = c("A", "C", "G", "T"),
  prob = rep(0.25, 4),
  shape1 = 1,
  shape2 = shape1,
  rate = NA,
  theta = NA,
  as = "DNAStringSet",
  ol = c("random", "first", "last"),
  ...
)
```

## Arguments

- n:

  The number of sequences to simulate

- width:

  Width of sequences to simulate

- pfm:

  List of Probability Weight/Frequency Matrices

- bg:

  Optional, pre-defined set of background sequences. Can be passed as an
  XStringSet or character vector. All sequences must be the same width

- nt:

  Nucleotides to include

- prob:

  Sampling probabilities for each nucleotide

- shape1, shape2:

  Passed to
  [rbetabinom.ab](https://rdrr.io/pkg/VGAM/man/betabinomUC.html)

- rate:

  The expected rate of motifs per sequence. Is equivalent to \\ \lambda
  \\ in [rpois](https://rdrr.io/r/stats/Poisson.html). If set to NULL or
  NA, all sequences will be simulated with a single motif, otherwise a
  Poisson distribution will be used

- theta:

  Overdispersion parameter passed to
  [rnegbin](https://rdrr.io/pkg/MASS/man/rnegbin.html). If set to NULL
  or NA the rate parameter will be passed to
  [rpois](https://rdrr.io/r/stats/Poisson.html). However if this value
  is set, the rate and theta parameters are passed to
  [rnegbin](https://rdrr.io/pkg/MASS/man/rnegbin.html) to simulate
  overdispersed counts

- as:

  ObjectClass to return objects as. Defaults to DNAStringSet, but other
  viable options may include 'character', 'CharacterList' or any other
  class from which a character vector may be coerced.

- ol:

  When randomly simulated positions overlap, choose one either at
  random, by the first occurring PFM in the list of PFMs, or by the
  last.

- ...:

  Not used

## Value

A DNAStringSet with mcols denoting the positions of all inserted motifs

## Details

Simulate a set of sequences with multiple motifs inserted using
different rates and distributions, as specified. All shape, rate and
theta parameters are recycled to match the length of the supplied motif
list, and can be supplied as vectors to tailor these parameters to each
provided element of the list of matrices

## Examples

``` r
data("ex_pfm")
## Simulate sequences including both ESR1 and ANDR, but with
## ESR1 being included at a higher rate
seq <- simMultiMotifs(10, 100, ex_pfm[1:2], rate = c(2, 1))
seq
#> DNAStringSet object of length 10:
#>      width seq
#>  [1]   100 GTGTCCTGGTGCTATAGGTCGCCCTGATCTGCT...AGGACAGAGTGACACACGCGCATCCCATCGTTA
#>  [2]   100 TCCAGTAGATGTTTACTGTTAAGTTTGTTGTAT...TCGAGATTGGCTTCGCTCAATCGTAATTAATAC
#>  [3]   100 ACTGTGCTAGTTAGATCCATTAGTCGTTTTAGC...TCTATGACATTTATTCTAGACTGACCCTTTTGG
#>  [4]   100 ACTTAGAAGCGGGTCTCCCTTTGGTTGTTTGCT...TAGTTAGCTCTATGACGCAGCGACGGCCCCACA
#>  [5]   100 CGAGTCATGCGCAGCATCGCGGAACACAGACTA...CCTGACCCACGGTTTTCGCTGGGAAGACCTGAG
#>  [6]   100 CAGGCCGCTGGCCGCGTGCTCTTTTCTAGGGCA...TCTTCGGTATGACTGTGACTGGGATTTCAGCAC
#>  [7]   100 CCATATCGCGGCGCATCTTGACAGACAGTCATA...TGGTTACACTGACCCGGGCCATATTGATACTGT
#>  [8]   100 GGGTCACATAGACCTCATTTGTTTGCATGTCCG...GTAACGTGTTAAGGGGGTTACGATTTCCTTTGG
#>  [9]   100 TACGATGCATGGGCGAAGATACGTATCGCCAGA...GCTGCGTAGCAATGTACCGTATAATTTGCAAAG
#> [10]   100 GGCCCAATGACGCGAGGTACCGCACGCACGGCC...GTCAAAGACCTAAGCGTAGTACAGAAACAAGTT
## The positions of the motifs are included in the mcols
mcols(seq)
#> DataFrame with 10 rows and 3 columns
#>             ESR1          ANDR  n_motifs
#>    <IntegerList> <IntegerList> <numeric>
#> 1       16,54,68                       3
#> 2             37            12         2
#> 3          63,80            65         3
#> 4       11,62,64         16,59         5
#> 5          42,61                       2
#> 6        2,28,47             3         4
#> 7       63,68,83                       3
#> 8              1            10         2
#> 9             38            80         2
#> 10            48         29,35         3

```
