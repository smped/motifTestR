# Simulate sequences using optional TFBMs

Simulate a set of fixed-width sequences using optional TFBMs

## Usage

``` r
simSeq(
  n,
  width,
  pfm = NULL,
  nt = c("A", "C", "G", "T"),
  prob = rep(0.25, 4),
  shape1 = 1,
  shape2 = shape1,
  rate = NULL,
  theta = NULL,
  as = "DNAStringSet",
  ...
)
```

## Arguments

- n:

  The number of sequences to simulate

- width:

  Width of sequences to simulate

- pfm:

  Probability Weight/Frequency Matrix

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
  or NA, the rate parameter will be passed to
  [rpois](https://rdrr.io/r/stats/Poisson.html). However if this value
  is set, the rate and theta parameters are passed to
  [rnegbin](https://rdrr.io/pkg/MASS/man/rnegbin.html) to simulate
  overdispersed counts

- as:

  ObjectClass to return objects as. Defaults to DNAStringSet, but other
  viable options may include 'character', 'CharacterList' or any other
  class from which a character vector may be coerced.

- ...:

  Not used

## Value

By default a DNAStringSet will be returned. If possible, the position of
any randomly sampled motifs will be included in the mcols element of the
returned object.

## Details

Using the nucleotide and probabilities provided as set of sequences can
be simulated. By default, this will effectively be a set of 'background'
sequences, with letters effectively chosen at random.

If a PWM/PFM is supplied, the shape parameters are first passed to
[rbetabinom.ab](https://rdrr.io/pkg/VGAM/man/betabinomUC.html) to
determine the random positions the motif will be placed, with the
default parameters representing a discrete uniform distribution.

The sequences to have a motif inserted will be selected, along with the
number of motifs, using the rate and theta parameters. If both are NULL,
every sequence will have a single motif inserted. If the rate is \> 0
and theta is NULL, sequences will be selected to have motifs inserted
using a poisson distribution. If theta is also provided, sequences will
be selected to contain motifs using a negative binomial distribution,
noting that smaller values of theta lead to higher over-dispersion

Once positions and sequences for the TFBM have been selected,
nucleotides will be randomly sampled using the probabilities provided in
the PWM and these motifs will be placed at the randomly sampled
positions

## Examples

``` r
## Randomly generate 10x50nt sequences without any TFBMs present
simSeq(10, 50)
#> DNAStringSet object of length 10:
#>      width seq
#>  [1]    50 GAAGGACAATAGGAATCGCAAGCCTTCGGGACGACGCTTCGAGTGGAAGG
#>  [2]    50 CGCCCGCAAGATGAGGAGTTAGGACGCAACGTTAATTGGCTACTTTGACT
#>  [3]    50 TGCTCTCTTAGGGAGATTTCTCTTTGGATTGGGCCCCGATCATCGTGTCC
#>  [4]    50 CCGAATGATAACATAATTGACTTAATGGTCCAATAGCTGACAAAGGGTCC
#>  [5]    50 ATCTCTCGAGCACCATCTCCGCGCATTTTCATACGCTCACCCAGAAATTA
#>  [6]    50 GAAGCCAATACGGTGTCGCGCAGGGGTGTTCAGTATCAACTGCTGTTGCT
#>  [7]    50 ATGCCATCTATGACCTACTGGAGATAGTGGATGAAAACAGCATTTTGAAG
#>  [8]    50 CACTGTGTTTTTTTCGTATTCCTATGCCCAGCATGCCTGTGAAGAGCCGT
#>  [9]    50 ACCAATTCCCATCATAAGATACGGTGGTATCGCAACTTCATGAACTACAG
#> [10]    50 CCTCGGTGTGCCTATCCGTCGGAATGTCCTCACTAGGTGATGTAGGTTGT

## Now place a motif at random positions
data('ex_pfm')
sim_seq <- simSeq(10, width = 20, pfm = ex_pfm$ESR1)
sim_seq
#> DNAStringSet object of length 10:
#>      width seq
#>  [1]    20 TATGTCAGCCTGACCCTTCA
#>  [2]    20 AGAGAGGTCACATGCACATA
#>  [3]    20 GGCAGAGGTCAAAATGGCCC
#>  [4]    20 AGGTCAAATTGCTAATACGG
#>  [5]    20 TTGATCAACCTGACCCTGCC
#>  [6]    20 AAAAGAGGTCAGACTGACCC
#>  [7]    20 GGGGCGCCCTGACCCGTAAC
#>  [8]    20 GTGGGGGCGCCGTGACTTTG
#>  [9]    20 TAAAGGTCAGCCCGAACTAC
#> [10]    20 ATGTTCGGTGAACATGACCC
## The position of the motif within each sequence is included in the mcols
mcols(sim_seq)
#> DataFrame with 10 rows and 2 columns
#>          pos  n_motifs
#>    <integer> <integer>
#> 1          2         1
#> 2          5         1
#> 3          6         1
#> 4          1         1
#> 5          2         1
#> 6          6         1
#> 7          1         1
#> 8          4         1
#> 9          4         1
#> 10         6         1

## Use this to extract the random motifs from the random sequences
library(IRanges)
i <- mcols(sim_seq)$pos + cumsum(width(sim_seq)) - width(sim_seq)
Views(unlist(sim_seq), start = i, width = 10)
#> Views on a 200-letter DNAString subject
#> subject: TATGTCAGCCTGACCCTTCAAGAGAGGTCACATG...TCAGCCCGAACTACATGTTCGGTGAACATGACCC
#> views:
#>        start end width
#>    [1]     2  11    10 [ATGTCAGCCT]
#>    [2]    25  34    10 [AGGTCACATG]
#>    [3]    46  55    10 [AGGTCAAAAT]
#>    [4]    61  70    10 [AGGTCAAATT]
#>    [5]    82  91    10 [TGATCAACCT]
#>    [6]   106 115    10 [AGGTCAGACT]
#>    [7]   121 130    10 [GGGGCGCCCT]
#>    [8]   144 153    10 [GGGGCGCCGT]
#>    [9]   164 173    10 [AGGTCAGCCC]
#>   [10]   186 195    10 [CGGTGAACAT]

```
