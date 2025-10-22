# Find matches from a PWM cluster within an XStringSet

Find matches from a PWM cluster within a set of sequences

## Usage

``` r
getClusterMatches(
  cl,
  stringset,
  rc = TRUE,
  min_score = "80%",
  best_only = FALSE,
  break_ties = c("all", "random", "first", "last", "central"),
  mc.cores = 1,
  ...
)

countClusterMatches(
  cl,
  stringset,
  rc = TRUE,
  min_score = "80%",
  mc.cores = 1,
  ...
)
```

## Arguments

- cl:

  A list of Position Weight Matrices, universalmotifs, with each element
  representing clusters of related matrices

- stringset:

  An XStringSet

- rc:

  logical(1) Also find matches using the reverse complement of PWMs in
  the cluster

- min_score:

  The minimum score to return a match

- best_only:

  logical(1) Only return the best match

- break_ties:

  Method for breaking ties when only returning the best match Ignored
  when all matches are returned (the default)

- mc.cores:

  Passed to [mclapply](https://rdrr.io/r/parallel/mclapply.html)

- ...:

  Passed to [matchPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html)

## Value

Output from getClusterMatches will be a list of DataFrames with columns:
`seq`, `score`, `direction`, `start`, `end`, `from_centre`, `seq_width`,
`motif` and `match`

The first three columns describe the sequence with matches, the score of
the match and whether the match was found using the forward or reverse
PWM. The columns `start`, `end` and `width` describe the where the match
was found in the sequence, whilst `from_centre` defines the distance
between the centre of the match and the centre of the sequence being
queried. The motif column denotes which individual motif was found to
match in this position, again noting that when matches overlap, only the
one with the highest relative score is returned. The final column
contains the matching fragment of the sequence as an `XStringSet`.

Output from countClusterMatches will be a simple integer vector the same
length as the number of clusters

## Details

This function extends
[getPwmMatches](https://smped.github.io/motifTestR/reference/getPwmMatches.md)
by returning a single set of results for set of clustered motifs. This
can help remove some of the redundancy in results returned for highly
similar PWMs, such as those in the GATA3 family.

Taking a set of sequences as an XStringSet, find all matches above the
supplied score (i.e. threshold) for a list of Position Weight Matrices
(PWMs), which have been clustered together as highly-related motifs. By
default, matches are performed using the PWMs as provided and the
reverse complement, however this can easily be disabled by setting
`rc = FALSE`.

The function relies heavily on
[matchPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html) and
[Views](https://rdrr.io/pkg/IRanges/man/Views-class.html) for speed.

Where overlapping matches are found for the PWMs within a cluster, only
a single match is returned. The motif with the highest relative score
(score / maxScore(PWM)) is selected.

When choosing to return the best match (`best_only = TRUE`), only the
match with the highest relative score is returned for each sequence.
Should there be tied scores, the best match can be chosen as either the
first, last, most central, all tied matches, or choosing one at random
(the default).

## Examples

``` r
# Load example PFMs
data("ex_pfm")
# Cluster using default settings
cl_ids <- clusterMotifs(ex_pfm)
ex_cl <- split(ex_pfm, cl_ids)
# Add optional names
names(ex_cl) <- vapply(ex_cl, \(x) paste(names(x), collapse = "/"), character(1))

# Load example sequences
data("ar_er_seq")
# Get all matches for each cluster
getClusterMatches(ex_cl, ar_er_seq)
#> $ESR1
#> DataFrame with 22 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          29   18.0880         F       193       207           0       400
#> 2          34   20.8412         R       321       335         128       400
#> 3          60   17.8088         R       154       168         -39       400
#> 4          62   17.5548         R       206       220          13       400
#> 5          98   20.2850         F        13        27        -180       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 18        478   18.9927         R       134       148         -59       400
#> 19        517   19.0738         F       223       237          30       400
#> 20        552   18.4739         F       232       246          39       400
#> 21        575   17.7611         R         4        18        -189       400
#> 22        646   17.5586         R       209       223          16       400
#>           motif           match
#>     <character>  <DNAStringSet>
#> 1          ESR1 AGGTCACCCTGGCCC
#> 2          ESR1 AGGTCACCGTGACCC
#> 3          ESR1 AGGTGACCCTGACCT
#> 4          ESR1 GGGTCACACTGTCCT
#> 5          ESR1 AGGTCACAATGACCT
#> ...         ...             ...
#> 18         ESR1 AGGTCACCCTGACCG
#> 19         ESR1 GGGTCAGCATGACCT
#> 20         ESR1 AGGACACACTGACCT
#> 21         ESR1 AGGTCACCCTAACCT
#> 22         ESR1 AGGTTAGCCTGACCT
#> 
#> $`ANDR/FOXA1`
#> DataFrame with 121 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          12   14.0694         F       199       210         4.5       400
#> 2          16   14.4605         F       177       188       -17.5       400
#> 3          17   15.1441         F       341       352       146.5       400
#> 4          18   14.4417         F       301       312       106.5       400
#> 5          21   14.6369         R       206       217        11.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 117       817   14.9185         R       291       302        96.5       400
#> 118       826   14.5836         F       261       272        66.5       400
#> 119       833   23.0102         R       167       184       -24.5       400
#> 120       844   14.3461         F        40        51      -154.5       400
#> 121       844   14.0222         F        72        83      -122.5       400
#>           motif              match
#>     <character>     <DNAStringSet>
#> 1         FOXA1       TGTTTGCTTTTG
#> 2         FOXA1       TGTTTACTTTCC
#> 3         FOXA1       TGTTTATTTAGG
#> 4         FOXA1       TGTTTATTCTGG
#> 5         FOXA1       TGTTTACTCAAC
#> ...         ...                ...
#> 117       FOXA1       TGTTTACACAGT
#> 118       FOXA1       TATTTACTTTAG
#> 119        ANDR TGTTCTTTTTTGTTTGTT
#> 120       FOXA1       TGTTTACTTTCT
#> 121       FOXA1       TGTTTGCTCTGC
#> 
#> $ZN143
#> DataFrame with 21 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          30   26.8427         F       166       187       -23.5       400
#> 2          30   26.8023         F       217       238        27.5       400
#> 3          67   29.0591         F       210       231        20.5       400
#> 4         118   29.0063         R       205       226        15.5       400
#> 5         182   28.0840         F       225       246        35.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 17        750   28.4138         R       206       227        16.5       400
#> 18        829   26.6710         R       151       172       -38.5       400
#> 19        836   28.9222         R       166       187       -23.5       400
#> 20        837   30.0534         F       216       237        26.5       400
#> 21        837   28.0840         F       352       373       162.5       400
#>           motif                  match
#>     <character>         <DNAStringSet>
#> 1         ZN143 AGCCTGCCGGGAGATGTAGTTC
#> 2         ZN143 GGCACGCCGGGAAATGTAGTTC
#> 3         ZN143 GGCATGCTGGGATTTGTAGTCT
#> 4         ZN143 TGCCTCCTGGGAAATGTAGTCC
#> 5         ZN143 TGCATGCTGGGAACTGTAGTCT
#> ...         ...                    ...
#> 17        ZN143 GGCATGCCGGGAGTTGTAGTCC
#> 18        ZN143 TGCCCGCTGGGAACTGTAGTCC
#> 19        ZN143 TGCATGCTGGGATTTGTAGTCC
#> 20        ZN143 TGCATGCTGGGAGTTGTAGTCT
#> 21        ZN143 TGCATGCTGGGAACTGTAGTCT
#> 
#> $ZN281
#> DataFrame with 14 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1         109   19.6553         R       369       383         176       400
#> 2         118   20.1871         F        60        74        -133       400
#> 3         122   19.3263         F        95       109         -98       400
#> 4         171   18.5467         R        84        98        -109       400
#> 5         192   19.1012         F       260       274          67       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 10        456   21.2625         R       343       357         150       400
#> 11        507   18.6243         R       175       189         -18       400
#> 12        668   20.0793         R       168       182         -25       400
#> 13        763   22.6040         R       274       288          81       400
#> 14        764   18.8379         R       310       324         117       400
#>           motif           match
#>     <character>  <DNAStringSet>
#> 1         ZN281 AGTTGGGGGAGGGGC
#> 2         ZN281 GGCGGGGGGAGGGGA
#> 3         ZN281 GAATGGGGGAGGGGC
#> 4         ZN281 GGATGGGGGAAGGGG
#> 5         ZN281 GGGAGGGGGCGGGGG
#> ...         ...             ...
#> 10        ZN281 CGGTGGGGGAGGGGG
#> 11        ZN281 GGGAGGGGGAGGGAG
#> 12        ZN281 GGGTGGGGGTGGGGG
#> 13        ZN281 GGGTGGGGGAGGGGG
#> 14        ZN281 AGTGGGGGGAGGGGA
#> 
# Or Just count them
countClusterMatches(ex_cl, ar_er_seq)
#>       ESR1 ANDR/FOXA1      ZN143      ZN281 
#>         22        121         21         14 
# Compare this to individual counts
countPwmMatches(ex_pfm, ar_er_seq)
#>  ESR1  ANDR FOXA1 ZN143 ZN281 
#>    22     8   113    21    14 
```
