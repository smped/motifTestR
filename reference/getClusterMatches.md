# Find matches from a PWM cluster within an XStringSet

Find matches from a PWM cluster within a set of sequences

## Usage

``` r
getClusterMatches(
  cl,
  stringset,
  rc = TRUE,
  min_score = "50%",
  best_only = FALSE,
  break_ties = c("all", "random", "first", "last", "central"),
  mc.cores = 1,
  ...
)

countClusterMatches(
  cl,
  stringset,
  rc = TRUE,
  min_score = "50%",
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
#> DataFrame with 190 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           1   17.3522         R       216       230          23       400
#> 2           2   11.9459         R       187       201          -6       400
#> 3          10   15.7958         R       176       190         -17       400
#> 4          24   13.5711         F       132       146         -61       400
#> 5          29   18.0880         F       193       207           0       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 186       824   11.1652         R        63        77        -130       400
#> 187       826   11.2094         R       196       210           3       400
#> 188       831   14.8580         R       377       391         184       400
#> 189       832   11.1978         R       212       226          19       400
#> 190       849   16.8796         F       313       327         120       400
#>           motif           match
#>     <character>  <DNAStringSet>
#> 1          ESR1 TGGTCACAGTGACCT
#> 2          ESR1 AGCCCAGAGTGACCT
#> 3          ESR1 GGGTCATCCTGTCCC
#> 4          ESR1 AGGCCACAGGGACCT
#> 5          ESR1 AGGTCACCCTGGCCC
#> ...         ...             ...
#> 186        ESR1 GGGTCGACCTGATCC
#> 187        ESR1 AGGTCAGAATGCTCA
#> 188        ESR1 AAGTCAGACTGTCCT
#> 189        ESR1 AGAACAAATTGACCT
#> 190        ESR1 AGGTCAGAATGACCG
#> 
#> $`ANDR/FOXA1`
#> DataFrame with 1088 rows and 9 columns
#>            seq     score direction     start       end from_centre seq_width
#>      <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1            5   10.3544         R       297       308       102.5       400
#> 2            7   13.7470         R       203       214         8.5       400
#> 3            9   11.9835         F       161       172       -33.5       400
#> 4           12   14.0694         F       199       210         4.5       400
#> 5           12   11.0321         F       342       353       147.5       400
#> ...        ...       ...       ...       ...       ...         ...       ...
#> 1084       845   9.43295         F       243       254        48.5       400
#> 1085       845   9.80593         R       259       270        64.5       400
#> 1086       846  10.58692         F       371       382       176.5       400
#> 1087       847  12.81220         R        19        36      -172.5       400
#> 1088       847  11.97986         R       208       219        13.5       400
#>            motif              match
#>      <character>     <DNAStringSet>
#> 1          FOXA1       TATTTGCACAGA
#> 2          FOXA1       TGTTTATTCTGT
#> 3          FOXA1       TATTTACAGAGC
#> 4          FOXA1       TGTTTGCTTTTG
#> 5          FOXA1       TGTTTATTGTTC
#> ...          ...                ...
#> 1084       FOXA1       TATTGACATTAA
#> 1085       FOXA1       TGTTGACTAAGT
#> 1086       FOXA1       TGTTTGCAATGG
#> 1087        ANDR TTTTTTTTTTTTTTTGCA
#> 1088       FOXA1       TGTTTATCTTTG
#> 
#> $ZN143
#> DataFrame with 76 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           3   23.2423         F       200       221        10.5       400
#> 2           3   25.2848         R       267       288        77.5       400
#> 3           3   24.3993         F       360       381       170.5       400
#> 4           6   18.0118         R       138       159       -51.5       400
#> 5           6   16.5389         F       178       199       -11.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 72        836   28.9222         R       166       187       -23.5       400
#> 73        837   30.0534         F       216       237        26.5       400
#> 74        837   21.8957         F       276       297        86.5       400
#> 75        837   28.0840         F       352       373       162.5       400
#> 76        848   22.0774         R       159       180       -30.5       400
#>           motif                  match
#>     <character>         <DNAStringSet>
#> 1         ZN143 CGCCCCCTGGGACTTGTAGTCT
#> 2         ZN143 GGGCCGCCGGGAGTTGTAGTTT
#> 3         ZN143 AGCGCCCTGGGAAATGTAGTCC
#> 4         ZN143 GGCCTGCCGGGCCTGGTAGTTC
#> 5         ZN143 CGCCTGCCGGTAGCTGTAGTCC
#> ...         ...                    ...
#> 72        ZN143 TGCATGCTGGGATTTGTAGTCC
#> 73        ZN143 TGCATGCTGGGAGTTGTAGTCT
#> 74        ZN143 GGCATGCAGGGAGTTGTAGTCG
#> 75        ZN143 TGCATGCTGGGAACTGTAGTCT
#> 76        ZN143 GGCACTGTGGGACTCGTAGTCT
#> 
#> $ZN281
#> DataFrame with 213 rows and 9 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           1   11.9978         R       160       174         -33       400
#> 2           2   12.4604         F       378       392         185       400
#> 3           9   11.7078         R        88       102        -105       400
#> 4          11   18.2091         R        35        49        -158       400
#> 5          11   13.0063         R        71        85        -122       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 209       815   12.9624         F        32        46        -161       400
#> 210       816   13.0790         F       307       321         114       400
#> 211       840   16.4984         F       366       380         173       400
#> 212       840   12.9600         F       372       386         179       400
#> 213       846   14.9216         R       284       298          91       400
#>           motif           match
#>     <character>  <DNAStringSet>
#> 1         ZN281 GGGGTGGGGCGGGGC
#> 2         ZN281 GGCAGGGGGTGGGCC
#> 3         ZN281 AGGTGTGGGAGGAGG
#> 4         ZN281 CGCGGGGGGAGGGGC
#> 5         ZN281 GAGCGGGGGAGGTGC
#> ...         ...             ...
#> 209       ZN281 GAGTGTGGGATGGGC
#> 210       ZN281 TGGGGGTGGAGGGGC
#> 211       ZN281 AGTAGGGGGTGGGGG
#> 212       ZN281 GGGTGGGGGAGAGAC
#> 213       ZN281 TAATGGGGGAGGGAA
#> 
# Or Just count them
countClusterMatches(ex_cl, ar_er_seq)
#>       ESR1 ANDR/FOXA1      ZN143      ZN281 
#>        199       1088         76        213 
# Compare this to individual counts
countPwmMatches(ex_pfm, ar_er_seq)
#>  ESR1  ANDR FOXA1 ZN143 ZN281 
#>   199   290  1041    76   213 
```
