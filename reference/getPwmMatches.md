# Find all PWM matches within an XStringSet

Find all PWM matches within a set of sequences

## Usage

``` r
getPwmMatches(
  pwm,
  stringset,
  rc = TRUE,
  min_score = "50%",
  best_only = FALSE,
  break_ties = c("all", "random", "first", "last", "central"),
  mc.cores = 1,
  ...
)
```

## Arguments

- pwm:

  A Position Weight Matrix, list of PWMs or universalmotif list

- stringset:

  An XStringSet

- rc:

  logical(1) Also find matches using the reverse complement of pwm

- min_score:

  The minimum score to return a match

- best_only:

  logical(1) Only return the best match

- break_ties:

  Method for breaking ties when only returning the best match Ignored
  when all matches are returned (the default)

- mc.cores:

  Passed to [mclapply](https://rdrr.io/r/parallel/mclapply.html) if
  passing multiple PWMs

- ...:

  Passed to [matchPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html)

## Value

A DataFrame with columns: `seq`, `score`, `direction`, `start`, `end`,
`from_centre`, `seq_width`, and `match`

The first three columns describe the sequence with matches, the score of
the match and whether the match was found using the forward or reverse
PWM. The columns `start`, `end` and `width` describe the where the match
was found in the sequence, whilst `from_centre` defines the distance
between the centre of the match and the centre of the sequence being
queried. The final column contains the matching fragment of the sequence
as an `XStringSet`.

When passing a list of PWMs, a list of the above DataFrames will be
returned.

## Details

Taking a set of sequences as an XStringSet, find all matches above the
supplied score (i.e. threshold) for a single Position Weight Matrix
(PWM), generally representing a transcription factor binding motif. By
default, matches are performed using the PWM as provided and the reverse
complement, however this can easily be disabled by setting `rc = FALSE`.

The function relies heavily on
[matchPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html) and
[Views](https://rdrr.io/pkg/IRanges/man/Views-class.html) for speed.

When choosing to return the best match (`best_only = TRUE`), only the
match with the highest score is returned for each sequence. Should there
be tied scores, the best match can be chosen as either the first, last,
most central, all tied matches, or choosing one at random (the default).

## Examples

``` r
## Load the example PWM
data("ex_pfm")
esr1 <- ex_pfm$ESR1

## Load the example Peaks
data("ar_er_seq")

## Return all matches
getPwmMatches(esr1, ar_er_seq)
#> DataFrame with 190 rows and 8 columns
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
#>               match
#>      <DNAStringSet>
#> 1   TGGTCACAGTGACCT
#> 2   AGCCCAGAGTGACCT
#> 3   GGGTCATCCTGTCCC
#> 4   AGGCCACAGGGACCT
#> 5   AGGTCACCCTGGCCC
#> ...             ...
#> 186 GGGTCGACCTGATCC
#> 187 AGGTCAGAATGCTCA
#> 188 AAGTCAGACTGTCCT
#> 189 AGAACAAATTGACCT
#> 190 AGGTCAGAATGACCG

## Just the best match
getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
#> DataFrame with 175 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           1   17.3522         R       216       230          23       400
#> 2           2   11.9459         R       187       201          -6       400
#> 3          10   15.7958         R       176       190         -17       400
#> 4          24   13.5711         F       132       146         -61       400
#> 5          29   18.0880         F       193       207           0       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 171       824   11.1652         R        63        77        -130       400
#> 172       826   11.2094         R       196       210           3       400
#> 173       831   14.8580         R       377       391         184       400
#> 174       832   11.1978         R       212       226          19       400
#> 175       849   16.8796         F       313       327         120       400
#>               match
#>      <DNAStringSet>
#> 1   TGGTCACAGTGACCT
#> 2   AGCCCAGAGTGACCT
#> 3   GGGTCATCCTGTCCC
#> 4   AGGCCACAGGGACCT
#> 5   AGGTCACCCTGGCCC
#> ...             ...
#> 171 GGGTCGACCTGATCC
#> 172 AGGTCAGAATGCTCA
#> 173 AAGTCAGACTGTCCT
#> 174 AGAACAAATTGACCT
#> 175 AGGTCAGAATGACCG

## Apply multiple PWMs as a list
getPwmMatches(ex_pfm, ar_er_seq, best_only = TRUE)
#> $ESR1
#> DataFrame with 175 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           1   17.3522         R       216       230          23       400
#> 2           2   11.9459         R       187       201          -6       400
#> 3          10   15.7958         R       176       190         -17       400
#> 4          24   13.5711         F       132       146         -61       400
#> 5          29   18.0880         F       193       207           0       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 171       824   11.1652         R        63        77        -130       400
#> 172       826   11.2094         R       196       210           3       400
#> 173       831   14.8580         R       377       391         184       400
#> 174       832   11.1978         R       212       226          19       400
#> 175       849   16.8796         F       313       327         120       400
#>               match
#>      <DNAStringSet>
#> 1   TGGTCACAGTGACCT
#> 2   AGCCCAGAGTGACCT
#> 3   GGGTCATCCTGTCCC
#> 4   AGGCCACAGGGACCT
#> 5   AGGTCACCCTGGCCC
#> ...             ...
#> 171 GGGTCGACCTGATCC
#> 172 AGGTCAGAATGCTCA
#> 173 AAGTCAGACTGTCCT
#> 174 AGAACAAATTGACCT
#> 175 AGGTCAGAATGACCG
#> 
#> $ANDR
#> DataFrame with 220 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          18   12.5485         F       200       217         8.5       400
#> 2          20   17.0683         R       176       193       -15.5       400
#> 3          21   14.9645         R       210       227        18.5       400
#> 4          22   13.9410         R       262       279        70.5       400
#> 5          26   12.7758         F       130       147       -61.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 216       833   23.0102         R       167       184       -24.5       400
#> 217       834   16.6073         F       189       206        -2.5       400
#> 218       840   14.0060         F       248       265        56.5       400
#> 219       841   16.1179         R       173       190       -18.5       400
#> 220       847   12.8122         R        19        36      -172.5       400
#>                  match
#>         <DNAStringSet>
#> 1   TGTGTTGAAATATTTACA
#> 2   TGTTCTAGATTATTTATA
#> 3   TCTCCTCTCTTGTTTACT
#> 4   TTTTATATTCTGTTTATA
#> 5   TTTTCCTACAAGTTTACT
#> ...                ...
#> 216 TGTTCTTTTTTGTTTGTT
#> 217 TGTTCTTTCGTGTTTGAC
#> 218 TGTGCTCTTCTCTTTGCA
#> 219 TCTGCTTTATTGTTTGTT
#> 220 TTTTTTTTTTTTTTTGCA
#> 
#> $FOXA1
#> DataFrame with 563 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           5   10.3544         R       297       308       102.5       400
#> 2           7   13.7470         R       203       214         8.5       400
#> 3           9   11.9835         F       161       172       -33.5       400
#> 4          12   14.0694         F       199       210         4.5       400
#> 5          14   11.6783         F       114       125       -80.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 559       843  10.09923         R        67        78      -127.5       400
#> 560       844  14.34613         F        40        51      -154.5       400
#> 561       845   9.43295         F       243       254        48.5       400
#> 562       846  10.58692         F       371       382       176.5       400
#> 563       847  11.97986         R       208       219        13.5       400
#>              match
#>     <DNAStringSet>
#> 1     TATTTGCACAGA
#> 2     TGTTTATTCTGT
#> 3     TATTTACAGAGC
#> 4     TGTTTGCTTTTG
#> 5     TGTTTGCAGAGC
#> ...            ...
#> 559   TGTTTGTCTTTG
#> 560   TGTTTACTTTCT
#> 561   TATTGACATTAA
#> 562   TGTTTGCAATGG
#> 563   TGTTTATCTTTG
#> 
#> $ZN143
#> DataFrame with 39 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           3   24.3993         F       360       381       170.5       400
#> 2           6   16.5389         F       178       199       -11.5       400
#> 3          11   19.9562         F       140       161       -49.5       400
#> 4          30   26.8427         F       166       187       -23.5       400
#> 5          67   29.0591         F       210       231        20.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 35        829   22.4081         F       206       227        16.5       400
#> 36        830   24.5149         F       143       164       -46.5       400
#> 37        836   28.9222         R       166       187       -23.5       400
#> 38        837   30.0534         F       216       237        26.5       400
#> 39        848   22.0774         R       159       180       -30.5       400
#>                      match
#>             <DNAStringSet>
#> 1   AGCGCCCTGGGAAATGTAGTCC
#> 2   CGCCTGCCGGTAGCTGTAGTCC
#> 3   AGCCTCATGGGGGTTGGAGTCC
#> 4   AGCCTGCCGGGAGATGTAGTTC
#> 5   GGCATGCTGGGATTTGTAGTCT
#> ...                    ...
#> 35  GGCATGCTAGGAGTTGTAGTGT
#> 36  TGGTTTCTGGGAATTGTAGTGT
#> 37  TGCATGCTGGGATTTGTAGTCC
#> 38  TGCATGCTGGGAGTTGTAGTCT
#> 39  GGCACTGTGGGACTCGTAGTCT
#> 
#> $ZN281
#> DataFrame with 139 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1           1   11.9978         R       160       174         -33       400
#> 2           2   12.4604         F       378       392         185       400
#> 3           9   11.7078         R        88       102        -105       400
#> 4          11   18.2091         R        35        49        -158       400
#> 5          15   16.4907         R       224       238          31       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 135       811   17.6496         R        69        83        -124       400
#> 136       815   16.6439         F         7        21        -186       400
#> 137       816   13.0790         F       307       321         114       400
#> 138       840   16.4984         F       366       380         173       400
#> 139       846   14.9216         R       284       298          91       400
#>               match
#>      <DNAStringSet>
#> 1   GGGGTGGGGCGGGGC
#> 2   GGCAGGGGGTGGGCC
#> 3   AGGTGTGGGAGGAGG
#> 4   CGCGGGGGGAGGGGC
#> 5   AGGTGGGGGTTGGGC
#> ...             ...
#> 135 CGGAGGGGGCGGGGC
#> 136 GGGTGGAGGTGGGGG
#> 137 TGGGGGTGGAGGGGC
#> 138 AGTAGGGGGTGGGGG
#> 139 TAATGGGGGAGGGAA
#> 
```
