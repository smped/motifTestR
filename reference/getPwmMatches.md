# Find all PWM matches within an XStringSet

Find all PWM matches within a set of sequences

## Usage

``` r
getPwmMatches(
  pwm,
  stringset,
  rc = TRUE,
  min_score = "80%",
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
#> DataFrame with 22 rows and 8 columns
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
#>               match
#>      <DNAStringSet>
#> 1   AGGTCACCCTGGCCC
#> 2   AGGTCACCGTGACCC
#> 3   AGGTGACCCTGACCT
#> 4   GGGTCACACTGTCCT
#> 5   AGGTCACAATGACCT
#> ...             ...
#> 18  AGGTCACCCTGACCG
#> 19  GGGTCAGCATGACCT
#> 20  AGGACACACTGACCT
#> 21  AGGTCACCCTAACCT
#> 22  AGGTTAGCCTGACCT

## Just the best match
getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
#> DataFrame with 22 rows and 8 columns
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
#>               match
#>      <DNAStringSet>
#> 1   AGGTCACCCTGGCCC
#> 2   AGGTCACCGTGACCC
#> 3   AGGTGACCCTGACCT
#> 4   GGGTCACACTGTCCT
#> 5   AGGTCACAATGACCT
#> ...             ...
#> 18  AGGTCACCCTGACCG
#> 19  GGGTCAGCATGACCT
#> 20  AGGACACACTGACCT
#> 21  AGGTCACCCTAACCT
#> 22  AGGTTAGCCTGACCT

## Apply multiple PWMs as a list
getPwmMatches(ex_pfm, ar_er_seq, best_only = TRUE)
#> $ESR1
#> DataFrame with 22 rows and 8 columns
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
#>               match
#>      <DNAStringSet>
#> 1   AGGTCACCCTGGCCC
#> 2   AGGTCACCGTGACCC
#> 3   AGGTGACCCTGACCT
#> 4   GGGTCACACTGTCCT
#> 5   AGGTCACAATGACCT
#> ...             ...
#> 18  AGGTCACCCTGACCG
#> 19  GGGTCAGCATGACCT
#> 20  AGGACACACTGACCT
#> 21  AGGTCACCCTAACCT
#> 22  AGGTTAGCCTGACCT
#> 
#> $ANDR
#> DataFrame with 8 rows and 8 columns
#>         seq     score direction     start       end from_centre seq_width
#>   <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1        27   18.9055         F        58        75      -133.5       400
#> 2       110   20.5599         F       235       252        43.5       400
#> 3       285   18.9230         F       205       222        13.5       400
#> 4       519   18.9718         R       264       281        72.5       400
#> 5       701   20.3572         F       329       346       137.5       400
#> 6       704   20.5870         F        68        85      -123.5       400
#> 7       708   20.9669         F       278       295        86.5       400
#> 8       833   23.0102         R       167       184       -24.5       400
#>                match
#>       <DNAStringSet>
#> 1 TGTTCTTTTTTGTTGATT
#> 2 TGTCCTTTTCTGTTTATT
#> 3 TGTTCCTCTCTGTTTACC
#> 4 TGTTCAGCTTTGTTTGCT
#> 5 TGTTCTTTTGTATTTGCT
#> 6 TGTTCTTCTATGTTTATT
#> 7 TGTTCTTTATTATTTGCT
#> 8 TGTTCTTTTTTGTTTGTT
#> 
#> $FOXA1
#> DataFrame with 107 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          12   14.0694         F       199       210         4.5       400
#> 2          16   14.4605         F       177       188       -17.5       400
#> 3          17   15.1441         F       341       352       146.5       400
#> 4          18   14.4417         F       301       312       106.5       400
#> 5          21   14.6369         R       206       217        11.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 103       793   14.3562         F       194       205        -0.5       400
#> 104       816   14.1175         R        63        74      -131.5       400
#> 105       817   14.9185         R       291       302        96.5       400
#> 106       826   14.5836         F       261       272        66.5       400
#> 107       844   14.3461         F        40        51      -154.5       400
#>              match
#>     <DNAStringSet>
#> 1     TGTTTGCTTTTG
#> 2     TGTTTACTTTCC
#> 3     TGTTTATTTAGG
#> 4     TGTTTATTCTGG
#> 5     TGTTTACTCAAC
#> ...            ...
#> 103   TGTTTACTTTAA
#> 104   TGTTTATTTTAG
#> 105   TGTTTACACAGT
#> 106   TATTTACTTTAG
#> 107   TGTTTACTTTCT
#> 
#> $ZN143
#> DataFrame with 15 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1          30   26.8427         F       166       187       -23.5       400
#> 2          67   29.0591         F       210       231        20.5       400
#> 3         118   29.0063         R       205       226        15.5       400
#> 4         182   28.0840         F       225       246        35.5       400
#> 5         330   27.3485         R       196       217         6.5       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 11        569   28.1714         R        10        31      -179.5       400
#> 12        750   28.4852         R       151       172       -38.5       400
#> 13        829   26.6710         R       151       172       -38.5       400
#> 14        836   28.9222         R       166       187       -23.5       400
#> 15        837   30.0534         F       216       237        26.5       400
#>                      match
#>             <DNAStringSet>
#> 1   AGCCTGCCGGGAGATGTAGTTC
#> 2   GGCATGCTGGGATTTGTAGTCT
#> 3   TGCCTCCTGGGAAATGTAGTCC
#> 4   TGCATGCTGGGAACTGTAGTCT
#> 5   AGCCTTGTGGGAGTTGTAGTTT
#> ...                    ...
#> 11  GGCATTTTGGGAGTTGTAGTTT
#> 12  CGCATGCTGGGAATTGTAGTTC
#> 13  TGCCCGCTGGGAACTGTAGTCC
#> 14  TGCATGCTGGGATTTGTAGTCC
#> 15  TGCATGCTGGGAGTTGTAGTCT
#> 
#> $ZN281
#> DataFrame with 13 rows and 8 columns
#>           seq     score direction     start       end from_centre seq_width
#>     <integer> <numeric>  <factor> <integer> <integer>   <numeric> <integer>
#> 1         109   19.6553         R       369       383         176       400
#> 2         118   20.1871         F        60        74        -133       400
#> 3         122   19.3263         F        95       109         -98       400
#> 4         171   18.5467         R        84        98        -109       400
#> 5         192   19.1012         F       260       274          67       400
#> ...       ...       ...       ...       ...       ...         ...       ...
#> 9         456   21.2625         R       343       357         150       400
#> 10        507   18.6243         R       175       189         -18       400
#> 11        668   20.0793         R       168       182         -25       400
#> 12        763   22.6040         R       274       288          81       400
#> 13        764   18.8379         R       310       324         117       400
#>               match
#>      <DNAStringSet>
#> 1   AGTTGGGGGAGGGGC
#> 2   GGCGGGGGGAGGGGA
#> 3   GAATGGGGGAGGGGC
#> 4   GGATGGGGGAAGGGG
#> 5   GGGAGGGGGCGGGGG
#> ...             ...
#> 9   CGGTGGGGGAGGGGG
#> 10  GGGAGGGGGAGGGAG
#> 11  GGGTGGGGGTGGGGG
#> 12  GGGTGGGGGAGGGGG
#> 13  AGTGGGGGGAGGGGA
#> 
```
