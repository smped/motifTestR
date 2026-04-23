# Form a set of random, matching ranges for bootstrapping or permuting

Form a set of ranges from y which (near) exactly match those in x for
use as a background set requiring matching

## Usage

``` r
makeRMRanges(x, y, ...)

# S4 method for class 'GRanges,GRanges'
makeRMRanges(
  x,
  y,
  exclude = GRanges(),
  n_iter = 1,
  n_total = NULL,
  replace = TRUE,
  ...,
  force_ol = TRUE
)

# S4 method for class 'GRangesList,GRangesList'
makeRMRanges(
  x,
  y,
  exclude = GRanges(),
  n_iter = 1,
  n_total = NULL,
  replace = TRUE,
  mc.cores = 1,
  ...,
  force_ol = TRUE,
  unlist = TRUE
)
```

## Arguments

- x:

  GRanges/GRangesList with ranges to be matched

- y:

  GRanges/GRangesList with ranges to select random matching ranges from

- ...:

  Not used

- exclude:

  GRanges of ranges to omit from testing

- n_iter:

  The number of times to repeat the random selection process

- n_total:

  Setting this value will over-ride anything set using n_iter. Can be
  vector of any length, corresponding to the length of x, when x is a
  GRangesList

- replace:

  logical(1) Sample with our without replacement when creating the set
  of random ranges.

- force_ol:

  logical(1) Enforce an overlap between every site in x and y

- mc.cores:

  Passsed to [mclapply](https://rdrr.io/r/parallel/mclapply.html)

- unlist:

  logical(1) Return as a sorted GRanges object, or leave as a
  GRangesList

## Value

A GRanges or GRangesList object

## Details

This function uses the width distribution of the 'test' ranges (i.e.
`x`) to randomly sample a set of ranges with matching width from the
ranges provided in `y`. The width distribution will clearly be exact
when a set of fixed-width ranges is passed to `x`, whilst random
sampling may yield some variability when matching ranges of variable
width.

When both x and y are GRanges objects, they are implicitly assumed to
both represent similar ranges, such as those overlapping a promoter or
enhancer. When passing two GRangesList objects, both objects are
expected to contain ranges annotated as belonging to key features, such
that the list elements in y must encompass all elements in x. For
example if `x` contains two elements named 'promoter' and 'intron', `y`
should also contain elements named 'promoter' and 'intron' and these
will be sampled as matching ranges for the same element in `x`. If
elements of `x` and `y` are not named, they are assumed to be in
matching order.

The default behaviour is to assume that randomly-generated ranges are
for iteration, and as such, ranges are randomly formed in multiples of
the number of 'test' ranges provided in `x`. The column `iteration` will
be added to the returned ranges. Placing any number into the `n_total`
argument will instead select a total number of ranges as specified here.
In this case, no `iteration` column will be included in the returned
ranges.

Sampling is assumed to be with replacement as this is most suitable for
bootstrapping and related procedures, although this can be disabled by
setting `replace = FALSE`

## Examples

``` r
## Load the example peaks
data("ar_er_peaks")
sq <- seqinfo(ar_er_peaks)
#> Loading required namespace: GenomeInfoDb
## Now sample size-matched ranges for two iterations from chr1
makeRMRanges(ar_er_peaks, GRanges(sq)[1], n_iter = 2)
#> GRanges object with 1698 ranges and 1 metadata column:
#>          seqnames              ranges strand | iteration
#>             <Rle>           <IRanges>  <Rle> | <integer>
#>      [1]     chr1         91036-91435      * |         2
#>      [2]     chr1       151570-151969      * |         2
#>      [3]     chr1       156793-157192      * |         1
#>      [4]     chr1       335925-336324      * |         2
#>      [5]     chr1       365758-366157      * |         1
#>      ...      ...                 ...    ... .       ...
#>   [1694]     chr1 248443254-248443653      * |         2
#>   [1695]     chr1 248493084-248493483      * |         2
#>   [1696]     chr1 248801396-248801795      * |         1
#>   [1697]     chr1 248869157-248869556      * |         2
#>   [1698]     chr1 249196752-249197151      * |         1
#>   -------
#>   seqinfo: 24 sequences from hg19 genome

## Or simply sample 100 ranges if not planning any iterative analyses
makeRMRanges(ar_er_peaks, GRanges(sq)[1], n_total = 100)
#> GRanges object with 100 ranges and 0 metadata columns:
#>         seqnames              ranges strand
#>            <Rle>           <IRanges>  <Rle>
#>     [1]     chr1     1594492-1594891      *
#>     [2]     chr1     6204139-6204538      *
#>     [3]     chr1     8926298-8926697      *
#>     [4]     chr1   10338571-10338970      *
#>     [5]     chr1   11417646-11418045      *
#>     ...      ...                 ...    ...
#>    [96]     chr1 237632083-237632482      *
#>    [97]     chr1 241318705-241319104      *
#>    [98]     chr1 243568836-243569235      *
#>    [99]     chr1 245708450-245708849      *
#>   [100]     chr1 247427291-247427690      *
#>   -------
#>   seqinfo: 24 sequences from hg19 genome

```
