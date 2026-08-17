# Descriptive Statistics for One or More Variables

Computes a compact set of descriptive statistics useful for data
screening and psychometric work (including skewness and excess
kurtosis), with an optional correlation matrix for the numeric
variables.

## Usage

``` r
descriptives(x, correlations = FALSE, listwise = FALSE)
```

## Arguments

- x:

  A `data.frame`, tibble, or `matrix` containing the variables to
  summarize. A `matrix` is coerced to a `data.frame`.

- correlations:

  Logical. If `TRUE`, also return a correlation matrix for the numeric
  variables (see Details).

- listwise:

  Logical. If `TRUE`, apply listwise deletion (drop every row containing
  any `NA`) to `x` *before* computing any statistics. If `FALSE` (the
  default), each variable is summarized using its own available
  (non-missing) observations.

## Value

A list with two elements, returned in a stable shape regardless of the
`correlations` argument:

- `descriptives`:

  A `data.frame` with one row per input variable and columns `variable`,
  `type`, `n`, `n_missing`, `prop_missing`, `mean`, `median`, `sd`,
  `min`, `max`, `q25`, `q75`, `skewness`, and `kurtosis` (excess
  kurtosis). For non-numeric variables, the numeric summary columns are
  `NA`.

- `correlations`:

  A \\p \times p\\ correlation matrix among the numeric variables, or
  `NULL` if `correlations = FALSE` or fewer than two numeric variables
  are available. This is returned as a plain matrix (not a table),
  because its natural structure is a symmetric two-dimensional array;
  that format reads well and plugs directly into downstream tools such
  as [`cov2cor`](https://rdrr.io/r/stats/cor.html) or factor analysis
  and SEM software.

## Details

**Skewness and kurtosis.** The reported values use the bias-corrected
formulas commonly referred to as SAS/SPSS Type 2: \$\$\mathrm{skewness}
= \frac{n}{(n-1)(n-2)} \sum\_{i=1}^{n} \left(\frac{x_i -
\bar{x}}{s}\right)^3,\$\$ \$\$\mathrm{kurtosis} =
\frac{n(n+1)}{(n-1)(n-2)(n-3)} \sum\_{i=1}^{n} \left(\frac{x_i -
\bar{x}}{s}\right)^4 - \frac{3(n-1)^2}{(n-2)(n-3)},\$\$ where \\s\\ is
the (divisor-\\n-1\\) sample standard deviation. Kurtosis is reported as
excess kurtosis, so a normal distribution has an expected value of 0. As
a rough guide to deciding whether normal-theory inference (e.g., maximum
likelihood in factor analysis or SEM) is defensible, values of
\\\|\mathrm{skewness}\| \> 2\\ or \\\|\mathrm{kurtosis}\| \> 7\\ are
frequently flagged as problematic.

**Variable-type handling.** Each column is classified as `"numeric"`,
`"integer"`, `"logical"`, `"factor"`, `"character"`, or whatever its
first class is otherwise. Columns of type `"numeric"`, `"integer"`, and
`"logical"` receive full distributional summaries (logicals are coerced
so that `mean` is the proportion of `TRUE`s). All other columns report
only `n`, `n_missing`, and `prop_missing`; their remaining columns are
`NA`.

**Correlations.** The correlation matrix, when requested, uses Pearson
correlations with `use = "pairwise.complete.obs"` so that each pairwise
correlation uses the maximum information available. Only numeric-type
variables are included.

## See also

[`cor`](https://rdrr.io/r/stats/cor.html),
[`quantile`](https://rdrr.io/r/stats/quantile.html)

Other descriptive statistics:
[`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md),
[`skewness()`](https://yelleknek.github.io/DMAR/reference/skewness.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Four cognitive tests from the Holzinger and Swineford (1939) study
# (301 children in two schools).
hs_tests <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                                    "t4_lozenges",
                                    "t6_paragraph_comprehension")]
descriptives(hs_tests)
#> $descriptives
#>                     variable    type   n n_missing prop_missing      mean
#> 1       t1_visual_perception integer 301         0            0 29.614618
#> 2                   t2_cubes integer 301         0            0 24.352159
#> 3                t4_lozenges integer 301         0            0 18.003322
#> 4 t6_paragraph_comprehension integer 301         0            0  9.182724
#>   median       sd min max q25 q75   skewness   kurtosis
#> 1     30 7.004593   4  51  25  34 -0.2569003  0.3553640
#> 2     24 4.709802   9  37  21  27  0.4747983  0.3808109
#> 3     17 9.047835   2  36  11  25  0.3872808 -0.8883744
#> 4      9 3.492349   0  19   7  11  0.2701735  0.1225896
#> 
#> $correlations
#> NULL
#> 

# Include a correlation matrix (useful during scale development).
descriptives(hs_tests, correlations = TRUE)
#> $descriptives
#>                     variable    type   n n_missing prop_missing      mean
#> 1       t1_visual_perception integer 301         0            0 29.614618
#> 2                   t2_cubes integer 301         0            0 24.352159
#> 3                t4_lozenges integer 301         0            0 18.003322
#> 4 t6_paragraph_comprehension integer 301         0            0  9.182724
#>   median       sd min max q25 q75   skewness   kurtosis
#> 1     30 7.004593   4  51  25  34 -0.2569003  0.3553640
#> 2     24 4.709802   9  37  21  27  0.4747983  0.3808109
#> 3     17 9.047835   2  36  11  25  0.3872808 -0.8883744
#> 4      9 3.492349   0  19   7  11  0.2701735  0.1225896
#> 
#> $correlations
#>                            t1_visual_perception  t2_cubes t4_lozenges
#> t1_visual_perception                  1.0000000 0.2973455   0.4406680
#> t2_cubes                              0.2973455 1.0000000   0.3398490
#> t4_lozenges                           0.4406680 0.3398490   1.0000000
#> t6_paragraph_comprehension            0.3727063 0.1529302   0.1586396
#>                            t6_paragraph_comprehension
#> t1_visual_perception                        0.3727063
#> t2_cubes                                    0.1529302
#> t4_lozenges                                 0.1586396
#> t6_paragraph_comprehension                  1.0000000
#> 

# Mixed-type data: numeric summaries for the test scores, and
# type = "factor" (with NA numeric columns) for school.
descriptives(holzinger_swineford[, c("school", "t1_visual_perception",
                                     "t2_cubes")])
#> $descriptives
#>               variable    type   n n_missing prop_missing     mean median
#> 1               school  factor 301         0            0       NA     NA
#> 2 t1_visual_perception integer 301         0            0 29.61462     30
#> 3             t2_cubes integer 301         0            0 24.35216     24
#>         sd min max q25 q75   skewness  kurtosis
#> 1       NA  NA  NA  NA  NA         NA        NA
#> 2 7.004593   4  51  25  34 -0.2569003 0.3553640
#> 3 4.709802   9  37  21  27  0.4747983 0.3808109
#> 
#> $correlations
#> NULL
#> 

# Data with missing values: the revised paper form board and flags tests
# were administered only in the Grant-White school, so 156 of the 301
# children have no score. Per-variable N and missingness are reported.
hs_partial <- holzinger_swineford[, c("t1_visual_perception",
                                      "t25_paper_form_board_r",
                                      "t26_flags")]
descriptives(hs_partial)
#> $descriptives
#>                 variable    type   n n_missing prop_missing     mean median
#> 1   t1_visual_perception integer 301         0    0.0000000 29.61462     30
#> 2 t25_paper_form_board_r integer 145       156    0.5182724 15.64828     16
#> 3              t26_flags integer 145       156    0.5182724 36.30345     37
#>         sd min max q25 q75   skewness   kurtosis
#> 1 7.004593   4  51  25  34 -0.2569003  0.3553640
#> 2 3.085655   5  23  14  18 -0.4433998  0.5627528
#> 3 8.339435  13  48  31  43 -0.6484687 -0.2552943
#> 
#> $correlations
#> NULL
#> 

# The same data with listwise deletion applied first.
descriptives(hs_partial, listwise = TRUE)
#> $descriptives
#>                 variable    type   n n_missing prop_missing     mean median
#> 1   t1_visual_perception integer 145         0            0 29.57931     30
#> 2 t25_paper_form_board_r integer 145         0            0 15.64828     16
#> 3              t26_flags integer 145         0            0 36.30345     37
#>         sd min max q25 q75   skewness    kurtosis
#> 1 6.913824  11  51  25  34 -0.1189539 -0.04569674
#> 2 3.085655   5  23  14  18 -0.4433998  0.56275284
#> 3 8.339435  13  48  31  43 -0.6484687 -0.25529427
#> 
#> $correlations
#> NULL
#> 
```
