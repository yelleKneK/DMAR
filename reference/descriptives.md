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
# Employee attitude survey data (30 respondents, 7 numeric rating items).
descriptives(attitude)
#> $descriptives
#>     variable    type  n n_missing prop_missing     mean median        sd min
#> 1     rating numeric 30         0            0 64.63333   65.5 12.172562  40
#> 2 complaints numeric 30         0            0 66.60000   65.0 13.314757  37
#> 3 privileges numeric 30         0            0 53.13333   51.5 12.235430  30
#> 4   learning numeric 30         0            0 56.36667   56.5 11.737013  34
#> 5     raises numeric 30         0            0 64.63333   63.5 10.397226  43
#> 6   critical numeric 30         0            0 74.76667   77.5  9.894908  49
#> 7    advance numeric 30         0            0 42.93333   41.0 10.288706  25
#>   max   q25   q75   skewness    kurtosis
#> 1  85 58.75 71.75 -0.3967148 -0.49460977
#> 2  90 58.50 77.00 -0.2387632 -0.38172518
#> 3  83 45.00 62.50  0.4202101 -0.04219127
#> 4  75 47.00 66.75 -0.0598894 -1.07638294
#> 5  88 58.25 71.00  0.2189518 -0.28201359
#> 6  92 69.25 80.00 -0.9596072  0.69175815
#> 7  72 35.00 47.75  0.9425594  1.07314413
#> 
#> $correlations
#> NULL
#> 

# Include a correlation matrix (useful during scale development).
descriptives(attitude, correlations = TRUE)
#> $descriptives
#>     variable    type  n n_missing prop_missing     mean median        sd min
#> 1     rating numeric 30         0            0 64.63333   65.5 12.172562  40
#> 2 complaints numeric 30         0            0 66.60000   65.0 13.314757  37
#> 3 privileges numeric 30         0            0 53.13333   51.5 12.235430  30
#> 4   learning numeric 30         0            0 56.36667   56.5 11.737013  34
#> 5     raises numeric 30         0            0 64.63333   63.5 10.397226  43
#> 6   critical numeric 30         0            0 74.76667   77.5  9.894908  49
#> 7    advance numeric 30         0            0 42.93333   41.0 10.288706  25
#>   max   q25   q75   skewness    kurtosis
#> 1  85 58.75 71.75 -0.3967148 -0.49460977
#> 2  90 58.50 77.00 -0.2387632 -0.38172518
#> 3  83 45.00 62.50  0.4202101 -0.04219127
#> 4  75 47.00 66.75 -0.0598894 -1.07638294
#> 5  88 58.25 71.00  0.2189518 -0.28201359
#> 6  92 69.25 80.00 -0.9596072  0.69175815
#> 7  72 35.00 47.75  0.9425594  1.07314413
#> 
#> $correlations
#>               rating complaints privileges  learning    raises  critical
#> rating     1.0000000  0.8254176  0.4261169 0.6236782 0.5901390 0.1564392
#> complaints 0.8254176  1.0000000  0.5582882 0.5967358 0.6691975 0.1877143
#> privileges 0.4261169  0.5582882  1.0000000 0.4933310 0.4454779 0.1472331
#> learning   0.6236782  0.5967358  0.4933310 1.0000000 0.6403144 0.1159652
#> raises     0.5901390  0.6691975  0.4454779 0.6403144 1.0000000 0.3768830
#> critical   0.1564392  0.1877143  0.1472331 0.1159652 0.3768830 1.0000000
#> advance    0.1550863  0.2245796  0.3432934 0.5316198 0.5741862 0.2833432
#>              advance
#> rating     0.1550863
#> complaints 0.2245796
#> privileges 0.3432934
#> learning   0.5316198
#> raises     0.5741862
#> critical   0.2833432
#> advance    1.0000000
#> 

# Mixed-type data: numeric summaries for the four measurements, and type="factor"
# (with NA numeric columns) for Species.
descriptives(iris)
#> $descriptives
#>       variable    type   n n_missing prop_missing     mean median        sd min
#> 1 Sepal.Length numeric 150         0            0 5.843333   5.80 0.8280661 4.3
#> 2  Sepal.Width numeric 150         0            0 3.057333   3.00 0.4358663 2.0
#> 3 Petal.Length numeric 150         0            0 3.758000   4.35 1.7652982 1.0
#> 4  Petal.Width numeric 150         0            0 1.199333   1.30 0.7622377 0.1
#> 5      Species  factor 150         0            0       NA     NA        NA  NA
#>   max q25 q75   skewness  kurtosis
#> 1 7.9 5.1 6.4  0.3149110 -0.552064
#> 2 4.4 2.8 3.3  0.3189657  0.228249
#> 3 6.9 1.6 5.1 -0.2748842 -1.402103
#> 4 2.5 0.3 1.8 -0.1029667 -1.340604
#> 5  NA  NA  NA         NA        NA
#> 
#> $correlations
#> NULL
#> 

# Data with missing values. Per-variable N and missingness are reported.
descriptives(airquality)
#> $descriptives
#>   variable    type   n n_missing prop_missing       mean median        sd  min
#> 1    Ozone integer 116        37   0.24183007  42.129310   31.5 32.987885  1.0
#> 2  Solar.R integer 146         7   0.04575163 185.931507  205.0 90.058422  7.0
#> 3     Wind numeric 153         0   0.00000000   9.957516    9.7  3.523001  1.7
#> 4     Temp integer 153         0   0.00000000  77.882353   79.0  9.465270 56.0
#> 5    Month integer 153         0   0.00000000   6.993464    7.0  1.416522  5.0
#> 6      Day integer 153         0   0.00000000  15.803922   16.0  8.864520  1.0
#>     max    q25    q75     skewness   kurtosis
#> 1 168.0  18.00  63.25  1.241796404  1.2903027
#> 2 334.0 115.75 258.75 -0.428044526 -0.9684668
#> 3  20.7   7.40  11.50  0.347817775  0.1114183
#> 4  97.0  72.00  85.00 -0.377884464 -0.4035054
#> 5   9.0   6.00   8.00 -0.002391498 -1.2975830
#> 6  31.0   8.00  23.00  0.002651853 -1.1988345
#> 
#> $correlations
#> NULL
#> 

# The same data with listwise deletion applied first.
descriptives(airquality, listwise = TRUE)
#> $descriptives
#>   variable    type   n n_missing prop_missing       mean median        sd  min
#> 1    Ozone integer 111         0            0  42.099099   31.0 33.275969  1.0
#> 2  Solar.R integer 111         0            0 184.801802  207.0 91.152302  7.0
#> 3     Wind numeric 111         0            0   9.939640    9.7  3.557713  2.3
#> 4     Temp integer 111         0            0  77.792793   79.0  9.529969 57.0
#> 5    Month integer 111         0            0   7.216216    7.0  1.473434  5.0
#> 6      Day integer 111         0            0  15.945946   16.0  8.707194  1.0
#>     max   q25   q75    skewness   kurtosis
#> 1 168.0  18.0  62.0  1.26526649  1.3165400
#> 2 334.0 113.5 255.5 -0.49293305 -0.9164162
#> 3  20.7   7.4  11.5  0.46190700  0.3504590
#> 4  97.0  71.0  84.5 -0.22819120 -0.6430780
#> 5   9.0   6.0   9.0 -0.29527316 -1.2482951
#> 6  31.0   9.0  22.5 -0.01300861 -1.0385216
#> 
#> $correlations
#> NULL
#> 
```
