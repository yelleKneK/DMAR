# Simulate Data From a One-Way Fixed-Effects ANOVA Model

Generates random data appropriate for a one-way fixed-effects analysis
of variance. Each group's observations are drawn from a normal
distribution with that group's population mean and a common (or
per-group) standard deviation. Per-group sample sizes may be equal or
unequal.

## Usage

``` r
simulate_anova_data(mu, sigma, a, n, seed = NULL)
```

## Arguments

- mu:

  A numeric vector of length `a` giving the population mean in each
  group.

- sigma:

  Within-group population standard deviation. Either a single number
  (homoscedastic; common across groups) or a numeric vector of length
  `a` (heteroscedastic; one SD per group).

- a:

  The number of fixed levels of the grouping factor (per the convention
  used throughout DMAR for fixed-factor designs).

- n:

  A single number (equal sample size per group) or a numeric vector of
  length `a` giving the sample size in each group.

- seed:

  Optional integer random seed for reproducibility (default `NULL`;
  supply an integer such as `113` for reproducible output).

## Value

A long-format `data.frame` with one row per simulated subject and two
columns:

- `group`:

  A factor with `a` levels.

- `y`:

  Numeric simulated outcome.

## Details

The fixed-effects ANOVA model assumes group-specific means and a common
within-group variance. Setting `sigma` to a vector relaxes the
homoscedasticity assumption; in that case the simulated data violate the
standard ANOVA assumption (a useful feature for studying robustness or
the performance of Welch-style alternatives).

## See also

[`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_regression_data`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md),
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md)

Other data simulators:
[`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_ancova_factorial_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
[`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Three-group ANOVA, equal n per group.
set.seed(113)
d <- simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3, n = 30)
aggregate(y ~ group, data = d, FUN = mean)
#>   group        y
#> 1     1 51.43086
#> 2     2 56.55109
#> 3     3 59.40916

# Same design with unequal n per group.
simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3,
                    n = c(40, 30, 20), seed = 113)
#>    group        y
#> 1      1 51.06684
#> 2      1 61.00177
#> 3      1 55.98972
#> 4      1 39.64920
#> 5      1 45.52983
#> 6      1 36.14041
#> 7      1 40.17117
#> 8      1 50.49970
#> 9      1 47.14297
#> 10     1 47.77378
#> 11     1 65.30042
#> 12     1 43.84951
#> 13     1 57.26283
#> 14     1 46.10848
#> 15     1 52.53280
#> 16     1 44.44187
#> 17     1 55.94404
#> 18     1 54.80560
#> 19     1 49.00888
#> 20     1 58.95260
#> 21     1 52.26633
#> 22     1 63.98192
#> 23     1 59.57886
#> 24     1 47.31220
#> 25     1 57.17865
#> 26     1 64.29169
#> 27     1 54.74173
#> 28     1 34.31343
#> 29     1 58.00016
#> 30     1 48.08852
#> 31     1 38.53479
#> 32     1 52.50063
#> 33     1 52.56658
#> 34     1 53.72451
#> 35     1 45.62643
#> 36     1 59.64269
#> 37     1 37.32579
#> 38     1 54.28905
#> 39     1 66.04419
#> 40     1 45.86928
#> 41     2 55.18918
#> 42     2 74.91562
#> 43     2 54.94984
#> 44     2 55.11040
#> 45     2 47.78754
#> 46     2 67.58346
#> 47     2 54.95464
#> 48     2 47.34617
#> 49     2 63.73710
#> 50     2 62.62330
#> 51     2 46.34935
#> 52     2 59.40388
#> 53     2 49.73349
#> 54     2 60.69886
#> 55     2 61.50479
#> 56     2 62.50074
#> 57     2 50.04821
#> 58     2 46.99661
#> 59     2 67.62389
#> 60     2 51.35184
#> 61     2 65.83470
#> 62     2 50.75471
#> 63     2 60.89853
#> 64     2 35.38485
#> 65     2 55.06033
#> 66     2 57.77432
#> 67     2 58.16157
#> 68     2 32.49224
#> 69     2 52.05944
#> 70     2 69.16834
#> 71     3 58.12576
#> 72     3 55.68143
#> 73     3 46.36344
#> 74     3 56.35509
#> 75     3 56.49787
#> 76     3 57.14659
#> 77     3 71.34016
#> 78     3 60.81402
#> 79     3 72.91256
#> 80     3 64.44470
#> 81     3 75.55203
#> 82     3 49.14457
#> 83     3 47.34800
#> 84     3 45.03782
#> 85     3 55.03748
#> 86     3 64.96932
#> 87     3 61.64875
#> 88     3 63.69860
#> 89     3 70.41855
#> 90     3 62.14909

# Heteroscedastic case: each group has its own SD.
simulate_anova_data(mu = c(50, 55, 60), sigma = c(5, 8, 12),
                    a = 3, n = 30, seed = 113)
#>    group        y
#> 1      1 50.66677
#> 2      1 56.87611
#> 3      1 53.74358
#> 4      1 43.53075
#> 5      1 47.20615
#> 6      1 41.33775
#> 7      1 43.85698
#> 8      1 50.31231
#> 9      1 48.21436
#> 10     1 48.60861
#> 11     1 59.56276
#> 12     1 46.15594
#> 13     1 54.53927
#> 14     1 47.56780
#> 15     1 51.58300
#> 16     1 46.52617
#> 17     1 53.71503
#> 18     1 53.00350
#> 19     1 49.38055
#> 20     1 55.59537
#> 21     1 51.41645
#> 22     1 58.73870
#> 23     1 55.98679
#> 24     1 48.32013
#> 25     1 54.48666
#> 26     1 58.93231
#> 27     1 52.96358
#> 28     1 40.19589
#> 29     1 55.00010
#> 30     1 48.80533
#> 31     2 43.53479
#> 32     2 57.50063
#> 33     2 57.56658
#> 34     2 58.72451
#> 35     2 50.62643
#> 36     2 64.64269
#> 37     2 42.32579
#> 38     2 59.28905
#> 39     2 71.04419
#> 40     2 50.86928
#> 41     2 55.18918
#> 42     2 74.91562
#> 43     2 54.94984
#> 44     2 55.11040
#> 45     2 47.78754
#> 46     2 67.58346
#> 47     2 54.95464
#> 48     2 47.34617
#> 49     2 63.73710
#> 50     2 62.62330
#> 51     2 46.34935
#> 52     2 59.40388
#> 53     2 49.73349
#> 54     2 60.69886
#> 55     2 61.50479
#> 56     2 62.50074
#> 57     2 50.04821
#> 58     2 46.99661
#> 59     2 67.62389
#> 60     2 51.35184
#> 61     3 76.25205
#> 62     3 53.63206
#> 63     3 68.84779
#> 64     3 30.57728
#> 65     3 60.09049
#> 66     3 64.16147
#> 67     3 64.74236
#> 68     3 26.23836
#> 69     3 55.58915
#> 70     3 81.25251
#> 71     3 57.18863
#> 72     3 53.52214
#> 73     3 39.54515
#> 74     3 54.53264
#> 75     3 54.74680
#> 76     3 55.71989
#> 77     3 77.01024
#> 78     3 61.22103
#> 79     3 79.36884
#> 80     3 66.66705
#> 81     3 83.32804
#> 82     3 43.71686
#> 83     3 41.02200
#> 84     3 37.55672
#> 85     3 52.55622
#> 86     3 67.45398
#> 87     3 62.47313
#> 88     3 65.54790
#> 89     3 75.62783
#> 90     3 63.22363
```
