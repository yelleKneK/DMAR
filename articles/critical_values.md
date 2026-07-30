# Reproducing the Textbook Critical-Value Tables

The tables of critical values printed in the back of a statistics
textbook, for the *t*, *F*, studentized range, studentized maximum
modulus, Dunnett, and Bryant–Paulson distributions, are not data to be
trusted on faith. Each one is the quantile of a known distribution, so
each one can be recomputed. This vignette shows that DMAR’s
critical-value family reproduces those tables, walks through the
parameter conventions that make the reproduction line up, and closes
with a case where recomputation caught a typographical error in a
printed table.

DMAR’s critical-value functions share a common surface (`alpha`, `df`, a
count of groups or comparisons, and a tidy one-row return):

| Function | Distribution | Appendix table it reproduces |
|----|----|----|
| [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md), [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md) | Student’s *t*, standard normal | A.1 |
| [`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md) | studentized range | A.4 |
| [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) | studentized maximum modulus | A.5 |
| [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md) | many-to-one (Dunnett) | A.6, A.7 |
| [`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md) | ANCOVA generalized studentized range | A.8 |
| [`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md) | Scheffe (an *F* transformation) | derived from A.2 |

The base-R quantile functions
[`qt()`](https://rdrr.io/r/stats/TDist.html),
[`qf()`](https://rdrr.io/r/stats/Fdist.html),
[`qtukey()`](https://rdrr.io/r/stats/Tukey.html), and
[`qchisq()`](https://rdrr.io/r/stats/Chisquare.html) cover the *t*, *F*,
studentized range, and chi square tables directly; DMAR adds the
multiple-comparison distributions that base R does not ship.

## The *t* Table (A.1)

Appendix Table A.1 is one-tailed *t* critical values at the
Bonferroni-spaced levels .05, .05/2, .05/4, .05/6, .05/8, .05/10. Each
cell is `qt(1 - alpha, df)`:

``` r
alpha1 <- 0.05 / c(1, 2, 4, 6, 8, 10)
tab <- sapply(alpha1, function(a) round(qt(1 - a, df = c(5, 10, 20, 60)), 2))
dimnames(tab) <- list(df = c(5, 10, 20, 60), alpha1 = signif(alpha1, 3))
tab
#>     alpha1
#> df   0.05 0.025 0.0125 0.00833 0.00625 0.005
#>   5  2.02  2.57   3.16    3.53    3.81  4.03
#>   10 1.81  2.23   2.63    2.87    3.04  3.17
#>   20 1.72  2.09   2.42    2.61    2.74  2.85
#>   60 1.67  2.00   2.30    2.46    2.58  2.66
```

The printed value at df = 1, alpha = .05 is 6.31 = 6.31; the infinite-df
row is the normal quantile
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md). A
subtlety worth noting: the fourth column’s header prints “.0083”, but
the exact level is 1/6 of .05 = .008333…, and it is the exact fraction
that reproduces the printed cells. A literal .0083 would move a dozen
cells by 0.01.

## The Studentized Range (A.4) and Maximum Modulus (A.5)

Tukey’s studentized range is
[`qtukey()`](https://rdrr.io/r/stats/Tukey.html), wrapped by
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md):

``` r
# Six groups, 14 error df, alpha = .05: Appendix A.4 prints 4.64.
qtukey(0.95, nmeans = 6, df = 14)
#> [1] 4.638538
```

The studentized maximum modulus
([`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md))
needs a word about its index. Its argument is `n_comparisons`, the
number *m* of simultaneous statistics. When the maximum modulus is used
for all pairwise comparisons among *a* groups, that number is *m =
a(a-1)/2*, not *a*. So the “Number of Groups” column of a printed SMM
table maps to `n_comparisons = a * (a - 1) / 2`:

``` r
# Appendix A.5, alpha = .05, "groups" a = 4, 5, 6 in the large-sample limit.
a <- c(4, 5, 6)
data.frame(groups = a, n_comparisons = a * (a - 1) / 2,
           cv = round(sapply(a * (a - 1) / 2, function(m)
             cv_smm(.05, df = Inf, n_comparisons = m, verbose = FALSE)$value), 2),
           printed = c(2.63, 2.80, 2.93))
#>   groups n_comparisons   cv printed
#> 1      4             6 2.63    2.63
#> 2      5            10 2.80    2.80
#> 3      6            15 2.93    2.93
```

([`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md)
evaluates the one-dimensional integral that defines the maximum modulus
and inverts it deterministically, so the value carries no Monte Carlo
error and needs no seed.)

## The Dunnett Tables (A.6, A.7)

[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
compares each of several treatments to a single control. Its
`n_comparisons` argument is the number of treatments, that is, *a - 1*
where *a* counts all groups including the control:

``` r
# Two-sided, a = 4 groups (3 treatments vs control), 10 error df: A.6 prints 2.76.
cv_dunnett(alpha_level = .05, df = 10, n_comparisons = 3, alternative = "not_equal")
```

| term     | value | area_less | area_greater |
|:---------|:------|:----------|:-------------|
| upper_cv | 2.76  | 0.95      | 0.05         |

With one treatment the many-to-one family collapses to a single
two-sided *t* test, a useful check on the argument mapping:

``` r
c(dunnett = cv_dunnett(.05, df = 10, n_comparisons = 1,
                       alternative = "not_equal")$value,
  t       = qt(1 - .05/2, df = 10))
#>  dunnett        t 
#> 2.228139 2.228139
```

Appendix A.6 and A.7 reproduce Dunnett’s own tables. Recomputed with the
exact integral above, the printed values agree throughout except for a
scattering of cells at the .01 level that differ by 0.01 in the last
digit, a property of the historical tabulation rather than a
transcription error. (Reaching for a Monte Carlo multivariate- quantile
instead of the deterministic integral would introduce sampling error of
its own, several hundredths at the .01 level, which is why
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
integrates rather than simulates.)

## Bryant–Paulson (A.8): Where Recomputation Caught an Error

The Bryant–Paulson generalized studentized range is the
analysis-of-covariance analogue of Tukey’s range, for comparing
covariate-**adjusted** means when the covariate is random. It is the one
distribution here that most software, including base R, does not compute
natively;
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
and the underlying
[`qbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
supply it. (For its use in a full ANCOVA analysis see the
`cv_bryant_paulson` and `bryant_paulson_simulation` vignettes.)

It reproduces the published anchor exactly:

``` r
# Bryant & Paulson (1976) Table 1a: q_{.05; p=1, a=6, nu=14} = 4.83.
cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1,
                  verbose = FALSE)$value
#> [1] 4.829856
```

Recomputing an entire block of Appendix Table A.8 (error df = 20, three
covariates, alpha = .01) against the printed values shows a single cell
out of step:

``` r
a_vals <- c(2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20)
recomputed <- sapply(a_vals, function(k)
  cv_bryant_paulson(.01, df = 20, groups = k, covariates = 3,
                    verbose = FALSE)$value)
printed <- c(4.35, 5.03, 5.45, 5.75, 6.09, 6.19, 6.36, 6.63, 6.85, 7.19, 7.45)
data.frame(a = a_vals, printed = printed, recomputed = round(recomputed, 2),
           diff = round(recomputed, 2) - printed)
#>     a printed recomputed diff
#> 1   2    4.35       4.35  0.0
#> 2   3    5.03       5.03  0.0
#> 3   4    5.45       5.45  0.0
#> 4   5    5.75       5.75  0.0
#> 5   6    6.09       5.99 -0.1
#> 6   7    6.19       6.19  0.0
#> 7   8    6.36       6.36  0.0
#> 8  10    6.63       6.63  0.0
#> 9  12    6.85       6.85  0.0
#> 10 16    7.19       7.19  0.0
#> 11 20    7.45       7.45  0.0
```

Every cell agrees except *a* = 6, where the table prints 6.09 and the
distribution gives 5.99. Two arguments settle which is right. First, the
critical value is concave in the number of means, so the increments
across the row must decrease; 5.99 keeps them decreasing while 6.09
forces a reversal:

``` r
rbind(printed_increments    = diff(printed),
      recomputed_increments = round(diff(recomputed), 2))
#>                       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> printed_increments    0.68 0.42  0.3 0.34  0.1 0.17 0.27 0.22 0.34  0.26
#> recomputed_increments 0.68 0.42  0.3 0.24  0.2 0.17 0.27 0.22 0.34  0.26
```

Second, a direct Monte Carlo of the statistic, drawn from its definition
as the range of *a* standard normals divided by an independent scale,
agrees with the computed 5.99, not with 6.09:

``` r
G <- 2e5; a <- 6; nu <- 20; p <- 3
R   <- apply(matrix(rnorm(G * a), ncol = a), 1, function(z) max(z) - min(z))
S   <- sqrt(rchisq(G, nu) / nu)
delta <- rbeta(G, (nu + 1) / 2, p / 2)                 # covariate-shrinkage factor
Q   <- R / (S * sqrt(delta))                           # Bryant-Paulson statistic
c(simulated_99th_pctile = round(quantile(Q, 0.99, names = FALSE), 2),
  cv_bryant_paulson     = round(cv_bryant_paulson(.01, df = 20, groups = 6,
                                covariates = 3, verbose = FALSE)$value, 2))
#> simulated_99th_pctile     cv_bryant_paulson 
#>                  6.01                  5.99
```

The printed 6.09 is a typographical error for 5.99. Recomputation is not
just a way to avoid typing a table by hand; it is a way to check the
tables one has.

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey’s method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Dunnett, C. W. (1964). New tables for multiple comparisons with a
control. *Biometrics, 20*, 482–491.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.
