# Kuder-Richardson Formula 20 (KR-20) With a Confidence Interval

Estimates Kuder-Richardson formula 20 (Kuder & Richardson, 1937) for a
homogeneous composite scored on dichotomous (0/1) items, and returns a
confidence interval for the population coefficient.

## Usage

``` r
reliability_kr20(
  data,
  ci_method = c("feldt", "bonett", "fisher", "hakstian_whalen", "ml", "ml_logistic",
    "adf", "adf_logistic", "bootstrap_se", "bootstrap_se_logistic", "percentile", "bca",
    "none"),
  conf_level = 0.95,
  B = 10000,
  seed = NULL,
  ncpus = 1L,
  progress = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of 0/1 item scores (rows are
  respondents, columns are items). Rows with any missing values are
  listwise-deleted. Non-binary values trigger an error, and so do fewer
  than two complete rows, from which no variance can be computed. A
  column that is not numeric, and an item that takes one value in every
  complete row (one no respondent passes, or every respondent does), are
  refused by name: such an item adds a constant to the total score,
  whose reliability it leaves as it is, while KR-20 would count it among
  the \\J\\ items, at \\J(J-2)/(J-1)^2\\ times the KR-20 of the other
  items, the coefficient of the same total. So is a total score that
  takes one value in every complete row (on two rows, as many items
  passed in one row as in the other), whose reliability is not defined.
  For incomplete data,
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  with `missing = "fiml"` computes the same coefficient (alpha of 0/1
  items is KR-20) from the incomplete rows as well. SPSS user-defined
  missing values, which `haven::read_sav(user_na = TRUE)` keeps as
  values of a `haven_labelled_spss` column, are missing, as SPSS treats
  them.

- ci_method:

  Method for constructing the confidence interval; see
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  for the full list. Defaults to `"feldt"`.

- conf_level:

  Confidence level for the interval (1 - Type I error rate). Defaults to
  `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected,
  a single whole number. Defaults to `10000`. The fewest accepted
  depends on the interval method and the level: at the default
  `conf_level`, 400 for `"percentile"` (10 replications beyond each
  limit), 2000 for `"bca"` (50 beyond each nominal limit, because the
  adjusted limits sit farther into the tails), and 100 for
  `"bootstrap_se"` and `"bootstrap_se_logistic"` (a standard deviation
  of the replications); at 0.99 the first two are 2000 and 10000. A
  smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). The count is held to the number requested.
  Replications set aside afterward are counted and reported in a warning
  of class `dmar_bootstrap_replications_dropped`, which names the
  function and says to request more when those that return a value are
  fewer than the number required; when fewer than half of the `B`
  replications return a value, the interval is withheld (limits and
  bootstrap standard error `NA`) under a warning of class
  `dmar_bootstrap_insufficient`. See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  which also describes `options(dmar.bootstrap_exploratory = TRUE)` for
  a first look with a smaller `B`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, under which the resamples are drawn from, and advance, the
  current random number stream; supply an integer for a reproducible
  interval, which is set for the call only, the caller's random number
  state being restored when it returns.

- ncpus:

  The number of processes over which the bootstrap replications are
  spread. The default `1` runs them in this R session; a larger number
  uses forked workers (boot's `parallel = "multicore"`) or, on Windows,
  which cannot fork, a socket cluster (`"snow"`), whose workers load the
  installed DMAR. boot draws every resample before any replication runs,
  so a given `seed` returns the same interval, counts, and replicates
  whatever `ncpus` is. A session that forks workers and also starts
  background processes through callr or processx needs the environment
  variable `PROCESSX_NOTIFY_OLD_SIGCHLD` set before processx is loaded
  (see
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)).
  Ignored unless a bootstrap method is selected.

- progress:

  `NULL` (the default), or a function of two arguments that the
  bootstrap calls as it runs, so that an application can show how far a
  long bootstrap has gone and let its user stop it (a Shiny app, for
  example). It is called as `progress(done, total)`, with the number of
  replications attempted so far and `B`: after each replication when
  `ncpus = 1`, and with more processes after each block of replications
  they complete together, about a fiftieth of `B`
  (`ncpus * ceiling(B / (50 * ncpus))` replications, the last block
  possibly fewer), so that `done` advances by the block. A BCa interval
  takes its acceleration from a regression on the resampling frequencies
  when at least as many replications return a value as there are cases,
  and from the jackknife, one leave-one-out refit per case, when fewer
  do or when that regression fails or gives influence values that are
  not all finite; the jackknife runs at most once in a call, and not at
  all when no interval is formed or every returned replicate has the
  same value. While it runs, the callback is called again after each
  refit (with more processes, after each block of
  `ncpus * ceiling(n / (50 * ncpus))` refits, `n` the number of cases
  resampled), with `done` and `total` both `B`, so it can be called more
  than `B` times, and `done` equal to `total` does not mean that the
  bootstrap has finished. If the callback returns `FALSE`, the bootstrap
  stops before the next replication or block with an error of class
  `dmar_bootstrap_cancelled`, whose elements `done` and `total` say how
  far it got (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md));
  any other value is ignored. An error inside the callback stops the
  call with an error of class `dmar_progress_error`, whose message names
  the callback and which carries `done`, `total`, and, as `parent`, the
  callback's own condition with its class intact, so that an application
  can tell a failing callback from a failing analysis. The callback runs
  outside the random number stream, so a given `seed` returns the same
  result with or without it, and a callback that draws random numbers
  leaves the caller's stream as it found it. With `NULL` nothing is
  called. Ignored unless a bootstrap method is selected.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample KR-20), `"se"` (standard error on the coefficient scale, `NA`
for methods that do not produce one; for the transformation-based
intervals `"fisher"`, `"bonett"`, and `"hakstian_whalen"` it is the
delta method back-transform of the transformation-scale standard error),
`"se_transformed"` (only for those transformation-based intervals: the
standard error on the transformation scale, with the scale named in the
attribute `se_transform_scale`: `"fisher_z"`, `"log(1-alpha)"`, or
`"cube_root"`), `"ci_lower"` and `"ci_upper"` (each clamped into \[0,
1\], so the limits stay in order), `"conf_level"`, `"N"` (effective
sample size after listwise deletion), `"N_complete"` (the complete
cases; equal to `"N"` here, and carried so the whole reliability family
returns one shape), and `"n_items"` (number of items). Attributes
`coefficient` (`"kr20"`) and `ci_method` record the computation;
bootstrap calls also record `B`, the replication counts in the
`bootstrap` attribute (see
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md);
KR-20 fits no model, so no replication is improper), and the number of
replications the interval uses in the attribute `B_used`, with
`n_dropped`, `B_required`, and `bootstrap_status` (see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).

## Details

Kuder and Richardson's (1937) formula 20 for a \\J\\-item composite of
binary items is \$\$KR\_{20} = \frac{J}{J-1}\left(1 - \frac{\sum\_{j}
p\_{j} q\_{j}}{s\_{Y}^{2}}\right),\$\$ where \\p\_{j}\\ is the
proportion of respondents endorsing item *j* (i.e., scoring 1), \\q\_{j}
= 1 - p\_{j}\\, and \\s\_{Y}^{2}\\ is the variance of the composite
score. For dichotomous items \\p\_{j} q\_{j}\\ is the item variance, so
KR-20 is algebraically identical to coefficient \\\alpha\\ (Guttman,
1945; Cronbach, 1951) computed from the item covariance matrix. KR-20
predates coefficient \\\alpha\\ by 14 years and Feldt's (1965)
*F*-distribution interval was derived specifically for the sampling
distribution of KR-20.

Because KR-20 is a special case of coefficient \\\alpha\\ (Guttman,
1945; Cronbach, 1951), the same considerations apply: KR-20 equals the
population reliability of the composite under essential
\\\tau\\-equivalence (i.e., equal factor loadings) and serves as a lower
bound otherwise. For dichotomous items see also
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
(categorical omega), which treats the items via a probit-link
single-factor model and does not assume equal loadings.

KR-20, like \\\alpha\\, is at most one but can fall below zero, when the
item covariances are negative on balance; on dichotomous items that can
arise when an item is scored in the opposite direction to the others (a
reversed key). A reliability coefficient cannot be negative, so a
negative estimate is returned as computed with a warning that says so,
and each confidence limit is clamped into \[0, 1\]; the interval around
an estimate far enough below zero is then \[0, 0\], which says that the
data support no positive reliability for the items as keyed.

Available confidence interval methods (set via `ci_method`) are the same
as for
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md);
see that function's *Details* for full descriptions. The default for
KR-20 is `"feldt"`, the *F*-distribution interval originally developed
for KR-20.

From two complete rows KR-20 is computed with Feldt's interval, which is
defined for any *N* of at least 2 (every item must differ between the
two rows, since an item they answer alike does not vary, and the two
totals must differ, since a total that does not vary has no
reliability). Feldt's *F* reference is exact under multivariate
normality and parallel items, which 0/1 items cannot meet, and for them
it is the approximation Feldt (1965) derived. The intervals that are not
defined there (`"bonett"`, `"fisher"`, `"ml"`, `"adf"`, their
`_logistic` variants, and the bootstrap) are refused, as
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
explains, and `"hakstian_whalen"`, which is defined, has an upper limit
of one there whatever the data, at any confidence level above about .90.

The menu also includes the nonparametric bootstrap intervals
`"percentile"`, `"bca"`, `"bootstrap_se"`, and
`"bootstrap_se_logistic"`, which resample the rows of `data` with
replacement `B` times and recompute KR-20 on each replication (Efron &
Tibshirani, 1993). They are worth the cost when the closed forms are
least trustworthy, which for dichotomous items means highly unbalanced
item difficulties or a sample size too small for the normal-theory
derivations behind `"feldt"` and `"bonett"`. No bootstrap runs unless
`ci_method` asks for one; when it does, the default is `B = 10000`
replications, and supplying `seed` makes the interval reproducible.
Replications on which KR-20 cannot be computed are set aside and
counted, among them a replication in which an item takes a single value,
refused by name as the data would be: a resample leaves out a given
respondent about 37 percent of the time (\\(1 - 1/N)^N\\, near
\\e^{-1}\\), so an item passed by one respondent loses its variation in
about a third of the resamples. The interval is computed from the
replications that return a value; the counts travel in the `bootstrap`
attribute described under
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
print under the table, and are stated in a warning whenever a
replication was set aside.

**Comparison with other packages.** The psych package computes the same
quantity via `alpha` (since \\\alpha\\ on 0/1 data *is* KR-20).
`reliability_kr20` restricts input to raw 0/1 data so the
dichotomous-items assumption cannot be quietly violated, presents the
historical formula in the documentation, and accompanies the point
estimate with a confidence interval drawn from the methods compared in
Kelley and Pornprasertmanit (2016).

## References

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297–334.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Feldt, L. S. (1965). The approximate sampling distribution of
Kuder-Richardson reliability coefficient twenty. *Psychometrika, 30*,
357–370.

Guttman, L. (1945). A basis for analyzing test-retest reliability.
*Psychometrika, 10*(4), 255–282.

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Kuder, G. F., & Richardson, M. W. (1937). The theory of the estimation
of test reliability. *Psychometrika, 2*, 151–160.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
`alpha`.

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Ten dichotomous items with a single underlying ability.
N <- 300
J <- 10
ability <- rnorm(N)
loadings <- rep(0.6, J)
latent <- outer(ability, loadings) +
          matrix(rnorm(N * J, sd = sqrt(1 - 0.6^2)), N, J)
items <- (latent > 0) * 1
colnames(items) <- paste0("y", seq_len(J))

reliability_kr20(data = items)
#>  term       value
#>  estimate   0.774
#>  se         <NA> 
#>  ci_lower   0.734
#>  ci_upper   0.810
#>  conf_level 0.950
#>  N          300  
#>  N_complete 300  
#>  n_items    10   
reliability_kr20(data = items, ci_method = "bonett")
#>  term           value
#>  estimate       0.774
#>  se             0.020
#>  se_transformed 0.086
#>  ci_lower       0.732
#>  ci_upper       0.809
#>  conf_level     0.950
#>  N              300  
#>  N_complete     300  
#>  n_items        10   

# A bootstrap interval recomputes KR-20 on each of B resamples of the
# rows. B = 500 keeps the example quick; a reported interval deserves
# the default B = 10000, and seed makes the interval reproducible.
reliability_kr20(data = items, ci_method = "percentile", B = 500,
                 seed = 113)
#>  term       value
#>  estimate   0.774
#>  se         0.019
#>  ci_lower   0.739
#>  ci_upper   0.810
#>  conf_level 0.950
#>  N          300  
#>  N_complete 300  
#>  n_items    10   
#> 
#> Bootstrap: 500 replications attempted, 500 returned a value
```
