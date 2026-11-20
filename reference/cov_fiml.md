# Maximum Likelihood Covariance Matrix and Means of Incomplete Data

For continuous variables some of whose values are missing, `cov_fiml()`
returns the maximum likelihood estimates of their covariance matrix and
mean vector under the multivariate normal model, computed from every
observed value rather than from the complete cases alone. These are the
estimates of the saturated (unrestricted) model fitted by full
information maximum likelihood (FIML), the model in which every
variance, covariance, and mean is free. They are computed by the
expectation maximization (EM) algorithm of Dempster, Laird, and Rubin
(1977) and are often called the EM covariance matrix and the EM means
(Enders & Peugh, 2004). The function requires lavaan, whose EM
implementation does the computing.

The usual alternatives each give something up. Listwise deletion keeps
only the complete cases, and its estimates are consistent in general
only when values are missing completely at random (MCAR). Pairwise
deletion estimates each covariance from a different subset of cases, so
the matrix need not be positive definite and has no single sample size
(Marsh, 1998). The EM estimates use every case and are consistent when
values are missing at random (MAR) given the variables in the model
(Schafer & Graham, 2002). They are therefore the matrix to hand to a
method that starts from a covariance or correlation matrix, such as an
exploratory factor analysis or a parallel analysis of incomplete items.

## Usage

``` r
cov_fiml(data, aux = NULL)
```

## Arguments

- data:

  A `data.frame` or numeric matrix of continuous variables, with missing
  values coded `NA`. Every column that `aux` does not name is an
  analysis variable, a row and column of the returned `S`. Every column
  must be numeric; a factor, character, or logical column is refused, as
  is a column with fewer than two distinct observed values.

- aux:

  Optional character vector naming auxiliary columns of `data`:
  variables that are not wanted in `S` but carry information about the
  missing values, because they predict which values are missing, the
  missing values themselves, or both. They enter the saturated model
  with the analysis variables, so their information reaches `S` and
  `mu`, and are then left out of what is returned (the saturated
  correlates treatment of Graham, 2003, the one the reliability family's
  `aux` argument uses).

## Value

A list with the same six elements whatever the data. A matrix is what
every use of the result takes
([`stats::cov2cor()`](https://rdrr.io/r/stats/cor.html),
`stats::factanal(covmat = )`, a function's `S` argument), so the result
holds the matrix itself, as
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) and
[`covmat_from_cfa`](https://yelleknek.github.io/DMAR/reference/covmat_from_cfa.md)
return theirs, rather than a table of its elements:

- `S`:

  The maximum likelihood covariance matrix of the analysis variables
  (divisor `N`), with their names as row and column names, in the order
  of the columns of `data`.

- `mu`:

  The maximum likelihood mean vector of the analysis variables, named,
  in the same order.

- `N`:

  The number of cases used: those with at least one observed value among
  the analysis variables and the auxiliaries.

- `N_complete`:

  The number of cases observed on every analysis variable.

- `patterns`:

  The number of distinct patterns of observed and missing values among
  the cases used, over every variable in the saturated model (the
  analysis variables and the auxiliaries);
  [`missing_patterns`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)
  tabulates them.

- `converged`:

  `TRUE` when the EM algorithm met its convergence criterion.

## Details

**What is estimated.** With complete data the estimates are the sample
moments with the maximum likelihood divisor: `S` equals
`cov(data) * (N - 1) / N` and `mu` equals `colMeans(data)`. With missing
values, the EM algorithm alternates two steps until the estimates stop
changing. The expectation step fills in the sums and cross products each
case would contribute, replacing a missing value by its regression
prediction from the observed values of its case under the current
estimates and adding the residual covariance to the cross products; the
maximization step recomputes the moments from those sums. No iteration
lowers the likelihood of the observed data. Meeting the stopping
criterion establishes numerical stability of the iterations, not a
finite global maximum or uniqueness. The estimates are those of the
unrestricted (h1) model of every FIML fit in lavaan:
`lavaan::lavInspect(fit, "sampstat.h1")` returns them for a model fitted
with `missing = "ml"` to the same variables. `cov_fiml()` computes them
with lavaan's EM
([`lavaan::lavCor()`](https://rdrr.io/pkg/lavaan/man/lavCor.html) with
`missing = "ml"`) and its default convergence criterion, which stops the
iterations when no mean or covariance changes by more than \\10^{-5}\\.
That is a criterion on the absolute change, and in the units of
variables whose standard deviations are near .001 it is met after an
iteration or two, far from the maximum. `cov_fiml()` therefore runs the
algorithm on the variables centered at their means and divided by their
standard deviations, each computed from the variable's observed values,
and transforms the estimates back, so that they do not depend on the
units: multiplying a variable by a constant multiplies its mean, and its
row and column of `S`, by that constant. They agree with lavaan's own
estimates of the same moments to the precision of the criterion, about
seven significant digits. When the algorithm stops at its iteration
limit instead, `converged` is `FALSE`, a warning says so, and `S` and
`mu` are the last iterate rather than the maximum likelihood estimates.

**When the estimates are consistent.** The estimates are consistent when
the values are MAR given the variables in the saturated model: the
probability that a value is missing may depend on the observed values of
the variables included, but not on the missing values themselves
(Schafer & Graham, 2002). MAR is defined relative to the variables
analyzed (Enders, 2001). If whether an item is missing depends on a
variable that is left out, MAR need not hold for the items alone, and
naming that variable in `aux` brings it back into the model without
making it part of `S`. Including such variables (the inclusive strategy)
costs little and can remove a bias that leaving them out would leave in
place (Collins, Schafer, & Kam, 2001). Under missingness that depends on
the missing values themselves (not at random), no choice of auxiliaries
guarantees consistency. Normality enters as the working likelihood:
Yuan, Yang-Wallentin, and Bentler (2012) summarize results showing that
the normal theory estimates can remain consistent when the population
distribution is unknown, and their simulations found the normal theory
variance estimates biased at smaller sample sizes when the data were not
normal, a bias that shrank as the sample grew.

**Cases used.** A case on which every variable, analysis and auxiliary,
is missing carries no information and is set aside, as lavaan would set
it aside. Every other case is used, including a case whose analysis
variables are all missing but whose auxiliaries are observed. Under MAR
given an auxiliary, such cases carry the auxiliary values of exactly the
people whose analysis variables are missing, which is how the model
learns how those people differ, and setting them aside would pull the
estimates back toward those of the cases that responded. `N` counts the
cases used and `N_complete` the cases complete on the analysis
variables, all that listwise deletion would keep.

**Refusals.** Two variables that are never observed together have a
covariance the data say nothing about: the likelihood is flat in that
direction, and an EM estimate of it reflects the starting values and the
other covariances rather than the data. The same holds, in substance,
for two variables observed together only in rows where one of them takes
a single value (or in one row only): those rows say nothing about how
the two vary together, and the estimate of their covariance rests on the
form of the normal model, or is not determined at all when one of them
is observed nowhere else. The function refuses either kind of pair and
names it;
[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
reports the number of cases behind every pair.

**Using the matrix downstream.** `N` is the number of cases used, not an
effective sample size, and the matrix carries less information than `N`
complete cases would. A method that takes a covariance matrix and a
sample size and treats them as complete data gets consistent point
estimates from `S`, but its standard errors and test statistics are too
optimistic, because they ignore the information the missing values took
away (Savalei & Bentler, 2009). Enders and Peugh (2004) studied which
sample size to give such a method when fitting a confirmatory factor
model: the smallest number of cases behind any covariance gave Type I
error rates near the nominal level when every indicator had missing
values (the smallest number behind any variance did when only some
indicators had them), and the harmonic mean of the numbers behind the
variances gave confidence intervals with coverage near the nominal
level.
[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
reports the number of cases behind every variance and covariance. The
two-stage approach of Savalei and Bentler (2009) corrects the standard
errors and the test statistic of a model fitted to the EM moments
(lavaan's `missing = "two.stage"`), and a model fitted to the raw data
by FIML needs no correction. The fraction of missing information
([`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md))
says what the missing values cost each parameter of such a fit. For
point estimates, such as the eigenvalues of a factor analysis or of a
parallel analysis, the EM matrix is the consistent basis.

**Standard errors.** The function returns point estimates only. The
sampling variability that matters is that of whatever is computed from
the matrix, and the analysis that computes it is where its interval
belongs.

## References

Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*(4), 330–351.
[doi:10.1037/1082-989X.6.4.330](https://doi.org/10.1037/1082-989X.6.4.330)

Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood
from incomplete data via the EM algorithm. *Journal of the Royal
Statistical Society: Series B (Methodological), 39*(1), 1–38.

Enders, C. K. (2001). A primer on maximum likelihood algorithms
available for use with missing data. *Structural Equation Modeling,
8*(1), 128–141.

Enders, C. K., & Peugh, J. L. (2004). Using an EM covariance matrix to
estimate structural equation models with missing data: Choosing an
adjusted sample size to improve the accuracy of inferences. *Structural
Equation Modeling, 11*(1), 1–19.

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Marsh, H. W. (1998). Pairwise deletion for missing data in structural
equation models: Nonpositive definite matrices, parameter estimates,
goodness of fit, and adjusted sample sizes. *Structural Equation
Modeling, 5*(1), 22–36.
[doi:10.1080/10705519809540087](https://doi.org/10.1080/10705519809540087)

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*(3), 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the
state of the art. *Psychological Methods, 7*(2), 147–177.
[doi:10.1037/1082-989X.7.2.147](https://doi.org/10.1037/1082-989X.7.2.147)

Yuan, K.-H., Yang-Wallentin, F., & Bentler, P. M. (2012). ML versus MI
for missing data with violation of distribution conditions.
*Sociological Methods & Research, 41*(4), 598–629.
[doi:10.1177/0049124112460373](https://doi.org/10.1177/0049124112460373)

## See also

[`missing_patterns`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)
and
[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
to describe the missing values first;
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
for what they cost the parameters of a FIML fit;
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
and
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
whose `missing = "fiml"` paths also bring auxiliary variables into a
saturated model.

Other missing data:
[`df_barnard_rubin()`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md),
[`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
[`fraction_missing_information()`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
[`missing_coverage()`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
[`missing_patterns()`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The six spatial tests of the Holzinger and Swineford data. The revised
# paper form board (t25) and flags (t26) tests were given only at
# Grant-White, so the 156 Pasteur pupils are missing both of them.
spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
em <- cov_fiml(holzinger_swineford[spatial])
em$N
#> [1] 301
em$N_complete
#> [1] 145
em$patterns
#> [1] 2
round(stats::cov2cor(em$S), 2)
#>                        t1_visual_perception t2_cubes t3_paper_form_board
#> t1_visual_perception                   1.00     0.30                0.37
#> t2_cubes                               0.30     1.00                0.24
#> t3_paper_form_board                    0.37     0.24                1.00
#> t4_lozenges                            0.44     0.34                0.31
#> t25_paper_form_board_r                 0.40     0.33                0.41
#> t26_flags                              0.47     0.21                0.31
#>                        t4_lozenges t25_paper_form_board_r t26_flags
#> t1_visual_perception          0.44                   0.40      0.47
#> t2_cubes                      0.34                   0.33      0.21
#> t3_paper_form_board           0.31                   0.41      0.31
#> t4_lozenges                   1.00                   0.31      0.49
#> t25_paper_form_board_r        0.31                   1.00      0.30
#> t26_flags                     0.49                   0.30      1.00

# Listwise deletion would keep the 145 Grant-White pupils alone and
# estimate every correlation from them.
round(stats::cor(holzinger_swineford[spatial], use = "complete.obs"), 2)
#>                        t1_visual_perception t2_cubes t3_paper_form_board
#> t1_visual_perception                   1.00     0.33                0.37
#> t2_cubes                               0.33     1.00                0.19
#> t3_paper_form_board                    0.37     0.19                1.00
#> t4_lozenges                            0.45     0.42                0.37
#> t25_paper_form_board_r                 0.40     0.31                0.40
#> t26_flags                              0.47     0.23                0.33
#>                        t4_lozenges t25_paper_form_board_r t26_flags
#> t1_visual_perception          0.45                   0.40      0.47
#> t2_cubes                      0.42                   0.31      0.23
#> t3_paper_form_board           0.37                   0.40      0.33
#> t4_lozenges                   1.00                   0.33      0.47
#> t25_paper_form_board_r        0.33                   1.00      0.30
#> t26_flags                     0.47                   0.30      1.00

# Values missing at random given a variable that is not analyzed: word
# meaning (t9) goes missing more often for pupils who scored high on
# general information (t5). The two visual tests say little about those
# pupils, so without t5 the mean of t9 stays near that of the pupils who
# took it; named as an auxiliary, t5 recovers most of the difference.
set.seed(113)
d <- holzinger_swineford[c("t1_visual_perception", "t4_lozenges",
                           "t9_word_meaning", "t5_general_information")]
z <- as.numeric(scale(d$t5_general_information))
p_missing <- stats::plogis(-0.5 + 2 * z)
d$t9_word_meaning[stats::runif(nrow(d)) < p_missing] <- NA
without_aux <- cov_fiml(d[1:3])
with_aux <- cov_fiml(d, aux = "t5_general_information")
rbind(without_aux = without_aux$mu, with_aux = with_aux$mu,
      full_data = colMeans(holzinger_swineford[names(d)[1:3]]))
#>             t1_visual_perception t4_lozenges t9_word_meaning
#> without_aux             29.61462    18.00332        12.42404
#> with_aux                29.61462    18.00332        14.66152
#> full_data               29.61462    18.00332        15.29900
```
