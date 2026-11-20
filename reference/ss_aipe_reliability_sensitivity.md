# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Reliability Coefficient

Quantifies how much misspecification of the population reliability
coefficient distorts an AIPE-based sample size plan for the reliability
of a composite score. On each replication the function simulates an *n*
\\\times\\ *J* data matrix of subjects by items from a single-factor
model of parallel tests whose population reliability of the sum score
equals `true_reliability`, fits the requested estimator (alpha or omega)
via the corresponding `reliability_*` function with the supplied
`ci_method`, and records the realized reliability estimate and its
confidence interval.

**Population model.** Each item has a single common-factor loading and
uncorrelated unique error. With per-item variance normalized to 1, the
loading and unique variance are chosen so the reliability of the sum
score equals `true_reliability`: \$\$\lambda^2 \\=\\ \frac{\rho}{J(1 -
\rho) + \rho}, \qquad \psi^2 \\=\\ 1 - \lambda^2,\$\$ where \\\rho =
\\`true_reliability` and \\J\\ is the number of items, `n_items`. Item
scores are \\y\_{ij} = \lambda T_i + e\_{ij}\\, with \\T_i \sim N(0,
1)\\ and \\e\_{ij} \sim N(0, \psi^2)\\.

**Missing data.** Two arguments carry a plan into incomplete data, and
they stand to each other as `supposed_reliability` stands to
`true_reliability`. `fmi` is the planning value: the fraction of missing
information the plan assumes for the coefficient, passed to
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
which enlarges the sample size accordingly. `missing_prop` is the data
generating mechanism: each simulated item score is deleted independently
with that probability (missing completely at random), and every
replication is analyzed by full information maximum likelihood. A plan
made with `fmi` is then judged by the widths the incomplete data
actually deliver, and a plan that ignored the missing data (`fmi = 0`)
is seen to fall short. Planning with `fmi` while simulating complete
data shows the opposite misspecification, widths narrower than planned.

## Usage

``` r
ss_aipe_reliability_sensitivity(
  true_reliability = NULL,
  supposed_reliability = NULL,
  n_items,
  width,
  specified_N = NULL,
  estimator = c("alpha", "omega"),
  ci_method = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  fmi = 0,
  missing_prop = 0,
  B = NULL,
  seed = NULL
)
```

## Arguments

- true_reliability:

  Population reliability coefficient (in \\\[0, 1)\\).

- supposed_reliability:

  Reliability supposed when the study is planned, a planning assumption
  rather than an estimate from data; the function passes the factor
  loading and residual variance it implies to
  [`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
  as `true_factor_loading` and `true_residual_var`.

- n_items:

  Number of items in the composite (\\J\\ in the formulas): a whole
  number of at least 2, and at least 3 for `estimator = "omega"`.

- width:

  Desired full width of the CI on reliability.

- specified_N:

  Sample size to evaluate (incompatible with `supposed_reliability`).

- estimator:

  One of `"alpha"` (default; coefficient alpha via
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md))
  or `"omega"` (composite reliability via
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).
  For a population of parallel tests the two coincide; differences in
  sample estimates reflect estimator-specific finite-sample bias and CI
  behavior.

- ci_method:

  CI method passed to the estimator. The default `NULL` is `"ml"` for
  alpha and omega, with complete or incomplete data: the normal theory
  Wald interval the plan is made for (see *Details*). Any other interval
  the estimator offers can be named.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to the planner, whose Monte
  Carlo search then runs at its defaults (`G_initial = 500`, `G = 5000`)
  and draws from the current random number stream, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes the plan and the simulation reproducible together. The plan is
  the one the planner reports for the parallel model with the normal
  theory interval, which the default `ci_method` evaluates;
  `specified_N` evaluates any other size.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample reliability estimate, the two confidence limits, the
  interval width, and two indicators of whether the interval missed
  `true_reliability` below or above, followed, when `missing_prop` is
  above 0, by the fraction of missing information the replication
  realized for the coefficient): nothing is written when `filename` is
  `NULL` (the default), a new file with a header row is created
  otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`.

- fmi:

  The fraction of missing information the plan assumes for the
  coefficient, passed to
  [`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
  when `supposed_reliability` resolves the sample size: a single number
  in \\\[0, 1)\\, with the default `0` for a complete data plan. With
  `specified_N` there is no plan for it to act on, and a value above 0
  is an error.

- missing_prop:

  The probability with which each simulated item score is missing,
  completely at random: a single number in \\\[0, 1)\\. The default `0`
  simulates complete data; above 0, every replication is analyzed by
  full information maximum likelihood (see *Details*).

- B:

  Number of bootstrap replications inside each simulated study when
  `ci_method` is a bootstrap method. The default, `NULL`, is the larger
  of 200 and one and a quarter times the number the interval method and
  the level require (see
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)):
  500 for `"percentile"`, 2500 for `"bca"`, and 200 for the two
  bootstrap standard error methods at the 95 percent level. The quarter
  is room for replications that return no value (a resample without a
  one-factor solution): the required number (400, 2000, 100) is held to
  the request, so a simulated study that loses a few replications still
  reports its interval, and with the quarter its limits still rest on
  about the required number. A simulated study in which more than half
  of the replications return no value reports no interval. Every
  simulated study runs `B` replications, so the run time is proportional
  to `G` times `B`. Not used by the other interval methods.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized
reliability and CI width, the proportion of intervals at or below
`width`, tail-specific and overall non-coverage of `true_reliability`,
and the input echoes `total_N`, `n_items`, `true_reliability`,
`supposed_reliability`, `width`, `conf_level`, and `assurance` (present
only when an assurance was supplied). The proportion of intervals at or
below `width`, the widths, and the noncoverage proportions are over the
replications that returned an interval. The `simulation` attribute
records the numbers of replications attempted, with an interval, and
without one (see the section on the record of the simulation), and,
here, the message of the first failure (`first_failure`). When `fmi` or
`missing_prop` is above 0, the rows are unchanged and the missing data
inputs travel as attributes: `fmi` (`NA` with `specified_N`),
`missing_prop`, and, when `fmi` is above 0, `necessary_N_complete`, the
planner's complete data sample size. When `missing_prop` is above 0, two
more attributes record the fraction of missing information the simulated
data realized for the coefficient: `mean_realized_fmi`, its mean over
the replications in which it was computed, and `n_realized_fmi`, the
number of those replications. With `fmi` and `missing_prop` both at 0 no
attribute is added. With a bootstrap `ci_method` the result also carries
the attributes `B`, `B_used`, `n_dropped`, `B_required`, and
`bootstrap_status` of
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md);
they record the request each simulated study ran (`B_used` equals `B`
and `n_dropped` is 0, because replications dropped inside a simulated
study are not totaled), and `bootstrap_status` is `"exploratory"` when
`B` is below `B_required`.

## Details

**The interval evaluated.** An AIPE plan is made for a particular
confidence interval procedure, so a check of the plan evaluates that
procedure. The plan here is the one
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
makes for the parallel model with the normal theory interval: its closed
form takes the sampling variance of alpha to be the compound symmetry
form \\2 (1 - \rho)^2 J / ((J - 1)(N - 1))\\ of van Zyl, Neudecker, and
Nel (2000), and its assurance search forms, on each simulated sample,
the Wald interval with that variance at the sample alpha. The default
`ci_method = "ml"` evaluates the interval of
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)`(ci_method = "ml")`,
the Wald interval with the general variance of van Zyl et al. at the
sample covariance matrix, which at a population of parallel items is the
same variance; for omega it evaluates the interval of
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)`(denominator = "model_implied", ci_method = "ml")`,
the delta method Wald interval of the model implied omega, whose
variance at a population of parallel items is that same compound
symmetry form (it is the variance
[`ss_aipe_reliability()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
plans with for the congeneric model, and the two agree to machine
precision there). The plan does not describe other intervals, among them
Bonett's (2002) for alpha and `"mlr"` for omega. With four items, a true
and planned reliability of .80, a width of .15, and an assurance of .80,
the plan of 99 cases gave alpha intervals no wider than .15 in 81.3
percent of 2,000 replications with `"ml"` and in 76.75 percent with
Bonett's interval (`set.seed(113)` before each call); for omega the two
intervals came closer (83.3 percent of 1,000 with `"ml"` and 81.9 with
`"mlr"`). Any other interval the estimator offers can still be named in
`ci_method`, which asks how the planned *N* serves an interval other
than the one it was planned for.

**The analysis under missing data.** With `missing_prop` above 0,
coefficient alpha is
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)`(missing = "fiml")`,
the classical coefficient computed from the full information maximum
likelihood estimate of the item covariance matrix, whose `"ml"`
interval, the default here as with complete data, takes its standard
error from the FIML information matrix of the saturated model.
Coefficient omega is
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)`(missing = "fiml")`,
whose `"ml"` interval takes its standard error from the FIML information
matrix of the single-factor model. These are the Wald intervals the
planner's `fmi` adjustment plans for, with the standard error the
missing data inflate. The complete data closed forms (Bonett, Feldt,
Fisher, Hakstian-Whalen) and the asymptotic distribution-free intervals
do not apply to partially observed rows, so, with `missing_prop` above
0, `ci_method` must be `"ml"` or `"ml_logistic"` for alpha and one of
`"ml"`, `"mlr"`, `"ml_logistic"`, or `"mlr_logistic"` for omega. A row
that loses every item score carries no information and is dropped, as
lavaan drops it.

**Choosing `fmi` for a given `missing_prop`.** The two are not the same
number. The FMI of the coefficient depends on the number of items and
how strongly they are correlated, because information lost on one item
is partly recovered from the others; under data missing completely at
random it is bounded by the largest rates of univariate and pairwise
missingness (Savalei & Rhemtulla, 2012). For the population of the
second example (four parallel items, reliability .80, width .15),
`missing_prop = .2` gives coefficient alpha and coefficient omega an FMI
of .184 relative to every case, the footing of the planner's `fmi` (from
one FIML fit to 200,000 simulated cases; it is .182 relative to the
cases the fit keeps, since a share \\.2^4 = .0016\\ of them loses every
item). In 1,000 replications with the default `ci_method = "ml"`
(`set.seed(113)` before each call), the complete data plan (\\N = 74\\)
gave alpha intervals that averaged .153 wide with complete data and .170
with `missing_prop = .2`, wider by close to the width inflation factor
\\1 / \sqrt{1 - .184} = 1.107\\, while the plan made with `fmi = .184`
(\\N = 91\\) averaged .152 with `missing_prop = .2`. The help page of
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
describes where a planning FMI comes from in practice.

**The fraction of missing information the data realized.** With
`missing_prop` above 0, every replication that returns an interval also
estimates the FMI of its coefficient from the same simulated data, as
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
does with its default unstructured estimator: the FIML fit behind the
coefficient is refit at its estimates with the information of the
observed data and with the information complete data would have carried,
and the FMI is one minus the ratio of the coefficient's complete data
sampling variance to its sampling variance with the data as observed.
For alpha the fit is the saturated model of the items whose covariance
matrix the classical coefficient is computed from, and the three
estimators of Chen and Savalei (2021) coincide for a saturated model (p.
9); for omega it is the single-factor model, and the FMI is that of its
model implied coefficient. The FMI is taken relative to all *n*
simulated cases, the *N* of the plan. A case that loses every item score
is set aside by the fit, and the FMI relative to the \\n\_{used}\\ cases
the fit used, \\\gamma\_{used}\\ (the FMI that
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
reports), leaves out the information those cases would have carried; it
is folded in as \\1 - (n\_{used} / n)(1 - \gamma\_{used})\\. The
difference matters with few items and much missing data: a share \\p^J\\
of the cases loses every item, with *p* the value of `missing_prop`,
which is .125 for three items at `missing_prop = .5`. The mean over the
replications is the `mean_realized_fmi` attribute and is printed under
the table, so the `fmi` a plan assumed can be read beside the FMI that
the simulated missing data cost the coefficient. A single sample FMI is
a poor guide: in the simulations of Chen and Savalei (2021, p. 9) the
sample FMI of a regression coefficient was largely unbiased but varied
from sample to sample, with root mean squared errors of about .10 or
more at *N* = 50 and 100. The mean over many replications is the stable
figure, and it is the value to carry into `fmi` for a study whose
missing data the simulation mimics. A replication whose FMI cannot be
computed (a fit that fails, or an improper solution that lavaan will not
refit at its estimates; see
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md))
is left out of the mean, and the `n_realized_fmi` attribute counts the
replications that entered it. The computation adds a model fit and two
refits to every replication, about doubling the time a replication
takes.

**Bootstrap intervals.** When `ci_method` is a bootstrap method
(`"bootstrap_se"`, `"bootstrap_se_logistic"`, `"percentile"`, or
`"bca"`), every one of the `G` replications runs its own bootstrap of
`B` resamples, so the simulation computes the coefficient \\G \times B\\
times. The default `B` is the larger of 200 and one and a quarter times
the number of replications the interval method and the level require
(see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)):
at the 95 percent level, 500 for `"percentile"`, 2500 for `"bca"`, and
200 for the two bootstrap standard error methods. The quarter leaves
room for replications that return no value: the number required is held
to the request, so a simulated study that loses a few replications still
reports its interval, and with the quarter its limits still rest on
about the number required. That default is far below the `B = 10000` of
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
and
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
and keeps the product manageable. It also makes each replication's
limits less precise than the reported interval's would be: a percentile
limit at `B = 500` rests on about 12 resamples in each tail, and the
Monte Carlo error of the limits adds to the spread of the realized
widths. A check of a bootstrap interval that will be reported with a
large `B` is closer to that interval with a larger `B` here, at a cost
that grows in proportion. A `B` below the required number is refused
unless `options(dmar.bootstrap_exploratory = TRUE)` is set; the run is
then marked by one warning of class `dmar_bootstrap_exploratory` for the
call, not one for each simulated study.

**The omega coefficient.** Omega is computed with
`denominator = "model_implied"`, the coefficient whose delta method
interval the planner evaluates and the one the Wald intervals are
defined for. For the parallel population simulated here the model
implied and observed denominators estimate the same value.

**Coverage at a true reliability of zero.** Coverage and the misses in
each tail are counted on the interval the estimation function reports,
whose limits are each clamped into \[0, 1\] (see
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)),
since that is the interval a plan is made for and a study reports.
Clamping can only raise the coverage of a value in \[0, 1\], and at the
boundary it does: with `true_reliability = 0`, an alpha interval lying
wholly below zero is reported as \[0, 0\], which contains zero, so
`pct_ci_miss_high` is zero by construction and `total_type_I_error` is
the rate of misses low. A miss high there would place the reliability
below zero, outside the range a reliability coefficient can take. With
four items, \\N = 40\\, and Bonett's interval (400 replications,
`set.seed(1)`) the intervals miss .0075 of the time, all of it low. For
any true reliability above zero, an interval lying wholly below zero
misses high.

**Replications that fail.** A replication that returns no interval (a
fit that did not converge, which is more common under full information
maximum likelihood with small samples and much missing data, a fit whose
likelihood has no maximum, which the estimation functions refuse with an
error of class `dmar_no_maximum` when an error variance runs off below
zero while the likelihood is flat or still rising, or an interval the
method could not form) is set aside: it counts in no summary, and one
warning states how many were set aside and the first reason. When no
replication returns an interval, the function stops with that reason
rather than return a table of missing values.

## The record of the simulation

Every sensitivity function of the package summarizes its replications
over one denominator: the replications that produced a confidence
interval (for a power sensitivity function, a test). The proportion of
intervals that met the width target (`pct_ci_less_width`), the
proportions of intervals that missed the true value, and the summaries
of the widths are all taken over those replications. A replication
without an interval is in none of these rows. It is counted, and the
returned table carries the counts in its attribute `simulation`, a list
read with `attr(x, "simulation")`:

- `attempted`:

  the number of replications asked for (`G`).

- `returned`:

  the number that produced an interval (a test).

- `failed`:

  the number that did not (`attempted - returned`).

- `denominator`:

  the number of replications the rows named above are taken over, which
  is `returned`.

- `width_target_probability`:

  the proportion of the replications with an interval whose interval met
  the width target, the value of the `pct_ci_less_width` row; `NA` for a
  function or a call without a width target, and when no replication
  produced an interval.

- `width_target_mc_ci`:

  the exact binomial (Clopper-Pearson) limits for that proportion, at
  the confidence level of the table. They describe the error of the
  simulation, which shrinks as `G` grows, and nothing else.

A function that records more about its run adds elements to the list,
and its help page names them. The line "Replications" beneath the
printed table gives `G`, and the three counts when a replication was
without an interval. With many replications without an interval the rows
describe only the studies in which the interval could be formed, so
`failed` is to be reported beside them.

## References

Bonett, D. G. (2002). Sample size requirements for testing and
estimating coefficient alpha. *Journal of Educational and Behavioral
Statistics, 27*(4), 335–340.
[doi:10.3102/10769986027004335](https://doi.org/10.3102/10769986027004335)

Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
missing information from full information maximum likelihood. *Frontiers
in Psychology, 12*, 667802.
[doi:10.3389/fpsyg.2021.667802](https://doi.org/10.3389/fpsyg.2021.667802)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
fraction of missing information from full information maximum
likelihood. *Structural Equation Modeling, 19*(3), 477–494.
[doi:10.1080/10705511.2012.687669](https://doi.org/10.1080/10705511.2012.687669)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution
of the maximum likelihood estimator of Cronbach's alpha. *Psychometrika,
65*(3), 271–280.
[doi:10.1007/BF02296146](https://doi.org/10.1007/BF02296146)

## See also

[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Reduced Monte Carlo sweep (small G) so the example runs quickly;
# raise G for a production sensitivity analysis.
set.seed(113)
ss_aipe_reliability_sensitivity(
  true_reliability     = 0.80,
  supposed_reliability = 0.80,
  n_items = 4, width = 0.15,
  estimator = "alpha",
  G = 20, print_iter = FALSE
)
#>  term                 value
#>  mean_reliability     0.782
#>  median_reliability   0.792
#>  sd_reliability       0.041
#>  mean_ci_width        0.163
#>  median_ci_width      0.155
#>  sd_ci_width          0.031
#>  pct_ci_less_width    0.450
#>  pct_ci_miss_low      0.000
#>  pct_ci_miss_high     0.050
#>  total_type_I_error   0.050
#>  total_N              74   
#>  n_items              4    
#>  true_reliability     0.800
#>  supposed_reliability 0.800
#>  width                0.150
#>  conf_level           0.950
#> 
#> Replications: 20
#> 
#> Confidence level: 95%

# The same plan made for an FMI of .18 and checked against data with each
# item score missing completely at random with probability .2, analyzed
# by full information maximum likelihood. The second line under the table
# gives the FMI those data realized for the coefficient, to read beside
# the .18 the plan assumed. Every replication fits a model, so G = 3 only
# shows the call; a reported check deserves G = 1000.
set.seed(113)
ss_aipe_reliability_sensitivity(
  true_reliability     = 0.80,
  supposed_reliability = 0.80,
  n_items = 4, width = 0.15,
  estimator = "alpha",
  G = 3, fmi = 0.18, missing_prop = 0.2
)
#>  term                 value
#>  mean_reliability     0.819
#>  median_reliability   0.830
#>  sd_reliability       0.030
#>  mean_ci_width        0.132
#>  median_ci_width      0.124
#>  sd_ci_width          0.023
#>  pct_ci_less_width    0.667
#>  pct_ci_miss_low      0.000
#>  pct_ci_miss_high     0.000
#>  total_type_I_error   0.000
#>  total_N              90   
#>  n_items              4    
#>  true_reliability     0.800
#>  supposed_reliability 0.800
#>  width                0.150
#>  conf_level           0.950
#> 
#> Replications: 3
#> 
#> Confidence level: 95%
#> 
#> Planned for a fraction of missing information of 0.18; simulated with each item score missing completely at random with probability 0.2 and analyzed by full information maximum likelihood.
#> The fraction of missing information the simulated data realized for the coefficient averaged 0.156 over the 3 replications in which it was computed.
```
