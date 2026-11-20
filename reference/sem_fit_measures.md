# Fit Measures for a Structural Equation Model

Reports in one table the measures researchers use to judge the global
fit of a structural equation model: the chi square test of exact fit,
the root mean square error of approximation (RMSEA) with its confidence
interval and the tests of close and not close fit, the incremental
indices CFI, TLI, NFI, IFI, and RNI, the SRMR, and gamma hat. The
function takes either a fitted lavaan model, whose own statistics,
sample size convention, and corrected (robust) versions it reads, or the
summary statistics a paper reports, so that a published fit table can be
reproduced and extended with an interval for the RMSEA and the tests it
implies. A fitted model requires lavaan (in Suggests).

## Usage

``` r
sem_fit_measures(
  fit = NULL,
  chisq = NULL,
  df = NULL,
  N = NULL,
  baseline_chisq = NULL,
  baseline_df = NULL,
  p = NULL,
  n_groups = 1,
  conf_level = 0.9,
  close_fit = 0.05,
  not_close_fit = 0.05,
  version = c("auto", "standard", "scaled", "robust"),
  sample_size_convention = c("N-1", "N")
)
```

## Arguments

- fit:

  A fitted lavaan model. It must have converged to an admissible
  solution, as Details describes (lavaan's post-estimation check, made
  on the model's own parameters when the fit has auxiliary variables).
  Supply either `fit` or the summary statistics below, not both.

- chisq, df:

  For summary statistic input, the model's chi square test statistic
  \\T\\ and its degrees of freedom: the likelihood ratio statistic of a
  normal theory maximum likelihood (ML) fit. Summary input gives the
  standard measures only.

- N:

  For summary statistic input, the total sample size, summed over
  groups: a whole number greater than `n_groups`.

- baseline_chisq, baseline_df:

  For summary statistic input, the chi square statistic and degrees of
  freedom of the baseline (independence) model fitted to the same data.
  The incremental indices need them and are `NA` without them. The
  baseline model is nested in the target model, so `baseline_df` must be
  positive and at least `df`.

- p:

  For summary statistic input, the number of observed variables (not the
  number of parameters). Gamma hat and adjusted gamma hat need it and
  are `NA` without it; they are defined here for a single group and a
  covariance structure (see Details), so they are `NA`, with a warning,
  when `df` exceeds the \\p(p + 1)/2\\ covariance moments, and `df` can
  be at most \\p(p + 3)/2\\. Below that bound summary input cannot show
  that a model restricts the means, and supplying `p` states that it
  does not.

- n_groups:

  For summary statistic input, the number of independent groups of a
  multiple group model, whose `chisq`, `df`, and `N` are totals over the
  groups. Default 1. A fitted model supplies its own number of groups.

- conf_level:

  Confidence level of the RMSEA interval and of the gamma hat and
  adjusted gamma hat intervals, which come from the same noncentrality
  limits. Defaults to 0.90, the conventional level of the RMSEA
  interval, whose two limits decide the tests of close and not close fit
  at .05 (see Details); `conf_level = 0.95` gives 95 percent intervals.

- close_fit:

  The RMSEA value \\\epsilon_0\\ under the null hypothesis of close fit,
  \\H_0: \epsilon \le \epsilon_0\\. Default 0.05.

- not_close_fit:

  The RMSEA value \\\epsilon_1\\ under the null hypothesis of not close
  fit, \\H_0: \epsilon \ge \epsilon_1\\. Default 0.05.

- version:

  Which version of the measures to report for a fitted model: `"auto"`
  (the default), `"standard"`, `"scaled"`, or `"robust"`. See Details.
  Summary statistic input accepts `"auto"` and `"standard"`, which
  coincide.

- sample_size_convention:

  For summary statistic input, the multiplier \\a\\ in \\T = a\hat{F}\\
  with which `chisq` was formed: `"N-1"` (the default), each group's
  sample size less one, so \\a = N - G\\ for \\G\\ groups, or `"N"`, so
  \\a = N\\. lavaan's default normal likelihood and Mplus form \\T\\
  with \\N\\, so their statistics take `"N"`. A fitted model supplies
  its own convention.

## Value

A `dmar_tbl` (class `dmar_sem_fit`) with columns `term` and `value` and,
always in this order, the rows `chi_square`, `df`, `p_chi_square` (the
test of exact fit); `rmsea`, `rmsea_ci_lower`, `rmsea_ci_upper`;
`close_fit_bound`, `p_close_fit`, `not_close_fit_bound`,
`p_not_close_fit`; `cfi`, `tli`, `nfi`, `ifi`, `rni`; `srmr`;
`gamma_hat`, `gamma_hat_ci_lower`, `gamma_hat_ci_upper`,
`adjusted_gamma_hat`, `adjusted_gamma_hat_ci_lower`,
`adjusted_gamma_hat_ci_upper`; `N` and `n_groups`. A row the input or
the version does not compute is `NA`. The `sem_fit` attribute records
the source, the version requested and reported and why, the statistic in
the `chi_square` row, the correction behind the robust measures, the
sample size convention and multiplier, why rows of a ULS fit of
continuous data are `NA` (`uls`; `NA` otherwise), the boundary
conventions that applied, and the lines of the printed footer.

## Details

### Exact Fit and the RMSEA

Write \\T\\ for the model's chi square statistic on \\df\\ degrees of
freedom, \\a\\ for the sample size multiplier with which it was formed
(\\T = a\hat{F}\\, where \\\hat{F}\\ is the minimized fit function),
\\G\\ for the number of groups, and \\\alpha\\ for 1 minus `conf_level`.
The *p*-value of the test of exact fit (row `p_chi_square`) is the upper
tail of the central chi square distribution at \\T\\. Under a model that
is wrong by an amount small relative to sampling error, \\T\\ is
approximately noncentral chi square with noncentrality \\\lambda = a
\cdot df \cdot \epsilon^2 / G\\, where \\\epsilon\\ is the population
RMSEA, the misfit per degree of freedom on the scale of the fit function
(Steiger & Lind, 1980; Browne & Cudeck, 1993). The point estimate is
\\\hat{\epsilon} = \sqrt{G \max(T - df, 0) / (a \cdot df)}\\, and the
confidence interval inverts the noncentral chi square distribution at
\\T\\
([`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)
does the inversion) and maps its noncentrality limits to the RMSEA
metric in the same way. The factor \\G\\ is Steiger's (1998, Eqs. 23 and
24) multiple group RMSEA, which divides the discrepancy by the average
degrees of freedom per group; applying the single group formula to a
multiple group statistic gives a value too small by a factor of
\\\sqrt{G}\\ (Steiger, 1998, Eq. 25).

When \\T\\ falls at or below the upper \\\alpha/2\\ quantile of the
central chi square distribution, the lower limit is 0. When it falls at
or below the lower \\\alpha/2\\ quantile, no noncentrality is consistent
with \\T\\ at that level, and the upper limit is reported as 0, as
lavaan reports it and as Yuan, Chan, Marcoulides, and Bentler (2016, Eq.
4) define the corresponding T-size RMSEA, so that the upper limit and
the test of not close fit always agree. This happens with probability
\\\alpha/2\\ for a model that fits exactly, and more often when the
normal theory statistic is deflated, for example by data with light
tails.

The RMSEA is conventionally reported with a 90 percent interval, the one
Steiger and Lind (1980) formed from the .05 and .95 percentile points of
the noncentrality, whose details Browne and Cudeck (1993) give and which
lavaan's summary prints. Its two limits are the two one-sided tests at
level .05 described below: a \\100(1 - 2\alpha)\\ percent interval
settles both tests at level \\\alpha\\, so the lower limit exceeds
\\\epsilon_0\\ exactly when close fit is rejected, and the upper limit
falls below \\\epsilon_1\\ exactly when not close fit is rejected
(MacCallum, Browne, & Sugawara, 1996, p. 136 and Table 1, with one bound
playing both roles). A lower limit of zero likewise says that the test
of exact fit does not reject at .05. The default `conf_level = 0.90`
therefore gives the conventional interval, which reproduces lavaan's
printed limits, and every RMSEA interval in DMAR defaults to it.
`conf_level = 0.95`, the level of the package's other intervals, gives a
wider interval whose limits decide the two tests at .025 rather than at
.05. The printed table states the level beneath it.

A statistic in the millions, as from a large model fitted to a very
large sample, puts the noncentrality where R's noncentral chi square
distribution function, which lavaan uses, loses its accuracy (above
about \\10^5\\) and then fails (above about \\2 \times 10^6\\, where it
once returned an interval collapsed to a point below its estimate). The
limits and the tests of close and not close fit are then computed here,
from the `chi_square` and `df` rows, with the Poisson mixture described
in
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md);
for a robust version whose corrected statistic DMAR cannot rebuild
(FIML, ordered data, or two-stage estimation) they are `NA`, with a
warning, once the statistic exceeds \\10^6\\.

### Tests of Close and Not Close Fit

The test of close fit (Browne & Cudeck, 1993) asks whether the data can
reject \\H_0: \epsilon \le \epsilon_0\\; its *p*-value is
\\P\\\chi^2\_{df}(\lambda_0) \ge T\\\\ with \\\lambda_0 = a \cdot df
\cdot \epsilon_0^2 / G\\ (row `p_close_fit`). Failing to reject close
fit does not establish it. The test of not close fit (MacCallum, Browne,
& Sugawara, 1996) reverses the roles: it asks whether the data can
reject \\H_0: \epsilon \ge \epsilon_1\\, so that a small *p*-value,
\\P\\\chi^2\_{df}(\lambda_1) \le T\\\\ with \\\lambda_1 = a \cdot df
\cdot \epsilon_1^2 / G\\ (row `p_not_close_fit`), is evidence that the
population RMSEA is below \\\epsilon_1\\. It is the equivalence test of
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
which reports it with its interval and verdict. MacCallum et al. (1996)
used 0.05 for both bounds, the defaults here; lavaan's default bound for
not close fit is 0.08, which `not_close_fit = 0.08` reproduces. These
values are conventions for framing a test, not cutoffs that certify a
model.

### Incremental Indices

The incremental indices compare the model with a baseline model, by
default the independence model, whose statistic and degrees of freedom
are \\T_b\\ and \\df_b\\. The comparative fit index is \\CFI = 1 -
\max(T - df, 0) / \max(T - df, T_b - df_b, 0)\\, equal to 1 when both
are zero (Bentler, 1990). Without the truncation at zero it is the
relative noncentrality index, \\RNI = 1 - (T - df)/(T_b - df_b)\\
(Bentler, 1990; McDonald & Marsh, 1990), which can fall outside \\\[0,
1\]\\. The Tucker-Lewis index (Tucker & Lewis, 1973), also called the
nonnormed fit index (Bentler & Bonett, 1980), is \\TLI = (T_b/df_b -
T/df)/(T_b/df_b - 1)\\; it rewards parsimony and can also exceed 1. The
normed fit index is \\NFI = (T_b - T)/T_b\\ (Bentler & Bonett, 1980),
and the incremental fit index is \\IFI = (T_b - T)/(T_b - df)\\ (Bollen,
1989, Eq. 3). No index is clipped to \\\[0, 1\]\\, and an index whose
denominator is zero or negative is `NA`. Brosseau-Liard and Savalei
(2014, Eqs. 3 to 10) give the two population values these sample indices
estimate: the NFI, CFI, RNI, and IFI all converge to the population CFI,
\\1 - F_H/F_B\\ for the population ML discrepancies of the model and the
baseline (the NFI poorly, since it does not subtract the degrees of
freedom), and the TLI to the population TLI, \\1 -
(F_H/F_B)(df_B/df_H)\\. None of them has a confidence interval here: the
usual asymptotic theory assumes small misfit, which the independence
model violates.

A fitted model with no degrees of freedom reproduces the moments
exactly, so its statistic is zero in every version, and so is the model
term \\T - c \cdot df\\ of the robust indices, whatever the scaling
constant \\c\\, which is 0/0 there (lavaan's robust and scaled indices
of such a model come out `NA` or not by rounding). Its CFI is therefore
1, and its RNI, NFI, and IFI are their values at a statistic of zero,
those of lavaan's standard version, in every version that has them; its
TLI, like the RMSEA, divides by the degrees of freedom and is `NA`. A
model with no degrees of freedom need not fit exactly: an inequality
constraint that binds, such as an error variance held at zero by
`bounds = "pos.var"`, fixes its parameter as an equality would, though
lavaan does not count it among the degrees of freedom. Such a model's
statistic is not zero (32.8 for one factor of three variables correlated
.8, .8, and .5, \\N\\ = 200), and its CFI, RNI, NFI, and IFI are formed
from it (.925, .925, .926, and .926), with the standard NFI computed as
\\1 - T/T_b\\ where lavaan reports 1 at zero degrees of freedom; the TLI
and the RMSEA rows are `NA` as before, and the warning and the footer
say that the model does not fit exactly.

### SRMR and Gamma Hat

The standardized root mean square residual (row `srmr`) is the root mean
square of the covariance residuals, each divided by the product of the
two sample standard deviations (lavaan's `srmr_bentler_nomean`). It is
read from a fitted model and is `NA` for summary input. It is the same
quantity in every version, since no test statistic enters it. lavaan's
own `srmr` also averages in the mean residuals when the model has a mean
structure, which makes it smaller by about \\\sqrt{t/(t + p)}\\, for \\t
= p(p + 1)/2\\ covariance moments of \\p\\ variables, when the means are
saturated, as under FIML; this table keeps the covariance definition
throughout. With auxiliary variables (lavaan's `aux` under FIML) the
SRMR is over the model's own variables, the auxiliary variables set
aside as in the admissibility check below. lavaan's measures average
over them as well, and an auxiliary variable, a saturated correlate,
adds residuals of about zero that lower the SRMR of the same model: by
about 8 percent with one auxiliary variable and 15 percent with two for
a three-factor model of nine items. The test statistic, and the measures
built from it, compare the model with the saturated model with the
auxiliary variables in both, so they need no such change. For the
unbiased estimate of the SRMR, its interval, and its tests see
[`sem_residual_fit`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md).

Gamma hat rescales the estimated noncentrality to a goodness of fit
metric, \\\hat{\gamma} = p/\[p + 2(T - df)/a\]\\ for \\p\\ observed
variables, and adjusted gamma hat penalizes it for the degrees of
freedom used, \\1 - \[p(p + 1)/(2 df)\](1 - \hat{\gamma})\\, counting
the \\p(p + 1)/2\\ covariance moments of one group. They estimate the
population values of the goodness of fit index (GFI) and the adjusted
GFI that LISREL reports, which Steiger (1989) showed to be these
functions of the population noncentrality for a model invariant under a
constant scaling factor (a result he restates in Steiger, 1990, p. 178).
Neither is truncated, so both exceed 1 when \\T \< df\\. Each is a
monotone function of the noncentrality, so the noncentrality limits of
the RMSEA interval give their limits, \\p/(p + 2\lambda/a)\\ at each
limit \\\lambda\\ and the adjusted index of that (Steiger & Fouladi,
1997, pp. 251–252), in the rows `gamma_hat_ci_lower` through
`adjusted_gamma_hat_ci_upper`. The limits are at most 1, the upper one
is 1 when the lower RMSEA limit is 0, and a point estimate above 1 lies
above its interval.

This gamma hat is the index lavaan 0.7-3 reports as `gfi_lrt`, with the
same limits; lavaan caps the point estimate at 1 and reports it only
when it is requested by name. The `gfi` lavaan prints by default is the
same function of Browne's residual based (reweighted least squares)
statistic in place of \\T\\, which differs from gamma hat in a finite
sample (0.959 against 0.957 for the three-factor model of the examples),
and lavaan's `gfi_lisrel` and `agfi_lisrel` are the sample GFI and
adjusted GFI of LISREL, not estimates of their population values. lavaan
also reports robust versions, and so does this table: for a robust ML
fit of complete continuous data with a Satorra-Bentler, Yuan-Bentler,
scaled and shifted, or mean and variance adjusted statistic (MLM, MLR,
MLMV, MLMVS), the robust version estimates the noncentrality as
\\s(T_s - df_s)\\, for the corrected statistic \\T_s\\, its scaling
constant \\s\\, and its degrees of freedom \\df_s\\, which is \\T - c
\cdot df\\ in each case (see the robust version below), the estimate
behind the robust RMSEA of Brosseau-Liard, Savalei, and Li (2012), and
takes its limits from the robust RMSEA interval; for the first three
this is lavaan's `gfi_lrt.robust` and its limits.

Steiger's result is for covariance structures. Neither index is computed
for a fitted model that restricts the means, whose statistic and degrees
of freedom include the misfit and the restrictions of the means while
\\p\\ and \\p(p + 1)/2\\ count the covariances: a growth model, or
intercepts constrained equal (by a shared label, with or without
lavaan's `ceq.simple = TRUE`), constrained to a constant, entered in a
linear constraint, tied to a loading, or held at the bound of an
inequality constraint (which lavaan's degrees of freedom do not count).
Whether a fit restricts its means is read from the derivatives of its
model implied moments, with the constraints that bind at the estimates
imposed (the equality constraints, and the inequality constraints lavaan
treats as active when it computes the standard errors): the means are
free when the parameter directions that leave every covariance unchanged
move the means in as many independent directions as there are modeled
means, so that the means add no degree of freedom (the means of
exogenous variables under `fixed.x = TRUE` are fixed at their sample
values and are not counted). A mean structure whose means are free, as
lavaan's default parameterization leaves them or as one indicator's
intercept fixed at zero with the latent mean free leaves them, changes
neither index. Summary input carries no mean structure to read, so it
shows restricted means only when `df` exceeds \\p(p + 1)/2\\, which only
restricted means can produce, and the indices are then `NA`; below that
bound they are computed from the summary statistics on the understanding
that the model is a covariance structure (the statistic of a linear
growth model of four waves, 12.2 on 5 df with \\N\\ = 400, gives 0.991
and 0.982, which the fitted model does not). The two indices are
computed for the standard and robust versions just described of a single
group, single level ML fit of continuous variables without FIML or
sampling weights, and for summary input when `p` is given, and are `NA`
otherwise.

### Standard, Scaled, and Robust Versions

A fitted model can carry more than one test statistic, and lavaan
computes the fit measures from each.

- `"standard"`:

  The measures above computed from the uncorrected statistic. They are
  the right ones for a normal theory ML fit of complete continuous data.

- `"scaled"`:

  The same formulas with a scaled (or mean and variance adjusted)
  statistic \\T_s\\ on its degrees of freedom \\df_s\\ put in place of
  \\T\\ and \\df\\: the RMSEA and its interval invert the noncentral chi
  square at \\T_s\\ on \\df_s\\, the `chi_square` and `df` rows of the
  same table (Brosseau-Liard, Savalei, & Li, 2012, Eqs. 11 and 13, where
  \\df_s = df\\ for a scaled statistic; Savalei, 2018, Table 1), and the
  tests of close and not close fit, which those sources do not treat,
  follow the same substitution so that they keep their duality with the
  interval. This is the substitution EQS and Mplus report, and the one
  lavaan reported when Brosseau-Liard et al. (2012, p. 911) wrote, and
  the one a published table of scaled measures was computed with. It is
  not lavaan's default for a Satorra-Bentler or Yuan-Bentler statistic
  (MLM, MLR) or a mean and variance adjusted one (MLMVS): lavaan 0.7-3
  refers the unscaled \\T\\ to the noncentral chi square on \\c \cdot
  df\\ degrees of freedom, the trace of \\U\Gamma\\ for scaling constant
  \\c\\, which gives the same point estimate but a different interval
  and different close fit *p*-values from those the rows beside them
  imply. A fit made with lavaan's `rmsea.scaled.mplus = TRUE` gives the
  values reported here, as lavaan's default does for a scaled and
  shifted statistic (MLMV, WLSMV). The scaled RMSEA and CFI do not
  consistently estimate the population values the standard formulas
  define (Brosseau-Liard et al., 2012; Brosseau-Liard & Savalei, 2014;
  Savalei, 2018, 2021); they are offered to reproduce a published table.

- `"robust"`:

  lavaan's corrected measures, which estimate the population values of
  the standard measures: for a robust ML fit of complete continuous data
  (MLM, MLR), the nonnormality corrections of Brosseau-Liard, Savalei,
  and Li (2012) for the RMSEA and of Brosseau-Liard and Savalei (2014)
  for the CFI, TLI, and RNI, with the correction of Savalei (2018) for a
  scaled and shifted statistic (MLMV). The point estimates do not depend
  on which corrected statistic accompanies them: each subtracts from the
  ML statistic \\T\\ its asymptotic expectation \\c \cdot df\\, with
  \\c\\ the trace of \\U\Gamma\\ over \\df\\ (Brosseau-Liard et al.,
  2012, Appendix; Savalei, 2018, Eq. 6). lavaan computes no robust
  measures beside the older mean and variance adjusted statistic (MLMVS,
  `test = "mean.var.adjusted"`, the statistic Savalei, 2018, Appendix B,
  writes \\T\_{MVold}\\ and no longer recommends; her mean and variance
  adjusted statistic is the scaled and shifted one), whose scaling
  constant \\s\\ and degrees of freedom \\df_s\\ give \\c = s \cdot df_s
  / df\\, so they are computed here, for the model and for its baseline,
  with the interval and the tests of close and not close fit from
  \\T/c\\ on \\df\\ (Brosseau-Liard et al., 2012, Eq. 10); they equal
  lavaan's robust measures for the same model fitted by MLM, and the
  footer says that DMAR computed them (for a saturated model, whose
  \\c\\ is 0/0 and whose robust CFI and RNI are 1 at a statistic of
  zero, nothing is computed from \\c\\, and the footer names the
  corrections alone, as for MLM). For ordered categorical data (WLSMV),
  the correction of Savalei (2021); and for a FIML fit
  (`missing = "ml"`), the FIML corrected indices of Zhang and Savalei
  (2023). The last are not corrections for nonnormality: they evaluate
  the complete data ML fit function at the FIML estimates, with the
  saturated FIML estimates in place of the sample moments, and so
  approximate the RMSEA and CFI the model would have had with complete
  data, whereas the uncorrected FIML indices change with the amount of
  missing data. lavaan uses their version 3, which does not assume
  normality, when the fit has a scaled statistic (MLR), and their
  version 6, which does, otherwise. The corrected values rest on the
  corrected statistic, which under FIML differs from the chi square
  shown in the table. After two-stage estimation they are the corrected
  indices of Zhang and Savalei (2023), which correct the complete data
  indices of the second stage for the missing data: after robust
  two-stage estimation (`missing = "robust.two.stage"`) their TS V2, the
  version they recommend, and after normal theory two-stage estimation
  (`"two.stage"`) its normal theory analog, which takes the asymptotic
  covariance matrix of the stage one moments to be the inverse of their
  observed information where TS V2 takes the sandwich of the observed
  and first order information, and which they did not study (see
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)).
  Under FIML and two-stage estimation alike, Zhang and Savalei define
  only the RMSEA and the CFI; the robust TLI and RNI, the RMSEA
  interval, and the tests of close and not close fit are lavaan's
  extension of the correction. NFI and IFI have no robust version, in
  lavaan or here, and are `NA`; gamma hat has the one described above.

The default `"auto"` reports the robust version for a fit that carries a
corrected test statistic (MLM, MLR, WLSMV, a sampling weighted fit with
lavaan's default test, and the like) or uses FIML, and the standard
version otherwise. A fit to ordered data by DWLS or ULS has only the
standard statistic unless a corrected test is requested (it is the
default test of `estimator = "DWLS"` and `"ULS"`, which `"WLSMV"` and
`"ULSMV"` extend with a corrected one), and lavaan gives that statistic
no reference distribution, so `"auto"` reports its standard measures
with a warning that the statistic is not chi square distributed and that
a corrected test is needed for the robust measures. The same holds for a
DWLS fit of continuous data whose test was set to `"standard"` (lavaan's
default test for continuous data adds Browne's residual based statistic,
which is chi square under normality and is then the one reported; the
warning says to leave out `test = "standard"` or to use a corrected
test), for a ULS fit of continuous data whose test was set so (see *ULS
Fits of Continuous Data*), and for a sampling weighted fit whose test
was set to `"standard"`, whose statistic is a weighted pseudo likelihood
ratio statistic. Full WLS is the exception among the estimators for
ordered data: its weight matrix is the inverse of the asymptotic
covariance matrix of the thresholds and polychoric correlations, its
statistic is asymptotically chi square (Browne, 1984; Muthén, 1984),
lavaan reports its *p*-value, and `"auto"` reports its standard measures
without a warning (a request for the robust version says that they are
the ones to report). The reason in the footer names the estimator. The
measures of a ULS fit of continuous data are reported only where they do
not depend on the units of the variables (see *ULS Fits of Continuous
Data*). Browne's residual based statistic, which lavaan computes beside
the standard one, is an alternative test and not a correction, and does
not count. When lavaan provides a corrected statistic but no robust
measures for it (and DMAR cannot compute them), `"auto"` reports, with a
warning, the standard version for a continuous ML fit, whose standard
measures estimate the population values of the ML measures as the robust
ones do and the scaled ones do not (Savalei, 2018), though their
interval and tests assume normal data, and the scaled version for any
other fit. The `chi_square`, `df`, and `p_chi_square` rows show the
scaled statistic under the scaled and robust versions when the fit has
one, and the standard statistic otherwise. The `sem_fit` attribute and
the printed footer name the version reported, why, and the statistic
used.

### ULS Fits of Continuous Data

The corrections of Brosseau-Liard and colleagues target the population
discrepancy of the ML fit function. Their use with another fit function
needs separate justification (Savalei, 2018). For continuous ULS,
`version = "robust"` is therefore refused for the entire robust family,
including CFI, TLI, and RNI. Invariance to multiplying every variable by
the same constant is not enough to establish the population
interpretation of an incremental index. This is a conservative
supported-method policy, not a claim that Savalei proved every
continuous-ULS index invalid. Ordered-data corrections with their own
derivation are handled separately.

With ULSM or ULSMV, `version = "auto"` selects eligible scaled indices.
These put a corrected test statistic into the index formulas; they are
descriptive approximations, do not estimate ML-population fit, and have
no established license to use conventional ML cutoffs. When the fit
restricts means, incompatible units in the mean and covariance portions
of the ULS discrepancy can make these indices change even under a common
rescaling. Those rows are withheld, with the existing mean-structure
qualification.

The standard statistic \\N F\_{ULS}\\ is in squared-covariance units.
Its RMSEA and incremental indices that subtract degrees of freedom
change with those units and are `NA`; the standard NFI, a ratio of two
such statistics, remains available when eligible. This does not apply to
a separately justified residual-based statistic. A fit by
`estimator = "ULS"` alone is tested in lavaan 0.7-3 by Browne's normal
theory residual based statistic, which weights the residuals by the
inverse of their normal theory covariance matrix and has no units as
long as the estimates at which it is evaluated do not change with the
units either. They do not when the model leaves the means free, or
restricts them through parameters the covariances do not share, as the
equal intercepts above do, and the standard measures formed from that
statistic are then reported as they are. When the model restricts the
means through parameters the covariances share, as the loadings enter
the means through the latent means of a multiple group model whose
loadings and intercepts are equal (the scalar rung of
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)),
the ULS fit function trades the mean residuals against the covariance
residuals, which are in other units, so the estimates change with the
units, and so do Browne's statistic and every measure formed from it: on
the scalar rung of the first example of
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
the CFI is .8874 for the raw scores and .8861 for the scores divided by
10, and the RMSEA .0680 and .0684. The RMSEA, its interval, the tests of
close and not close fit, and the CFI, TLI, RNI, NFI, and IFI are then
`NA`, with the same warning, element, and line; the `chi_square` row,
whose chi square reference distribution holds at any root-\\n\\
consistent estimates (Browne, 1984, Proposition 4), and the SRMR, which
describe the fit at hand, are reported. Whether the two parts share a
parameter is read from the same derivatives as whether the means are
restricted: they share one when the rank of the derivatives of every
model implied moment falls short of the rank of the covariance rows plus
that of the mean rows. With ordered data the robust measures of a fit
with `"ULSMV"` are the categorical correction of Savalei (2021), formed
from polychoric correlations, which have no units, and none of this
applies. Maximum likelihood is scale invariant; `estimator = "MLR"` is
the estimator for nonnormal continuous data.

### Two-Stage Fits

A fit by two-stage maximum likelihood (Yuan & Bentler, 2000; Savalei &
Bentler, 2009) fits the model to the EM estimates of the means and
covariances as if they were the moments of complete data. lavaan's
standard statistic for it is that naive one, which is not chi square
distributed, so its standard measures (TS V0 in the terms of Zhang &
Savalei, 2023) are refused, and so are its scaled measures, which put an
approximate statistic into the standard formulas; the default `"auto"`
reports the robust version. The `chi_square` row is then the residual
based statistic of Savalei and Bentler (2009, Eq. 8; Eq. 14 with
auxiliaries) after normal theory two-stage estimation
(`missing = "two.stage"`), computed by DMAR from lavaan's moments and
the observed information of the saturated model, and lavaan's scaled
statistic after robust two-stage estimation (Savalei & Falk, 2014): the
statistics
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
and
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
report for the same fits. The indices are not computed from either
statistic but from the naive one and its expectation, as Zhang and
Savalei (2023) define them: TS V2 after robust two-stage estimation and
its normal theory analog after normal theory two-stage estimation.

A case observed on an auxiliary variable (lavaan's `aux`) but on none of
the model's variables is read by the first stage, whose EM estimates of
the means and covariances it moves, and is left out of the sample size
lavaan 0.7-3 reports (`lavaan::lavInspect(fit, "nobs")`), which counts
it as empty. The naive statistic, and with it lavaan's scaled statistic,
the corrected indices, the RMSEA multiplier, and the `N` row, then rest
on fewer cases than the moments, while the residual based statistic of a
normal theory two-stage fit uses every case, so the rows of the table do
not describe one set of cases (and lavaan's two-stage standard errors
are too large by the square root of the ratio of the two counts).
`sem_fit_measures()` reports the rows of a fit as lavaan computed them,
and warns when a fit has such cases and counts them.
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)`(missing = "two_stage")`
and
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
keep such cases and compute their fits again with the sample size set to
every case the first stage read, so that the fits they return, and the
tables of this function for them, count every case; FIML
(`missing = "ml"`) with the same auxiliary variables keeps them as well.
Setting the cases aside instead would condition the analysis on an
observed variable of the model, which biases the estimates when the
items went missing together at random given the auxiliary variables.

### The Sample Size Convention

The multiplier \\a\\ must be the one with which \\T\\ was formed. A
fitted model supplies it: lavaan forms \\T = N\hat{F}\\ for ML under its
default normal likelihood, and \\T = (N - G)\hat{F}\\ under
`likelihood = "wishart"` and for the least squares estimators, and the
table follows the fit, so its values equal lavaan's (for the scaled
version, those lavaan gives under `rmsea.scaled.mplus = TRUE`; see
above). For summary input the default `"N-1"` is the convention of
Browne and Cudeck (1993), MacCallum et al. (1996), Kelley and Lai
(2011), and DMAR's
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md),
and
[`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md).
A statistic reported by lavaan under its default likelihood or by Mplus
was formed with \\N\\ and takes `sample_size_convention = "N"`. The
wrong multiplier changes the RMSEA, its limits, the close fit tests, and
gamma hat by an amount of order \\1/N\\.

### What Is Refused

A fitted model is refused when it did not converge (an error of class
`dmar_not_converged`) and when its solution is inadmissible (a refusal
of class `dmar_inadmissible`): fit measures of such a solution describe
estimates outside the parameter space and do not speak to the specified
model. A solution is inadmissible when lavaan's post-estimation check
finds a negative variance estimate (a Heywood case among them) or a
covariance matrix that is not positive definite: that of the latent
variables, that of the residuals of the observed variables (theta), or
that of the residuals of the structural part of the model (psi). The
message names the variable or the matrix, and the condition carries them
(see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).
With auxiliary variables (lavaan's `aux` under FIML) the check is made
on the model's own parameters, with the auxiliary variables set aside:
lavaan models an auxiliary variable as a saturated correlate, whose
covariances with the items run through their residuals, and its own
check of theta, which includes those covariances, fails once an
auxiliary variable is related to a factor strongly enough, however
proper the model's own solution. For one factor whose items have
uncorrelated residuals, and an auxiliary variable related to the items
only through the factor, that happens when the squared correlation
between the auxiliary variable and the factor exceeds \\(1 - H) / H\\,
with *H* the maximal reliability of the factor's items (coefficient *H*;
see
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)),
which is a correlation above .33 when *H* is .9. Auxiliary variables are
recognized when they are given through lavaan's `aux` argument;
saturated correlates written into the model by hand, as
`semTools::cfa.auxiliary()` writes them, stay in theta and are judged
with it, so refit with `aux` to have them set aside. The same check
serves
[`sem_residual_fit`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`ecvi`](https://yelleknek.github.io/DMAR/reference/ecvi.md), and
[`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
judges its fit and its replications by it. Resolve the improper solution
first, by changing the model or by judging that a negative variance is
within sampling error of zero and fixing it, or pass the reported
statistics as summary input if you accept that responsibility.
Multilevel fits are refused, since one set of indices mixes the two
levels; lavaan's level specific fit measures
(`lavaan::fitMeasures(fit, level = ...)`) are the tool for them. A fit
whose cases were removed by listwise deletion is evaluated on the
complete cases with a warning, since that assumes the data are missing
completely at random; refit it with `missing = "ml"` to use every case.

For a fit by normal theory two-stage estimation
(`missing = "two.stage"`) the residual based statistic, and with it the
table, is refused when the fit has sampling weights or clustered cases,
conditions on an observed exogenous variable (`fixed.x = TRUE`, the
default of [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html)
for a model that regresses on one, or `conditional.x = TRUE`), has
inequality constraints or a parameter on a bound, or has auxiliary
variables that the installed lavaan does not keep where DMAR reads them.
Each message says what to do instead; robust two-stage estimation, whose
scaled statistic lavaan computes, is the alternative throughout. A fit
that conditions on an exogenous variable is refused with an error of
class `dmar_two_stage_fixed_x`, which carries the exogenous variables,
and the advice to refit with `fixed.x = FALSE`, and with
`conditional.x = FALSE` too, since lavaan does not take the first beside
`conditional.x = TRUE`. The refit leaves the estimates as they are when
the exogenous variables are complete, and otherwise also keeps the cases
lavaan deletes for a missing exogenous variable when it conditions on
them (robust two-stage estimation deletes them too), so that its
estimates differ; to keep the fit's estimates, refit on the rows the fit
used, those on which every exogenous variable is observed, which
`lavaan::lavInspect(fit, "case.idx")` lists.

### Reading the Measures

No row of the table certifies a model, and the table reports no cutoff.
The conventions in use are descriptions of a population RMSEA: MacCallum
et al. (1996, p. 134) credit the view that values below .05 indicate
close fit to Steiger (1989) and others, describe .05 to .08 as fair and
above .10 as poor after Browne and Cudeck (1993), and call .08 to .10
mediocre themselves. Chen, Curran, Bollen, Kirby, and Paxton (2008)
found no empirical support for .05, or any other value, as a universal
cutoff: a fixed bound rejects too many correct models in samples of 100
or fewer, too few misspecified ones in samples of 800 or more, and the
test of close fit with a fixed bound performs inconsistently across
sample sizes, degrees of freedom, and specifications. They recommend
reading the RMSEA beside the chi square statistic and the other indices
rather than as a single accept or reject decision; the interval and the
two tests above are how this table states the uncertainty, and the
researcher chooses the bound.

The cutoffs in common use for the other indices, near .95 for the CFI
and TLI, .06 for the RMSEA, and .08 for the SRMR, applied in pairs, come
from the simulations of Hu and Bentler (1999), who themselves wrote that
no single cutoff works equally well across conditions and that the
trade-off between the two kinds of error at samples of 250 or fewer
cannot be judged from the sample. Marsh, Hau, and Wen (2004) showed that
the rationale behind those simulations treated a minor misspecification
as a model to be rejected, so that the rules reward rejecting models
that fit acceptably and cannot hold across sample sizes; they warn
against reading the rules of thumb as golden rules. What a value of an
index means depends on the model: on the size of the loadings and unique
variances, which lower the chi square, RMSEA, and SRMR of a misspecified
model when the loadings are low (Heene, Hilbert, Draxler, Ziegler, &
Buehner, 2011); on the number of observed variables, which inflates the
statistic and every index formed from it (Moshagen, 2012); and on the
reliability of the factors, the number of items and factors, and the
degrees of freedom, which is why McNeish and Wolf (2023) derive cutoffs
for the model at hand by simulation, at several sizes of
misspecification, so that an index reads as an effect size rather than
as a verdict. This table reports the indices and their uncertainty; the
bound a model is held to is a judgment about the model, stated in
advance.

With few degrees of freedom the RMSEA is unstable (Kenny, Kaniskan, &
McCoach, 2015): for correctly specified models with 1 to 3 degrees of
freedom and \\N\\ = 400, the sample RMSEA exceeded .05 in 13 to 15
percent of samples, against none at 50 degrees of freedom, and its
interval at \\N\\ = 50 and 1 degree of freedom was more than sixteen
times as wide as at \\N\\ = 1,000 and 100 degrees of freedom. They
advise against computing the RMSEA when the degrees of freedom are small
and suggest diagnosing poor fit by fitting the models that free the
omitted parameters instead
([`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
compares them). The interval in this table makes the instability visible
rather than hiding it.

The incremental indices measure fit against the baseline model, and
their metric is the baseline's misfit: under maximum likelihood the
independence baseline's statistic is \\-(N - 1)\log\|R\|\\ for the
sample correlation matrix \\R\\ (van Laar & Braeken, 2022). Two
correctly specified models with the same chi square and RMSEA
distributions have different CFI distributions when their \\\|R\|\\
differ, and misspecified models can rank differently on the CFI and on
the chi square, so a CFI is comparable across data sets only when the
baselines are; van Laar and Braeken recommend reporting \\\|R\|\\ and
the baseline statistic and reading an incremental index as the share of
the baseline's misfit the model removes. A fitted model's baseline
statistic and degrees of freedom are
`lavaan::fitMeasures(fit, c("baseline.chisq", "baseline.df"))`; summary
input takes them as given.

The raw-data H1 screen is sufficient and nonexhaustive. In models with
fixed or conditional observed predictors, it screens response variables
only, since predictor moments are not free. It does not assess
conditional design rank or exact response-on-predictor relations. No
detected obstacle is not a proof of likelihood existence or uniqueness.

## References

Bentler, P. M. (1990). Comparative fit indexes in structural models.
*Psychological Bulletin, 107*(2), 238–246.
[doi:10.1037/0033-2909.107.2.238](https://doi.org/10.1037/0033-2909.107.2.238)

Bentler, P. M., & Bonett, D. G. (1980). Significance tests and goodness
of fit in the analysis of covariance structures. *Psychological
Bulletin, 88*, 588–606.

Bollen, K. A. (1989). A new incremental fit index for general structural
equation models. *Sociological Methods & Research, 17*(3), 303–316.
[doi:10.1177/0049124189017003004](https://doi.org/10.1177/0049124189017003004)

Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
indices for nonnormality. *Multivariate Behavioral Research, 49*(5),
460–470.
[doi:10.1080/00273171.2014.933697](https://doi.org/10.1080/00273171.2014.933697)

Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An investigation of
the sample performance of two nonnormality corrections for RMSEA.
*Multivariate Behavioral Research, 47*(6), 904–930.
[doi:10.1080/00273171.2012.715252](https://doi.org/10.1080/00273171.2012.715252)

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*(1), 62–83.

Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model
fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation
models* (pp. 136–162). Sage.

Chen, F., Curran, P. J., Bollen, K. A., Kirby, J., & Paxton, P. (2008).
An empirical evaluation of the use of fixed cutoff points in RMSEA test
statistic in structural equation models. *Sociological Methods &
Research, 36*(4), 462–494.
[doi:10.1177/0049124108314720](https://doi.org/10.1177/0049124108314720)

Heene, M., Hilbert, S., Draxler, C., Ziegler, M., & Buehner, M. (2011).
Masking misfit in confirmatory factor analysis by increasing unique
variances: A cautionary note on the usefulness of cutoff values of fit
indices. *Psychological Methods, 16*(3), 319–336.
[doi:10.1037/a0024917](https://doi.org/10.1037/a0024917)

Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in
covariance structure analysis: Conventional criteria versus new
alternatives. *Structural Equation Modeling, 6*(1), 1–55.
[doi:10.1080/10705519909540118](https://doi.org/10.1080/10705519909540118)

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*(1), 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015). The performance of
RMSEA in models with small degrees of freedom. *Sociological Methods &
Research, 44*(3), 486–507.
[doi:10.1177/0049124114543236](https://doi.org/10.1177/0049124114543236)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Marsh, H. W., Hau, K.-T., & Wen, Z. (2004). In search of golden rules:
Comment on hypothesis-testing approaches to setting cutoff values for
fit indexes and dangers in overgeneralizing Hu and Bentler's (1999)
findings. *Structural Equation Modeling, 11*(3), 320–341.
[doi:10.1207/s15328007sem1103_2](https://doi.org/10.1207/s15328007sem1103_2)

McDonald, R. P., & Marsh, H. W. (1990). Choosing a multivariate model:
Noncentrality and goodness of fit. *Psychological Bulletin, 107*,
247–255.

McNeish, D., & Wolf, M. G. (2023). Dynamic fit index cutoffs for
confirmatory factor analysis models. *Psychological Methods, 28*(1),
61–88. [doi:10.1037/met0000425](https://doi.org/10.1037/met0000425)

Moshagen, M. (2012). The model size effect in SEM: Inflated
goodness-of-fit statistics are due to the size of the covariance matrix.
*Structural Equation Modeling, 19*(1), 86–98.
[doi:10.1080/10705511.2012.634724](https://doi.org/10.1080/10705511.2012.634724)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Savalei, V. (2018). On the computation of the RMSEA and CFI from the
mean-and-variance corrected test statistic with nonnormal data in SEM.
*Multivariate Behavioral Research, 53*(3), 419–429.
[doi:10.1080/00273171.2018.1455142](https://doi.org/10.1080/00273171.2018.1455142)

Savalei, V. (2021). Improving fit indices in structural equation
modeling with categorical data. *Multivariate Behavioral Research,
56*(3), 390–407.
[doi:10.1080/00273171.2020.1717922](https://doi.org/10.1080/00273171.2020.1717922)

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*, 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
robust full information maximum likelihood with incomplete nonnormal
data. *Structural Equation Modeling, 21*(2), 280–302.
[doi:10.1080/10705511.2014.882692](https://doi.org/10.1080/10705511.2014.882692)

Steiger, J. H. (1989). *EzPATH: A supplementary module for SYSTAT and
SYGRAPH*. SYSTAT.

Steiger, J. H. (1990). Structural model evaluation and modification: An
interval estimation approach. *Multivariate Behavioral Research, 25*(2),
173–180.

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
fit index. *Structural Equation Modeling, 5*(4), 411–419.
[doi:10.1080/10705519809540115](https://doi.org/10.1080/10705519809540115)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

Steiger, J. H., & Lind, J. C. (1980). *Statistically-based tests for the
number of common factors*. Paper presented at the annual Spring meeting
of the Psychometric Society, Iowa City, IA.

Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
likelihood factor analysis. *Psychometrika, 38*(1), 1–10.

van Laar, S., & Braeken, J. (2022). Caught off base: A note on the
interpretation of incremental fit indices. *Structural Equation
Modeling, 29*(6), 935–943.
[doi:10.1080/10705511.2022.2050730](https://doi.org/10.1080/10705511.2022.2050730)

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

Yuan, K.-H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016).
Assessing structural equation models by equivalence testing with
adjusted fit indexes. *Structural Equation Modeling, 23*(3), 319–330.
[doi:10.1080/10705511.2015.1065414](https://doi.org/10.1080/10705511.2015.1065414)

Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
following FIML and TS estimation with missing data. *Psychological
Methods, 28*(2), 263–283.
[doi:10.1037/met0000445](https://doi.org/10.1037/met0000445)

## See also

[`sem_residual_fit`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md)
for the unbiased SRMR with its interval and tests,
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md)
for the test of not close fit with its verdict,
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md) for
an interval from a reported RMSEA,
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
and
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
for comparing models, and the vignette
[`vignette("sem_fit_measures", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/sem_fit_measures.md).

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
[`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
[`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The statistics lavaan reports for the three-factor model of nine of
# the holzinger_swineford tests (N = 301), with its baseline.
# lavaan forms the statistic with N under its default likelihood.
sem_fit_measures(chisq = 85.3055, df = 24, N = 301,
                 baseline_chisq = 918.8516, baseline_df = 36, p = 9,
                 sample_size_convention = "N")
#>  term                        value   
#>  chi_square                  85.305  
#>  df                          24      
#>  p_chi_square                < 0.0001
#>  rmsea                       0.092   
#>  rmsea_ci_lower              0.071   
#>  rmsea_ci_upper              0.114   
#>  close_fit_bound             0.050   
#>  p_close_fit                 0.0007  
#>  not_close_fit_bound         0.050   
#>  p_not_close_fit             0.9993  
#>  cfi                         0.931   
#>  tli                         0.896   
#>  nfi                         0.907   
#>  ifi                         0.931   
#>  rni                         0.931   
#>  srmr                        <NA>    
#>  gamma_hat                   0.957   
#>  gamma_hat_ci_lower          0.936   
#>  gamma_hat_ci_upper          0.974   
#>  adjusted_gamma_hat          0.919   
#>  adjusted_gamma_hat_ci_lower 0.879   
#>  adjusted_gamma_hat_ci_upper 0.950   
#>  N                           301     
#>  n_groups                    1       
#> 
#> Confidence level: 90%
#> 
#> Version: standard (summary statistics give the standard measures only)
#> Sample size: N = 301 in 1 group; multiplier a = 301 (sample_size_convention = "N")

# The same model fitted here. The interval is the conventional 90
# percent one, the default, whose limits are the two one-sided .05 tests
# of close and not close fit.
fit <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
   speed   =~ t10_addition + t12_counting_groups_of_dots +
              t13_straight_and_curved_capitals",
  data = holzinger_swineford, std.lv = TRUE)
sem_fit_measures(fit)
#>  term                        value   
#>  chi_square                  85.306  
#>  df                          24      
#>  p_chi_square                < 0.0001
#>  rmsea                       0.092   
#>  rmsea_ci_lower              0.071   
#>  rmsea_ci_upper              0.114   
#>  close_fit_bound             0.050   
#>  p_close_fit                 0.0007  
#>  not_close_fit_bound         0.050   
#>  p_not_close_fit             0.9993  
#>  cfi                         0.931   
#>  tli                         0.896   
#>  nfi                         0.907   
#>  ifi                         0.931   
#>  rni                         0.931   
#>  srmr                        0.065   
#>  gamma_hat                   0.957   
#>  gamma_hat_ci_lower          0.936   
#>  gamma_hat_ci_upper          0.974   
#>  adjusted_gamma_hat          0.919   
#>  adjusted_gamma_hat_ci_lower 0.879   
#>  adjusted_gamma_hat_ci_upper 0.950   
#>  N                           301     
#>  n_groups                    1       
#> 
#> Confidence level: 90%
#> 
#> Version: standard (auto: a normal theory ML fit of complete data with no corrected statistic)
#> Chi square: standard
#> Sample size: N = 301 in 1 group; multiplier a = 301 (lavaan's normal likelihood forms T with N)

# The 95 percent interval, wider, decides the two tests at .025.
sem_fit_measures(fit, conf_level = 0.95)
#>  term                        value   
#>  chi_square                  85.306  
#>  df                          24      
#>  p_chi_square                < 0.0001
#>  rmsea                       0.092   
#>  rmsea_ci_lower              0.067   
#>  rmsea_ci_upper              0.118   
#>  close_fit_bound             0.050   
#>  p_close_fit                 0.0007  
#>  not_close_fit_bound         0.050   
#>  p_not_close_fit             0.9993  
#>  cfi                         0.931   
#>  tli                         0.896   
#>  nfi                         0.907   
#>  ifi                         0.931   
#>  rni                         0.931   
#>  srmr                        0.065   
#>  gamma_hat                   0.957   
#>  gamma_hat_ci_lower          0.931   
#>  gamma_hat_ci_upper          0.976   
#>  adjusted_gamma_hat          0.919   
#>  adjusted_gamma_hat_ci_lower 0.871   
#>  adjusted_gamma_hat_ci_upper 0.956   
#>  N                           301     
#>  n_groups                    1       
#> 
#> Confidence level: 95%
#> 
#> Version: standard (auto: a normal theory ML fit of complete data with no corrected statistic)
#> Chi square: standard
#> Sample size: N = 301 in 1 group; multiplier a = 301 (lavaan's normal likelihood forms T with N)

# A robust ML fit reports the robust measures by default.
fit_mlr <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
   speed   =~ t10_addition + t12_counting_groups_of_dots +
              t13_straight_and_curved_capitals",
  data = holzinger_swineford, std.lv = TRUE, estimator = "MLR")
sem_fit_measures(fit_mlr)
#>  term                        value   
#>  chi_square                  87.132  
#>  df                          24      
#>  p_chi_square                < 0.0001
#>  rmsea                       0.092   
#>  rmsea_ci_lower              0.072   
#>  rmsea_ci_upper              0.114   
#>  close_fit_bound             0.050   
#>  p_close_fit                 0.0005  
#>  not_close_fit_bound         0.050   
#>  p_not_close_fit             0.9995  
#>  cfi                         0.930   
#>  tli                         0.895   
#>  nfi                         <NA>    
#>  ifi                         <NA>    
#>  rni                         0.930   
#>  srmr                        0.065   
#>  gamma_hat                   0.956   
#>  gamma_hat_ci_lower          0.935   
#>  gamma_hat_ci_upper          0.973   
#>  adjusted_gamma_hat          0.918   
#>  adjusted_gamma_hat_ci_lower 0.879   
#>  adjusted_gamma_hat_ci_upper 0.949   
#>  N                           301     
#>  n_groups                    1       
#> 
#> Confidence level: 90%
#> 
#> Version: robust (auto: the fit carries the corrected statistic yuan.bentler.mplus)
#> Chi square: yuan.bentler.mplus (scaled); robust measures: nonnormality corrections of Brosseau-Liard, Savalei, and Li (2012) and Brosseau-Liard and Savalei (2014)
#> Sample size: N = 301 in 1 group; multiplier a = 301 (lavaan's normal likelihood forms T with N)

# The model fitted in each of the two schools: the multiple group RMSEA
# needs the number of groups (Steiger, 1998).
sem_fit_measures(chisq = 115.8513, df = 48, N = 301, n_groups = 2,
                 sample_size_convention = "N")
#>  term                        value   
#>  chi_square                  115.851 
#>  df                          48      
#>  p_chi_square                < 0.0001
#>  rmsea                       0.097   
#>  rmsea_ci_lower              0.075   
#>  rmsea_ci_upper              0.120   
#>  close_fit_bound             0.050   
#>  p_close_fit                 0.0006  
#>  not_close_fit_bound         0.050   
#>  p_not_close_fit             0.9994  
#>  cfi                         <NA>    
#>  tli                         <NA>    
#>  nfi                         <NA>    
#>  ifi                         <NA>    
#>  rni                         <NA>    
#>  srmr                        <NA>    
#>  gamma_hat                   <NA>    
#>  gamma_hat_ci_lower          <NA>    
#>  gamma_hat_ci_upper          <NA>    
#>  adjusted_gamma_hat          <NA>    
#>  adjusted_gamma_hat_ci_lower <NA>    
#>  adjusted_gamma_hat_ci_upper <NA>    
#>  N                           301     
#>  n_groups                    2       
#> 
#> Confidence level: 90%
#> 
#> Version: standard (summary statistics give the standard measures only)
#> Sample size: N = 301 in 2 groups; multiplier a = 301 (sample_size_convention = "N")
```
