# Mediation Analysis via Model-Based Constrained Optimization

Tests hypotheses about mediation effects with the model-based
constrained optimization (MBCO) procedure of Tofighi and Kelley (2020):
a likelihood ratio test that compares the full mediation model to a null
model fit subject to the nonlinear constraint that the effect of
interest (an indirect effect, a total effect, or any smooth function of
path coefficients) equals zero. The model is specified in lavaan syntax
but is fit with OpenMx, whose optimizers support the nonlinear equality
constraints the procedure requires (see Details). The function
enumerates the mediation effects implied by the model (the total effect,
the direct effect, the total indirect effect, and every specific
indirect pathway from `x` to `y`), reports each with a confidence
interval, and tests each with the MBCO likelihood ratio statistic and
its *p*-value. Models may include observed or latent variables, one or
many mediators, in parallel or in sequence, and raw data or summary
statistics may be supplied. With a grouping variable the model becomes a
multiple-group SEM: every effect is estimated within each group and the
between-group difference of each effect is tested, which is moderated
mediation with a categorical moderator. With a numeric `moderator`,
interactions stated in the model are probed: conditional effects at
chosen moderator values, the index of moderated mediation, and a joint
constancy test, all with MBCO likelihood ratio inference (see Details).

## Usage

``` r
mediation_mbco(
  model,
  data = NULL,
  S = NULL,
  M = NULL,
  N = NULL,
  group = NULL,
  x = NULL,
  y = NULL,
  moderator = NULL,
  probe_values = NULL,
  hypotheses = NULL,
  ci_method = c("profile_likelihood", "monte_carlo", "wald"),
  conf_level = 0.95,
  B = 10000,
  optimizer = c("SLSQP", "CSOLNP", "NPSOL"),
  seed = NULL
)
```

## Arguments

- model:

  A character string with the model in lavaan syntax (regressions `~`,
  latent variable definitions `=~`, residual covariances `~~`,
  optionally labeled coefficients and defined quantities via `:=`), or a
  fitted lavaan object whose parameter table and data are then reused.
  Parameter labels must be valid R names (e.g., `b1`, not `1b`).

- data:

  A `data.frame` of raw data containing the observed variables. Rows
  with missing values are retained and handled by full information
  maximum likelihood. Supply either `data` or `S` (with `N`), not both.

- S:

  A covariance matrix of the observed variables (with dimnames naming
  the variables), used with `N` and optionally `M` when raw data are not
  available. Complete-data maximum likelihood depends on the data only
  through the sample means and covariance matrix, so an analysis from
  summary statistics reproduces the raw-data analysis exactly (see
  Details). For a multiple-group analysis, a named list of such
  matrices, one per group.

- M:

  Optional named numeric vector of variable means accompanying `S` (for
  multiple groups, a named list of such vectors). When omitted, means of
  zero are used; intercepts are then uninteresting but every effect,
  test, and confidence interval is unaffected.

- N:

  Total sample size. Required with `S`; ignored for raw data. For a
  multiple-group analysis, a vector with one sample size per group in
  `S` (matched by name when named).

- group:

  Optional name (single character string) of the grouping variable in
  `data` for a multiple-group analysis. With summary input, supply `S`,
  `M`, and `N` as lists instead and leave `group` alone (it then merely
  names the constructed grouping column). See the Details section on
  multiple groups, including how lavaan's labeling rules decide which
  parameters differ by group.

- x, y:

  Names (single character strings) of the focal predictor and the
  outcome between which effects are traced. When both are omitted the
  function looks for a unique variable with no incoming regression path
  (the source) and a unique variable with no outgoing regression path
  (the outcome); if either is ambiguous an error asks for `x` and `y`
  explicitly. When the model defines quantities via `:=` or `hypotheses`
  are supplied, `x` and `y` may be omitted entirely and only those
  quantities are estimated and tested.

- moderator:

  Optional name (single character string) of a numeric moderator
  variable, turning on probing: every pathway effect that involves an
  interaction with the moderator is estimated at each probe value, the
  change in the effect per unit moderator is reported (for a pathway
  moderated in one place, the index of moderated mediation), and a
  constancy test asks whether the pathway effect depends on the
  moderator at all. State the interaction in the model syntax with a `:`
  term (e.g., `m ~ x + w + x:w`) or include a precomputed product column
  among the predictors; both are recognized. A categorical moderator
  goes through `group` instead.

- probe_values:

  Numeric vector of moderator values at which moderated effects are
  estimated; names, when given, label the rows. Defaults to the
  moderator's mean and one standard deviation either side (labeled
  `low`, `mean`, `high`), or to the two observed values of a two-valued
  moderator.

- hypotheses:

  Optional named character vector of additional quantities to estimate
  and test against zero, written as R expressions in the parameter
  labels and effect names (e.g.,
  `c(difference = "indirect_via_imagery - indirect_via_repetition")`
  tests the equality of two specific indirect effects). A named *list*
  may mix such scalar expressions with character vectors of several
  expressions; a multi-expression element is tested jointly (all
  constrained to zero at once, with as many degrees of freedom as
  expressions), and its row reports the test alone, with no scalar
  estimate.

- ci_method:

  Confidence interval method for every reported effect:
  `"profile_likelihood"` (default; Neale & Miller, 1997, computed by
  OpenMx), `"monte_carlo"` (simulate the coefficients from their joint
  normal approximation and form percentile intervals of the effect;
  MacKinnon, Lockwood, & Williams, 2004), or `"wald"` (delta method
  symmetric interval; reported for comparison, not recommended for
  indirect effects, whose sampling distributions are skewed).

- conf_level:

  Confidence level. Defaults to 0.95.

- B:

  Number of Monte Carlo replications when `ci_method = "monte_carlo"`.
  Defaults to 10000.

- optimizer:

  Optimizer used by OpenMx for all fits: `"SLSQP"` (default),
  `"CSOLNP"`, or `"NPSOL"`. The CRAN build of OpenMx ships SLSQP and
  CSOLNP; NPSOL, the optimizer used in Tofighi and Kelley (2020), is
  available only in the build distributed from the OpenMx website.

- seed:

  Optional integer seed for the Monte Carlo interval, used locally (the
  caller's random number generator state is restored on exit). Default
  `NULL` leaves the random number generator state alone.

## Value

A `data.frame` (classes `dmar_mediation_mbco`, `dmar_tbl`) with one row
per effect and columns `pathway` (the traced pathway or defining
expression), `term`, `estimate`, `se` (delta method), `ci_lower`,
`ci_upper`, `lrt` (the MBCO likelihood ratio statistic), `df`,
`p_value`, `delta_aic`, and `delta_bic`. A joint test row (a
several-constraint moderation or hypothesis test) reports the test
columns with `df` equal to the number of constraints; its `estimate`,
`se`, and interval are `NA` because no single number summarizes several
constraints. Attributes: `"conf_level"`, `"ci_method"`, `"optimizer"`,
`"x"`, `"y"`, `"N"`, `"deviance"`, `"aic"`, `"bic"`, `"n_par"` (full
model fit information), `"groups"` (the group labels, reference first;
multiple-group fits only), `"R2"` (named vector, endogenous variables
under the full model, per group when grouped), `"delta_R2"` (matrix of
full minus null \\R^2\\, variables by tested effects), `"moderation"`
(for a moderated analysis: the moderator name, the probe values, and
each moderated effect's moderator polynomial, which is what
[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
draws), and `"mx_model"` (the fitted OpenMx full model, an escape hatch
for further OpenMx work). Use `tidy()` and `glance()` for broom-style
views.

## Details

**The MBCO procedure.** A hypothesis about mediation is recast as a
model comparison (Maxwell, Delaney, & Kelley, 2027). The full model
estimates all path coefficients freely. The null model is the same model
estimated subject to the constraint that the effect of interest is zero,
for example \\H_0\\: \beta_1 \beta_2 = 0\\ for the indirect effect of
\\X\\ on \\Y\\ through one mediator. Writing \\D =
-2\\\mathrm{log\\likelihood}\\ for the deviance of each fitted model,
the test statistic is \$\$\mathrm{LRT}\_{\mathrm{MBCO}} =
D\_{\mathrm{null}} - D\_{\mathrm{full}},\$\$ which has a large sample
chi square distribution with degrees of freedom equal to the number of
constraints (here 1) when the null hypothesis is true (Wilks, 1938). The
resulting *p*-value is a continuous measure of the compatibility of the
data with the null model, not merely a reject or fail-to-reject
decision, and in the simulations of Tofighi and Kelley (2020) the test
controls the Type I error rate more robustly than the tests based on
confidence intervals in common use, especially when the indirect effect
is truly zero.

**Why the constraint is nonlinear.** Null hypotheses about single
parameters are linear constraints: fixing \\\beta_3 = 0\\ restricts the
parameter space to a flat hyperplane, something every structural
equation modeling program does by fixing the parameter and refitting.
The null hypothesis of no mediation is different in kind. The constraint
\\\beta_1 \beta_2 = 0\\ involves a product of parameters, so it is a
nonlinear function of the parameter vector, and its solution set is not
a hyperplane but the union of two hyperplanes: the set where \\\beta_1 =
0\\ (with \\\beta_2\\ free) and the set where \\\beta_2 = 0\\ (with
\\\beta_1\\ free). Mediation is absent in infinitely many ways, and no
single fixed parameter expresses them all: fixing \\\beta_1 = 0\\ alone
tests a more restrictive hypothesis than \\H_0\\: \beta_1 \beta_2 = 0\\.
Maximizing the likelihood over such a constraint set requires an
optimizer that handles general nonlinear equality constraints
(sequential quadratic programming methods such as SLSQP and NPSOL, or
CSOLNP); among R packages for structural equation modeling, OpenMx
provides them, which is why the model is fit there even though it is
specified in lavaan syntax. The same geometry explains why Wald-type
(delta method) tests of a product behave badly near the null: at
\\\beta_1 = \beta_2 = 0\\ the two hyperplanes intersect, the gradient of
\\\beta_1 \beta_2\\ vanishes, and the usual normal approximation for the
product breaks down. The MBCO procedure sidesteps that approximation by
comparing maximized likelihoods directly.

**Local solutions and the search strategy.** Because the null set is a
union of surfaces, the constrained deviance surface generally has one
local minimum per branch (one where the mediator does not respond to
\\X\\, one where the outcome does not respond to the mediator, and so on
along longer chains). A constrained optimizer started from the
full-model estimates can converge to whichever branch is nearest rather
than to the branch that fits best. The likelihood ratio test is defined
by the globally best-fitting null model, so `mediation_mbco()` refits
each null model from several starting configurations (the full-model
estimates, and the estimates with each coefficient entering the
constrained effect set to zero in turn), verifies that each candidate
solution actually satisfies the constraint, and keeps the feasible
solution with the smallest deviance. In the memory example below this
matters: for the single mediator model the best-fitting null model sets
the imagery-to-recall path to zero (\\\mathrm{LRT}\_{\mathrm{MBCO}} =
72.54\\), while the branch that sets the instruction-to-imagery path to
zero fits worse (\\\mathrm{LRT}\_{\mathrm{MBCO}} = 175.77\\, the value
reported for this example in Tofighi & Kelley, 2020, whose optimizer
stopped on that branch); either way the null model is overwhelmingly
incompatible with the data, so the substantive conclusion is the same.

**Effects estimated and tested.** With `x` and `y` resolved, the
function enumerates every directed pathway from `x` to `y` along
regression (`~`) paths. Each pathway through at least one intermediate
variable contributes a specific indirect effect, the product of its path
coefficients, reported as `indirect_via_` followed by the intermediate
variable names. The direct effect is the `x` to `y` coefficient when
that path is in the model. The total indirect effect is the sum of the
specific indirect effects (reported when there are two or more), and the
total effect is the direct effect plus the total indirect effect
(reported when the direct path is in the model). Quantities defined in
the model syntax via `:=` and any `hypotheses` are estimated and tested
the same way, so contrasts of indirect effects, proportions mediated, or
any other smooth function of parameters can be examined; each is tested
against zero, so write an equality of two effects as their difference.
Every reported row carries its estimate, a delta method standard error
(via [`mxSE`](https://rdrr.io/pkg/OpenMx/man/mxSE.html)), the
`ci_method` confidence interval, and the MBCO likelihood ratio test with
its degrees of freedom and *p*-value.

**Multiple groups (moderated mediation across groups).** With `group`
(or list-form `S`, `M`, `N`), the model is fit as a multiple-group SEM:
the same structure in every group, with group-specific parameters. Every
enumerated effect is then estimated within each group (its term carries
the group label as a suffix), and for each effect the difference from
the reference group (the first group label) is estimated and tested with
the same constrained-optimization machinery, since a difference of two
products is itself a nonlinear function of the parameters. That
difference test is moderated mediation with a categorical moderator:
"does the indirect effect differ across groups?" is exactly the
between-group contrast of the conditional indirect effects. lavaan's
labeling rules decide what varies: a single label such as `b1` on a path
applies to every group and therefore imposes cross-group equality (the
corresponding difference is identically zero and its row is dropped); to
let a path differ by group, leave it unlabeled or give per-group labels
with the vector form `c("b1_f", "b1_m")*x`. With more than two groups,
each non-reference group is compared with the reference group.

**Moderated mediation, probed.** With `moderator`, every regression
coefficient along a pathway becomes a linear function of the moderator
wherever the model contains the matching interaction (stated as `x:w` in
the syntax, or as a product column among the predictors). A pathway
effect, the product of its edge coefficients, is then a polynomial in
the moderator, and the function derives that polynomial symbolically.
Three kinds of rows follow for each moderated effect. First, the
conditional effect at each probe value, each with its own confidence
interval and MBCO test. Second, the polynomial's moderator coefficients:
for a pathway moderated in one place the single such coefficient is
exactly the index of moderated mediation (Hayes, 2015), here tested by
likelihood ratio rather than bootstrap; a pathway moderated in several
places gets one row per power of the moderator. Third, for a pathway
moderated in several places (so the conditional effect is curved in the
moderator), a joint constancy test of all moderator coefficients at
once, with as many degrees of freedom as constraints, asking whether the
pathway effect depends on the moderator at all; as a joint test its row
reports no scalar estimate. Unmoderated effects in the same model keep
their single rows. The null set of a conditional-effect constraint is
again a union of branches (one edge's conditional coefficient or
another's must vanish at the probed value), and the null fits are
started on each branch, and at the model with all interactions removed
for the moderation tests, so the reported statistics reflect the
best-fitting null models. For example:

    mediation_mbco("m ~ x + w + x:w
                    y ~ m + x + w",
                   data = d, x = "x", y = "y", moderator = "w")

[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
draws the same conditional effects as curves over the moderator's whole
range with a confidence band, the visual companion to the probe rows.
Bespoke conditional quantities can still be written directly with `:=`
definitions or `hypotheses` when the built-in probing does not cover
them.

**Moderated mediation and mediated moderation.** Moderated mediation
asks whether an indirect effect depends on a moderator \\W\\: the
conditional indirect effect varies with \\w\\ (Muller, Judd, & Yzerbyt,
2005; Preacher, Rucker, & Hayes, 2007). Mediated moderation asks whether
an observed \\X \times W\\ interaction on \\Y\\ is transmitted through
the mediator. The two share their algebra: with first-stage moderation
\\M = a_1 X + a_2 W + a_3 XW + e_M\\ and \\Y = b M + \ldots\\, the
quantity \\a_3 b\\ is at once the index of moderated mediation (Hayes,
2015) and the indirect effect of the product term through \\M\\; what
differs is the question and the reporting emphasis. This function
reports the conditional-indirect-effect framing: the probe rows and the
moderation rows above. For a categorical moderator, use `group`.

**Refusals, warnings, and judgment calls.** The function stops where a
computed number would be mislabeled: a nonrecursive (feedback)
regression structure, categorical-endogenous syntax (thresholds),
explicit `==` constraints, a model with no indirect pathway between `x`
and `y`. It warns where the data make trouble detectable: an endogenous
variable with only two distinct values (a binary mediator or outcome,
where the product of coefficients is not the causal indirect effect), an
interaction term among the predictors that no declared `moderator`
accounts for (declaring the moderator resolves the warning by probing
the moderation), an interaction included without its matching main
effect (the principle of marginality), null models that converge
imperfectly, profile bounds that fail. What it cannot check is left to
the analyst, and stated rather than assumed: the no omitted confounder
assumption, linearity, and the causal direction of the arrows come from
the design, not from the fit.

**Information criteria.** The columns `delta_aic` and `delta_bic` report
AIC and BIC for the null model minus the same criterion for the full
model; positive values favor the full model (the effect improves fit by
more than the parsimony penalty). A scalar equality constraint reduces
the effective number of free parameters by one, so the differences equal
\\\mathrm{LRT}\_{\mathrm{MBCO}} - 2\\ for the AIC and
\\\mathrm{LRT}\_{\mathrm{MBCO}} - \log N\\ for the BIC. (The OpenMx
summary counts the same number of estimated parameters in both models,
so its printed AIC difference equals the likelihood ratio statistic;
DMAR counts the constraint against the null model, consistent with the
degrees of freedom of the test.)

**Change in explained variance.** Following the reporting recommendation
of Tofighi and Kelley (2020), the function computes \\R^2\\ for every
endogenous variable under the full model (the `"R2"` attribute) and the
drop in each \\R^2\\ under every null model (the `"delta_R2"` attribute,
variables by tested effects), so the fit cost of removing an effect can
be read as a change in effect size and not only as a test statistic.

**Summary statistics input.** With complete data the multivariate normal
log likelihood depends on the data only through the sample means and the
sample covariance matrix. When `S` (and optionally `M`) is supplied, the
function therefore constructs an internal data set with exactly those
moments (via [`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) with
`empirical = TRUE`) and proceeds as with raw data; the estimates,
likelihood ratio tests, and confidence intervals are identical to what
the raw data would give, whatever internal data set realizes the
moments. `S` is treated as the unbiased (divisor \\N - 1\\) covariance
matrix, the form reported in articles. This is how a published mediation
analysis can be reproduced, and its hypotheses re-tested with the MBCO
procedure, from a table of descriptive statistics alone.

**Assumptions.** The causal reading of any mediation analysis rests on
the no omitted confounder assumption for the predictor to mediator and
mediator to outcome relations, in addition to the usual distributional
assumptions (residuals multivariate normal, linear relations, no
treatment by mediator interaction); randomizing \\X\\ supports the first
link but not the second (MacKinnon, 2008; Tofighi & Kelley, 2020). With
a binary randomized \\X\\ the variable enters the model as numeric 0/1,
its exogenous variance freely estimated; the normality assumption
applies to the residuals of the endogenous variables.

## References

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*(4), 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Tofighi, D., & Kelley, K. (2020). Indirect effects in sequential
mediation models: Evaluating methods for hypothesis testing and
confidence interval formation. *Multivariate Behavioral Research,
55*(2), 188–210.
[doi:10.1080/00273171.2019.1618545](https://doi.org/10.1080/00273171.2019.1618545)

Hayes, A. F. (2015). An index and test of linear moderated mediation.
*Multivariate Behavioral Research, 50*(1), 1–22.
[doi:10.1080/00273171.2014.962683](https://doi.org/10.1080/00273171.2014.962683)

MacKinnon, D. P. (2008). *Introduction to statistical mediation
analysis*. Erlbaum.

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
limits for the indirect effect: Distribution of the product and
resampling methods. *Multivariate Behavioral Research, 39*(1), 99–128.
[doi:10.1207/s15327906mbr3901_4](https://doi.org/10.1207/s15327906mbr3901_4)

MacKinnon, D. P., Valente, M. J., & Wurpts, I. C. (2018). Benchmark
validation of statistical models: Application to mediation analysis of
imagery and memory. *Psychological Methods, 23*(4), 654–671.
[doi:10.1037/met0000174](https://doi.org/10.1037/met0000174)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Muller, D., Judd, C. M., & Yzerbyt, V. Y. (2005). When moderation is
mediated and mediation is moderated. *Journal of Personality and Social
Psychology, 89*(6), 852–863.
[doi:10.1037/0022-3514.89.6.852](https://doi.org/10.1037/0022-3514.89.6.852)

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2016). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika, 81*(2), 535–549.
[doi:10.1007/s11336-014-9435-8](https://doi.org/10.1007/s11336-014-9435-8)

Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based
confidence intervals in genetic models. *Behavior Genetics, 27*(2),
113–120.
[doi:10.1023/A:1025681223921](https://doi.org/10.1023/A%3A1025681223921)

Preacher, K. J., Rucker, D. D., & Hayes, A. F. (2007). Addressing
moderated mediation hypotheses: Theory, methods, and prescriptions.
*Multivariate Behavioral Research, 42*(1), 185–227.
[doi:10.1080/00273170701341316](https://doi.org/10.1080/00273170701341316)

Wilks, S. S. (1938). The large-sample distribution of the likelihood
ratio for testing composite hypotheses. *Annals of Mathematical
Statistics, 9*(1), 60–62.
[doi:10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360)

## See also

[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
for the conditional-effect display of a moderated analysis;
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md) for
the regression-based simple mediation model with bootstrap intervals;
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)
for planning;
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)
for the delta method variance.

Other mediation:
[`mediate()`](https://yelleknek.github.io/DMAR/reference/mediate.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A quick simulated single-mediator example. Requires the OpenMx and
# lavaan packages to be installed.
set.seed(113)
n <- 200
x <- rnorm(n)
m <- 0.5 * x + rnorm(n, 0, sqrt(1 - 0.25))
y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
d <- data.frame(x = x, m = m, y = y)
mediation_mbco("m ~ x \n y ~ m + x", data = d, x = "x", y = "y",
               ci_method = "wald")
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of x on y; N = 200; optimizer: SLSQP
#>   Full model: deviance = 1551.498904, AIC = 1569.498904, BIC = 1599.183761
#>   R2 (full model): m = 0.212, y = 0.280
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway               term           estimate se     ci_lower ci_upper lrt  df
#>  x -> y (all pathways) total_effect   0.375    0.0633 0.251    0.499    32.4 1 
#>  x -> y                direct_effect  0.194    0.0656 0.0652   0.322    8.54 1 
#>  x -> m -> y           indirect_via_m 0.182    0.0391 0.105    0.258    33.2 1 
#>  p_value  delta_aic delta_bic
#>  < 0.0001 30.4      27.1     
#>  0.0035   6.54      3.24     
#>  < 0.0001 31.2      27.9     
#> 
#> Confidence level: 95%

# \donttest{
# Replicate the memory experiment analyses of Tofighi and Kelley
# (2020) from the summary statistics in their Table 1 (data from
# MacKinnon, Valente, & Wurpts, 2018; N = 369). Instruction is 1 for
# imagery rehearsal instructions and 0 for repetition instructions.
vars <- c("instruction", "imagery", "repetition", "recall")
sds  <- c(0.50, 2.96, 2.84, 3.40)
R_tk <- matrix(c(1.00,  .62, -.67,  .32,
                  .62, 1.00, -.56,  .51,
                 -.67, -.56, 1.00, -.28,
                  .32,  .51, -.28, 1.00), 4, 4,
               dimnames = list(vars, vars))
S_tk <- outer(sds, sds) * R_tk
M_tk <- c(instruction = 0.51, imagery = 5.66, repetition = 6.08,
          recall = 12.07)

# Single-mediator model: instruction -> imagery -> recall.
single <- "
  imagery ~ b1*instruction
  recall  ~ b2*imagery + b3*instruction
"
mediation_mbco(single, S = S_tk, M = M_tk, N = 369,
               x = "instruction", y = "recall",
               ci_method = "monte_carlo", seed = 113)
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of instruction on recall; N = 369; optimizer: SLSQP
#>   Full model: deviance = 4040.806567, AIC = 4058.806567, BIC = 4094.003737
#>   R2 (full model): imagery = 0.384, recall = 0.260
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway                              term                 estimate se   
#>  instruction -> recall (all pathways) total_effect         2.18     0.335
#>  instruction -> recall                direct_effect        0.042    0.388
#>  instruction -> imagery -> recall     indirect_via_imagery 2.13     0.279
#>  ci_lower ci_upper lrt    df p_value  delta_aic delta_bic
#>  1.51     2.84     39.9   1  < 0.0001 37.9      34       
#>  -0.724   0.804    0.0117 1  0.9139   -1.99     -5.9     
#>  1.61     2.7      71.3   1  < 0.0001 69.3      65.4     
#> 
#> Confidence level: 95%
# The indirect effect is about 2.1 words (the paper reports 2.121,
# SE = 0.276, 95% Monte Carlo CI [1.600, 2.682]).

# Parallel two-mediator model, with the contrast of the two specific
# indirect effects (the paper's Research Questions 2 and 3).
parallel <- "
  imagery    ~ b1*instruction
  repetition ~ b3*instruction
  recall     ~ b2*imagery + b4*repetition + b5*instruction
  imagery ~~ repetition
"
mediation_mbco(parallel, S = S_tk, M = M_tk, N = 369,
               x = "instruction", y = "recall",
               hypotheses = c(imagery_minus_repetition =
                 "indirect_via_imagery - indirect_via_repetition"),
               ci_method = "monte_carlo", seed = 113)
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of instruction on recall; N = 369; optimizer: SLSQP
#>   Full model: deviance = 5613.915747, AIC = 5641.915747, BIC = 5696.6669
#>   R2 (full model): imagery = 0.384, repetition = 0.449, recall = 0.260
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway                                        term                    
#>  instruction -> recall (all pathways)           total_effect            
#>  instruction -> recall                          direct_effect           
#>  instruction -> recall (all indirect pathways)  total_indirect          
#>  instruction -> imagery -> recall               indirect_via_imagery    
#>  instruction -> repetition -> recall            indirect_via_repetition 
#>  indirect_via_imagery - indirect_via_repetition imagery_minus_repetition
#>  estimate se    ci_lower ci_upper lrt    df p_value  delta_aic delta_bic
#>  2.18     0.335 1.51     2.84     39.9   1  < 0.0001 37.9      34       
#>  0.0943   0.447 -0.799   0.978    0.0445 1  0.8329   -1.96     -5.87    
#>  2.08     0.356 1.4      2.8      35.4   1  < 0.0001 33.4      29.5     
#>  2.15     0.286 1.6      2.74     68.1   1  < 0.0001 66.1      62.2     
#>  -0.0669  0.284 -0.61    0.488    0.0556 1  0.8136   -1.94     -5.86    
#>  2.22     0.444 1.34     3.09     25.7   1  < 0.0001 23.7      19.8     
#> 
#> Confidence level: 95%
# The indirect effect through repetition is near zero (the paper
# reports LRT = 0.083, p = .773), while the contrast shows the
# imagery pathway is larger (the paper reports LRT = 25.828,
# difference = 2.222, SE = 0.445).

# Moderated mediation across groups: the same model in two groups,
# with the between-group difference of every effect tested. Leaving
# the paths unlabeled lets them differ by group.
set.seed(113)
n <- 150
two_groups <- rbind(
  within(data.frame(x = rnorm(n), condition = "treatment"), {
    m <- 0.6 * x + rnorm(n)
    y <- 0.5 * m + 0.2 * x + rnorm(n)
  }),
  within(data.frame(x = rnorm(n), condition = "control"), {
    m <- 0.2 * x + rnorm(n)
    y <- 0.1 * m + 0.2 * x + rnorm(n)
  }))
mediation_mbco("m ~ x \n y ~ m + x", data = two_groups,
               group = "condition", x = "x", y = "y",
               ci_method = "wald")
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of x on y; N = 300; optimizer: SLSQP
#>   Groups (reference first): treatment, control
#>   Full model: deviance = 2514.474242, AIC = 2550.474242, BIC = 2617.142327
#>   R2 (full model): m (treatment) = 0.231, y (treatment) = 0.356, m (control) = 0.000, y (control) = 0.105
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway                                    
#>  x -> y (all pathways) [treatment]          
#>  x -> y (all pathways) [control]            
#>  x -> y (all pathways) [control - treatment]
#>  x -> y [treatment]                         
#>  x -> y [control]                           
#>  x -> y [control - treatment]               
#>  x -> m -> y [treatment]                    
#>  x -> m -> y [control]                      
#>  x -> m -> y [control - treatment]          
#>  term                                   estimate se     ci_lower ci_upper
#>  total_effect_treatment                 0.518    0.0859 0.349    0.686   
#>  total_effect_control                   0.232    0.0701 0.0945   0.369   
#>  total_effect_control_minus_treatment   -0.286   0.111  -0.503   -0.0682 
#>  direct_effect_treatment                0.26     0.0876 0.0879   0.431   
#>  direct_effect_control                  0.235    0.0687 0.1      0.369   
#>  direct_effect_control_minus_treatment  -0.0249  0.111  -0.243   0.193   
#>  indirect_via_m_treatment               0.258    0.057  0.146    0.37    
#>  indirect_via_m_control                 -0.00274 0.014  -0.0302  0.0247  
#>  indirect_via_m_control_minus_treatment -0.261   0.0587 -0.376   -0.146  
#>  lrt    df p_value  delta_aic delta_bic
#>  32.5   1  < 0.0001 30.5      26.8     
#>  10.6   1  0.0012   8.56      4.86     
#>  6.56   1  0.0104   4.56      0.853    
#>  8.53   1  0.0035   6.53      2.83     
#>  11.2   1  0.0008   9.23      5.53     
#>  0.05   1  0.8231   -1.95     -5.65    
#>  33.4   1  < 0.0001 31.4      27.7     
#>  0.0363 1  0.8490   -1.96     -5.67    
#>  29.8   1  < 0.0001 27.8      24.1     
#> 
#> Confidence level: 95%

# Probing moderated mediation with a continuous moderator: the
# effect of x on m depends on w (the x:w term). The output has the
# conditional indirect effect at the moderator's mean and one
# standard deviation either side, the index of moderated mediation
# (the change in the indirect effect per unit w), and its MBCO
# likelihood ratio test.
set.seed(113)
n <- 300
x <- rnorm(n)
w <- rnorm(n)
m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
d_mod <- data.frame(x = x, w = w, m = m, y = y)
mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d_mod,
               x = "x", y = "y", moderator = "w",
               ci_method = "wald")
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of x on y; N = 300; optimizer: SLSQP
#>   Full model: deviance = 4194.186376, AIC = 4232.186376, BIC = 4302.558243
#>   R2 (full model): m = 0.403, y = 0.421
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway                                  term                      estimate
#>  x -> y (all pathways) at w = -0.8997     total_effect_at_low       0.31    
#>  x -> y (all pathways) at w = 0.1036      total_effect_at_mean      0.537   
#>  x -> y (all pathways) at w = 1.107       total_effect_at_high      0.764   
#>  x -> y (all pathways), change per unit w total_effect_moderation   0.226   
#>  x -> y                                   direct_effect             0.219   
#>  x -> m -> y at w = -0.8997               indirect_via_m_at_low     0.0905  
#>  x -> m -> y at w = 0.1036                indirect_via_m_at_mean    0.318   
#>  x -> m -> y at w = 1.107                 indirect_via_m_at_high    0.545   
#>  x -> m -> y, change per unit w           indirect_via_m_moderation 0.226   
#>  se     ci_lower ci_upper lrt  df p_value  delta_aic delta_bic
#>  0.0732 0.166    0.453    16   1  < 0.0001 14        10.3     
#>  0.0619 0.416    0.658    67.1 1  < 0.0001 65.1      61.4     
#>  0.071  0.625    0.903    113  1  < 0.0001 111       107      
#>  0.0369 0.154    0.299    52.2 1  < 0.0001 50.2      46.5     
#>  0.0625 0.0967   0.342    12   1  0.0005   10        6.35     
#>  0.0446 0.00307  0.178    4.24 1  0.0394   2.24      -1.46    
#>  0.0427 0.234    0.401    94.5 1  < 0.0001 92.5      88.8     
#>  0.0662 0.415    0.674    137  1  < 0.0001 135       131      
#>  0.0369 0.154    0.299    52.2 1  < 0.0001 50.2      46.5     
#> 
#> Confidence level: 95%
# }
```
