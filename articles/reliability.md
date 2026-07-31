# Reliability in DMAR: Choosing the Coefficient, Choosing the Interval

``` r

library(DMAR)
```

Reporting the reliability of a composite score involves two separable
decisions, following the framework of Kelley and Pornprasertmanit
(2016): choosing the coefficient and choosing the interval.

- **Choosing the coefficient.** The coefficient reported is a claim
  about the measurement model and about which composite is scored. Each
  coefficient below estimates the proportion of composite variance
  attributable to the common factor, under different commitments about
  the model. Choosing a coefficient is choosing which commitments the
  analysis makes.
- **Choosing the interval.** Given the coefficient, the confidence
  interval method is an empirical performance question: which procedures
  cover the population value at their nominal rate across realistic
  conditions. The Monte Carlo studies of Kelley and
  Pornprasertmanit (2016) addressed that question, and the DMAR defaults
  encode their recommendations.

The map for choosing the coefficient:

| Function | The composite | The claim being made |
|----|----|----|
| [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md) | unit-weighted sum | one factor with *equal loadings* (essential tau-equivalence) |
| [`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md) | sum of 0/1 items | alpha’s claim, for binary items |
| [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md) (robust omega, the default) | unit-weighted sum | one congeneric factor for the *numerator only*; the composite variance is estimated from the data |
| `reliability_omega(denominator = "model_implied")` | unit-weighted sum | one congeneric factor; the fitted model reproduces the composite variance |
| [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md) | sum of observed ordinal categories | a probit threshold model links items to one factor |
| [`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md) | *optimally weighted* composite | one congeneric factor; a different composite, not a different assumption |

The vignette walks the table top to bottom on one running example, then
turns to the interval.

## A Running Example

A six-item congeneric scale with unequal loadings, which is the normal
state of real items:

``` r

set.seed(113)
N <- 300
lambda <- c(0.4, 0.5, 0.6, 0.7, 0.75, 0.8)
eta <- rnorm(N)
items <- sweep(matrix(rep(eta, 6), N, 6), 2, lambda, `*`) +
  matrix(rnorm(N * 6), N, 6) %*% diag(sqrt(1 - lambda^2))
colnames(items) <- paste0("y", seq_len(6))
```

## Alpha Assumes Equal Loadings

Coefficient alpha equals the population reliability when the items are
essentially tau-equivalent, that is, when every item loads equally on
the common factor. When loadings differ, alpha underestimates. On the
running example the loadings range from .4 to .8, so alpha sits a little
below omega:

``` r

res_alpha <- reliability_alpha(data = items, ci_method = "none")
res_omega <- reliability_omega(data = items,
                               denominator = "model_implied",
                               ci_method = "none")
c(alpha = res_alpha$value[res_alpha$term == "estimate"],
  omega = res_omega$value[res_omega$term == "estimate"])
#>     alpha     omega 
#> 0.7976157 0.8067366
```

The gap here is modest because the loadings, though unequal, are all
positive and of similar magnitude. The gap becomes large when the
tau-equivalence assumption fails badly. The first five `mtcars`
variables, used here as a five-variable example, are not positively
keyed and several covary negatively; the negative covariances lower
alpha but not omega, whose signed loadings recover the common variance
regardless of keying:

``` r

S_cars <- cov(mtcars[, 1:5])
a_cars <- reliability_alpha(S = S_cars, N = 32, ci_method = "none")
o_cars <- reliability_omega(S = S_cars, N = 32,
                            denominator = "model_implied",
                            ci_method = "none")
c(alpha = a_cars$value[a_cars$term == "estimate"],
  omega = o_cars$value[o_cars$term == "estimate"])
#>     alpha     omega 
#> 0.4671591 0.9050374
```

Alpha here answers under an assumption these data violate. The pattern
generalizes: each coefficient is correct under its stated model, and the
choice among them is a choice among models.

## The Denominator Choice in Omega

Coefficient omega requires an estimate of the composite’s total variance
in its denominator, and
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
offers two. The default, `denominator = "observed"` (robust omega),
estimates the composite variance directly from the data, which is
consistent for the actual variance of the scale score whether or not the
one-factor model is right. The alternative,
`denominator = "model_implied"`, uses the total variance reproduced by
the fitted single-factor model, which is correct exactly when the
one-factor model reproduces the composite variance.

The coefficient was built for the following situation. Model misfit is
treated as a set of minor common factors, visible as residual
correlations: the scale is then not perfectly unidimensional, but a
single-factor approximation is retained. The coefficient isolates the
variance attributable to the general factor only, excluding the minor
factors, and expresses it relative to the observed variance of the
unweighted composite. Kelley and Pornprasertmanit (2016) named it
hierarchical omega after that hierarchical factor logic; upon
reflection, we would have named it observed omega or robust omega, and
DMAR uses *robust omega*. Two qualifications keep the name accurate. The
robustness is to misspecification of the total variance only, because
the numerator remains model based under either denominator; and robust
omega is distinct both from outlier-robust estimation of alpha and omega
(Zhang & Yuan, 2016) and from the robust maximum-likelihood standard
errors available through `ci_method`.

When the one-factor model holds, the two agree, as on the running
example:

``` r

o_mi <- reliability_omega(data = items,
                          denominator = "model_implied",
                          ci_method = "none")
o_ob <- reliability_omega(data = items, ci_method = "none")
c(model_implied = o_mi$value[o_mi$term == "estimate"],
  observed      = o_ob$value[o_ob$term == "estimate"])
#> model_implied      observed 
#>     0.8067366     0.8056993
```

Now contaminate the scale with a nuisance doublet, a shared specific
factor on the first two items, the sort of minor structure real item
sets carry:

``` r

doublet <- rnorm(N)
items_d <- items
items_d[, 1] <- items[, 1] + 0.6 * doublet
items_d[, 2] <- items[, 2] + 0.6 * doublet

o_mi_d <- reliability_omega(data = items_d,
                            denominator = "model_implied",
                            ci_method = "none")
o_ob_d <- reliability_omega(data = items_d, ci_method = "none")
c(model_implied = o_mi_d$value[o_mi_d$term == "estimate"],
  observed      = o_ob_d$value[o_ob_d$term == "estimate"])
#> model_implied      observed 
#>     0.7846721     0.7693093
```

The two now diverge, and only robust omega still answers the question a
reliability coefficient is supposed to answer, the proportion of the
variance of the composite you actually computed that is attributable to
the common factor. The divergence itself is diagnostic: a nontrivial gap
between the two denominators says the one-factor model is not
reproducing the composite variance, and the right response to a large
gap is to model the structure, not to choose a denominator. Note also
what the observed denominator does not fix: the fitted loadings absorb
part of the doublet under either setting, so checking the one-factor
model (for example with
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md))
remains part of reporting any omega.

## Categorical Omega for Ordered Items

When items are ordered categories (Likert responses, symptom counts),
the composite people score is a sum of those observed categories, and
the category thresholds make the item-factor relationship nonlinear.
Categorical omega (Green & Yang, 2009) fits the probit threshold model
and returns the reliability of the categorical sum score on its own
metric. Cutting the running example’s items into four skewed categories:

The threshold configuration decides how much this matters. First the
benign case, every item cut at the same thresholds:

``` r

items_same <- apply(items, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1, -0.2, 0.6, Inf))))
colnames(items_same) <- colnames(items)

oc_same <- reliability_omega_categorical(data = items_same, ci_method = "none")
o_same  <- reliability_omega(data = items_same, ci_method = "none")
c(categorical_omega     = oc_same$value[oc_same$term == "estimate"],
  treated_as_continuous = o_same$value[o_same$term == "estimate"])
#>     categorical_omega treated_as_continuous 
#>             0.7820753             0.7755024
```

With similar threshold patterns across items the two approaches nearly
agree, which is precisely the condition under which Kelley and
Pornprasertmanit (2016, Study 3) found the continuous treatment held its
interval coverage. In real item sets the threshold patterns typically
differ across items: some items skew one way, some the other. Cutting
the same underlying responses with threshold patterns that differ across
items:

``` r

breaks_hi <- c(-Inf,  0.5,  1.2,  1.9, Inf)
breaks_lo <- c(-Inf, -1.9, -1.2, -0.5, Inf)
breaks_by_item <- list(breaks_hi, breaks_hi, breaks_hi,
                       breaks_lo, breaks_lo, breaks_lo)
items_cat <- sapply(seq_len(6), function(j)
  as.integer(cut(items[, j], breaks = breaks_by_item[[j]])))
colnames(items_cat) <- colnames(items)

oc <- reliability_omega_categorical(data = items_cat, ci_method = "none")
o_as_cont <- reliability_omega(data = items_cat, ci_method = "none")
c(categorical_omega     = oc$value[oc$term == "estimate"],
  treated_as_continuous = o_as_cont$value[o_as_cont$term == "estimate"])
#>     categorical_omega treated_as_continuous 
#>             0.7190046             0.6755575
```

The two now diverge, here by about 0.04, with the continuous treatment
understating on this configuration. The direction and size of the
continuous treatment’s error depend on the threshold configuration,
which is exactly what makes it untrustworthy as a general practice; in
Study 3 its interval coverage failed as threshold patterns diverged
across items. Categorical omega is constructed on the metric of the
summed categories themselves, so it is the coefficient to report for
ordinal items, which is why
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
routes integer-coded, few-category items to
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
automatically.

Categorical omega and robust omega share a design principle in Kelley
and Pornprasertmanit (2016): in both, the total variance in the
denominator is not taken from the fitted factor model. Robust omega
estimates it from the sample covariances; categorical omega assembles it
from the saturated polychoric correlations and the thresholds. In each
case the fitted single factor supplies only the numerator, so the
coefficient stays anchored to the composite that is actually scored even
when the factor model is an approximation.

## Coefficient H: A Different Composite

Coefficient H (Hancock & Mueller, 2001) makes the same congeneric claim
as omega but changes the composite: it is the reliability of the
*optimally weighted* combination of the items, the best possible
composite the items can form, rather than the unit-weighted sum. The two
coincide when loadings are equal and diverge as loadings spread, because
optimal weights exploit the stronger items:

``` r

std <- cfa_1(S = cov(items), N = N, output = "standardized")
lam <- std[std$op == "=~", ]
reliability_H(loadings = lam$est.std, se_loadings = lam$se)
```

| term          | value    |
|:--------------|:---------|
| reliability_H | 0.847    |
| lower_limit   | 0.819    |
| upper_limit   | 0.872    |
| var_H         | 0.000182 |

Confidence level: 95%

Report omega when people will score the scale by summing; report H when
scores will come from the factor model itself (or to see the ceiling the
item set could reach with optimal weighting). Reporting both, with one
sentence on which composite each describes, is often the most
informative choice under unequal loadings. Reliability coefficients
matter because researchers score and use composites; an analysis
conducted entirely within the SEM framework, relating latent variables
directly, has less need of them.

KR-20 completes the family as the binary-item special case of alpha,
provided for the classical test theory literature that names it
separately;
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md)
and
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
agree on 0/1 data.

## Choosing the Interval

Choosing the interval is an empirical question, and it is the subject of
Kelley and Pornprasertmanit (2016): across coefficients, sample sizes,
loading patterns, distributions, and (for categorical items) threshold
patterns, which interval procedures cover at their nominal rate. The
DMAR defaults are their recommendations:

| Coefficient | Default behavior | Recommended interval and basis |
|----|----|----|
| alpha | Bonett (2002) interval | well behaved among the closed forms under normality |
| omega, model implied | Wald interval with the robust ML standard error (`"mlr"`) | the best performing closed form when item distributions deviate from normality |
| robust omega | point estimate; a message names the bootstrap call | percentile or BCa bootstrap, the pairing recommended when unidimensionality is approximate; no closed form exists for the observed denominator |
| categorical omega | point estimate; a message names the bootstrap call | BCa bootstrap, which held coverage across threshold patterns, category counts, and sample sizes |

A bootstrap is never run unless it is requested. When one is requested,
`B = 10000` replications is the default.

Every function exposes its full menu through `ci_method` (Feldt, Fisher,
Bonett, Hakstian-Whalen, ML, robust ML, ADF, their logistic transformed
variants, and the bootstrap family, as applicable), so a reviewer’s
requested method is available; the default is the one with the
simulation evidence behind it. Two practical notes. The bootstrap
methods need raw data, not a covariance matrix, and take a `seed`
argument for reproducibility. And a bootstrap interval on the
categorical coefficient refits an ordinal factor model per replicate, so
it is the slowest interval in the family; the default `B = 10000` is a
deliberate accuracy choice, worth the wait for a final analysis and
worth reducing while exploring.

``` r

# The default: robust omega, point estimate, with the message naming
# the bootstrap call.
reliability_omega(data = items)
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
```

| term        | value |
|:------------|:------|
| estimate    | 0.806 |
| se          | NA    |
| lower_limit | NA    |
| upper_limit | NA    |
| conf_level  | 0.95  |
| N           | 300   |
| N_complete  | 300   |
| J           | 6     |

``` r


# The recommended interval, requested explicitly (B reduced here to
# keep the vignette fast; the default is B = 10000).
reliability_omega(data = items, ci_method = "percentile",
                  B = 500, seed = 113)
```

| term        | value  |
|:------------|:-------|
| estimate    | 0.806  |
| se          | 0.0169 |
| lower_limit | 0.771  |
| upper_limit | 0.838  |
| conf_level  | 0.95   |
| N           | 300    |
| N_complete  | 300    |
| J           | 6      |

``` r


# The model implied denominator has a closed form and reports its
# robust ML Wald interval by default.
reliability_omega(data = items, denominator = "model_implied")
```

| term        | value  |
|:------------|:-------|
| estimate    | 0.807  |
| se          | 0.0165 |
| lower_limit | 0.774  |
| upper_limit | 0.839  |
| conf_level  | 0.95   |
| N           | 300    |
| N_complete  | 300    |
| J           | 6      |

## Missing Data: FIML and Auxiliary Variables

The [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
vignette makes the case at length: listwise deletion is unbiased only
when values are missing completely at random, and it is inefficient
always. The same argument applies to a reliability analysis, so
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
and
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
(and the
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
wrapper) take the same position
[`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) does. Two
arguments govern the treatment:

- `missing = c("listwise", "fiml")`. The default remains listwise
  deletion, so no previously computed result changes and the
  missing-data treatment is always a visible, deliberate choice.
  `"fiml"` keeps every case with at least one observed item and
  estimates by full information maximum likelihood, consistent and
  efficient under the weaker missing at random (MAR) assumption.
- `aux`, a character vector naming *auxiliary variables*: columns of
  `data` that are not part of the composite but are correlated with the
  items or with the reasons values are missing. They enter as saturated
  correlates (Graham, 2003): correlated freely with each other and with
  every item’s residual, never loading on the factor and never entering
  the composite. The measurement model is undisturbed while FIML uses
  their information, and missingness that depends on the auxiliary
  becomes MAR once the auxiliary is conditioned on (Collins, Schafer, &
  Kam, 2001). Supplying `aux` implies `missing = "fiml"`.

``` r

# Missingness on y2 that depends on an auxiliary z (MAR given z):
# listwise deletion is biased here, FIML with z is not.
z <- eta + rnorm(N, sd = 0.5)
d <- data.frame(items, z = z)
d$y2[runif(N) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA

# Listwise: the analysis quietly drops the incomplete rows.
reliability_alpha(data = d[, paste0("y", 1:6)])
```

| term        | value |
|:------------|:------|
| estimate    | 0.762 |
| se          | 0.106 |
| lower_limit | 0.707 |
| upper_limit | 0.807 |
| conf_level  | 0.95  |
| N           | 215   |
| N_complete  | 215   |
| J           | 6     |

``` r


# FIML with the auxiliary: every case with at least one observed item
# contributes, and z informs the estimation.
reliability_alpha(data = d, aux = "z")
```

| term        | value  |
|:------------|:-------|
| estimate    | 0.796  |
| se          | 0.0188 |
| lower_limit | 0.759  |
| upper_limit | 0.833  |
| conf_level  | 0.95   |
| N           | 300    |
| N_complete  | 215    |
| J           | 6      |

The returned table always reports both `N` (the cases used) and
`N_complete` (the complete cases), so the cost of listwise deletion is
legible at a glance. The classical (analytic) alpha under FIML is still
the classical coefficient: the same formula, applied to the FIML
estimate of the item covariance matrix rather than the complete-case
one. Interval methods that cannot be made correct under FIML (the
complete-data closed forms, ADF, the profile likelihood) are refused
with an explanation rather than silently reverting to listwise deletion;
the delta method (`"ml"`, `"mlr"`) and bootstrap intervals remain
available. Users of MBESS may recognize the `aux` idea from
[`ci.reliability()`](https://rdrr.io/pkg/MBESS/man/ci.reliability.html),
whose implementation no longer runs on current semTools; DMAR implements
the saturated correlates model directly.

## One Entry Point and the Broom Verbs

[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
dispatches on `type` (or auto-detects: ordered few-category items go to
categorical omega, otherwise omega) and forwards `denominator`; every
family member returns the same tidy shape with `tidy()` and `glance()`
methods:

``` r

res <- reliability(data = items, type = "omega",
                   denominator = "model_implied")
generics::tidy(res)
#>    term  estimate         se  ci_lower  ci_upper
#> 1 omega 0.8067366 0.01652757 0.7743432 0.8391301
generics::glance(res)
#>   coefficient  estimate         se  ci_lower  ci_upper conf_level nobs n_items
#> 1       omega 0.8067366 0.01652757 0.7743432 0.8391301       0.95  300       6
#>   ci_method
#> 1       mlr
```

## Planning the Study

Reliability estimation is also a design problem: a reliability
coefficient reported without a narrow interval is a weak claim.
[`ss_aipe_reliability()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
plans the sample size so the interval for alpha or omega achieves a
target width (with an `assurance` probability if desired), which closes
the loop between this family and DMAR’s accuracy in parameter estimation
(AIPE) tradition.

## What Reliability Coefficients Cannot Do

Reliability is a property of a test score in a particular population,
not a property of the items in the abstract. A high alpha or omega does
*not* mean:

- the items measure the construct of interest (that is *validity*, not
  reliability);
- the scale will achieve the same reliability in a different population;
- averaging across many items has fixed the problem of poor individual
  item quality.

Single-administration coefficients also book item-specific systematic
variance as error, so all of the coefficients here are lower bounds on
classical test-retest reliability. Use the family for what it is: an
internal-consistency analysis conditional on a substantive theory of
what the items mean.

## See Also

- [`?reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md),
  [`?reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`?reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  [`?reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  [`?reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
  [`?reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md).
- [`?ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
  for sample size planning.
- [`?cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) for
  the one-factor CFA underlying the omega family, and
  [`?average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
  and [`?htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for
  the validity companions.

## References

Bonett, D. G. (2002). Sample size requirements for testing and
estimating coefficient alpha. *Journal of Educational and Behavioral
Statistics, 27*, 335–340.

Collins, L. M., Schafer, J. L., & Kam, C. M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*, 330–351.

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*, 80–100.

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sorbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Lawrence
Erlbaum Associates.
