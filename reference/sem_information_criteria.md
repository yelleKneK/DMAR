# Information Criteria for Structural Equation Models

Computes the information criteria used to choose among structural
equation models fitted to the same data: AIC (Akaike, 1974), BIC
(Schwarz, 1978), the sample size adjusted BIC (SABIC; Sclove, 1987),
Haughton's (1988) BIC (HBIC), and the information matrix based BIC
(IBIC) and scaled unit information prior BIC (SPBIC) that Bollen, Ray,
Zavisca, and Harden (2012) proposed and Bollen, Harden, Ray, and Zavisca
(2014) evaluated for structural equation models. Each criterion is minus
twice the maximized log likelihood plus a penalty for the free
parameters, smaller values are preferred, and only the difference
between two models fitted to the same data under the same likelihood
carries meaning. The function takes either a fitted lavaan model, whose
log likelihood, parameter count, sample size, and (for IBIC and SPBIC)
information matrices it reads, or those quantities as summary input, so
that a published comparison can be reproduced.
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
compares a set of fitted models by one criterion. A fitted model
requires lavaan (in Suggests).

## Usage

``` r
sem_information_criteria(
  fit = NULL,
  log_likelihood = NULL,
  npar = NULL,
  N = NULL,
  measures = c("AIC", "BIC", "SABIC", "HBIC"),
  parameter_estimate = NULL,
  prior_mean = NULL,
  zero_parameters = NULL,
  information_observed = NULL,
  information_expected = NULL,
  information_scale = c("total", "per_observation"),
  spbic_dimension = NULL
)
```

## Arguments

- fit:

  A fitted lavaan model: a maximum likelihood fit of continuous
  variables under the normal likelihood (estimator ML or its robust
  variants MLM, MLMV, and MLR) that converged to an admissible solution,
  as
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
  judges it. See Details for what is refused. Supply either `fit` or the
  summary input below, not both.

- log_likelihood:

  For summary input, the maximized log likelihood \\\ell\\, including
  its constants (lavaan's `logl`).

- npar:

  For summary input, the number *d* of independent free parameters of
  the likelihood whose maximum is `log_likelihood`, a positive whole
  number, counted the same way for every model to be compared. A fitted
  model supplies lavaan's count, `attr(logLik(fit), "df")`, which leaves
  out the saturated means of a model without a mean structure and, under
  `fixed.x = TRUE`, the moments of the exogenous observed variables; a
  difference between two models with the same settings is unaffected.

- N:

  For summary input, the sample size *N* in the penalties of the BIC
  family, greater than 1. A fitted model supplies its total over groups.

- measures:

  The criteria to compute: any of `"AIC"`, `"BIC"`, `"SABIC"`, `"HBIC"`,
  `"IBIC"`, and `"SPBIC"`, in either case. The default computes the four
  that need only \\\ell\\, *d*, and *N*. Every row of the table is
  reported whatever is requested, with `NA` for a criterion not
  computed.

- parameter_estimate:

  For SPBIC from summary input, the vector \\\hat{\theta}\\ of the *d*
  free parameter estimates, in the order of the rows of
  `information_observed`.

- prior_mean:

  For SPBIC, the vector \\\theta^\*\\ of prior means, in the order of
  `parameter_estimate`. Supply this or `zero_parameters`. With a fitted
  model and neither, the rule of Bollen et al. (2014, p. 6) applies:
  free loadings and regression coefficients are centered at zero and
  every other free parameter at its estimate (see Details).

- zero_parameters:

  For SPBIC, the positions or unique names of the entries of
  `parameter_estimate` to center at zero, every other entry being
  centered at its estimate. A saturated model needs an explicit choice;
  Bollen et al. (2014, p. 6) center its covariances at zero and its
  variances (and means) at their estimates.

- information_observed:

  For SPBIC from summary input, the observed information matrix
  \\I_O(\hat{\theta})\\: the negative Hessian of the log likelihood at
  the estimates, in the free parameters, in the order of
  `parameter_estimate`.

- information_expected:

  For IBIC from summary input, the expected (Fisher) information matrix
  \\I_E(\hat{\theta})\\ at the estimates, in the free parameters. The
  inverse of a robust (sandwich) covariance matrix of the estimates is
  not an information matrix and is no substitute.

- information_scale:

  For summary input, whether the supplied information matrices are for
  the whole sample (`"total"`, the default; the inverse of the
  covariance matrix of the estimates under ML) or per observation
  (`"per_observation"`, what
  [`lavaan::lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html)
  reports), in which case they are multiplied by `N`. A fitted model
  supplies per observation matrices, which are multiplied by *N*, and
  the argument is refused with one.

- spbic_dimension:

  Optional: the dimension used in the penalty of SPBIC, a positive whole
  number no greater than `npar`. By default it is `npar`, as in the
  general expression. Bollen et al.'s (2014) Table 6 used smaller counts
  for SPBIC, printed in parentheses (15 for the saturated model rather
  than 20), and states no rule for them; supply such a count only to
  reproduce that table, never to tune a ranking. It changes neither
  `parameter_estimate` nor the information matrices.

## Value

A `dmar_tbl` (class `dmar_sem_information`) with columns `term` and
`value` and, always in this order, the rows `AIC`, `BIC`, `SABIC`,
`HBIC`, `IBIC`, and `SPBIC` (`NA` for a criterion not requested),
`log_likelihood`, `npar`, and `N`. The `sem_information` attribute
records the source; for a fitted model the likelihood, the missing data
method, the exogenous observed variables it conditions on, the mean
structure, the number of groups, and, when IBIC or SPBIC is computed,
the identification of each latent variable; the scale and provenance of
the information matrices; the log determinant of the whole sample
expected information for IBIC; for SPBIC the estimates, the prior means
and the rule that set them, *Q*, the penalty dimension and its rule, the
penalty, and whether the boundary branch applied; for a two-stage fit
whose first stage read cases observed on an auxiliary variable alone,
`two_stage_rows`, the cases read and lavaan's count of them in each
group, how many such cases there are, whether the fit was computed again
with all of them (`corrected`), and whether its estimates depend on the
counts, so that they were computed again too (`reestimated`); and the
lines of the printed footer.

## Details

### The Criteria

Write \\\ell\\ for the maximized log likelihood, *d* for the number of
free parameters, and *N* for the sample size. Each criterion is
\\-2\ell\\ plus a penalty: \\AIC = -2\ell + 2d\\ (Akaike, 1974), \\BIC =
-2\ell + d\log N\\ (Schwarz, 1978), \\SABIC = -2\ell + d\log\\(N +
2)/24\\\\ (Sclove, 1987; lavaan's `bic2`), and \\HBIC = -2\ell +
d\log\\N/(2\pi)\\\\ (Haughton, 1988).

Bollen et al. (2012) derive the BIC family from the Laplace
approximation to the log of a model's marginal likelihood, \$\$\log
P(Y\|M) = \ell + \log P(\hat{\theta}\|M) + \frac{d}{2}\log(2\pi) -
\frac{d}{2}\log N - \frac{1}{2}\log\|\bar{I}\_E(\hat{\theta})\| +
O(N^{-1/2})\$\$ (Bollen et al., 2014, Eq. 9), where
\\P(\hat{\theta}\|M)\\ is the prior density at the estimates and
\\\bar{I}\_E = I_E/N\\ the expected information per observation.
Multiplied by \\-2\\, BIC keeps the first and fourth terms, HBIC the
third as well, and IBIC every term but the prior: \\IBIC = -2\ell +
d\log\\N/(2\pi)\\ + \log\|\bar{I}\_E(\hat{\theta})\|\\ (Bollen et al.,
2012, Eq. 20; 2014, Eq. 16), which equals \\-2\ell - d\log(2\pi) +
\log\|I_E(\hat{\theta})\|\\ for the whole sample information.
Recomputing the worked example of Bollen et al. (2012, Table 1) from its
printed inputs shows that the log determinant it adds to
\\d\log\\N/(2\pi)\\\\ is that of the whole sample information (its
elements are about \\N\\ times the per observation ones the log
likelihood implies), which counts \\\log N\\ twice relative to their Eq.
20; this function follows Eq. 20, which reproduces the IBIC row of
Bollen et al. (2014, Table 6).

SPBIC replaces BIC's unit information prior with a normal prior centered
at \\\theta^\*\\ (`prior_mean`) whose variance, a multiple of the
inverse unit information, is chosen to maximize the marginal likelihood,
an empirical Bayes prior. With \\Q = (\hat{\theta} - \theta^\*)'
I_O(\hat{\theta}) (\hat{\theta} - \theta^\*)\\ for the estimates
\\\hat{\theta}\\ (`parameter_estimate`) and the whole sample observed
information \\I_O\\, \\SPBIC = -2\ell + d\\1 + \log(Q/d)\\\\ when \\Q \>
d\\ (Bollen et al., 2012, Eq. 14; 2014, Eq. 18). When \\Q \le d\\ the
maximizing prior variance is zero, the prior is a point mass at
\\\theta^\*\\, and \\SPBIC = -2\ell + Q\\ (Bollen et al., 2012, Eq. 15);
the two branches meet at \\Q = d\\. Bollen et al. (2014, note 7) never
met that branch and caution that its behavior is little known, and the
function warns when it applies.

SABIC replaces *N* in BIC's penalty with \\(N + 2)/24\\, the penalty of
Rissanen's (1978) shortest data description criterion as Sclove (1987,
criteria 4 and 5) gives it, which comes from the length of a coded
description of the data rather than from the Laplace approximation; Lin,
Huang, and Weng (2017) write it ABIC, and Mplus and lavaan's `bic2` call
it the sample size adjusted BIC. Its penalty per parameter is smaller
than BIC's, and it is not positive when \\N \le 22\\, where the function
warns. Lin et al. (2017) found that it improved on BIC in selecting path
models mainly with small samples or small misspecifications, and that
HBIC and SPBIC did best overall. AIC estimates the relative expected
Kullback-Leibler discrepancy of a model rather than its marginal
likelihood (Burnham & Anderson, 2004). AICc is not offered: its
correction is derived for the univariate linear model and is not
justified for structural equation models in general.

### Units and Identification

AIC, BIC, SABIC, and HBIC depend on a model only through \\\ell\\ and
*d*, which do not change when an observed variable is rescaled or a
latent variable is identified another way. IBIC and SPBIC also depend on
the information matrix in the free parameters as they are parameterized.

IBIC keeps \\\log\|\bar{I}\_E\|\\ but drops the prior density of the
Laplace approximation, and that is the term whose Jacobian would cancel
a change of parameterization. When the units of a free parameter are
multiplied by \\c^k\\, \\\log\|\bar{I}\_E\|\\ changes by \\-2k\log c\\.
Rescaling an observed variable by *c* therefore changes the IBIC of one
model relative to another by \\-2k\log c\\ for each free parameter the
first has and the second lacks whose units carry that variable's (\\k =
1\\ for a loading or for a covariance with an unscaled variable, \\k =
2\\ for its variance), and the ranking of two models fitted to the same
data can reverse when one variable's units change. Identifying the
latent variables by fixing their variances (`std.lv = TRUE`) rather than
by a marker indicator changes IBIC differences too, since the map
between the two parameterizations is nonlinear and its Jacobian differs
between models. An IBIC comparison is defined only for a stated choice
of units and identification. Standardizing every observed variable
before fitting every model, or keeping their original metrics, fixes the
units, but it is a convention rather than an invariant answer: report
it, and read IBIC beside BIC, HBIC, and SPBIC.

SPBIC's *Q* is invariant to a linear reparameterization when the prior
mean is transformed with the estimates, so rescaling an observed
variable does not change it under the default prior. That prior,
however, centers the free loadings and regression coefficients at zero
and the other free parameters at their estimates, and which parameters
are free loadings depends on the identification: under `std.lv = TRUE`
every loading is centered at zero, while under a marker indicator, the
convention of the example of Bollen et al. (2014), the scale sits in the
factor variance, which stays at its estimate. A mean structure adds the
saturated means to *d*, and their information to \\\bar{I}\_E\\. Models
compared by IBIC or SPBIC must therefore share their identification and
mean structure;
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
requires every latent variable of every model to be identified the same
way (each by a marker indicator, or each by a fixed variance) and the
same mean structure, and refuses a set that does not, and for a fitted
model the `sem_information` attribute records both.

### What a Fitted Model Supplies

\\\ell\\ is lavaan's `logLik(fit)`, *d* its `df` attribute, and *N* the
total sample size. In a multiple group model \\\ell\\ is the sum over
the groups and *N* the total, as in lavaan. Under `fixed.x = TRUE` (the
default of [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html))
or `conditional.x = TRUE`, \\\ell\\ is the log likelihood conditional on
the exogenous observed variables, and criteria compare only models with
the same exogenous variables. Under FIML (`missing = "ml"`) \\\ell\\ is
the observed data log likelihood. Under `missing = "two.stage"` or
`"robust.two.stage"` it is the complete data normal log likelihood at
the stage one EM estimates of the means and covariances, so that AIC and
BIC are the two-stage criteria of Lai (2021, Eqs. 32 and 35), up to a
constant shared by all models fitted to the same data; Lai shows that,
unlike the FIML criteria, they select in large samples the model
complete data would select when the values are missing at random, and
that in samples of the usual size they miss it more often than complete
data would. They are comparable only with other two-stage fits of the
same variables with the same auxiliary variables in the first stage, and
they are the `AIC` and `BIC` rows
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports after two-stage estimation. For IBIC and SPBIC, the expected and
observed information matrices are lavaan's per observation
`lavInspect(fit, "information.expected")` and
`lavInspect(fit, "information.observed")`, multiplied by *N*; they are
the matrices the criteria are defined on under lavaan's defaults
`h1.information = "structured"` (the expected information at the fitted
moments) and `observed.information = "hessian"` (the negative Hessian),
and other settings are refused.

A two-stage fit can have cases observed on an auxiliary variable
(lavaan's `aux`) but on none of the model's variables. lavaan's first
stage reads them, so they move the EM estimates, but its sample size
(`lavaan::lavInspect(fit, "nobs")`) counts them as empty and leaves them
out, and it forms the log likelihood as that sample size times the
complete data log likelihood per case at the EM estimates. Its log
likelihood, and the *N* of the penalties of BIC, SABIC, and HBIC, then
rest on fewer cases than the moments. Both are computed here from the
fit computed again with every case the first stage read: in one group,
and in several groups that no constraint holds together, the estimates
do not depend on the count, and the log likelihood of each group is its
value per case times the cases its first stage read, with *N* their
total. A line under the table says so, and the `sem_information`
attribute records the counts. On six Holzinger and Swineford tests with
the whole block missing for 77 of 301 pupils, at random given general
information (the auxiliary), lavaan's *N* for the two-factor model is
224 and its log likelihood -4040.25, which is 224/301 of the -5429.09
reported here with *N* = 301 (BIC 10966.61 here, 8183.32 from lavaan's
values). In several groups held together by a constraint, such as equal
loadings or intercepts, the second stage weights each group by its
count, so when the groups' counts grow in different ratios the estimates
move as well, and the log likelihood is the one at the estimates with
every case read (those of semTools' `twostage()`), which the line says.
On the same tests by school with equal loadings and intercepts (130 and
94 cases with an observed test, 156 and 145 read), it is -5393.34 with
BIC 10957.89, against -5393.52 and 10958.25 at lavaan's estimates and
-4020.79 and 8203.92 from lavaan's values.
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)`(missing = "two_stage")`
and
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
return fits that count the cases already, and give the same criteria.
FIML (`missing = "ml"`) keeps the cases in its likelihood and in *N* as
well.

### What Is Refused

A fitted model is refused when it did not converge or has an improper
solution, when it has ordered categorical indicators or was fitted by an
estimator other than ML (its fit function is then not a log likelihood),
when it uses the Wishart likelihood, and when it was fitted with
sampling weights: lavaan then maximizes a weighted pseudo log likelihood
(Asparouhov, 2005), for which neither the penalties nor the information
matrices have their likelihood justification. For a two-level model the
sample size in the BIC family's penalty is ambiguous (Bollen et al.,
2012, p. 295), so only AIC is computed from the fit, and the others need
summary input with the *N* the user chooses. Automatic IBIC and SPBIC
further need a single group, single level ML fit (not MLM or MLR) of
complete data whose free parameters are unconstrained, with no
exploratory (rotated) blocks (their information matrices are singular)
and with the exogenous observed variables modeled (`fixed.x = FALSE`).
Each refusal says why and what to do; for another parameterization,
supply the quantities as summary input.

### Relative to the Saturated Model

Bollen et al. (2014) report each criterion as the difference between the
hypothesized model and the saturated model, so that BIC becomes \\T -
df\log N\\ for the likelihood ratio statistic \\T\\ on \\df\\ degrees of
freedom. Subtracting the saturated model's criterion from the
hypothesized model's gives that scale. A single model's criterion cannot
be reconstructed from \\T\\ and \\df\\ alone.

## References

Akaike, H. (1974). A new look at the statistical model identification.
*IEEE Transactions on Automatic Control, 19*(6), 716–723.
[doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705)

Asparouhov, T. (2005). Sampling weights in latent variable modeling.
*Structural Equation Modeling, 12*(3), 411–434.

Bollen, K. A., Harden, J. J., Ray, S., & Zavisca, J. (2014). BIC and
alternative Bayesian information criteria in the selection of structural
equation models. *Structural Equation Modeling, 21*(1), 1–19.
[doi:10.1080/10705511.2014.856691](https://doi.org/10.1080/10705511.2014.856691)

Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A
comparison of Bayes factor approximation methods including two new
methods. *Sociological Methods & Research, 41*(2), 294–324.
[doi:10.1177/0049124112452393](https://doi.org/10.1177/0049124112452393)

Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
Understanding AIC and BIC in model selection. *Sociological Methods &
Research, 33*(2), 261–304.
[doi:10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644)

Haughton, D. M. A. (1988). On the choice of a model to fit data from an
exponential family. *The Annals of Statistics, 16*(1), 342–355.
[doi:10.1214/aos/1176350709](https://doi.org/10.1214/aos/1176350709)

Lai, K. (2021). Using information criteria under missing data: Full
information maximum likelihood versus two-stage estimation. *Structural
Equation Modeling, 28*(2), 278–291.
[doi:10.1080/10705511.2020.1780925](https://doi.org/10.1080/10705511.2020.1780925)

Lin, L.-C., Huang, P.-H., & Weng, L.-J. (2017). Selecting path models in
SEM: A comparison of model selection criteria. *Structural Equation
Modeling, 24*(6), 855–869.
[doi:10.1080/10705511.2017.1363652](https://doi.org/10.1080/10705511.2017.1363652)

Rissanen, J. (1978). Modeling by shortest data description. *Automatica,
14*(5), 465–471.
[doi:10.1016/0005-1098(78)90005-5](https://doi.org/10.1016/0005-1098%2878%2990005-5)

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics, 6*(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

Sclove, S. L. (1987). Application of model-selection criteria to some
problems in multivariate analysis. *Psychometrika, 52*(3), 333–343.
[doi:10.1007/BF02294360](https://doi.org/10.1007/BF02294360)

## See also

[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
to compare fitted models by one of these criteria,
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
for the fit measures of one model, and the vignette
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
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The three-factor model of nine of the holzinger_swineford
# tests, identified by a marker indicator (lavaan's default).
fit <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
   speed   =~ t10_addition + t12_counting_groups_of_dots +
              t13_straight_and_curved_capitals",
  data = holzinger_swineford)
sem_information_criteria(fit)
#>  term           value    
#>  AIC            19198.131
#>  BIC            19275.980
#>  SABIC          19209.380
#>  HBIC           19237.384
#>  IBIC           <NA>     
#>  SPBIC          <NA>     
#>  log_likelihood -9578.065
#>  npar           21       
#>  N              301      
#> 
#> Smaller is preferred; only a difference between models fitted to the same data under the same likelihood is meaningful.
#> Not requested: IBIC, SPBIC (see 'measures').
#> Log likelihood: the normal log likelihood of the complete data.
sem_information_criteria(fit, measures = c("BIC", "IBIC", "SPBIC"))
#>  term           value    
#>  AIC            <NA>     
#>  BIC            19275.980
#>  SABIC          <NA>     
#>  HBIC           <NA>     
#>  IBIC           19101.759
#>  SPBIC          19259.359
#>  log_likelihood -9578.065
#>  npar           21       
#>  N              301      
#> 
#> Smaller is preferred; only a difference between models fitted to the same data under the same likelihood is meaningful.
#> Not requested: AIC, SABIC, HBIC (see 'measures').
#> Log likelihood: the normal log likelihood of the complete data.
#> IBIC and SPBIC depend on how the latent variables are identified (here marker indicator) and on the mean structure, IBIC also on the units; compare only models that share them.
#> SPBIC prior: Bollen et al. (2014, p. 6): free loadings and regression coefficients centered at zero, other free parameters at their estimates; Q = 1053.78, dimension 21.

# Summary input: the autoregressive model of Bollen et al. (2014,
# Table 6), whose BIC relative to the saturated model is T - df log N.
ar <- sem_information_criteria(log_likelihood = -534.27 / 2, npar = 14,
                               N = 3912)
saturated <- sem_information_criteria(log_likelihood = 0, npar = 20,
                                      N = 3912)
ar$value[ar$term == "BIC"] - saturated$value[saturated$term == "BIC"]
#> [1] 484.6392

# SPBIC from summary input, with an explicit prior mean.
sem_information_criteria(log_likelihood = -120, npar = 2, N = 200,
                         measures = "SPBIC",
                         parameter_estimate = c(0.4, 1),
                         prior_mean = c(0, 1),
                         information_observed = diag(c(100, 50)))
#>  term           value   
#>  AIC            <NA>    
#>  BIC            <NA>    
#>  SABIC          <NA>    
#>  HBIC           <NA>    
#>  IBIC           <NA>    
#>  SPBIC          246.159 
#>  log_likelihood -120.000
#>  npar           2       
#>  N              200     
#> 
#> Smaller is preferred; only a difference between models fitted to the same data under the same likelihood is meaningful.
#> Not requested: AIC, BIC, SABIC, HBIC, IBIC (see 'measures').
#> SPBIC depends on how the latent variables are identified and on the mean structure; compare only models that share them.
#> SPBIC prior: explicit prior mean; Q = 16, dimension 2.
```
