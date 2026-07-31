# Design, Measurement, and Analysis in R

A modern R package for design, measurement, and analysis in
human-centered research, with special strength in effect sizes,
confidence intervals, sample size planning, reliability and agreement,
mediation analysis, equivalence testing, meta-analysis, experimental and
quasi-experimental designs, repeated measures, and model
comparison-based inference. DMAR (pronounced “Dee-Mar”) is heavily
methodological in nature, drawing on the psychometric and statistical
traditions, and is aligned with the methodological and applied research
program of the author. Much of the package traces to the author's
methodological work and collaborations, including sample size planning
via accuracy in parameter estimation (AIPE; Kelley & Maxwell, 2003;
Kelley & Rausch, 2006; Maxwell, Kelley, & Rausch, 2008), the definition
and communication of effect sizes (Kelley & Preacher, 2012; Preacher &
Kelley, 2011), and the model comparison perspective of Maxwell, Delaney,
and Kelley (2027). It aims to be methodologically sound and particularly
well-suited to research in which the independent or dependent variables
involve the person, across psychology, sociology, education, management,
marketing, and information systems.

## Details

The package makes accessible to researchers a variety of methods that
are easy to use, including from sample estimates or results reported in
published articles: effect size estimation, confidence intervals for
effect sizes, sample size planning, multivariate methods, factor
analysis, and certain latent variable models. Particular strengths
include sample size planning under several complementary frameworks:
accuracy in parameter estimation (AIPE), power analysis, minimum risk,
and equivalence. Most exported functions return a tidy `data.frame` with
a `term` column and a numeric `value` column (some carry additional
typed columns per term), and ggplot2 is used for the plotting functions.
A few functions return the shape their task calls for instead: model
fits such as `mlmr` return a richer list-like object with `coef` /
`vcov` / `confint` methods, `descriptives` returns a list of summary
tables, and a small number of scalar utilities such as `skewness` and
`kurtosis` return a bare numeric. The interface is consistent, modern,
and opinionated, and is designed for clarity and reproducibility.

DMAR builds heavily on the MBESS package (Kelley, 2007a, 2007b), which
has been on CRAN for about two decades. MBESS was originally framed for
the behavioral, educational, and social sciences, but its use has grown
well beyond that scope; DMAR is a more modern, more general, and greatly
expanded reimagining that reflects how the methods are now applied. The
notational and programming changes were too extensive to ship under the
MBESS name without breaking compatibility, so the work ships under a new
name; MBESS itself remains stable.

**Function families.** A user-facing tour:

- Effect sizes:

  `smd`, `smd_c`, `smd_trimmed`, `eta_squared`, `eta_squared_partial`,
  `eta_squared_generalized`, `omega_squared`, `omega_squared_partial`,
  `cohen_f`, `cles`, `cliff_delta`, `vargha_delaney_A`,
  `proportion_of_superiority`, `probability_of_superiority_paired`,
  `lin_ccc`.

- Confidence intervals on effect sizes:

  `ci_smd`, `ci_smd_c`, `ci_R2`, `ci_r`, `ci_rc`, `ci_src`,
  `ci_eta_squared` (and partial / generalized variants),
  `ci_omega_squared`, `ci_pvaf`, `ci_snr`, `ci_srsnr`, `ci_mahalanobis`,
  `ci_eigenvalue`, `ci_cv`, `ci_sm`, `ci_reg_coef`, `ci_cc`, `ci_rmsea`.

- Maximum likelihood regression:

  `mlmr` (univariate full information maximum likelihood (FIML),
  lm-like), `mlmr_mv` (multivariate FIML).

- ANOVA and ANCOVA:

  `ancova`, `anova_within_two_way`, `mixed_anova`, `manova_split_plot`,
  `simple_effects_AB`, `contrast_test`, `pairwise_within`,
  `mauchly_test`, `obrien_test`.

- Reliability and agreement:

  `reliability`, `reliability_alpha`, `reliability_omega` (with a model
  implied or observed total-variance denominator, and a
  `reliability_omega_categorical` for ordered items),
  `reliability_kr20`, `reliability_H`, `cohen_kappa`, `fleiss_kappa`,
  `krippendorff_alpha`, `gwet_ac`, `loa`.

- Mediation:

  `mediate` (the simple mediation model with bootstrap, Monte Carlo, and
  Sobel intervals), `mediation_mbco` (likelihood ratio tests of
  arbitrary mediation effects by model-based constrained optimization,
  with multiple groups and moderated mediation probing), and
  `plot_mediation_mbco` (conditional effect curves with confidence
  bands).

- Confirmatory factor and SEM tools:

  `cfa_1`, `cov_sem`, `covmat_from_cfm`, `compare_cov_structures`.

- Sample size planning (AIPE):

  `ss_aipe_smd`, `ss_aipe_R2`, `ss_aipe_reg_coef`, `ss_aipe_partial_r`,
  `ss_aipe_omega_squared`, `ss_aipe_icc`, `ss_aipe_cv`, `ss_aipe_pcm`,
  `ss_aipe_rmsea`, the cluster-randomized planners `ss_aipe_crd_*`, plus
  their Monte Carlo sensitivity companions `ss_aipe_*_sensitivity`.

- Sample size planning (power):

  `ss_power_smd`, `ss_power_R2`, `ss_power_r`, `ss_power_reg_coef`,
  `ss_power_sem`, `ss_power_c`, `ss_power_c_ancova`,
  `ss_power_contrast`, `ss_power_factorial_anova`,
  `ss_power_split_plot_anova`, `ss_power_mixed_effects`,
  `ss_power_one_way_anova`, `ss_power_pcm`, `ss_power_rm_anova`,
  `ss_power_sc`.

- Critical values and tests:

  `cv_t`, `cv_z`, `cv_smm`, `cv_scheffe`, `cv_tukey_hsd`, `cv_dunnett`,
  `dunnett_ci`, `tukey_kramer_ci`, `scheffe_ci`, `welch_t`,
  `summary_t_test`, `correlations_test`, `power_fisher_exact`,
  `randomization_test_paired`, `tost_smd`, `tost_r`,
  `power_equivalence_md`.

- Parameterization conversions:

  `convert_R2_f` / `convert_f_R2`, `convert_R2_lambda` /
  `convert_lambda_R2`, `convert_delta_lambda` / `convert_lambda_delta`,
  `convert_r_z` / `convert_z_r`, `convert_cor_cov`.

- Visualization:

  `plot_smd`, `plot_ci`, `plot_R2`, `plot_trajectories`,
  `plot_trajectories_fitted`.

- Multilevel and clustering:

  `icc`, `icc_lmer`, `variance_components_mls`, `deft` (Kish design
  effect), `ss_aipe_crd_*`.

- Data sets:

  `bessel_errors` (Bessel's 1818 grouped distribution of Bradley's
  astronomical observation errors), `diagnosis_agreement` (Cohen's 1968
  weighted kappa illustration), `drinks_trial` (Smith, Meyers, and
  Delaney's 1998 Community Reinforcement Approach drinking trial),
  `holzinger_swineford` (the 1939 factor analysis study, with `HS_Data`
  alias), `prime_time_achievement` (the Indiana Prime Time third grade
  achievement evaluation, with `Prime_Time` alias), `pygmalion`
  (Rosenthal and Jacobson's 1968 teacher-expectancy data),
  `teacher_expectancy` (Raudenbush's 1984 meta-analysis of 18
  teacher-expectancy experiments), and `test_market` (Bryant and
  Bruvold's 1980 controlled test-market experiment for ANCOVA with a
  random covariate).

**Feedback.** Bug reports, feature requests, and suggestions for new
methods are welcomed by email to Ken Kelley <kkelley@nd.edu> (please put
“DMAR” in the subject line). See <https://kenkelley.org> for Ken
Kelley's web site, <https://kenkelley.org/publications/> for related
publications, and <https://github.com/yelleKneK/DMAR> for the project's
GitHub page.

## References

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*(2), 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

## See also

Useful links:

- <https://kenkelley.org>

- <https://yelleknek.github.io/DMAR/>

- <https://github.com/yelleKneK/DMAR>

- Report bugs at <https://github.com/yelleKneK/DMAR/issues>

## Author

Ken Kelley <kkelley@nd.edu>
