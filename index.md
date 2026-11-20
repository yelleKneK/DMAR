# DMAR

**DMAR** (Design, Measurement, and Analysis in R, pronounced “Dee-Mar”)
is a package for human-centered research: effect sizes, confidence
intervals for effect sizes, sample size planning, reliability and
agreement, mediation analysis (from the simple mediation model to
likelihood ratio tests of arbitrary indirect effects, with moderated
mediation probing), equivalence testing, meta-analysis, experimental and
quasi-experimental designs, repeated measures, and model
comparison-based inference. It draws on the psychometric and statistical
traditions and reflects the methodological and applied research program
of its author, Ken Kelley (University of Notre Dame). It aims to be
methodologically sound and useful across psychology, sociology,
education, management, marketing, and information systems, especially
where the independent or dependent variables involve the person.

## Installation

``` r

# The released version, from CRAN:
install.packages("DMAR")

# The development version, from GitHub:
# install.packages("remotes")
remotes::install_github("yelleKneK/DMAR")
```

Each CRAN release is tagged on GitHub (`v1.0.0` for the first); between
releases the GitHub copy carries the changes of the next one, listed in
NEWS.

R \>= 4.0.0 is required. Most estimation, inference, and planning
functions return a tidy `data.frame(term, value)` that composes
naturally with `dplyr` / `ggplot2` pipelines and with
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html).
Heavier dependencies (e.g., `lavaan`, `lme4`, `OpenMx`) are in
`Suggests` and gated by
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html); install
them only when you need the functions that use them.

## Function Families

DMAR is organized into these function families:

- **Effect size estimates.**
  [`smd()`](https://yelleknek.github.io/DMAR/reference/smd.md),
  [`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
  [`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)
  (sometimes called Cohen’s U3),
  [`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
  [`vargha_delaney_A()`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
  [`cohen_f()`](https://yelleknek.github.io/DMAR/reference/cohen_f.md),
  [`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
  [`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
  [`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md).
- **Confidence intervals for effect sizes.**
  [`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
  [`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
  [`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  [`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
  [`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
  [`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
  [`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
  [`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
  [`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
  [`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
  [`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
  [`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
  [`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
  [`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md)
  and
  [`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md)
  (Fieller’s set),
  [`ci_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md)
  (two proportions on the same units),
  [`ci_crossed()`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
  (units crossed with raters), and
  [`ci_two_phase()`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
  (a costly measurement with an inexpensive proxy).
- **Sample size planning under AIPE** (accuracy in parameter
  estimation). The `ss_aipe_*()` family covers SMD, R², ω², CV, partial
  / semipartial r, regression coefficients (standardized and
  unstandardized), ANCOVA contrasts, RMSEA, SEM paths, polynomial
  change, mixed-effects fixed effects, ICC, reliability, Cliff’s delta,
  a single mean
  ([`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)),
  crossed designs of units and raters
  ([`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)),
  and two-phase sampling
  ([`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
  [`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)).
  Most planners have a Monte Carlo `*_sensitivity()` sibling that
  quantifies the impact of misspecified planning values;
  [`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)
  also evaluates the designs of
  [`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md).
- **Sample size planning under power.**
  [`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
  [`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
  [`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
  [`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
  [`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
  [`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
  [`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md),
  [`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
  [`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
  [`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
  [`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md).
- **Minimum risk / sequential estimation.**
  [`mr_smd()`](https://yelleknek.github.io/DMAR/reference/mr_smd.md),
  [`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.md),
  [`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)
  (with
  [`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)),
  and
  [`ss_seq_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_seq_two_phase.md)
  (a mean under two-phase sampling).
- **Equivalence testing.**
  [`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
  [`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
  (with nonsuperiority),
  [`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
  [`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
  [`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md).
- **Reliability.**
  [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  [`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
  [`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
  [`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
  and
  [`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
  (omega hierarchical and its relatives),
  [`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
  the umbrella
  [`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
  and
  [`ss_aipe_reliability()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md).
- **Factor analysis and structural equation models.**
  [`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  (with [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
  and [`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md)),
  [`fa_exploratory()`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
  (including target rotation),
  [`factor_count()`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
  [`parallel_analysis()`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md),
  [`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
  [`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
  [`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
  [`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
  [`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  [`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
  [`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
  [`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
  [`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
  [`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
  [`sem_sam()`](https://yelleknek.github.io/DMAR/reference/sem_sam.md),
  [`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md), and
  the common method variance procedures
  ([`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
  [`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
  [`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md)).
- **Missing data.** Full information maximum likelihood (the default for
  continuous items in the factor analysis functions and in
  [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  and
  [`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md))
  and two-stage maximum likelihood (in
  [`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  [`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
  [`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  and
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)),
  with
  [`fraction_missing_information()`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
  [`missing_patterns()`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md),
  [`missing_coverage()`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
  [`cov_fiml()`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
  [`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
  and
  [`df_barnard_rubin()`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md)
  for what the missing values cost.
- **ANOVA / ANCOVA / contrasts.**
  [`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md),
  [`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
  [`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
  [`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
  [`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
  [`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
  [`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
  [`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
  [`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
  [`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
  [`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
  [`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
  [`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
  [`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md).
- **Multiple-comparison procedures.**
  [`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
  [`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
  [`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
  [`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
  [`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
  and their critical-value helpers
  [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
  [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
  [`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
  [`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
  [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
  [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md).
- **Agreement and reliability of raters.**
  [`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
  [`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
  [`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
  [`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
  [`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
  [`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
  [`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
  [`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md).
- **Multilevel and multivariate methods.**
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md),
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md),
  [`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
  [`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
  [`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
  [`covmat_from_cfa()`](https://yelleknek.github.io/DMAR/reference/covmat_from_cfa.md),
  [`cov_sem()`](https://yelleknek.github.io/DMAR/reference/cov_sem.md),
  [`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md).
- **Longitudinal designs.**
  [`ss_aipe_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md),
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
  [`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md),
  [`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md).
- **Design utilities.**
  [`design_effect()`](https://yelleknek.github.io/DMAR/reference/design_effect.md)
  (Kish’s design effect and DEFT),
  [`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
  [`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
  [`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
  [`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
  [`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
  [`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md).
- **Parameterization conversions** (`convert_*`). Invertible maps
  between equivalent metrics (`convert_r_Z` / `convert_Z_r`,
  `convert_R2_f` / `convert_f_R2`, `convert_lambda_R2` /
  `convert_R2_lambda`, `convert_delta_lambda` / `convert_lambda_delta`,
  `convert_d_r` / `convert_r_d`, `convert_d_or` / `convert_or_d`,
  `convert_F_chisq` / `convert_chisq_F`, `convert_cor_cov`).
- **Plots.**
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  and
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)
  show an effect size with its CI and, given the sample sizes, the
  sample size;
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
  and
  [`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)
  show observed and fitted individual trajectories.
- **Variance utilities.**
  [`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
  [`var_R2()`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
  [`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
  [`var_partial_r()`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
  [`var_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
  [`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
  [`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
  [`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
  [`var_icc()`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
  [`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md).

## Canonical Examples

### 1. Standardized Mean Difference With a Noncentral *t* Confidence Interval

``` r

library(DMAR)

# From summary statistics:
smd(mean_1 = 60, mean_2 = 55, s_1 = 10, s_2 = 10, n_1 = 30, n_2 = 30)
```

| term | value |
|:-----|:------|
| smd  | 0.500 |

``` r


# CI on the standardized mean difference (Cohen's d) via the noncentral t:
ci_smd(smd = 0.5, n_1 = 30, n_2 = 30, conf_level = 0.95)
```

| term     | value  |
|:---------|:-------|
| ci_lower | -0.016 |
| smd      | 0.500  |
| ci_upper | 1.012  |

Confidence level: 95%

### 2. R² With a Noncentral *F* Confidence Interval and AIPE Sample Size

``` r

ci_R2(R2 = 0.30, N = 100, p = 4, conf_level = 0.95)
```

| term     | value | prob_less | prob_greater |
|:---------|:------|:----------|:-------------|
| ci_lower | 0.130 | 0.025     | 0.975        |
| R2       | 0.300 | NA        | NA           |
| ci_upper | 0.430 | 0.975     | 0.025        |

Confidence level: 95%

``` r


# What N is needed so the 95% CI on R^2 has full width <= 0.20?
ss_aipe_R2(true_R2 = 0.30, width = 0.20, p = 4)
```

| term        | value |
|:------------|:------|
| necessary_N | 228   |

Confidence level: 95%

### 3. AIPE Planning With a Sensitivity Check

``` r

# Plan N for a half-width of .10 on the SMD at a population SMD of 0.40.
ss_aipe_smd(true_smd = 0.40, width = 0.20)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 784   |
| true_smd              | 0.400 |
| width                 | 0.200 |

Confidence level: 95%

``` r


# How does that plan hold up if the population SMD is actually 0.25?
set.seed(113)
ss_aipe_smd_sensitivity(
  true_smd      = 0.25,
  supposed_smd  = 0.40,
  width         = 0.20,
  G               = 200,
  print_iter      = FALSE
)
```

| term                | value |
|:--------------------|:------|
| mean_smd            | 0.253 |
| median_smd          | 0.249 |
| sd_smd              | 0.052 |
| mean_ci_width       | 0.199 |
| median_ci_width     | 0.199 |
| sd_ci_width         | 0.000 |
| mean_ci_width_lower | 0.099 |
| mean_ci_width_upper | 0.099 |
| pct_ci_less_width   | 1.000 |
| pct_ci_miss_low     | 0.030 |
| pct_ci_miss_high    | 0.025 |
| total_type_I_error  | 0.055 |
| n_per_group         | 784   |
| total_N             | 1568  |
| true_smd            | 0.250 |
| supposed_smd        | 0.400 |
| width               | 0.200 |
| conf_level          | 0.950 |

Replications: 200

Confidence level: 95%

### 4. Welch’s *t* Test With a Tidy Return

``` r

set.seed(113)
x <- rnorm(20, mean = 100, sd = 15)
y <- rnorm(20, mean = 110, sd = 25)
welch_t(x, y, conf_level = 0.95)
```

| term            | value   |
|:----------------|:--------|
| mean_difference | -16.871 |
| t_value         | -2.410  |
| df              | 28.349  |
| p_value         | 0.0227  |
| ci_lower        | -31.204 |
| ci_upper        | -2.537  |
| mean_x          | 100.297 |
| mean_y          | 117.168 |
| sd_x            | 14.289  |
| sd_y            | 27.860  |
| n_x             | 20      |
| n_y             | 20      |

Confidence level: 95%

### 5. Reliability With a Confidence Interval

``` r

# Cronbach's alpha plus the Bonett (2002) CI from a covariance matrix.
S <- matrix(
  c(1.0, 0.6, 0.5, 0.5,
    0.6, 1.0, 0.6, 0.5,
    0.5, 0.6, 1.0, 0.6,
    0.5, 0.5, 0.6, 1.0),
  nrow = 4)
reliability_alpha(S = S, N = 250, ci_method = "bonett")
```

| term           | value |
|:---------------|:------|
| estimate       | 0.830 |
| se             | 0.018 |
| se_transformed | 0.104 |
| ci_lower       | 0.792 |
| ci_upper       | 0.861 |
| conf_level     | 0.950 |
| N              | 250   |
| N_complete     | 250   |
| n_items        | 4     |

### 6. Kish’s Design Effect in a Clustered Sample

``` r

# 25 schools with attendance ranging from 0 (closed) to 30 (full class),
# ICC = .10. How much clustering information do we actually have?
sizes <- c(0, 0, 1, 1, 1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 22,
           22, 25, 26, 28, 28, 30, 30, 30, 30, 30)
design_effect(cluster_sizes = sizes, icc = 0.10)
```

| term                   | value   |
|:-----------------------|:--------|
| design_effect          | 3.327   |
| deft                   | 1.824   |
| effective_n            | 119.328 |
| N                      | 397     |
| n_clusters_total       | 25      |
| n_clusters_with_data   | 23      |
| n_clusters_empty       | 2       |
| n_clusters_singletons  | 3       |
| n_clusters_informative | 20      |
| m_bar                  | 17.261  |
| m_kish                 | 24.270  |
| icc                    | 0.100   |

## Numerical Validation

DMAR ships a dedicated numerical-correctness test suite
(`tests/testthat/test-numerical-correctness.R` and companion files) that
pins its output against authoritative reference implementations and
published values, so results are reproducible and checkable:

| Quantity | Reference | Agreement |
|----|----|----|
| Confidence intervals for R², SMD; noncentral *t* / *F* limits; AIPE sample sizes | `MBESS` | 1e-5 to 1e-6 (integer *N*) |
| Within-subjects sphericity epsilons; factorial ANOVA / ANCOVA | `car` | machine precision |
| CFA, coefficient omega, RMSEA | `lavaan`, `semTools` | machine precision |
| [`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md), [`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md) | `irr`, `irrCAC` | machine precision (validated designs) |
| [`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md) (Nakagawa & Schielzeth), [`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md) (Rights & Sterba) | `r2mlm` | machine precision (validated model classes) |
| Multiple-comparison critical values and CIs | `multcomp`, published tables | published values |
| Dunnett’s T3; the paired-proportion score interval; Fieller’s ratio of means; bifactor and second-order omegas | `PMCMRplus`, `PropCIs`, `mratios`, `semTools` | pinned values, compared live by the maintainer’s release tooling |

The machine-precision rows hold for the model and design classes the
test suite pins, including the boundary cases the tests cover. The
mixed-effects R-squared functions match `r2mlm` and the direct Johnson
quadratic form across random-intercept, correlated random-slope, split
(`||`) random-slope, multi-slope, and nested two-level models; on split
random-slope models `performance`/`insight` omits the random-slope
variance, so DMAR agrees with `r2mlm` and the Johnson form there and not
with `performance`. The agreement procedures match `irr` and `irrCAC` on
the ordinary rating designs and on the boundary cases the tests cover:
the zero-distance case of ratio-metric Krippendorff’s alpha, all-missing
units in
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
and incomplete category sets in
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md).

## Learn More

- **Package website:** <https://yelleknek.github.io/DMAR/>
- **Start here:**
  [`vignette("DMAR")`](https://yelleknek.github.io/DMAR/articles/DMAR.md)
  (overview),
  [`vignette("dmar_output")`](https://yelleknek.github.io/DMAR/articles/dmar_output.md)
  (the tidy `dmar_tbl` output and
  `tidy()`/`glance()`/[`as_kable()`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  methods),
  [`vignette("mbess_to_dmar")`](https://yelleknek.github.io/DMAR/articles/mbess_to_dmar.md)
  (MBESS → DMAR migration).
- **Full function reference** and topical articles are on the website
  under *Reference* and *Articles*.

## Relation to MBESS

DMAR is a more general and greatly expanded reimagining of the **MBESS**
package (Kelley, 2007a, *Journal of Statistical Software*, 20(8), 1–24;
2007b, *Behavior Research Methods*, 39(4), 979–984), which has been on
CRAN for about two decades. MBESS was named for the behavioral,
educational, and social sciences, but the methods it implements apply
across the empirical sciences, and DMAR is framed for that broader
range. The notational and programming changes (uniform `snake_case`,
tidy `data.frame(term, value)` returns, native `ggplot2` plotting,
broader coverage of within-subjects, multivariate, mixed-effects, and
sequential-estimation methods) were too extensive to ship under the
MBESS name without breaking compatibility, so the work ships under a new
name. MBESS remains on CRAN, and new methods are added to DMAR. See
[`vignette("mbess_to_dmar")`](https://yelleknek.github.io/DMAR/articles/mbess_to_dmar.md)
for the MBESS → DMAR migration mapping.

## Relation to Maxwell, Delaney, and Kelley (MDK, 4th ed.)

Many of the methods in DMAR are discussed in Maxwell, Delaney, and
Kelley (2027), *Designing Experiments and Analyzing Data: A Model
Comparison Perspective* (4th ed., Routledge). Functions cite the
specific MDK chapter in their `@references` block where the underlying
method is discussed in detail. DMAR is broader than MDK: it covers
methods (e.g., minimum risk sequential estimation, equivalence testing,
mediation effect size CIs, design effects for cluster sampling) that are
outside the scope of a textbook on experimental design. Conversely, MDK
covers conceptual material that falls outside the scope of an R package.
The two are complementary.

## Related Packages

`BUCSS` (Bias and Uncertainty Corrected Sample Size; Kelley, Anderson,
and Maxwell) plans the sample size of a future study for statistical
power from the summary statistics of a prior study, correcting the prior
study’s effect for publication bias and for uncertainty. Ken Kelley is a
coauthor of `BUCSS`. Its planners cover many of the designs of DMAR’s
`ss_power_*()` family; DMAR’s power planners take the planning effect as
given, and DMAR’s `ss_aipe_*()` planners plan for the width of a
confidence interval rather than for power.

## Citation

If you use DMAR in published work, please cite

> Kelley, K. (2026). *DMAR: Design, Measurement, and Analysis in R*. R
> package version 1.1.0. <https://doi.org/10.5281/zenodo.22960335>

and, where appropriate, the original MBESS references for the
methodological lineage:

> Kelley, K. (2007). Confidence Intervals for Standardized Effect Sizes:
> Theory, Application, and Implementation. *Journal of Statistical
> Software*, 20(8), 1–24. <doi:10.18637/jss.v020.i08>

> Kelley, K. (2007). Methods for the Behavioral, Educational, and Social
> Sciences: An R package. *Behavior Research Methods*, 39(4), 979–984.
> <doi:10.3758/BF03192993>

`citation("DMAR")` reproduces the canonical BibTeX entry from R.

## Feedback

Bug reports, feature requests, and suggestions for new methods are
welcomed by email to Ken Kelley <kkelley@nd.edu> (please put “DMAR” in
the subject line). See <https://kenkelley.org> and
<https://kenkelley.org/publications/> for related publications.
