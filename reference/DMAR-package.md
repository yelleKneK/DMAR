# Design, Measurement, and Analysis in R

An R package for design, measurement, and analysis, covering effect
sizes, confidence intervals, sample size planning, reliability and
agreement, mediation analysis, equivalence testing, meta-analysis,
experimental and quasi-experimental designs, repeated measures, and
model comparison-based inference. DMAR (pronounced “Dee-Mar”) is
methodological in nature, drawing on the psychometric and statistical
traditions, and follows the methodological and applied research program
of the author. Many aspects of the package trace to the author's
methodological work and collaborations, including sample size planning
via accuracy in parameter estimation (AIPE; Kelley & Maxwell, 2003;
Kelley & Rausch, 2006; Maxwell, Kelley, & Rausch, 2008), the definition
and communication of effect sizes (Kelley & Preacher, 2012; Preacher &
Kelley, 2011), and the model comparison perspective of Maxwell, Delaney,
and Kelley (2027). Its methods and examples are chosen for research in
which the independent or dependent variables involve the person, across
psychology, sociology, education, behavioral economics, management,
marketing, and information systems.

## Details

The package covers effect size estimation, confidence intervals for
effect sizes, sample size planning, multivariate methods, factor
analysis, and certain latent variable models, and many of its functions
accept the sample estimates or results reported in published articles as
well as raw data. Sample size planning is offered under several
complementary frameworks: accuracy in parameter estimation (AIPE), power
analysis, minimum risk, and equivalence. Most exported functions return
a tidy `data.frame` with a `term` column and a numeric `value` column
(some carry additional typed columns per term), and ggplot2 is used for
the plotting functions. A few functions return the shape their task
calls for instead: model fits such as `mlmr` return a richer list-like
object with `coef` / `vcov` / `confint` methods, `descriptives` returns
a list of summary tables, and a small number of scalar utilities such as
`skewness` and `kurtosis` return a bare numeric.

DMAR builds heavily on the MBESS package (Kelley, 2007a, 2007b), which
has been on CRAN for more than two decades. MBESS was originally framed
for the behavioral, educational, and social sciences; the methods of
both packages are not specific to those fields.

**Function families.** Every exported function is listed once below, in
the family the reference index of the package site also uses
(`ci_bifactor_indices` appears in two).

- Effect sizes:

  [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
  [`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
  [`smd_trimmed`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md),
  [`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
  [`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
  [`cohen_f`](https://yelleknek.github.io/DMAR/reference/cohen_f.md),
  [`cohen_h`](https://yelleknek.github.io/DMAR/reference/cohen_h.md),
  [`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
  [`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
  [`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
  [`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
  [`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
  [`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
  [`responder_analysis`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
  [`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
  [`signal_to_noise_R2`](https://yelleknek.github.io/DMAR/reference/signal_to_noise_R2.md),
  [`expected_partial_r`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
  [`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
  [`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md),
  [`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
  [`unbiased_R2`](https://yelleknek.github.io/DMAR/reference/unbiased_R2.md),
  [`sd_unbiased`](https://yelleknek.github.io/DMAR/reference/sd_unbiased.md),
  [`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md),
  [`skewness`](https://yelleknek.github.io/DMAR/reference/skewness.md),
  [`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md).

- Confidence intervals:

  [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
  [`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
  [`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
  [`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
  [`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md),
  [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
  [`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
  [`ci_dunnett_t3`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
  [`ci_eigenvalue`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
  [`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
  [`ci_mahalanobis`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
  [`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md),
  [`ci_proportion_paired`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md),
  [`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  [`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
  [`ci_R`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
  [`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
  [`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
  [`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
  [`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
  [`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
  [`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
  [`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
  [`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
  [`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
  [`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
  [`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
  [`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
  [`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
  [`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md).

- Sample size planning, accuracy in parameter estimation (AIPE):

  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
  [`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md),
  [`ss_aipe_c_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova_sensitivity.md),
  [`ss_aipe_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
  [`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
  [`ss_aipe_cliff_delta_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
  [`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
  [`ss_aipe_crd_both_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crd_both_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crd_es_both_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_es_both_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_es_n_clusters_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_es_n_clusters_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_es_n_individuals_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_es_n_individuals_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md),
  [`ss_aipe_crd_n_clusters_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crd_n_clusters_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crd_n_individuals_fixed_budget`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crd_n_individuals_fixed_width`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md),
  [`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
  [`ss_aipe_crossed_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
  [`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
  [`ss_aipe_equivalence_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
  [`ss_aipe_equivalence_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
  [`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
  [`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
  [`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
  [`ss_aipe_icc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
  [`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
  [`ss_aipe_indirect_effect_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
  [`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
  [`ss_aipe_mean_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
  [`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
  [`ss_aipe_mixed_effects_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
  [`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
  [`ss_aipe_omega_squared_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
  [`ss_aipe_partial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
  [`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md),
  [`ss_aipe_pcm_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
  [`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
  [`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
  [`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
  [`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
  [`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md),
  [`ss_aipe_rc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc_sensitivity.md),
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
  [`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
  [`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
  [`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
  [`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md),
  [`ss_aipe_rmsea_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea_sensitivity.md),
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
  [`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md),
  [`ss_aipe_sc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md),
  [`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md),
  [`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md),
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
  [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
  [`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md),
  [`ss_aipe_sm_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm_sensitivity.md),
  [`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
  [`ss_aipe_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md),
  [`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md),
  [`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md),
  [`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
  [`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
  [`ss_aipe_two_phase_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md).

- Sample size planning, power:

  [`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
  [`ss_power_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
  [`ss_power_composite_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
  [`ss_power_composite_ancova_2group`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
  [`ss_power_composite_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
  [`ss_power_composite_factorial_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
  [`ss_power_composite_factorial_ancova_het`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
  [`ss_power_composite_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
  [`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
  [`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
  [`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
  [`ss_power_factorial_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
  [`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
  [`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
  [`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
  [`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
  [`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
  [`ss_power_r`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
  [`ss_power_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
  [`ss_power_rc`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
  [`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
  [`ss_power_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
  [`ss_power_rm_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
  [`ss_power_sc`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
  [`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
  [`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
  [`ss_power_split_plot_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md),
  [`power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
  [`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
  [`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
  [`power_density_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
  [`power_fisher_exact`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md).

- Design analysis and design utilities:

  [`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
  [`design_effect`](https://yelleknek.github.io/DMAR/reference/design_effect.md),
  [`effects_coding`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
  [`helmert_coding`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
  [`orthogonal_polynomial`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md),
  [`is_orthogonal_set`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md).

- Minimum risk and sequential estimation:

  [`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md),
  [`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md),
  [`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md),
  [`ss_seq_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md),
  [`ss_seq_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_seq_two_phase.md).

- Equivalence testing:

  [`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
  [`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
  [`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
  [`equivalence_proportion_paired`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
  [`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md).

- Reliability and agreement:

  [`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md),
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
  [`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  [`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
  [`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
  [`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
  [`R2_mixed_effects_decomposition`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
  [`regions_of_significance`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
  [`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
  [`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
  [`krippendorff_alpha`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
  [`gwet_ac`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
  [`lin_ccc`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
  [`limits_of_agreement`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
  [`loa`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
  [`content_validity_index`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md).

- ANOVA, ANCOVA, and contrasts:

  [`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md),
  [`factorial_anova`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
  [`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
  [`anova_within_two_way`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
  [`mixed_anova`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
  [`manova_split_plot`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
  [`pairwise_within`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
  [`simple_effects_AB`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
  [`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
  [`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
  [`adjusted_means`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
  [`dunn_test`](https://yelleknek.github.io/DMAR/reference/dunn_test.md),
  [`obrien_test`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
  [`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
  [`epsilon_corrections`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
  [`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md),
  [`summary_t_test`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
  [`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
  [`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
  [`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md).

- Meta-analysis and evidence synthesis:

  [`meta_contrast`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md),
  [`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
  [`meta_r`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
  [`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md),
  [`combine_p`](https://yelleknek.github.io/DMAR/reference/combine_p.md).

- Bayesian analyses:

  [`bayes_independent_t`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md),
  [`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
  [`bayes_paired_t`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md).

- Mediation:

  [`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md),
  [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md).

- Longitudinal change:

  [`analysis_of_change`](https://yelleknek.github.io/DMAR/reference/analysis_of_change.md).

- Maximum likelihood regression (FIML):

  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md),
  [`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md).

- Missing data:

  [`missing_patterns`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md),
  [`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
  [`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
  [`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
  [`df_barnard_rubin`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md).

- Multivariate and latent variable methods:

  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  [`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  [`cfa_k`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
  [`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md),
  [`efa_k`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md),
  [`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
  [`factor_congruence`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md),
  [`sem_sam`](https://yelleknek.github.io/DMAR/reference/sem_sam.md),
  [`compare_cov_structures`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
  [`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md),
  [`covmat_from_cfa`](https://yelleknek.github.io/DMAR/reference/covmat_from_cfa.md),
  [`procrustes_phi`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
  [`variance_components_mls`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md),
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md),
  [`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
  [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
  [`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
  [`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
  [`common_method_adjusted_r`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
  [`common_method_single_factor`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
  [`sem_residual_fit`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
  [`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
  [`ecvi`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
  [`simple_structure`](https://yelleknek.github.io/DMAR/reference/simple_structure.md),
  [`measurement_alignment`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
  [`dmacs`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
  [`irt_grm`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
  [`irt_information`](https://yelleknek.github.io/DMAR/reference/irt_information.md).

- Critical values:

  [`cv_bonferroni_f`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
  [`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
  [`cv_chisq`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
  [`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
  [`cv_f`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
  [`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
  [`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
  [`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
  [`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
  [`cv_z`](https://yelleknek.github.io/DMAR/reference/cv_z.md),
  [`cv`](https://yelleknek.github.io/DMAR/reference/cv.md),
  [`dbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md),
  [`pbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md),
  [`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md).

- Parameterization conversions:

  [`convert_chisq_F`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
  [`convert_cor_cov`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
  [`convert_d_or`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
  [`convert_d_r`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
  [`convert_delta_lambda`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
  [`convert_F_chisq`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
  [`convert_f_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
  [`convert_lambda_delta`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
  [`convert_lambda_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
  [`convert_or_d`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
  [`convert_r_d`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
  [`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
  [`convert_R2_f`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
  [`convert_R2_lambda`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
  [`convert_z_normal`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md),
  [`convert_Z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md).

- Variance utilities:

  [`var_alpha`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
  [`var_cv`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
  [`var_ete`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
  [`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
  [`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
  [`var_omega_squared`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
  [`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
  [`var_r`](https://yelleknek.github.io/DMAR/reference/var_r.md),
  [`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
  [`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
  [`var_smd`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
  [`var_smd_trimmed`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md).

- Noncentral distributions:

  [`moments_nc_chisq`](https://yelleknek.github.io/DMAR/reference/moments_nc_chisq.md),
  [`moments_nc_F`](https://yelleknek.github.io/DMAR/reference/moments_nc_F.md),
  [`moments_nc_t`](https://yelleknek.github.io/DMAR/reference/moments_nc_t.md),
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
  [`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
  [`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md).

- Data simulation:

  [`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
  [`simulate_ancova_factorial_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
  [`simulate_anova_data`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
  [`simulate_longitudinal_gompertz`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_gompertz.md),
  [`simulate_longitudinal_logistic`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
  [`simulate_longitudinal_negative_exponential`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md),
  [`simulate_longitudinal_polynomial`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
  [`simulate_longitudinal_richards`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md),
  [`simulate_regression_data`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md).

- Plots:

  [`plot_cfa_k`](https://yelleknek.github.io/DMAR/reference/plot_fa_confirmatory.md),
  [`plot_ci`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  [`plot_equivalence`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
  [`plot_fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/plot_fa_confirmatory.md),
  [`plot_forest`](https://yelleknek.github.io/DMAR/reference/plot_forest.md),
  [`plot_irt_information`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
  [`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
  [`plot_R2`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
  [`plot_randomization_test`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md),
  [`plot_regions_of_significance`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
  [`plot_smd`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_trajectories`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
  [`plot_trajectories_fitted`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).

- Display helpers:

  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md),
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md),
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  [`format_p`](https://yelleknek.github.io/DMAR/reference/format_p.md),
  [`as_kable`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  [`print_anova`](https://yelleknek.github.io/DMAR/reference/print_anova.md),
  [`print_summary`](https://yelleknek.github.io/DMAR/reference/print_summary.md),
  [`now`](https://yelleknek.github.io/DMAR/reference/now.md).

- Broom tidiers and model methods:

  `tidy()` and `glance()` methods for the tables above (see
  [`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)),
  and [`anova()`](https://rdrr.io/r/stats/anova.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), and
  [`confint()`](https://rdrr.io/r/stats/confint.html) methods for `mlmr`
  and `mlmr_mv` fits.

- Data sets:

  [`bessel_errors`](https://yelleknek.github.io/DMAR/reference/bessel_errors.md)
  (Bessel's grouped distribution of Bradley's astronomical observation
  errors),
  [`depression_bdi`](https://yelleknek.github.io/DMAR/reference/depression_bdi.md)
  (the three-group depression study with a pretest covariate of Maxwell,
  Delaney, and Kelley),
  [`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md)
  (Cohen's weighted kappa illustration),
  [`drinks_trial`](https://yelleknek.github.io/DMAR/reference/drinks_trial.md)
  (Smith, Meyers, and Delaney's Community Reinforcement Approach
  drinking trial),
  [`holzinger_swineford`](https://yelleknek.github.io/DMAR/reference/holzinger_swineford.md)
  (the 1939 factor analysis study),
  [`prime_time_achievement`](https://yelleknek.github.io/DMAR/reference/prime_time_achievement.md)
  (the Indiana Prime Time third grade achievement evaluation),
  [`pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md)
  (Rosenthal and Jacobson's teacher-expectancy data),
  [`teacher_expectancy`](https://yelleknek.github.io/DMAR/reference/teacher_expectancy.md)
  (Raudenbush's meta-analysis of 18 teacher-expectancy experiments), and
  [`test_market`](https://yelleknek.github.io/DMAR/reference/test_market.md)
  (Bryant and Bruvold's controlled test-market experiment for ANCOVA
  with a random covariate). Each data set's page gives its source.

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
