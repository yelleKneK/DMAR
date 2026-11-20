#' Design, Measurement, and Analysis in R
#'
#' An R package for design, measurement, and analysis, covering
#' effect sizes, confidence intervals, sample size planning,
#' reliability and agreement, mediation analysis, equivalence testing,
#' meta-analysis, experimental and quasi-experimental designs, repeated
#' measures, and model comparison-based inference. DMAR (pronounced
#' \dQuote{Dee-Mar}) is methodological in nature, drawing on the
#' psychometric and statistical traditions, and follows the
#' methodological and applied research program of the author.
#' Many aspects of the package trace to the author's methodological work
#' and collaborations, including sample size planning via accuracy in
#' parameter estimation (AIPE; Kelley & Maxwell, 2003;
#' Kelley & Rausch, 2006; Maxwell, Kelley, &
#' Rausch, 2008), the definition and communication of effect sizes
#' (Kelley & Preacher, 2012; Preacher & Kelley, 2011), and the model
#' comparison perspective of Maxwell, Delaney, and Kelley (2027). Its
#' methods and examples are chosen for research in which the
#' independent or dependent variables involve the person, across
#' psychology, sociology, education, behavioral economics, management,
#' marketing, and information systems.
#'
#' The package covers effect size estimation, confidence intervals for
#' effect sizes, sample size planning, multivariate methods, factor
#' analysis, and certain latent variable models, and many of its
#' functions accept the sample estimates or results reported in
#' published articles as well as raw data. Sample size planning is
#' offered under several complementary frameworks: accuracy in
#' parameter estimation (AIPE), power analysis, minimum risk, and
#' equivalence. Most exported
#' functions return a tidy \code{data.frame} with a \code{term}
#' column and a numeric \code{value} column (some carry additional
#' typed columns per term), and \pkg{ggplot2} is used for the plotting
#' functions. A few functions return the shape their task calls for
#' instead: model fits such as \code{mlmr} return a richer list-like
#' object with \code{coef} / \code{vcov} / \code{confint} methods,
#' \code{descriptives} returns a list of summary tables, and a small
#' number of scalar utilities such as \code{skewness} and
#' \code{kurtosis} return a bare numeric.
#'
#' DMAR builds heavily on the \pkg{MBESS} package (Kelley, 2007a,
#' 2007b), which has been on CRAN for more than two decades. \pkg{MBESS}
#' was originally framed for the behavioral, educational, and social
#' sciences; the methods of both packages are not specific to those
#' fields.
#'
#' \strong{Function families.} Every exported function is listed once
#' below, in the family the reference index of the package site also
#' uses (\code{ci_bifactor_indices} appears in two).
#' \describe{
#'   \item{Effect sizes}{\code{\link{smd}}, \code{\link{smd_c}},
#'     \code{\link{smd_trimmed}}, \code{\link{eta_squared}},
#'     \code{\link{eta_squared_generalized}}, \code{\link{omega_squared}},
#'     \code{\link{cohen_f}},
#'     \code{\link{cohen_h}}, \code{\link{cles}},
#'     \code{\link{cliff_delta}}, \code{\link{vargha_delaney_A}},
#'     \code{\link{proportion_of_superiority}},
#'     \code{\link{probability_of_superiority_paired}},
#'     \code{\link{nnt_from_smd}}, \code{\link{responder_analysis}},
#'     \code{\link{correction_for_attenuation}},
#'     \code{\link{signal_to_noise_R2}}, \code{\link{expected_partial_r}},
#'     \code{\link{expected_r}}, \code{\link{expected_R2}},
#'     \code{\link{expected_smd}}, \code{\link{unbiased_R2}},
#'     \code{\link{sd_unbiased}}, \code{\link{descriptives}},
#'     \code{\link{skewness}}, \code{\link{kurtosis}}.}
#'   \item{Confidence intervals}{\code{\link{ci_bifactor_indices}},
#'     \code{\link{ci_c}}, \code{\link{ci_c_ancova}},
#'     \code{\link{ci_c_ancova_bp}}, \code{\link{ci_crossed}},
#'     \code{\link{ci_cv}}, \code{\link{ci_dunnett}},
#'     \code{\link{ci_dunnett_t3}}, \code{\link{ci_eigenvalue}},
#'     \code{\link{ci_eta_squared}},
#'     \code{\link{ci_eta_squared_generalized}},
#'     \code{\link{ci_games_howell}}, \code{\link{ci_mahalanobis}},
#'     \code{\link{ci_omega_squared}}, \code{\link{ci_proportion}},
#'     \code{\link{ci_proportion_paired}}, \code{\link{ci_pvaf}},
#'     \code{\link{ci_r}}, \code{\link{ci_R}}, \code{\link{ci_R2}},
#'     \code{\link{ci_ratio_of_means}},
#'     \code{\link{ci_ratio_of_means_paired}}, \code{\link{ci_rc}},
#'     \code{\link{ci_reg_coef}}, \code{\link{ci_rmsea}},
#'     \code{\link{ci_sc}}, \code{\link{ci_sc_ancova}},
#'     \code{\link{ci_scheffe}}, \code{\link{ci_sm}},
#'     \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, \code{\link{ci_snr}},
#'     \code{\link{ci_src}}, \code{\link{ci_srsnr}},
#'     \code{\link{ci_tukey_kramer}}, \code{\link{ci_two_phase}}.}
#'   \item{Sample size planning, accuracy in parameter estimation
#'     (AIPE)}{\code{\link{ss_aipe_c}}, \code{\link{ss_aipe_c_ancova}},
#'     \code{\link{ss_aipe_c_ancova_sensitivity}},
#'     \code{\link{ss_aipe_c_sensitivity}},
#'     \code{\link{ss_aipe_cliff_delta}},
#'     \code{\link{ss_aipe_cliff_delta_sensitivity}},
#'     \code{\link{ss_aipe_composite_sem}},
#'     \code{\link{ss_aipe_crd_both_fixed_budget}},
#'     \code{\link{ss_aipe_crd_both_fixed_width}},
#'     \code{\link{ss_aipe_crd_es_both_fixed_budget}},
#'     \code{\link{ss_aipe_crd_es_both_fixed_width}},
#'     \code{\link{ss_aipe_crd_es_n_clusters_fixed_budget}},
#'     \code{\link{ss_aipe_crd_es_n_clusters_fixed_width}},
#'     \code{\link{ss_aipe_crd_es_n_individuals_fixed_budget}},
#'     \code{\link{ss_aipe_crd_es_n_individuals_fixed_width}},
#'     \code{\link{ss_aipe_crd_n_clusters_fixed_budget}},
#'     \code{\link{ss_aipe_crd_n_clusters_fixed_width}},
#'     \code{\link{ss_aipe_crd_n_individuals_fixed_budget}},
#'     \code{\link{ss_aipe_crd_n_individuals_fixed_width}},
#'     \code{\link{ss_aipe_crossed}},
#'     \code{\link{ss_aipe_crossed_sensitivity}},
#'     \code{\link{ss_aipe_cv}}, \code{\link{ss_aipe_cv_sensitivity}},
#'     \code{\link{ss_aipe_equivalence_r}},
#'     \code{\link{ss_aipe_equivalence_r_sensitivity}},
#'     \code{\link{ss_aipe_equivalence_smd}},
#'     \code{\link{ss_aipe_equivalence_smd_sensitivity}},
#'     \code{\link{ss_aipe_icc}}, \code{\link{ss_aipe_icc_sensitivity}},
#'     \code{\link{ss_aipe_indirect_effect}},
#'     \code{\link{ss_aipe_indirect_effect_sensitivity}},
#'     \code{\link{ss_aipe_mean}}, \code{\link{ss_aipe_mean_sensitivity}},
#'     \code{\link{ss_aipe_mixed_effects}},
#'     \code{\link{ss_aipe_mixed_effects_sensitivity}},
#'     \code{\link{ss_aipe_omega_squared}},
#'     \code{\link{ss_aipe_omega_squared_sensitivity}},
#'     \code{\link{ss_aipe_partial_r}},
#'     \code{\link{ss_aipe_partial_r_sensitivity}},
#'     \code{\link{ss_aipe_pcm}}, \code{\link{ss_aipe_pcm_sensitivity}},
#'     \code{\link{ss_aipe_r}}, \code{\link{ss_aipe_r_sensitivity}},
#'     \code{\link{ss_aipe_R2}}, \code{\link{ss_aipe_R2_sensitivity}},
#'     \code{\link{ss_aipe_rc}}, \code{\link{ss_aipe_rc_sensitivity}},
#'     \code{\link{ss_aipe_reg_coef}},
#'     \code{\link{ss_aipe_reg_coef_sensitivity}},
#'     \code{\link{ss_aipe_reliability}},
#'     \code{\link{ss_aipe_reliability_sensitivity}},
#'     \code{\link{ss_aipe_rmsea}},
#'     \code{\link{ss_aipe_rmsea_sensitivity}}, \code{\link{ss_aipe_sc}},
#'     \code{\link{ss_aipe_sc_ancova}},
#'     \code{\link{ss_aipe_sc_ancova_sensitivity}},
#'     \code{\link{ss_aipe_sc_sensitivity}},
#'     \code{\link{ss_aipe_sem_path}},
#'     \code{\link{ss_aipe_sem_path_sensitivity}},
#'     \code{\link{ss_aipe_semipartial_r}},
#'     \code{\link{ss_aipe_semipartial_r_sensitivity}},
#'     \code{\link{ss_aipe_sm}}, \code{\link{ss_aipe_sm_sensitivity}},
#'     \code{\link{ss_aipe_smd}}, \code{\link{ss_aipe_smd_sensitivity}},
#'     \code{\link{ss_aipe_src}}, \code{\link{ss_aipe_src_sensitivity}},
#'     \code{\link{ss_aipe_two_phase}},
#'     \code{\link{ss_aipe_two_phase_cost}},
#'     \code{\link{ss_aipe_two_phase_sensitivity}}.}
#'   \item{Sample size planning, power}{\code{\link{ss_power_c}},
#'     \code{\link{ss_power_c_ancova}},
#'     \code{\link{ss_power_composite_ancova}},
#'     \code{\link{ss_power_composite_ancova_2group}},
#'     \code{\link{ss_power_composite_anova}},
#'     \code{\link{ss_power_composite_factorial_ancova}},
#'     \code{\link{ss_power_composite_factorial_ancova_het}},
#'     \code{\link{ss_power_composite_factorial_anova}},
#'     \code{\link{ss_power_composite_sem}},
#'     \code{\link{ss_power_contrast}},
#'     \code{\link{ss_power_equivalence_c}},
#'     \code{\link{ss_power_factorial_ancova}},
#'     \code{\link{ss_power_factorial_anova}},
#'     \code{\link{ss_power_indirect_effect}},
#'     \code{\link{ss_power_mixed_effects}},
#'     \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_pcm}},
#'     \code{\link{ss_power_r}}, \code{\link{ss_power_R2}},
#'     \code{\link{ss_power_R2_change}},
#'     \code{\link{ss_power_R2_sensitivity}}, \code{\link{ss_power_rc}},
#'     \code{\link{ss_power_reg_coef}},
#'     \code{\link{ss_power_reg_coef_sensitivity}},
#'     \code{\link{ss_power_rm_anova}}, \code{\link{ss_power_sc}},
#'     \code{\link{ss_power_sem}}, \code{\link{ss_power_smd}},
#'     \code{\link{ss_power_split_plot_anova}},
#'     \code{\link{power_equivalence_c}},
#'     \code{\link{power_equivalence_md}},
#'     \code{\link{power_equivalence_md_plot}},
#'     \code{\link{power_density_equivalence_md}},
#'     \code{\link{power_fisher_exact}}.}
#'   \item{Design analysis and design
#'     utilities}{\code{\link{design_consequences}},
#'     \code{\link{design_effect}}, \code{\link{effects_coding}},
#'     \code{\link{helmert_coding}}, \code{\link{orthogonal_polynomial}},
#'     \code{\link{is_orthogonal_set}}.}
#'   \item{Minimum risk and sequential estimation}{\code{\link{mr_smd}},
#'     \code{\link{mr_cv}}, \code{\link{ss_seq_c}},
#'     \code{\link{ss_seq_c_sensitivity}},
#'     \code{\link{ss_seq_two_phase}}.}
#'   \item{Equivalence testing}{\code{\link{equivalence_smd}},
#'     \code{\link{equivalence_r}}, \code{\link{equivalence_c}},
#'     \code{\link{equivalence_proportion_paired}},
#'     \code{\link{equivalence_rmsea}}.}
#'   \item{Reliability and agreement}{\code{\link{reliability}},
#'     \code{\link{reliability_alpha}}, \code{\link{reliability_H}},
#'     \code{\link{reliability_kr20}}, \code{\link{reliability_omega}},
#'     \code{\link{reliability_omega_categorical}}, \code{\link{icc}},
#'     \code{\link{icc_lmer}}, \code{\link{R2_mixed_effects}},
#'     \code{\link{R2_mixed_effects_decomposition}},
#'     \code{\link{regions_of_significance}}, \code{\link{cohen_kappa}},
#'     \code{\link{fleiss_kappa}}, \code{\link{krippendorff_alpha}},
#'     \code{\link{gwet_ac}}, \code{\link{lin_ccc}},
#'     \code{\link{limits_of_agreement}}, \code{\link{loa}},
#'     \code{\link{content_validity_index}}.}
#'   \item{ANOVA, ANCOVA, and contrasts}{\code{\link{ancova}},
#'     \code{\link{factorial_anova}}, \code{\link{anova_within}},
#'     \code{\link{anova_within_two_way}}, \code{\link{mixed_anova}},
#'     \code{\link{manova_split_plot}}, \code{\link{pairwise_within}},
#'     \code{\link{simple_effects_AB}}, \code{\link{contrast_test}},
#'     \code{\link{contrast_adjusted}}, \code{\link{adjusted_means}},
#'     \code{\link{dunn_test}}, \code{\link{obrien_test}},
#'     \code{\link{mauchly_test}}, \code{\link{epsilon_corrections}},
#'     \code{\link{welch_t}}, \code{\link{summary_t_test}},
#'     \code{\link{correlations_test}}, \code{\link{randomization_test}},
#'     \code{\link{randomization_test_paired}}.}
#'   \item{Meta-analysis and evidence
#'     synthesis}{\code{\link{meta_contrast}}, \code{\link{meta_es}},
#'     \code{\link{meta_r}}, \code{\link{meta_smd}},
#'     \code{\link{combine_p}}.}
#'   \item{Bayesian analyses}{\code{\link{bayes_independent_t}},
#'     \code{\link{bayes_one_sample_t}}, \code{\link{bayes_paired_t}}.}
#'   \item{Mediation}{\code{\link{mediate}},
#'     \code{\link{mediation_mbco}}.}
#'   \item{Longitudinal change}{\code{\link{analysis_of_change}}.}
#'   \item{Maximum likelihood regression (FIML)}{\code{\link{mlmr}},
#'     \code{\link{mlmr_mv}}.}
#'   \item{Missing data}{\code{\link{missing_patterns}},
#'     \code{\link{missing_coverage}}, \code{\link{cov_fiml}},
#'     \code{\link{fraction_missing_information}},
#'     \code{\link{df_missing}}, \code{\link{df_barnard_rubin}}.}
#'   \item{Multivariate and latent variable methods}{\code{\link{cfa_1}},
#'     \code{\link{cfa_2}}, \code{\link{fa_confirmatory}},
#'     \code{\link{cfa_k}}, \code{\link{fa_bifactor}},
#'     \code{\link{fa_second_order}}, \code{\link{fa_exploratory}},
#'     \code{\link{efa_k}}, \code{\link{factor_count}},
#'     \code{\link{factor_congruence}}, \code{\link{parallel_analysis}},
#'     \code{\link{sem_sam}}, \code{\link{compare_cov_structures}},
#'     \code{\link{cov_sem}}, \code{\link{covmat_from_cfa}},
#'     \code{\link{procrustes_phi}},
#'     \code{\link{variance_components_mls}},
#'     \code{\link{measurement_invariance}}, \code{\link{htmt}},
#'     \code{\link{average_variance_extracted}},
#'     \code{\link{bifactor_indices}}, \code{\link{ci_bifactor_indices}},
#'     \code{\link{factor_scores}}, \code{\link{common_method_marker}},
#'     \code{\link{common_method_adjusted_r}},
#'     \code{\link{common_method_single_factor}},
#'     \code{\link{sem_fit_measures}},
#'     \code{\link{sem_exact_fit_bootstrap}},
#'     \code{\link{sem_residual_fit}},
#'     \code{\link{sem_information_criteria}}, \code{\link{compare_sem}},
#'     \code{\link{ecvi}}, \code{\link{simple_structure}},
#'     \code{\link{measurement_alignment}}, \code{\link{dmacs}},
#'     \code{\link{irt_grm}}, \code{\link{irt_information}}.}
#'   \item{Critical values}{\code{\link{cv_bonferroni_f}},
#'     \code{\link{cv_bryant_paulson}}, \code{\link{cv_chisq}},
#'     \code{\link{cv_dunnett}}, \code{\link{cv_f}},
#'     \code{\link{cv_scheffe}}, \code{\link{cv_smm}}, \code{\link{cv_t}},
#'     \code{\link{cv_tukey_hsd}}, \code{\link{cv_z}}, \code{\link{cv}},
#'     \code{\link{dbryant_paulson}}, \code{\link{pbryant_paulson}},
#'     \code{\link{qbryant_paulson}}.}
#'   \item{Parameterization conversions}{\code{\link{convert_chisq_F}},
#'     \code{\link{convert_cor_cov}}, \code{\link{convert_d_or}},
#'     \code{\link{convert_d_r}}, \code{\link{convert_delta_lambda}},
#'     \code{\link{convert_F_chisq}}, \code{\link{convert_f_R2}},
#'     \code{\link{convert_lambda_delta}},
#'     \code{\link{convert_lambda_R2}}, \code{\link{convert_or_d}},
#'     \code{\link{convert_r_d}}, \code{\link{convert_r_Z}},
#'     \code{\link{convert_R2_f}}, \code{\link{convert_R2_lambda}},
#'     \code{\link{convert_z_normal}}, \code{\link{convert_Z_r}}.}
#'   \item{Variance utilities}{\code{\link{var_alpha}},
#'     \code{\link{var_cv}}, \code{\link{var_ete}}, \code{\link{var_icc}},
#'     \code{\link{var_indirect_effect}}, \code{\link{var_omega_squared}},
#'     \code{\link{var_partial_r}}, \code{\link{var_r}},
#'     \code{\link{var_R2}}, \code{\link{var_semipartial_r}},
#'     \code{\link{var_smd}}, \code{\link{var_smd_trimmed}}.}
#'   \item{Noncentral distributions}{\code{\link{moments_nc_chisq}},
#'     \code{\link{moments_nc_F}}, \code{\link{moments_nc_t}},
#'     \code{\link{ci_nc_t}}, \code{\link{ci_nc_F}},
#'     \code{\link{ci_nc_chisq}}.}
#'   \item{Data simulation}{\code{\link{simulate_ancova_data}},
#'     \code{\link{simulate_ancova_factorial_data}},
#'     \code{\link{simulate_anova_data}},
#'     \code{\link{simulate_longitudinal_gompertz}},
#'     \code{\link{simulate_longitudinal_logistic}},
#'     \code{\link{simulate_longitudinal_negative_exponential}},
#'     \code{\link{simulate_longitudinal_polynomial}},
#'     \code{\link{simulate_longitudinal_richards}},
#'     \code{\link{simulate_regression_data}}.}
#'   \item{Plots}{\code{\link{plot_cfa_k}}, \code{\link{plot_ci}},
#'     \code{\link{plot_equivalence}}, \code{\link{plot_fa_confirmatory}},
#'     \code{\link{plot_forest}}, \code{\link{plot_irt_information}},
#'     \code{\link{plot_mediation_mbco}}, \code{\link{plot_R2}},
#'     \code{\link{plot_randomization_test}},
#'     \code{\link{plot_regions_of_significance}}, \code{\link{plot_smd}},
#'     \code{\link{plot_trajectories}},
#'     \code{\link{plot_trajectories_fitted}}.}
#'   \item{Display helpers}{\code{\link{dmar_tbl}},
#'     \code{\link{dmar_metadata}}, \code{\link{dmar_bootstrap}},
#'     \code{\link{format_p}}, \code{\link{as_kable}},
#'     \code{\link{results_sentence}}, \code{\link{print_anova}},
#'     \code{\link{print_summary}}, \code{\link{now}}.}
#'   \item{Broom tidiers and model methods}{\code{tidy()} and
#'     \code{glance()} methods for the tables above (see
#'     \code{\link{dmar_tidiers}}), and \code{anova()}, \code{coef()},
#'     \code{vcov()}, and \code{confint()} methods for \code{mlmr} and
#'     \code{mlmr_mv} fits.}
#'   \item{Data sets}{\code{\link{bessel_errors}} (Bessel's grouped
#'     distribution of Bradley's astronomical observation errors),
#'     \code{\link{depression_bdi}} (the three-group depression study
#'     with a pretest covariate of Maxwell, Delaney, and Kelley),
#'     \code{\link{diagnosis_agreement}} (Cohen's weighted kappa
#'     illustration), \code{\link{drinks_trial}} (Smith, Meyers, and
#'     Delaney's Community Reinforcement Approach drinking trial),
#'     \code{\link{holzinger_swineford}} (the 1939 factor analysis study),
#'     \code{\link{prime_time_achievement}}
#'     (the Indiana Prime Time third grade achievement evaluation),
#'     \code{\link{pygmalion}}
#'     (Rosenthal and Jacobson's teacher-expectancy data),
#'     \code{\link{teacher_expectancy}} (Raudenbush's meta-analysis of
#'     18 teacher-expectancy experiments), and \code{\link{test_market}}
#'     (Bryant and Bruvold's controlled test-market experiment
#'     for ANCOVA with a random covariate). Each data set's page gives
#'     its source.}
#' }
#'
#' \strong{Feedback.} Bug reports, feature requests, and suggestions
#' for new methods are welcomed by email to Ken Kelley
#' \email{kkelley@@nd.edu} (please put \dQuote{DMAR} in the subject
#' line). See \url{https://kenkelley.org} for Ken Kelley's web site,
#' \url{https://kenkelley.org/publications/} for related publications,
#' and \url{https://github.com/yelleKneK/DMAR} for the project's
#' GitHub page.
#'
#' @references
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#' sizes: Theory, application, and implementation. \emph{Journal of
#' Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#' social sciences: An R package. \emph{Behavior Research Methods,
#' 39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#' regression: Obtaining regression coefficients that are accurate, not
#' simply significant. \emph{Psychological Methods, 8}(3), 305--321.
#' \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#' \emph{Psychological Methods, 17}(2), 137--152. \doi{10.1037/a0028086}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#' standardized mean difference: Accuracy in parameter estimation via
#' narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#' 363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#' experiments and analyzing data: A model comparison perspective}
#' (4th ed.). Routledge.
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size
#' planning for statistical power and accuracy in parameter estimation.
#' \emph{Annual Review of Psychology, 59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for
#' mediation models: Quantitative strategies for communicating indirect
#' effects. \emph{Psychological Methods, 16}(2), 93--115.
#' \doi{10.1037/a0022658}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @aliases DMAR dmar
#' @keywords package
"_PACKAGE"

# `.data` is the ggplot2/rlang pronoun for tidy-data column references inside
# aes() mappings (e.g., aes(x = .data[[time_var]])); declaring it here keeps
# R CMD check's "no visible binding for global variable" check quiet without
# requiring rlang to be in Imports.
utils::globalVariables(".data")
