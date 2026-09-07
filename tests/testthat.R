library(testthat)
library(DMAR)

# On CRAN (NOT_CRAN unset) the check runs a curated subset of the suite:
# the mechanical hygiene detectors the review rounds promised run there,
# the published-value and oracle anchor files that pin the package's
# numerics across platforms (every closed-form ss_aipe_* planner among
# them), and per-family smoke. Everywhere else (devtools::test(), the
# release gate, development) NOT_CRAN is set and every file runs, Monte
# Carlo blocks included. The subset was curated 2026-09-07 in response
# to the incoming pretest's overall-checktime NOTE; it is maintained
# alongside tools/release_gate.R, and no file below may be removed
# without rechecking the promises in cran-comments.md.
cran_files <- c(
  "R2_mixed_effects_decomposition", "bayes_t", "cfa_k", "ci_R2",
  "ci_dunnett", "ci_games_howell", "ci_nc_F", "ci_nc_t", "ci_rmsea",
  "cles", "combine_p", "common_method_single_factor", "contrast_adjusted",
  "cv_smm", "dmacs", "dmar_tbl", "ecvi", "expected_partial_r", "expected_r",
  "factorial_anova", "fleiss_kappa", "hyperg_2F1", "icc", "irt_information",
  "lin_ccc", "meta_contrast", "meta_es", "meta_r", "meta_smd", "moments_nc_F",
  "moments_nc_t", "nnt_from_smd", "numerical-correctness", "obrien_test",
  "omega_squared", "power_equivalence_c", "power_equivalence_quadrature",
  "power_fisher_exact", "randomization_test", "rd_hygiene", "reliability_alpha",
  "simple_structure", "smd_trimmed", "ss_aipe_R2", "ss_aipe_c", "ss_aipe_c_ancova",
  "ss_aipe_c_sensitivity", "ss_aipe_cliff_delta", "ss_aipe_composite_sem",
  "ss_aipe_crd", "ss_aipe_crd_es", "ss_aipe_cv", "ss_aipe_equivalence_r",
  "ss_aipe_equivalence_smd", "ss_aipe_icc", "ss_aipe_indirect_effect",
  "ss_aipe_mixed_effects", "ss_aipe_omega_squared", "ss_aipe_partial_r",
  "ss_aipe_pcm", "ss_aipe_r", "ss_aipe_rc", "ss_aipe_reg_coef", "ss_aipe_reliability",
  "ss_aipe_rmsea", "ss_aipe_rmsea_sensitivity", "ss_aipe_sc", "ss_aipe_sc_ancova",
  "ss_aipe_sem_path", "ss_aipe_semipartial_r", "ss_aipe_sm", "ss_aipe_smd",
  "ss_power_contrast", "ss_power_pcm", "ss_power_reg_coef", "ss_sensitivity_filename",
  "var_ete", "var_indirect_effect", "vargha_delaney_A", "variance_components_mls",
  "welch_t"
)

if (nzchar(Sys.getenv("NOT_CRAN"))) {
  test_check("DMAR")
} else {
  test_check("DMAR",
             filter = paste0("^(", paste(cran_files, collapse = "|"), ")$"))
}
