library(testthat)
library(DMAR)

# On CRAN (NOT_CRAN unset) the check runs a curated subset of the suite:
# the mechanical detectors the review rounds promised (rd_hygiene and
# ss_sensitivity_filename, which run whole), and the published-value and
# oracle anchors that pin the package's numerics across platforms, one
# family to a file: the pinned MBESS values (numerical-correctness),
# metafor (meta_es), multcomp (ci_dunnett), car (factorial_anova),
# lavaan (ecvi), the Welch-Satterthwaite df of Maxwell, Delaney, and
# Kelley's Equation 5.14 (ci_games_howell), Kelley's (2007) example,
# planned at 101 and 142 by the exact noncentral t past the range of
# pt(), where the article reports 99 and 141 (ss_aipe_cv), the exact
# noncentral t against its oracle (ci_nc_t), SIEVER's recorded numbers
# (reliability_H), the published studentized maximum modulus (cv_smm),
# Yuan et al.'s (2016) RMSEA (equivalence_rmsea), Shrout and Fleiss's
# (1979) tables (icc), Kelley and Rausch's (2011) Table 1 (ss_aipe_pcm),
# and Raudenbush and Liu's (2001) benchmark (ss_power_pcm). Twenty-four
# files keep one block of a closed-form sample size planner (no Monte
# Carlo, no model fitting) on CRAN: the G*Power 3.1 manual
# (ss_power_R2_change), Kelley and Maxwell's (2003) plans for
# ss_aipe_src and ss_aipe_reg_coef (ss_aipe_rc), the pinned plans of Lai
# and Kelley's (2012) approximation (ss_aipe_sc_ancova), the exact
# probability of an ss_aipe_sc plan computed apart from the package, the
# exact noncentral t oracle past the range of pt() (ss_power_c,
# ss_power_smd, ss_power_indirect_effect,
# ss_power_composite_ancova_2group), and, in the other planner files, a
# closed form, scan, or integral written in the test file itself.
# Closed-form planners that no file below tests are tested in the full
# suite only, among them ss_aipe_c (the 1.0.0 search over a grid, 4.5 s)
# and ss_aipe_R2 (the MBESS plan, a search over N, 0.2 s), whose blocks
# cost more than 150 ms, and ss_aipe_reliability, ss_aipe_rc, and
# ss_power_equivalence_c; cran-comments.md names the fifteen that 1.0.0
# tested here. In each anchor file but numerical-correctness the anchor
# sits at or near the top (moved there where it was not), and a
# top-level skip_on_cran() after it skips the rest of the file on CRAN
# with a single skip: a skip_on_cran() inside a block costs about 8 ms
# under the check's reporter, more than most blocks take to run. Where
# NOT_CRAN is true (devtools::test() and the release gate's full suite
# set it to "true") every file runs, Monte Carlo blocks included.
# NOT_CRAN is read as testthat's skip_on_cran() reads it, with
# as.logical(), so the two agree: "false", or any value that is not
# true, selects the subset and skip_on_cran() skips. Through 1.0.x any
# nonempty value, "false" included, selected every file while
# skip_on_cran() still skipped. The subset was curated 2026-09-07 in
# response to the incoming pretest's overall-checktime NOTE and cut on
# 2026-11-13 from 110 files (4,678 expectations, 23 s) to 16 (266
# expectations, under 5 s), to which the planner anchors then added 24
# (40 files, 527 expectations, under 6 s); the per-family smoke files it
# dropped run in the full suite as before. It runs in one process, is
# maintained alongside tools/release_gate.R, and no file below may be
# removed without rechecking the promises in cran-comments.md.
cran_files <- c(
  "ci_dunnett", "ci_games_howell", "ci_nc_t", "cv_smm", "ecvi",
  "equivalence_rmsea", "factorial_anova", "icc", "meta_es",
  "numerical-correctness", "rd_hygiene", "reliability_H",
  "ss_aipe_crossed", "ss_aipe_cv", "ss_aipe_equivalence_r",
  "ss_aipe_equivalence_smd", "ss_aipe_icc", "ss_aipe_indirect_effect",
  "ss_aipe_mean", "ss_aipe_mixed_effects", "ss_aipe_partial_r",
  "ss_aipe_pcm", "ss_aipe_r", "ss_aipe_rc", "ss_aipe_sc",
  "ss_aipe_sc_ancova", "ss_aipe_semipartial_r", "ss_aipe_two_phase",
  "ss_aipe_two_phase_cost", "ss_power_R2_change", "ss_power_c",
  "ss_power_composite_ancova_2group", "ss_power_contrast",
  "ss_power_factorial_ancova", "ss_power_indirect_effect", "ss_power_pcm",
  "ss_power_reg_coef", "ss_power_smd", "ss_power_split_plot_anova",
  "ss_sensitivity_filename"
)

if (isTRUE(as.logical(Sys.getenv("NOT_CRAN")))) {
  test_check("DMAR")
} else {
  Sys.setenv(TESTTHAT_PARALLEL = "false")
  test_check("DMAR",
             filter = paste0("^(", paste(cran_files, collapse = "|"), ")$"))
}
