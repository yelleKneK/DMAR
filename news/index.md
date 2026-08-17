# Changelog

## DMAR 1.0.0

First public release of DMAR (pronounced “Dee-Mar,” for “Design,
Measurement, and Analysis in R”), a greatly expanded reimagining of the
MBESS package.

### An Adversarial Quality Control Pass

- A package-wide audit (numerical verification against complex-step
  derivatives and Monte Carlo, documentation contracts run rather than
  read, and mechanical convention gates) closed out with these fixes.
  [`analysis_of_change()`](https://yelleknek.github.io/DMAR/reference/analysis_of_change.md)
  now fits an intercept-only polynomial (`order = 0`) under both
  methods; it previously stopped inside
  [`stats::poly()`](https://rdrr.io/r/stats/poly.html), and its no-fit
  message for the polynomial now counts occasions instead of suggesting
  `start`, which the polynomial rejects.

- The equivalence sensitivity siblings speak their parents’ language:
  [`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md)
  takes `delta_lower` / `delta_upper` and
  [`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md)
  takes `rho_lower` / `rho_upper`, positive magnitudes with the upper
  bound required and the lower defaulting to symmetric, exactly as in
  [`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
  and
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md).
  The signed `equivalence_lower` / `equivalence_upper` spellings, and
  their silent default of plus or minus 0.20, are gone.

- Every AIPE planner that reports an expected interval width names the
  row `ci_width_expected` and reports the expected full width:
  [`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
  and
  [`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md)
  previously reported the half-width as `ci_half_width_expected` while
  [`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md)
  reported the full width, a trap for anyone extracting the row
  programmatically.

- The equivalence trio attaches the `conf_level` attribute the `ci_*`
  family always carried, and the
  [`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
  verdict label is spelled `"Noninferior only"`, matching the solid
  `noninferior` the rest of its page uses.
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md)
  now validates raw data (matching lengths, at least 4 complete pairs)
  the way its summary-statistics path always did.
  [`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)
  and
  [`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)
  insist on a single value, as documented; vector input previously
  recycled the term column into duplicated rows. A follow-up ruling
  extended the same single-value guard to the whole scalar conversion
  family
  ([`convert_R2_f()`](https://yelleknek.github.io/DMAR/reference/convert_R2.md)
  and its three siblings,
  [`convert_delta_lambda()`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md)
  and its inverse, and
  [`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)),
  with domain checks where a map’s algebra ends (an `R2` at or past 1, a
  nonpositive group size), so no conversion can recycle vector input
  into duplicated rows.

- Documentation corrections from the same audit: the four nonlinear
  simulators state that the first-order delta method behind
  `reliability` and `reliability_by_occasion` can drift from the
  realized variance ratio when random variances are large, and document
  the `schedule` attribute; the Richards example comment places the
  delta = 3 inflection at about 63% of total change, per its own
  formula; the mixed-effects examples on
  [`?analysis_of_change`](https://yelleknek.github.io/DMAR/reference/analysis_of_change.md)
  now run live.

### One Prefix for Every Confidence Interval

- The four multiple-comparison interval functions moved into the `ci_*`
  family: `ci_dunnett`, `ci_scheffe`, `ci_games_howell`, and
  `ci_tukey_kramer` replace `dunnett_ci`, `scheffe_ci`,
  `games_howell_ci`, and `tukey_kramer_ci`. The suffix forms were the
  only four exports naming a confidence interval outside the family
  prefix, and each function now pairs with its critical-value sibling
  (`cv_dunnett`, `cv_scheffe`, `cv_tukey_hsd`). The package is
  unreleased, so the old names are gone rather than aliased.

- The correlation intervals now carry the names of their estimands,
  matching the rest of the correlation family (`ss_aipe_r`,
  `ss_power_r`, `var_r`, `expected_r`, `equivalence_r`):
  [`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  is the confidence interval for the simple Pearson correlation (renamed
  from the MBESS-heritage `ci_cc()`, with the estimate row renamed from
  `est_cor` to `r`), and
  [`ci_R()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  is the interval for the population multiple correlation coefficient
  (its former lowercase alias spelling is gone). Cross-references across
  the package were audited against this distinction and every one now
  points at the interval it meant. Because the two natural file names
  differ only by case, which case-insensitive filesystems refuse, both
  functions live in `R/ci_correlation.R` and share the `ci_correlation`
  help page.

- Kish’s design effect function is
  [`design_effect()`](https://yelleknek.github.io/DMAR/reference/design_effect.md);
  the abbreviation `deft()` is gone as a function name, while the
  returned table keeps its `design_effect` and `deft` rows, deft being
  the standard term for the standard error inflation factor, the square
  root of the design effect.

- The remaining exported alias pairs were resolved to single names.
  Fisher’s Z is written with a capital Z throughout the package, because
  it is the variance-stabilizing transform of a correlation and not a
  z-score, so the converts are
  [`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)
  and
  [`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)
  only (the lowercase spellings are gone) and the prose and math were
  swept to match. The limits of agreement function is
  [`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
  with
  [`loa()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md)
  kept as its short alias and `bland_altman_loa()` gone (the method is
  named for what it computes; Bland and Altman are credited in the
  references).
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  and
  [`covmat_from_cfa()`](https://yelleknek.github.io/DMAR/reference/covmat_from_cfa.md)
  are the only names for those functions; the `reliability_omega_c()`
  and `covmat_from_cfm()` aliases are gone. The dataset alias bindings
  `HS_Data` and `Prime_Time` are gone as well: every dataset carries
  exactly one name, the documented snake_case one, and unlike the alias
  bindings the canonical names work with
  [`data()`](https://rdrr.io/r/utils/data.html).

### Fitting Change Models, Linear and Nonlinear

- The vocabulary of the longitudinal simulators names what the rows are:
  the identifier column identifies units (persons, animals, trees,
  classrooms), `n` counts units, and the former `group` column is
  `population`, one level per data generating parameter vector.

- Measurement schedules can now be unit-specific, in every simulator: in
  place of a shared `target_times` grid, `time_range = c(lower, upper)`
  draws each unit’s own measurement times uniformly between the bounds,
  with `occasions` fixing the number of times per unit or `c(min, max)`
  letting it vary, so designs such as age at testing rather than grade
  at testing are simulated directly. The polynomial simulator
  additionally requires every unit’s count of occasions to reach
  `P + 1`, so each simulated trajectory identifies the polynomial it
  came from.

- [`analysis_of_change()`](https://yelleknek.github.io/DMAR/reference/analysis_of_change.md)
  closes the loop the simulators open: it fits any of the four nonlinear
  change models, or a polynomial of any order, to longitudinal data. The
  default two-stage method fits one curve per unit, from that unit’s
  data alone, with data-driven starting values, and summarizes the
  unit-level parameters by their mean, their standard deviation and
  variance across units (the individual differences), and the standard
  error of the mean; `method = "mixed"` instead estimates the proper
  random-coefficients model simultaneously,
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) for the
  polynomial and
  [`nlme::nlme()`](https://rdrr.io/pkg/nlme/man/nlme.html) for the
  nonlinear curves started at the two-stage estimates, so `sd_units`
  becomes a variance component purged of estimation noise. The help page
  positions the function candidly against
  [`nlme::lmList()`](https://rdrr.io/pkg/nlme/man/lmList.html),
  [`nlme::nlsList()`](https://rdrr.io/pkg/nlme/man/nlsList.html), and
  the self-starting `SS*` curves: what it adds is the landmark
  parameterizations, the simulator match, and one tidy interface across
  linear and nonlinear change. With a single trajectory (or `id = NULL`)
  it reduces to one nonlinear least squares fit reported with asymptotic
  standard errors, so N = 1 is one unit’s change. Non-converging units
  are dropped with a single counted warning and the effective count
  travels as the `n_used` attribute; the full matrix of unit-level
  estimates rides along for plotting or as starting values for a
  simultaneous [`nlme::nlme()`](https://rdrr.io/pkg/nlme/man/nlme.html)
  fit, and the help page states plainly that the between-unit spread of
  estimates includes estimation noise, which a variance-component model
  separates.

### Nonlinear Growth Curve Simulators

- Four random-coefficients nonlinear change simulators join
  [`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
  sharing its interface (several populations of units, between-unit
  parameter variances and correlations, level-one error by variance or
  by target reliability, error correlation structures, assessment-time
  jitter) and its long-format return that feeds
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
  and nonlinear mixed-model fitters:
  [`simulate_longitudinal_negative_exponential()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md)
  (asymptotic regression),
  [`simulate_longitudinal_logistic()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
  [`simulate_longitudinal_gompertz()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_gompertz.md),
  and
  [`simulate_longitudinal_richards()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md),
  whose shape parameter delta subsumes the logistic (delta = 1) and the
  Gompertz (delta -\> 0) as special cases, relations the tests pin
  exactly. The parameterizations are those of Kelley (2005,
  dissertation; 2008, Methodology), in which every parameter is a
  landmark of the change process (floor, ceiling, moment of fastest
  change, curvature, inflection height) and the intercept-shifting zeta
  frees the lower asymptote from zero; the Richards application follows
  Guo, Cheng, and Kelley (2016). A new vignette, “Nonlinear Growth
  Curves and the Meaning of Their Parameters”, illustrates all four
  curves, the Richards unification, and the polynomial comparison: a
  ninth-order polynomial needs ten uninterpretable coefficients to track
  a four parameter Gompertz inside the data and still collapses the
  moment it extrapolates.

### One Word for the Equivalence Family

- The equivalence tests are named for the question they answer:
  [`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
  and
  [`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
  replace `tost_smd()`, `tost_r()`, and `tost_c()`, and the AIPE
  planners follow
  ([`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
  [`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
  with their sensitivity siblings). The whole family now speaks the same
  word as the planning and plotting surfaces that always did
  ([`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
  [`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
  [`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md)).
  The procedure is unchanged, Schuirmann’s two one-sided tests, and each
  help page says so by that name; the reported TOST -value keeps its
  standard `p_tost` term. The old names are gone rather than aliased.

### Sample Size Planning for Correlation Equivalence

- [`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md)
  plans the sample size for an equivalence question about a Pearson
  correlation: the smallest whose 100(1 - 2 alpha)% Fisher’s interval
  (the interval
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md)
  and
  [`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  invert) has expected width at or below the target, with an optional
  Monte Carlo assurance correction. At the conservative default planning
  value of 0 the answer has a closed form, and the tests anchor the
  search to it. Its Monte Carlo sibling
  [`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md)
  follows the family API and reports, alongside the family’s width and
  coverage summaries, the realized proportion of equivalence verdicts
  inside the chosen bounds. This fills the gap noted when the
  equivalence family was reviewed: the SMD had an AIPE planner for its
  equivalence interval and the correlation did not.

### The Examples Teach From the Package’s Own Data

- The help-page examples now use DMAR’s own datasets rather than the
  base-R teaching sets: Dunnett-style comparisons run on
  `depression_bdi` with its wait list control, the studentized range
  procedures on `test_market`’s six panels, the unequal-variance and
  rank methods on `drinks_trial`’s skewed outcome, the factorial and
  generalized eta squared examples on `pygmalion`’s manipulated
  treatment crossed with measured grade, and the correlation,
  regression, reliability, and multivariate examples on the
  `holzinger_swineford` battery, whose second-form tests carry real
  missingness that now powers the FIML demonstrations in
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)
  in place of artificially punched holes. Examples that fit lme4 models
  keep
  [`lme4::sleepstudy`](https://rdrr.io/pkg/lme4/man/sleepstudy.html),
  since those pages require lme4 regardless. The
  [`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md) example
  and tests likewise moved from lavaan’s copy of the 1939 data to the
  package’s own `holzinger_swineford`. Numeric claims in example
  comments were recomputed against the new output throughout.
- Spell checking is now wired in: DESCRIPTION declares `Language: en-US`
  and `inst/WORDLIST` carries the package’s technical vocabulary, so
  `devtools::spell_check()` runs clean.

### A Quieter, Faster Test Suite and ggplot2 4.0 Compatibility

- [`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md)
  no longer passes the `fatten` argument of
  [`ggplot2::geom_pointrange()`](https://ggplot2.tidyverse.org/reference/geom_linerange.html),
  which ggplot2 4.0.0 deprecates. The point size is set through the
  `size` aesthetic at the value that draws the same figure on every
  supported ggplot2 (\>= 3.4.0), so building the plot under ggplot2 4.0
  no longer raises a deprecation warning.
- [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  and
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  now validate their inputs before the courtesy message that a default
  call reports no confidence interval, so a call that is about to fail
  with an informative error no longer receives advice first.
- The test suite runs its files in parallel
  (`Config/testthat/parallel: true`, with the slowest files scheduled
  first through `Config/testthat/start-first`; the testthat floor moves
  to 3.2.0). The heaviest Monte Carlo checks were recalibrated to
  smaller replication counts whose assertions still hold with room to
  spare, duplicated planner calls across neighboring tests were
  consolidated into single calls, and tests no longer leak messages,
  notes, or printed tables into the run’s output. A full
  `devtools::test()` on the maintainer’s machine dropped from about 20
  minutes to about 5, with no warnings and nothing skipped.

### One Return Schema Across the Sensitivity Family

- Every `ss_aipe_*_sensitivity()` member now returns the family’s
  documented schema: `mean_ci_width` / `median_ci_width` / `sd_ci_width`
  for the realized interval widths, `pct_ci_less_w` for the proportion
  of intervals at or below the planning width, `pct_ci_miss_low` /
  `pct_ci_miss_high` / `total_type_I_error` for the empirical
  non-coverage (proportions on the 0 to 1 scale, so the total is the sum
  of its tails), `mean_X` / `median_X` / `sd_X` for the member’s own
  estimand, and input echoes named for their unit (`total_N` or
  `n_per_group` for the evaluated size, `true_X`, `estimated_X`,
  `width`, `conf_level`, and, when one was supplied, `assurance`). The
  core term vector lives once, as the internal registry constant
  `.SS_AIPE_SENS_CORE_TERMS` in `R/dmar_tidiers.R`, and a family-wide
  contract test asserts every member against it. Every member’s
  `@return` now lists exactly the rows the function returns.
- The term names that departed from the schema moved to it, loudly (the
  package is unreleased; no old name survives).
  [`ss_aipe_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md)
  drops the lowercase `_r2` suffixes and `pct_less_w` (`mean_r2` is
  `mean_R2`, `mean_ci_width_r2` is `mean_ci_width`, the realized-limit
  rows are `mean_lower_limit` and siblings, the one-sided widths are
  `mean_ci_width_lower` / `mean_ci_width_upper`); the four term-matching
  reads inside
  [`ss_aipe_R2()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)’s
  assurance search were updated with it.
  [`ss_aipe_sc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md),
  [`ss_aipe_sc_ancova_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md),
  and
  [`ss_aipe_sm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm_sensitivity.md)
  retire `mean_full_width` and `pct_Width_obs_narrower` /
  `pct_width_obs_narrower`;
  [`ss_aipe_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md)
  retires `pct_less_desired`;
  [`ss_aipe_c_ancova_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova_sensitivity.md)
  retires `width_narrower` and `mean_width_obs` (its standard error
  comparison is now `mean_se_ratio`);
  [`ss_aipe_rmsea_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea_sensitivity.md)
  retires `rmsea_pop`, `desired_width`, `mean_width`, and, most
  importantly, an output row named `assurance` that actually reported
  the realized width-attainment proportion, now `pct_ci_less_w` so the
  name no longer collides with the planning input;
  [`ss_aipe_sem_path_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md)
  retires `width_less_than_desired` and the `type_I_err*` trio.
- Four members scaled their non-coverage rows by 100 while the rest of
  the family reported proportions:
  [`ss_aipe_sm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm_sensitivity.md),
  [`ss_aipe_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md),
  [`ss_aipe_sc_ancova_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)
  (both divisor branches), and
  [`ss_aipe_c_ancova_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova_sensitivity.md)
  (which also scaled its width-attainment row). All are now proportions
  on the 0 to 1 scale, completing the sweep that earlier fixed
  [`ss_aipe_sc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md),
  [`ss_aipe_cv_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
  and
  [`ss_aipe_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md).
- The family contract test caught a real defect while it was being
  written: the `divisor = "s_anova"` branch of
  [`ss_aipe_sc_ancova_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)
  read the confidence limits out of
  [`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
  by row position (`[2, 2]` and `[4, 2]`) when that function returns
  three rows, so the branch treated the point estimate as the lower
  limit and NA as the upper: its realized widths were NA and its “Type I
  error” rows were nonsense. The limits are now read by term name, as
  the `s_ancova` branch always did, and the branch’s summaries are
  meaningful for the first time.
- Members that reported only width summaries now also report their
  estimand: `mean_psi` and siblings in the three contrast members,
  `mean_sm`, `mean_smd`, `mean_rmsea`, and `mean_path`, computed from
  the same replications the widths come from. The echo rows the schema
  calls for were added where missing
  ([`ss_aipe_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
  [`ss_aipe_cv_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
  [`ss_aipe_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
  and its `rc` / `src` wrappers, the contrast members, and
  [`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md),
  which now echoes `half_width`, `true_psi`, `true_sigma`,
  `alpha_level`, and `m0`).
- The power siblings joined the sweep:
  [`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md)
  and
  [`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md)
  echo their planning inputs (`p`, `true_R2` / `true_b_j`,
  `estimated_R2` / `estimated_b_j`, `desired_power`, NA when a size was
  specified directly, and `alpha_level`), and their realized-`R^2` rows
  follow the meaningful-capital rule (`mean_R2`, and in the omnibus
  member `mean_F` / `F_crit`). `tidy()` and `glance()` for the
  `dmar_ss_power_sensitivity` class are unchanged; the echoes ride along
  in `glance()`.

### The ICC Sample Size Planner Now Plans the Average-of-k Forms

- [`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
  accepted and documented the six Shrout-Fleiss `type` labels but
  planned every one of them on the single-rater scale, so average-of-k
  plans were materially undersized. The planning value and target width
  are now interpreted on the scale of the requested form: an
  average-of-k value is mapped to the single-rater scale through the
  inverse Spearman-Brown relation and each candidate confidence limit is
  mapped back (the convention
  [`var_icc()`](https://yelleknek.github.io/DMAR/reference/var_icc.md)
  uses), and `type` is now validated, so an unrecognized label fails
  instead of silently planning `ICC(1,1)`. Average-of-k recommendations
  move: at `rho = .70`, `k = 3`, and `width = .20`, planning for
  `ICC(1,k)` now recommends n = 110 (realized mean width .20 across
  2,000 Monte Carlo replications of the *F*-based interval) where the
  undersized plan recommended n = 69 (realized width .26); at
  `width = .10` the recommendation is n = 421, and at `rho = .90`,
  `width = .10` it is n = 52. Single-rater plans are unchanged. The
  planned form travels as the `icc_type` attribute on the returned
  table.
- [`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md)
  had the matching defect on the generator side: it treated `true_rho`
  as the single-rater ICC no matter the `type`. For the average-of-k
  forms it now simulates data whose population ICC at the average-of-k
  level equals `true_rho`, so the realized estimates, widths, and
  coverage refer to the form being planned; `type` is validated the same
  way.
- The claim in
  [`?ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
  that the assurance correction over-recommends sample size by 25 to 40
  subjects was backwards. At the page’s own condition (`rho = .70`,
  `k = 3`, `width = .20`, `assurance = .80`) the recommended n = 79
  delivers an empirical assurance of about .77, and n = 81 is the
  smallest sample size that reaches .80 (10,000 Monte Carlo
  replications). The Details section now says so and points to
  [`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md)
  for checking a strict assurance target.

### The Noncentral F Clamp Warning Speaks the Caller’s Language

- When an observed *F* falls below the `alpha_lower` critical value of
  the central *F*-distribution, the lower noncentrality limit is 0 by
  construction; that is a normal consequence of a small observed effect,
  not a failure.
  [`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
  still warns once in that case, but the message now states the
  consequence for the interval (the lower confidence limit is 0) instead
  of describing achieved tail probabilities, and every function that
  builds its interval by inverting the noncentral *F* through it
  ([`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
  [`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
  [`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  [`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) with
  fixed predictors,
  [`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  [`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
  [`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
  [`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)
  with the parametric method, and
  [`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md))
  restates it for its own effect size, for example “the lower confidence
  limit on the signal-to-noise ratio is 0”, at most once per call.
  Previously
  [`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md) and
  [`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md)
  surfaced the inner wording, which pointed users to a `prob_greater`
  column those functions do not return. The warning carries the
  condition class `dmar_ncf_clamp`, which the iterative callers
  ([`ss_aipe_R2()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
  [`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
  [`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
  [`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md))
  now match by class when deduplicating.
- A genuine failure of the inner root finding is no longer a bare
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) message: the error
  now names the function that was called, reports the *F*-statistic,
  degrees of freedom, and tail probabilities involved, and says what to
  try.

### design_consequences() Now Uses the Exact Noncentral *t* When df Is Finite

- With finite `df`, the significance lens of
  [`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
  (`power`, `type_s_error`, `exaggeration_ratio`) is now computed from
  the noncentral *t* distribution of the test statistic, the exact
  distribution when the standard error is estimated from the data and
  the same sampling model the precision lens already used; the
  exaggeration ratio integrates the truncated normal moments over the
  chi distribution of the estimated standard error. The previous code
  evaluated a location-shifted central *t*, the known-se approximation
  behind Gelman and Carlin’s `retrodesign()`, so the help page promised
  the noncentral *t* while the code delivered the approximation.
  Small-df results move: at *n* = 5 per group with a true effect of *d*
  = 1, power is now 0.2863 (was 0.2469; a two-million replication Monte
  Carlo of the design gives 0.2860 with simulation standard error
  0.0003, and base R’s `power.t.test(strict = TRUE)` agrees exactly),
  the Type S error is 0.00129 (was 0.00937), and the exaggeration ratio
  is 1.658 (was 1.915). The differences fade as `df` grows (negligible
  by about 60 per group), and the `df = Inf` normal case is unchanged.

### The Robust Standardized Mean Difference Scales Correctly at Every Trim

- [`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)
  applied the reciprocal of the Algina-Keselman-Penfield scaling
  constant at every trimming proportion other than the 0.20 default: the
  internal constant computed the Winsorized variance of a standard
  normal correctly and then returned `1 / sqrt(win_var)` where the
  definition calls for `sqrt(win_var)`, and a hard-coded 0.642 at
  `trim = 0.20` masked the error at the default while making the
  estimate discontinuous there (`trim = 0.1999` returned 2.43 times the
  value `trim = 0.20` returned). The constant is now `sqrt(win_var)` at
  every trim and the hard-coded branch is gone. At the default the
  constant moves from the rounded 0.642 to its exact value 0.6419398, so
  estimates at `trim = 0.20` change by less than one part in ten
  thousand; at any other trim the correction is substantial (previous
  estimates were 1.47 times too large at `trim = 0.10`, 4.99 times too
  large at `trim = 0.30`, and 17.99 times too large at `trim = 0.40`).
  The help page identity is corrected to match: 0.642 is SD(X_W) /
  SD(X), that is sqrt(Var(X_W) / Var(X)), for a standard normal
  Winsorized at 0.20, not Var(X) / Var(X_W).

### The Robust Standardized Mean Difference Interval Now Uses the Yuen-Welch Degrees of Freedom

- [`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)’s
  help page promised a confidence interval on the Yuen-Welch degrees of
  freedom, but the code inverted the noncentral *t* at `h_1 + h_2 - 2`
  degrees of freedom with noncentrality
  `d_R * sqrt(h_1 * h_2 / (h_1 + h_2))`; the Yuen-Welch value was
  computed, returned in the `df_yuen` row, and never used. The interval
  now follows the construction in Keselman, Algina, Lix, Wilcox, and
  Deering (2008): Yuen’s *t*-statistic on the trimmed-mean difference
  (their Equation 8) is referred to a noncentral *t* distribution with
  the Yuen-Welch approximate degrees of freedom (their Equation 9), and
  the noncentrality limits are rescaled to the `d_R` metric. Feeding the
  summary statistics of the paper’s worked example (their Tables 1
  and 3) through this construction reproduces the printed robust effect
  size intervals \[0.31, 3.37\] and \[0.10, 1.11\] (p. 119) to the
  precision the rounded published inputs support, while `h_1 + h_2 - 2`
  degrees of freedom give \[0.40, 3.31\] for the first: the paper’s
  intervals use the Yuen-Welch value. The former noncentrality also
  understated the observed statistic (by about 16% at the default trim
  under equal Winsorized variances), so intervals move even where the
  two degrees of freedom nearly agree: on the help page example the
  interval is now \[-0.88, 0.08\] where the old construction gave
  \[-0.97, 0.17\].

### A Consistent Confidence Interval for Lin’s CCC

- [`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md)
  now builds its confidence interval on Lin’s (1989) *z*-transformed
  standard error, with the correction noted in Lin (2000), and the
  former default `method = "king_chinchilli"` is removed. The removed
  variance was not a consistent estimator of the sampling variance of
  the *z*-transformed CCC: when the two means and variances are equal it
  inflates the correct asymptotic variance by exactly the square of (1 +
  rho^2) / (1 - rho^2), already about 2.8-fold at a CCC of .5 and
  several hundred-fold at .95 (and by a comparable factor otherwise), so
  its intervals were nearly vacuous. The help page’s own first example
  returned CCC = .928 with CI \[-0.977, 1.000\]; the same example now
  returns \[0.871, 0.960\], matching the independent `DescTools::CCC()`
  *z*-transform interval to ten decimals, and the Lin interval’s
  simulated coverage is .952 at nominal .95 (bivariate normal, rho = .5,
  *n* = 50, where the removed default covered .982 with intervals 61%
  wider). The removed formula was also not King and Chinchilli’s (2001)
  estimator, so keeping it under that name would have credited a wrong
  formula to real authors; a correct King-Chinchilli variance may return
  in a later release. Calls that request `method = "king_chinchilli"`
  now fail loudly.

### The Fleiss Kappa Test Uses the Corrected Null Variance

- [`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
  computed its *z* statistic and *p*-value from a null variance that
  treats the category marginals as known constants, not from the
  estimated-marginals null variance of Fleiss, Nee, and Landis (1979),
  the correction of the standard errors in Fleiss (1971), even though
  the help page credited that paper. The two expressions coincide at
  uniform marginals, which is how the error hid, and diverge as the
  marginals skew: at *N* = 1000, *m* = 5, and marginals (.85, .10, .05),
  the old null standard deviation was 3.8 times the empirical one (.0303
  versus .0080). Because the known-marginals expression is never smaller
  than the corrected one (they are equal only at uniform marginals),
  every affected *z* was understated and every affected *p*-value
  overstated; the test was conservative, never anticonservative. On the
  help page’s own Fleiss
  1971. Table 1 example the *z* statistic moves from 15.64 to 17.65, now
        matching `irr::kappam.fleiss` to ten decimals. The point
        estimate and the confidence interval do not move: kappa is
        unchanged, and the `se`, `lower_limit`, and `upper_limit`
        columns come from the Gwet (2008) linearization, which was
        correct all along.

### One-Sided Dunnett Intervals Now Contain Their Point Estimates

- `ci_dunnett(alternative = "less")` formed its simultaneous upper
  bounds with the signed critical value that
  [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
  reports for that alternative, a negative lower-tail quantile, which
  placed every bound at `diff - |d| * se`: below the point estimate it
  was supposed to bound, and in contradiction with the adjusted
  *p*-values printed beside it (on the PlantGrowth example both “less”
  rows reported an upper limit below zero, a rejection, beside adjusted
  *p*-values of 0.162 and 0.989, no rejection). Simultaneous coverage of
  the true differences measured 0.00115 against the nominal .95 in a
  balanced null design with four groups and *n* = 10 per group. The
  critical value’s magnitude is now applied on the side the alternative
  dictates, so the “less” bound is `diff + |d| * se`; coverage in the
  same design measures 0.94685, and the PlantGrowth “less” upper limits
  move from -0.928 and -0.063 to 0.186 and 1.051, agreeing with
  `multcomp::glht()`’s one-sided limits to about four decimal places
  (the remaining daylight is multcomp’s simulated-quantile error) and
  with the adjusted *p*-values on every row. The
  `alternative = "greater"` interval already used the positive critical
  value and is unchanged, as are the two-sided intervals and all
  adjusted *p*-values.

### anova_within_two_way() Declines to Estimate an Inestimable Epsilon

- When an effect’s numerator degrees of freedom exceed n - 1, the
  covariance matrix of the effect’s orthonormal contrasts is singular
  and the Greenhouse-Geisser formula returns an artifact bounded above
  by (n - 1)/df regardless of the population epsilon (even under exact
  sphericity, where the population value is 1); the Huynh-Feldt value
  derived from it is equally meaningless.
  [`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)
  previously reported these artifacts as estimates. The
  Greenhouse-Geisser and Huynh-Feldt rows for such an effect now carry
  NA in epsilon, the adjusted degrees of freedom, and the *p*-value,
  with one warning naming the condition;
  [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html) likewise
  declines to report the corrections when the effect’s error matrix is
  singular. The unadjusted and lower-bound rows are unchanged: the lower
  bound 1/df is Geisser and Greenhouse’s a priori bound, not an
  estimate, and remains valid however few the subjects (Maxwell,
  Delaney, & Kelley, 2027, Chapters 11 and 13).

### ci_r() Now Requires at Least Four Observations

- [`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  accepted any sample size and went quietly wrong below n = 4: at n = 3
  it returned the vacuous interval \[-1, 1\] (the Fisher’s Z variance
  1/(n - 3) is infinite there), and at n = 2 it returned NaN limits. The
  function now stops with an informative error for n \< 4, and the n = 3
  message explains why no usable interval exists at that boundary.
- The same unguarded 1/(N - 3) sat behind `ci_method = "fisher"` in
  [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  and
  [`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md);
  the Fisher interval there now stops for N \< 4 with the same
  explanation. The package’s other Fisher’s Z consumers
  ([`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
  [`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
  [`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
  [`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
  [`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
  [`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
  [`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md))
  already validated their sample sizes.

### The AIPE Sensitivity Family Reports Coverage as Proportions

- Three members of the `ss_aipe_*_sensitivity` family scaled some or all
  of their coverage and width-attainment rows by 100 while the rest of
  the family reported proportions.
  [`ss_aipe_sc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md)
  reported `type_I_error_upper` and `type_I_error_lower` as percentages;
  [`ss_aipe_cv_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)
  reported `pct_ci_less_w`, `pct_ci_miss_low`, and `pct_ci_miss_high` as
  percentages beside a `total_type_I_error` that was already a
  proportion, so its own rows did not add up; and
  [`ss_aipe_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
  (inherited by
  [`ss_aipe_rc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc_sensitivity.md)
  and
  [`ss_aipe_src_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md))
  scaled all four. These rows are now proportions on the 0 to 1 scale
  everywhere, so `total_type_I_error` equals the sum of
  `pct_ci_miss_low` and `pct_ci_miss_high`, and a realized Type I error
  compares directly to 1 - `conf_level` with no per-function rescaling.
  Term names are unchanged; only the scale moved.

### The MBCO Likelihood Ratio Is Now Branch-Deterministic

- [`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)’s
  null hypothesis for an indirect effect is a union of branches (the
  product is zero when any factor is), and the constrained search could
  converge to a worse-fitting branch on some platforms and OpenMx
  builds, inflating the likelihood ratio statistic. When every
  constrained algebra is a pure product of free parameters, the branches
  are now also fit directly as ordinary unconstrained models with one
  factor fixed to zero, and the reported statistic is defined by the
  best-fitting branch on every platform. Non-product constraints keep
  the multi-start constrained search.

### tidy() and glance() Speak DMAR’s Names, on Every Table

- **The broom verbs now return DMAR’s native column names**: `p_value`,
  `se`, `ci_lower`, `ci_upper`, `p_adjusted`, `conf_level`,
  `std_estimate`, `R2`, `adj_R2`, `df_residual`, and so on, replacing
  broom’s dotted `p.value` / `std.error` / `conf.low` across every
  tidier in the package. One naming system now covers every DMAR
  surface; a pipeline that feeds a tool expecting broom’s dotted schema
  renames the columns at that boundary.
- **The wide tables answer the verbs too.**
  [`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
  [`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
  [`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  and
  [`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md)
  gained bespoke `tidy()`/`glance()` pairs (items, ladder rungs, and
  groups as terms; the scale and model level summaries as the one-row
  glance), and a generic wide branch in the default methods serves
  [`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
  [`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
  and any future wide table (term from the label columns, estimate from
  the first numeric column, remaining columns passed through). The
  output vignette’s promise that the verbs answer everywhere is now
  literally true.

### One Bootstrap Vocabulary, Three New Bootstrap Intervals

- **`B` is the number of bootstrap replications everywhere.**
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)
  rename `boot_R`,
  [`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md)
  and
  [`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)
  rename `R`, and
  [`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md)
  renames `n_boot`; the effective-count row is now `B_used`. Old
  argument names fail loudly.
- **[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
  and
  [`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
  gain bootstrap intervals** (`ci_method = "percentile"` or `"bca"`,
  `B = 10000`, `seed`), redeeming the kappa page’s own advice that the
  Wald interval can have poor small-sample coverage (Blackman & Koval,
  2000; Zapf, Castell, Morawietz, & Karch, 2016). Subjects are
  resampled, `table` input is expanded to the equivalent pairs, and the
  `se`, `z_value`, and `p_value` columns keep their asymptotic
  definitions; only the interval changes.
- **[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
  gains a percentile bootstrap** (`ci_method = "percentile"`,
  `B = 1000`, refitting the model per replication) and now always
  returns `ci_lower` / `ci_upper` columns (`NA` under
  `ci_method = "none"`). No interval is possible from `loadings` alone,
  and the function says so.
- **[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md)
  no longer bootstraps by default** (`boot = FALSE`), matching the
  package-wide rule that no analysis runs a bootstrap unless asked; its
  verbal reliability benchmarks were removed in favor of reporting the
  coefficient with its interval.
- **Failed bootstrap replications are dropped, never fatal**:
  [`mediate()`](https://yelleknek.github.io/DMAR/reference/mediate.md)
  no longer errors when a degenerate resample returns no indirect
  effect, and
  [`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md) no
  longer aborts when a resample leaves a block’s average
  within-construct correlation nonpositive. Both drop the replication,
  warn once with the count, and stop only when fewer than 100
  replications survive.
- **No help page wraps its examples in `\donttest{}` or `\dontrun{}` any
  more.** Under `R CMD check --as-cran`, one `\donttest{}` block
  anywhere makes R run the whole example corpus twice, so the wrapper
  cost time rather than saving it. The examples that were slow are now
  fast. Replication counts in the cheap demonstrations were lowered,
  with a comment naming what a reported analysis deserves. **No example
  runs a bootstrap confidence interval**: those calls are carried as
  commented code, so the syntax is still on the page for a reader who
  wants it, and the surrounding prose explains the interval and when to
  ask for it. The randomization tests keep their resampling, since
  permuting the data is what those functions do and their examples cost
  a tenth of a second. Anything else expensive is commented out the same
  way, each passage introduced by a sentence saying what it does and why
  it is not run. The data behind every example are unchanged. The
  slowest help page now takes 0.38 seconds, and all 302 together take 8
  seconds.

### The CFA Family: One General Function, Two Convenience Wrappers

- **[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) is
  now the one factor special case of
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)**, a
  convenience wrapper rather than an independent implementation, and the
  new
  **[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md)**
  is the two factor sibling.
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) only
  requires the data (and, when the data hold more than the items, a
  vector of item names);
  [`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) takes
  the items of each factor as `factor_1` and `factor_2`. Both forward
  everything to
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md), so
  there is now a single fitting implementation, one output schema (with
  confidence interval columns), and one set of fit index choices; under
  a robust estimator every CFA surface now reports the robust index
  versions. The former
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
  extras moved or dissolved: composite reliability with the observed
  denominator is `reliability_omega(denominator = "observed")`, the
  `omega` output mode is the `omega_f1` row of the standard table, and
  the term names follow the
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)
  convention (`lambda_f1_<item>`, `psi_f1_<item>`, `omega_f1`).
- **`data` and `S` are now separate arguments across the CFA family.**
  Raw data are passed as `data` and a covariance matrix as `S` (with
  `N`); an argument never means both.
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)
  refuses a square symmetric matrix passed as `data`, naming the fix.
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) keeps
  its conveniences: unnamed input is auto-named `y1`, `y2`, …, and
  `items` defaults to every column.

### Measurement Invariance Beyond the Exact-Invariance Ladder

- **[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  now fits any measurement model, not only the one-factor case.** It
  accepts `model` as lavaan syntax or as a named list mapping factors to
  their items (the `items` argument remains as the one-factor
  convenience), along with `ordered`, `missing` for full information
  maximum likelihood, `group_partial` for the partial invariance case,
  and `parameterization`. With ordered indicators the ladder itself
  changes: thresholds carry the location information, so a
  **thresholds** rung is fitted between configural and metric, following
  Wu and Estabrook (2016) and Millsap and Yun-Tein (2004). Constraining
  loadings before thresholds, as the continuous ladder does, tests the
  wrong hypothesis for ordered items. The chi square difference test is
  the scaled one whenever the estimator is robust or the data are
  ordered, and the fitted lavaan objects are returned on a `"fits"`
  attribute so a caller can run score tests and partial invariance
  refits without paying for the ladder twice.

- **[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md)**
  reports the dMACS effect size of measurement noninvariance (Nye &
  Drasgow, 2011): the expected difference between two groups’
  measurement equations for an item, integrated over the focal group’s
  latent distribution and standardized by the pooled item standard
  deviation. A score test says a loading or intercept differs
  detectably; dMACS says whether the difference is large enough to
  matter. The defining integral has a closed form under a normal latent
  variable, so no numerical integration is used, and the two agree to
  machine precision in the tests. Accepts a fitted multiple group lavaan
  model or the parameters a paper reports. The help page states plainly
  that the index is meaningless from a configural fit, where the groups
  share no metric.

- **[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md)**
  implements the alignment method of Asparouhov and Muthen (2014). Exact
  invariance essentially never holds across many groups, which leaves
  the ladder stalled at configural and group comparison blocked.
  Alignment estimates the group factor means and variances that make the
  measurement parameters as nearly invariant as possible, minimizing a
  simplicity function whose fourth-root component loss tolerates a few
  large differences and punishes many small ones. Both the fixed and
  free identifications are available, the optimizer runs from several
  starts because the surface has local minima, and the number of
  distinct optima found is reported so a solution is never presented as
  unique when it is not.

### Item Response Theory, as the Categorical Factor Model

- **[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md)**
  fits Samejima’s (1969) graded response model. It does so the way the
  rest of the package works, by fitting the categorical factor analysis
  model with lavaan and converting the solution: discrimination
  `a_i = lambda_i / sqrt(1 - lambda_i^2)` and location
  `b_ik = tau_ik / lambda_i`. The two models are the same model in
  different parameterizations (Takane & de Leeuw, 1987), so this adds
  item response theory without a second estimation engine and without
  leaving the factor analytic tradition. Both the normal ogive and the
  logistic metric are available.

- **[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md)**
  and
  **[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md)**
  give the item and test information functions and the standard error of
  the latent trait, `SE(theta) = 1 / sqrt(I(theta))`. This is precision
  as a function of where the respondent sits on the trait, which is what
  a single reliability coefficient cannot express: a scale can have
  excellent omega and still measure poorly over the range a study cares
  about. Verified against the closed form for the dichotomous case,
  which the graded model must reproduce exactly.

### Content Validity

- **[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md)**
  computes the item level content validity index from a panel of expert
  relevance ratings, with the modified kappa that corrects it for chance
  agreement (Polit, Beck, & Owen, 2007), Lawshe’s (1975) content
  validity ratio, and the scale level S-CVI/Ave and S-CVI/UA. Each I-CVI
  carries an exact binomial confidence interval, because an index
  computed from five or six experts is a proportion with real
  uncertainty and reporting it as a bare point estimate overstates what
  a small panel establishes. This is evidence about the items, gathered
  before any data are collected.

### Mediation inference via model-based constrained optimization

- **[`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)**
  implements the model-based constrained optimization (MBCO) procedure
  of Tofighi and Kelley (2020, *Psychological Methods*): a likelihood
  ratio test of any smooth function of path coefficients (an indirect
  effect, a total effect, a contrast of two indirect effects) formed by
  refitting the mediation model subject to the nonlinear constraint that
  the function equals zero. The model is specified in lavaan syntax
  (observed or latent variables; parallel or sequential mediators) and
  fit in OpenMx, whose optimizers support the nonlinear equality
  constraints the null model requires. The function enumerates the
  total, direct, total indirect, and every specific indirect pathway
  from `x` to `y`, reports each with a delta method standard error and a
  profile likelihood, Monte Carlo, or Wald confidence interval, tests
  each with the MBCO likelihood ratio statistic and its *p*-value, and
  reports AIC and BIC differences and the change in each endogenous
  R-squared under every null model. Because the null set of a product
  constraint is a union of surfaces, each null model is refit from
  several starting configurations and the best feasible solution is
  kept, so the reported statistic reflects the globally best-fitting
  null model rather than a local branch. Accepts raw data or the summary
  statistics a paper reports (covariance matrix, means, and sample
  size), which reproduce the raw-data analysis exactly; the tests
  replicate the published empirical example from its Table 1 moments. A
  `group` argument turns the model into a multiple-group SEM in which
  every effect is estimated per group and its between-group difference
  is tested, which is moderated mediation with a categorical moderator.
  A `moderator` argument probes continuous moderation stated in the
  syntax with `:` interaction terms (or precomputed product columns):
  each moderated pathway effect becomes a symbolically derived
  polynomial in the moderator, reported as conditional effects at probe
  values (the mean and one standard deviation either side by default),
  the index of moderated mediation (Hayes, 2015) tested by likelihood
  ratio rather than bootstrap, and, when a pathway is moderated in
  several places, a joint constancy test whose null model imposes
  several nonlinear constraints at once. `hypotheses` accepts a named
  list whose multi-expression elements are likewise tested jointly.
  Constrained null models are started on every branch of the
  constraint’s null set (including exact conditional-coefficient branch
  starts for probed effects) and each solution is polished by a warm
  restart, so the reported statistic reflects the best-fitting feasible
  null model. Structural guardrails refuse pathway enumeration for
  nonrecursive (feedback) structures and warn on binary endogenous
  variables, on interaction terms no declared moderator accounts for,
  and on interactions missing their main effect (the principle of
  marginality).
- **[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)**
  draws the conditional effects of a moderated
  [`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
  analysis as curves over the moderator’s observed range, one per
  moderated pathway effect, with a pointwise Monte Carlo confidence
  band, the probed values marked, a dashed zero line whose band
  crossings estimate the Johnson-Neyman boundaries, and a rug of the
  observed moderator values so extrapolation is visible. The band’s
  polynomial is the same stored quantity the table probes, evaluated
  over the whole range; the help page is explicit that the band is
  pointwise and that boundary locations read from it are estimates, with
  the table’s constancy test as the formal companion. Returns a plain
  ggplot object in the Okabe-Ito palette.

### API and Documentation Consistency

- **[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md)**
  computes the variance of the estimated treatment effect at selected
  covariate values in a two-group ANCOVA with heterogeneity of
  regression and a random covariate (Li, McLouth, & Delaney, 2020), the
  reimplementation of `MBESS::var.ete()`; tested against the MBESS
  reference on every branch and dogfooded in the Pygmalion vignette.

- **`tidy()` and `glance()` answer on every DMAR result.** Two default
  methods on the `dmar_tbl` class are the floor under the package: any
  result table answers `tidy()` with the broom-shaped long view (`term`,
  `estimate`) and `glance()` with the one-row wide view, repeated terms
  disambiguated rather than dropped. The family methods sit ahead of the
  defaults and keep their richer views, and two families join them: the
  `ss_aipe_*` planners carry a new `dmar_ss_aipe` class whose `tidy()`
  reports the planned size beside the desired width it was planned
  against (`term`, `estimate`, `width`), through the same size registry
  the power planners use and a width registry beside it, and whose
  `glance()` widens the echoed inputs.
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
  the one power planner missing its family class, joins `dmar_ss_power`,
  so its `tidy()`/`glance()` no longer error, and
  [`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)
  joins its siblings in `dmar_ci_anova`.

- **One sample size vocabulary across every planner.** A planner’s
  answer row is now named for what it counts: `necessary_N` (total),
  `necessary_n_per_group`, `necessary_n_per_cell`, or
  `necessary_n_clusters`; a user-fixed size echoes as `specified_*`; and
  `total_N` appears only as the implied-total companion beside a
  per-unit answer, never as an answer’s name. Before this sweep the same
  word meant different things in different planners: `sample_size` was
  the total in
  [`ss_aipe_R2()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)
  and its seven siblings but the per-group size in
  [`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
  and its three,
  [`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
  said `sample_size_per_group`,
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
  carried the MBESS-era `ss_c` and `ss_t` pair (two rows for one number
  in a balanced design, now one branch-named row),
  [`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md)
  used one bare `n_per_cell` for both its planned and its user-supplied
  branch, and the cluster planners’ answers had no prefix at all.
  Sensitivity outputs echo the evaluated size under its bare unit name
  (`total_N`, `n_per_group`), and the six sensitivity functions whose
  `specified_N` argument actually meant a per-group size now call it
  `n_per_group`
  (`ss_aipe_c/sc/sc_ancova/smd/equivalence_smd/pcm_sensitivity`). Two
  term names in the cluster family that contained literal spaces are
  underscored. The tidy()/glance() size recognizer no longer accepts the
  legacy names, so a stray old producer fails a test instead of slipping
  through, and a 216-table characterization grid recorded before the
  sweep reproduces identically after it, so only names moved, never
  values. The package is unreleased; old names fail loudly.

- **The number of predictors is `p` everywhere.**
  [`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
  [`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
  [`ss_aipe_rc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md),
  [`ss_aipe_src()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md),
  [`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
  and
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)
  renamed `J` to `p`, and
  [`ci_R()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  renamed `K` to `p`, matching
  [`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
  [`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
  and the rest of the regression family. `J` remains only in the partial
  and semipartial correlation family, where it counts the variables
  partialed out (a different quantity; `J = 0` is the simple
  correlation). The package is unreleased, so the old names fail loudly
  rather than being aliased.

- **[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) and
  [`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md) are
  documented as the thin wrappers they are** around
  [`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
  the general engine, with titles that distinguish the unstandardized,
  standardized, and general cases.

- **Every function help-page title is now AP title case** with no
  trailing period, matching base R convention (previously the package
  mixed sentence case and title case).

### Ten Help Pages Now State Exactly What the Code Computes

A numerical audit compared every displayed formula and quoted magnitude
on these pages against the quantity the function computes. In each case
the computation was correct and the page was not, so no returned value
changes; the documentation now matches the code, and new tests recompute
each corrected expression independently and assert agreement with the
function output.

- **[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md)**
  displayed the expectation with every index one lower than the
  Olkin-Pratt formula under the n to n - J substitution the code makes;
  as printed, the formula can exceed 1, an impossible value for the
  expectation of a correlation. The page now prints E\[r\] = rho
  2F1(1/2, 1/2; (n - J + 1)/2; rho^2) Gamma((n - J)/2)^2 / (Gamma((n -
  J - 1)/2) Gamma((n - J + 1)/2)), which matches the function to machine
  precision. The bias prose also described E\[r\] - rho while the
  returned `bias` column is rho - E\[r\]; the sign convention now
  matches the column, and the quoted magnitudes are the true +0.00624
  (rho = 0.4, n = 30, J = 2) and +0.00888 (J = 10).
- **[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md)**
  quoted bias magnitudes with the wrong sign and size (-0.024 and
  -0.006); the true values are +0.021 (rho = 0.5, n = 10) and +0.0065 (n
  = 30), as the page’s own example already reported.
- **[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md)**
  printed the variance of the Bonett (2002) L-transform as 1 / (2 (n -
  2)), omitting the k / (k - 1) factor the code applies; the page now
  prints k / (2 (k - 1) (n - 2)). A reader hand-building the interval
  from the old page at k = 2 would have an SE too small by a factor of
  sqrt(2) and roughly .84 coverage instead of .95.
- **[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)**
  printed the lower-bound epsilon as 1 / (df - 1); the code computes 1 /
  df, the attainable infimum, with df the effect’s numerator degrees of
  freedom.
- **[`loa()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md)**
  printed a symmetric central t interval for the CIs on the limits of
  agreement; the code computes the Carkeet (2015) exact noncentral t
  interval, which is asymmetric about the sample LoA. The page now shows
  the noncentral form: quantiles of the noncentral t with noncentrality
  parameter k sqrt(n) for the upper LoA and its negative for the lower.
- **[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)**
  printed a symmetric interval with unsquared constants and a -MS_b MS_w
  / n cross term that appears nowhere in the method; the code computes
  the genuine Burdick-Graybill (1992, equations 2.4.1–2.4.5) form, and
  the page now prints V_L and V_U with the squared constants and the
  G_12 / H_12 cross terms.
- **[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md)**
  printed Var(beta-hat) with a trailing design-effect factor the code
  does not apply; for the cluster-mean-centered level-1 slope the
  variance is sigma2_y (1 - rho_I) / (N sigma2_x), and the page now
  explains why the design effect does not enter.
- **[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)**
  pointed paired designs to
  [`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
  which is the interval for Glass’s estimator (two independent groups,
  control group SD) and takes no correlation between paired
  measurements, so following the pointer reproduced essentially the
  independent groups interval. The page now states plainly that a
  paired-design SMD interval is not currently provided.
- **[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md)**
  attributed the alternative distribution to Wallenius; the code
  computes Fisher’s noncentral hypergeometric, the conditional
  distribution of one binomial count given the total of two independent
  binomials, which is what conditioning on the margins of the 2 x 2
  table produces. The prose and references now cite Fisher
  1935. and Fog’s (2008) sampling-methods paper.
- **[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md)**
  said the default returns a - 1 pairwise contrasts; it returns all a
  (a - 1) / 2 of them, as the sibling pages already stated.

### A Second Documentation Pass: Quoted Values Now Match What Runs

A continuation of the audit above, covering another two dozen help
pages. As before, the computations were correct and the pages were not;
no returned value changes, and where a corrected number is load-bearing
a test now recomputes it independently.

- **Quoted numbers now match what the examples print.**
  [`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)’s
  90 percent example described its upper limit as landing below 0.05;
  the limit is 0.052, just above the Browne and Cudeck close fit
  threshold, and the page now draws the conclusion that follows (close
  fit is not established, even though the point estimate sits below the
  threshold).
  [`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md)’s
  Hunter example quoted p = .2595, which matches neither the quoted F =
  1.29 nor the exact statistic; the page now quotes p = .260, as
  computed.
  [`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md)’s
  Edgington sum for the Raudenbush example prints 7, not “near 6.9”.
  [`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md)’s
  hypertension example described 24 subjects across group sizes of 4, 6,
  5, and 5; it says 20.
- **[`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
  separates the published memory-example likelihood ratios from what its
  example computes.** Run from the rounded Table 1 summary statistics,
  the two null branches give LRT = 71.31 (the best-fitting branch, the
  statistic reported) and 179.02; the 72.54 and 175.77 the page
  previously quoted are the full-precision-moments values, and 175.77 is
  now attributed to Tofighi and Kelley (2020) as the published value,
  from the worse-fitting branch on which their optimizer stopped. The
  page notes that the difference comes from running from the published
  (rounded) summary statistics.
- **[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)’s
  example no longer equates the approximation with the simulation
  benchmark.** The joint significance approximation returns
  `necessary_N` = 65 for a = b = .39 at power .80; raw-data simulation
  puts the power at 65 nearer .77 and reaches .80 near N = 70, in line
  with the somewhat larger requirement in Fritz and MacKinnon’s (2007)
  simulation-based table, and the example comment now says exactly that.
  The raw-data validation test’s relative tolerance is widened from 0.03
  to 0.06 with the reason recorded in place: the approximation sits
  about 3 percent above the simulation at those settings on any seed, so
  the old tolerance failed on some seeds for the approximation gap
  alone.
- **`prime_time_achievement`’s corporation summaries are recomputed from
  the derived `corp_id` key.** The bare `corp` column merges the two
  corporations that share code 2400, so the page overstated the largest
  corporation: corporations have 15 to 756 students (median 117), not 15
  to 808 (median 124), and the three-level null model decomposition is
  16.29 / 22.72 / 240.43 with a corporation ICC of 0.058 (previously
  16.50 / 22.70 / 240.43 and 0.059). The data set’s tests now anchor the
  corrected values.
- **`bessel_errors`’ expected counts imply a normal sigma near 0.22, not
  0.2.** A least squares fit of the half-normal bin probabilities to
  Bessel’s expected frequencies gives 0.216; at 0.2 the first bin alone
  is off by about eight observations.
- **The `which_width` contract of
  [`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
  [`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
  and
  [`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
  now says what the code does.** The pages described `"Lower"` /
  `"Upper"` as a one-sided half-width; both settings interpret `width`
  as half the full width and return the same sample size. The pages now
  say so, note that an asymmetric interval’s realized half-widths differ
  from each other and from half the full width, and state that a
  genuinely one-sided width target is not currently offered.
- **[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md)
  documents the tail arguments its backend accepts.** The page promised
  that supplying `alpha_lower` and `alpha_upper` recomputes
  `conf_level`; that call errors, because the tails pass straight
  through to
  [`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
  which refuses a non-NULL `conf_level` beside them. The page now
  directs users to set `conf_level = NULL` and supply both alphas, and a
  test pins the contract.
- **[`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
  and
  [`conf_limits_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)
  state the monotonicity each search actually exploits.** The pages
  claimed each tail probability is strictly decreasing in the
  noncentrality parameter; the lower-limit condition works on the upper
  tail, which is strictly increasing. Each page now gives the direction
  per tail and notes that both roots are located on the decreasing
  lower-tail scale.
- **Smaller contract corrections.**
  [`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
  described an expected half-width target of omega where the code
  targets omega / 2 (full width omega), and a 1-row return where the
  table has 5 rows; both now match, and the return’s stale term names
  are corrected. The stale `sample_size` return term on the
  [`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
  [`ss_aipe_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md),
  and
  [`ss_aipe_sm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md)
  pages now reads `necessary_n_per_group`, `necessary_n_per_group`, and
  `necessary_N`.
  [`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)’s
  `N` is documented as the total number of observations (for `aovlist`
  fits, one more than the sum of all effect and residual degrees of
  freedom), matching its Details.
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md)’s
  `adj_R2` is documented as using the complete cases (`N_complete`),
  which the code deliberately uses, rather than the lavaan reported N.
  [`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)’s
  `contrast_type` entry dropped inline width expressions that omitted
  the 1/sqrt(2) and sqrt(1/n) factors and now points to the correct
  display in Details.
  [`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md)’s
  note no longer contains a sentence directing users to pass the error
  variance where the argument (and the rest of the page) wants the error
  standard deviation. The multiplicative `diff_size` illustration on the
  `ss_aipe_crd_*` pages was a verbatim copy of the additive one; it now
  works the multiplicative form’s arithmetic (cluster size 25 with
  `diff_size = c(0.8, 1, 1.2)` gives sizes 20, 25, 30, recycled across
  clusters).

### Equivalence and noninferiority for linear contrasts

- **[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)**
  performs the two one-sided tests procedure and the companion
  noninferiority test for a linear contrast of group means against
  bounds stated in raw units of the response, with one pooled error
  term, a summary-statistic and a direct (estimate, SE, df) interface, a
  `benchmark` argument for comparisons against a known constant, and a
  five-way verdict (equivalent, superior, inferior, non-inferior only,
  inconclusive). Reproduces
  `emmeans::test(..., side = "equivalence" | "noninferiority")` to
  machine precision.
- **[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md)**
  computes the exact TOST power for a contrast by integration over the
  chi distribution of the estimated error standard deviation,
  generalizing
  [`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md)
  to arbitrary weights and unequal group sizes, and the noninferiority
  power in closed noncentral *t* form.
- **[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)**
  finds the smallest per-group sample size whose exact TOST (or
  noninferiority) power reaches a target, the declaration-probability
  counterpart of the width-targeting
  [`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md).
- **[`ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.md)**
  and
  **[`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)**
  implement the purely sequential fixed-width confidence interval
  procedure for a contrast of Chattopadhyay, Bandyopadhyay, Kelley, and
  Padalunkal (2025), with cost-optimal allocation across groups; the
  sensitivity sibling verifies first-order efficiency and near-nominal
  coverage by Monte Carlo.
- **[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md)**
  draws contrast estimates and their intervals against the equivalence
  region, colored by the
  [`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
  verdict.

### Mixed-effects R-squared

- **[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md)**
  reports the Nakagawa and Schielzeth (2013) marginal and conditional
  R-squared for a fitted `lme4` or `nlme` mixed model. It agrees to
  machine precision with `r2mlm` and with the direct Johnson quadratic
  form across the model classes the tests pin: random-intercept,
  correlated random-slope, split (`||`) random-slope, multi-slope, and
  nested two-level models. On split random-slope models
  `performance`/`insight` omits the random-slope variance from its
  marginal and conditional R-squared, so DMAR matches `r2mlm` and the
  Johnson form there and, by design, not `performance::r2_nakagawa`; the
  two agree on random-intercept and correlated-slope models.
- **[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md)**
  implements the Rights and Sterba
  2019. integrative framework of mixed-effects model R-squared measures:
        the full family of total, within-cluster, and between-cluster
        measures from a complete five-source decomposition of the
        outcome variance, matching the authors’ `r2mlm` reference
        implementation to machine precision across the same tested model
        classes.

### Output, tidiers, and contrasts

- **Publication-ready output helpers**:
  [`knit_print.dmar_tbl()`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  renders DMAR results as formatted tables in knitted documents,
  [`as_kable()`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  produces a [`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html)
  view, and
  [`results_sentence()`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  writes an APA-style “estimate, CI” sentence from any interval-carrying
  result.
- **`broom` support** (`tidy()` / `glance()`) added for the elementary
  tests
  ([`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md),
  [`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
  [`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md))
  and the simultaneous-comparison intervals
  ([`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
  [`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
  [`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md)).
- **The broom summary now spans the power-based sample size planners.**
  `tidy()` and `glance()` summarize every closed-form `ss_power_*`
  planner that reports one size and one power, from the effect size
  planners to the ANOVA, ANCOVA, contrast, cluster, and mediation
  designs, reporting the design’s planning unit (per group, per cell,
  per subject, or per cluster, or the total for the one-way ANOVA)
  beside its power. The lookup that resolves the size row learned the
  design-specific names (`n_per_cell`, `n_subjects`, `J_per_arm`,
  `ss_per_group`), so the meaningful names are kept rather than
  homogenized. The Monte Carlo sensitivity siblings
  [`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md)
  and
  [`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md)
  carry a `dmar_ss_power_sensitivity` class whose `tidy()` places the
  empirical and analytic power side by side.
- **[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md)**
  plans the per-group sample size for composite power in a two-group
  ANCOVA: the probability that the group effect, the covariate effect,
  and the group by covariate interaction are all significant in the same
  study, the quantity a design must be planned against when its
  conclusion needs more than one result to hold at once. Composite power
  is not the product of the marginal powers, because every test divides
  by the same error estimate and the tests are positively dependent even
  when the effects are orthogonal. A one dimensional integral over the
  chi square distribution of that estimate evaluates the composite
  deterministically, with no simulation, and the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draws
  the population effects the plan rests on. This opens the new composite
  power family.
- **[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md)
  and
  [`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md)**
  carry composite power to any balanced factorial design: name any set
  of main effects and interactions in `effects`, each with its Cohen’s
  *f* or partial eta squared, and the planner returns the per-cell
  sample size at which all of them are jointly significant. The effects
  are noncentral *F* tests sharing one error estimate, so the same
  shared-error integral applies; a single effect reproduces
  [`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md)
  (or
  [`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md)
  with a covariate) exactly, and the two-effect composite matches a
  direct simulation. The ANOVA function is the no-covariate case named
  directly and admits no covariate. With no covariate the composite is
  exact to quadrature precision. The effects can be stated as sizes
  (Cohen’s *f* or partial eta squared) or as a full array of population
  cell means with a common within-cell standard deviation, from which
  each effect’s *f* is read off the analysis of variance decomposition
  of the means. The
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draws
  the purported population values: the cell-mean pattern (with error
  bars of one within-cell SD) when means were given, or the effect sizes
  annotated with their marginal power otherwise, with the composite
  power in the subtitle either way.
- **[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md)**
  carries composite power to a factorial ANCOVA whose covariate slope
  differs across the cells. The average slope (the covariate main
  effect) and the factor by covariate slope heterogeneity are then
  testable effects that can join the composite alongside the factorial
  mean effects. The full model fits the means, the covariate, and every
  factor by covariate slope, so the residual has N minus twice the cells
  degrees of freedom; a “mean”, “covariate”, or “slope” effect is named
  in `effects` and sized by a Cohen’s *f* or read from population values
  (cell means and a covariate outcome correlation per cell, with a
  common within-cell SD). The one-factor two-level case reproduces
  [`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md)
  to machine precision, and a three-effect composite matches a direct
  simulation of the heterogeneous-slope model. Kept separate from the
  common-slope
  [`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md)
  for ease of use. The
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draws
  the population regression line in each cell, so heterogeneous slopes
  read as lines of different angle.
- **[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md)
  and
  [`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md)**
  are the general entry points to the composite power family, covering a
  one-way design with any number of groups as well as any factorial
  arrangement. The ANCOVA planner takes a `slopes` argument,
  `"homogeneous"` for one common covariate slope or `"heterogeneous"` to
  let the slope differ across cells and make the covariate and
  slope-heterogeneity effects testable; its heterogeneous one-way case
  is the *a*-group generalization of the two-group
  [`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
  and its two-level case reproduces it exactly. The ANOVA planner is the
  no-covariate design named directly, with a covariate-free interface
  that points to the ANCOVA version when a covariate is present. Both
  accept effect sizes or population values (cell means with a common
  within-cell SD) and forward to the factorial planners, so the same
  broom `tidy()`/`glance()` summaries and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) figures
  apply.
- **[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md)
  and
  [`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md)**
  carry composite sample size planning to structural equation models,
  for both design goals. The researcher states the population as a fully
  fixed lavaan model (or its
  [`cov_sem()`](https://yelleknek.github.io/DMAR/reference/cov_sem.md)
  covariance matrix), labels the parameters of interest in the free
  analysis model (structural paths, loadings, covariances, or
  `:=`-defined quantities such as an indirect effect), and the planner
  finds the smallest *N* at which every labeled parameter is
  statistically significant in the same study with the desired
  probability
  ([`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md)),
  or at which every confidence interval is sufficiently narrow, in
  expectation or with a stated assurance for the joint event
  ([`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md)).
  No closed form covers a set of dependent SEM estimates, so both
  planners run an a priori Monte Carlo simulation (Muthén & Muthén,
  2002; Maxwell, Kelley, & Rausch, 2008): data are drawn from the
  population covariance matrix, the analysis model is fit `G` times per
  candidate *N*, and a search seeded by the analytic Wald approximation
  brackets and bisects to the smallest integer meeting the goals; each
  reported power or proportion travels with its simulation standard
  error, and a `seed` argument makes a plan reproducible. The power
  planner joins the `dmar_ss_power` broom family; marginal
  `power_<label>` and `width_within_desired_<label>` rows show which
  parameter binds the design. Both planners handle a population mean
  structure, so latent growth curve targets such as the slope factor’s
  mean are planned the same way:
  [`cov_sem()`](https://yelleknek.github.io/DMAR/reference/cov_sem.md)
  now also returns `mu_theta`, the model implied means of the observed
  variables, a `mu` argument supplies means beside a hand-built `Sigma`,
  and the simulated data carry those means whenever the analysis model
  has a mean structure. The “Composite Sample Size Planning for SEM”
  vignette
  ([`vignette("composite_sem_planning")`](https://yelleknek.github.io/DMAR/articles/composite_sem_planning.md))
  walks through the full workflow for a mediation model with observed
  variables and for a linear latent growth curve, planning both
  composite power and joint accuracy.
- **[`cohen_h()`](https://yelleknek.github.io/DMAR/reference/cohen_h.md)**
  returns Cohen’s *h*, the effect size for the difference between two
  proportions on the arcsine (variance-stabilizing) scale, *h* = 2
  asin(sqrt(p1)) minus 2 asin(sqrt(p2)). It is signed and is the
  proportion analogue of the standardized mean difference
  ([`smd()`](https://yelleknek.github.io/DMAR/reference/smd.md)), so a
  given *h* carries the same detectability wherever the proportions sit,
  which a raw difference does not.
- **[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)**
  tests an arbitrary contrast among covariate-adjusted cell means in a
  factorial ANCOVA, matching `emmeans::contrast()`.
- More estimators now route through the tidy `dmar_tbl` display layer,
  and literature-synonym `@concept` tags (Cohen’s d, Cronbach’s alpha,
  Cohen’s U3, …) make the marquee estimators searchable by their
  eponyms.
- **[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md)
  and
  [`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
  echo their user-supplied planning inputs as rows of the returned
  table**, so the assumptions a design was planned under travel with the
  result as tidy data rather than a printed footer. The supposed effect
  is labeled `supposed_smd` to make clear it is a value the researcher
  posits (a minimally important effect or a value believed to be true in
  the population), not a sample estimate; power planning also echoes
  `desired_power`, `alpha_level`, and a numeric `tails` (2 or 1), and
  AIPE echoes `width` (and `assurance` when supplied). The `value`
  column stays numeric.

### Correctness

- `reliability_omega_h()` now places the observed total-variance
  denominator on the same maximum-likelihood (N-divisor) metric as the
  fitted loadings, matching
  `MBESS::ci.reliability(type = "hierarchical")` and
  `semTools::compRelSEM(obs.var = TRUE)`.

DMAR is the modern, more general reimplementation of MBESS that reflects
how the methods are now used across quantitative psychology, sociology,
education, management, marketing, and information systems. MBESS remains
stable on CRAN; DMAR is the recommended path forward for new users.

### New measurement functions

- **[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  and
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  handle missing data by full information maximum likelihood, with
  auxiliary variables.** The measurement family thereby catches up with
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md), whose
  vignette makes the case against listwise deletion that these functions
  previously ignored. Two new arguments, also passed through by
  [`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md):
  `missing = c("listwise", "fiml")`, defaulting to listwise deletion so
  no existing result changes, and `aux`, a character vector naming
  auxiliary columns of `data` that enter as saturated correlates
  (Graham, 2003): correlated freely with each other and with every
  item’s residual, never loading on the factor. Supplying `aux` implies
  `missing = "fiml"`. The analytic alpha applies the classical formula
  to the FIML estimate of the item covariance matrix; the model based
  estimators fit with lavaan’s `missing = "ml"`; robust omega’s observed
  total variance comes from the FIML covariance matrix. MBESS documents
  an `aux` argument on `ci.reliability()` that no longer runs (upstream
  API drift in semTools), so DMAR implements the saturated correlates
  model directly in lavaan syntax rather than depending on it. The
  returned table now always carries an `N_complete` row beside `N`,
  making the cost of listwise deletion visible at a glance, and the
  treatment is recorded in `missing` and `aux` attributes. Interval
  methods that cannot be made correct under FIML (the complete-data
  closed forms, ADF, the profile likelihood) are errors, never silent
  fallbacks; the bootstrap resamples partially observed rows and refits
  by FIML. Under MAR missingness driven by an auxiliary, the test suite
  shows FIML-with-auxiliary reducing the bias of listwise deletion
  roughly 38-fold across 200 replications, and a hand-specified
  saturated correlates model in lavaan reproduces the estimates exactly.

- **[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)**
  gains a `denominator` argument selecting how the total variance in the
  denominator of coefficient omega is estimated: `"observed"` (the
  default: robust omega, the variance of the composite estimated
  directly from the data, the coefficient Kelley and Pornprasertmanit,
  2016, call hierarchical omega and `MBESS::ci.reliability()` calls type
  `"hierarchical"`) or `"model_implied"` (the textbook form). The two
  definitions coincide in the population when the single-factor model is
  correctly specified; only robust omega retains its interpretation as
  the proportion of the variance of the composite actually computed when
  the model is misspecified, which is why it is the default. The help
  page states the properties of each choice and the distinction from the
  bifactor omega-hierarchical of Zinbarg, Revelle, Yovel, and Li (2005).
  The returned object records the choice in a `denominator` attribute.

- **[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)
  accepts ordered-categorical items** (`ordered = TRUE` or a vector of
  item names; raw data required, each factor all ordered or all
  continuous). The model is fit by WLSMV to polychoric correlations with
  thresholds in the theta parameterization, which keeps every defined
  measurement quantity available, with robust standard errors. Because a
  sum score of ordered items lives on the metric of the observed
  categories rather than the latent response metric, each ordered
  factor’s omega is the Green and Yang (2009) categorical sum score
  omega computed from the same fit; the substitution is announced in a
  message, recorded per factor in an `omega_metric` attribute, and the
  delta method interval columns are NA for those rows (use
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  for a bootstrap interval). A single-factor ordered
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)
  reproduces
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  to 1e-6 in the regression tests despite the different
  parameterizations, which pins the parameterization invariance of the
  computation.

- **[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  gains an `estimator` argument** carrying the two routes to coefficient
  alpha, which were briefly two functions. `estimator = "analytic"` (the
  default) is the classical closed-form equation applied to the observed
  covariance matrix, the number a hand calculation produces;
  `estimator = "model_implied"` is the reliability implied by the
  tau-equivalent (equal loadings) single-factor model fit by maximum
  likelihood, which brings a delta method standard error, a robust
  (Satorra-Bentler) variant, the profile likelihood interval, and a
  testable fit of the equal-loadings claim. Users of MBESS will know
  them as `ci.reliability(type = "alpha")` and `type = "alpha-cfa")`;
  the point estimates match those to 1e-6 in the regression tests. The
  value `"model_implied"` is deliberately the word
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  already uses for the analogous choice, so one vocabulary covers both
  coefficients, and the
  [`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
  wrapper gains a matching `estimator` argument beside its
  `denominator`.

  This replaces the short-lived `reliability_alpha_analytic()`, whose
  name was inherited MBESS jargon that described neither the estimand
  nor the contrast. Merging the two also removes a defect their
  separation created: `reliability_alpha(ci_method = "likelihood")` used
  to report the classical point estimate beside an interval profiling
  the model implied coefficient, two different quantities, so on a
  misspecified model the interval could exclude the estimate printed
  above it (on the `psych` `bfi` A scale, estimate 0.431 against an
  interval of \[0.591, 0.636\]). Each interval method now belongs to the
  estimator that can supply it, and asking for one the chosen estimator
  cannot give is an error naming the estimator to use instead.

- **Profile likelihood confidence intervals** join the closed forms:
  `ci_method = "likelihood"` in
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  (model implied denominator) and
  [`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  (profiling the tau-equivalent, CFA-based alpha). The interval is the
  set of population values not rejected by the likelihood ratio test,
  computed by refitting the single-factor model under a nonlinear
  constraint on the model implied reliability. It respects \[0, 1\], is
  not forced to be symmetric, works from raw data or a covariance
  matrix, and matches `MBESS::ci.reliability(interval.type = "ll")` to
  the fourth decimal in the regression tests.

- **No bootstrap runs unless the user requests one, anywhere in the
  reliability family.** Robust omega and categorical omega, whose
  confidence intervals are bootstrap based, report the point estimate by
  default with a message naming the exact call that produces the
  recommended interval (percentile or BCa for robust omega; BCa for
  categorical omega). The model implied omega keeps its closed-form
  robust ML interval as the default, and alpha and KR-20 keep their
  closed forms. When a bootstrap is requested, `B = 10000` replications
  is the default.

- **`reliability_omega_h()` was removed.** Its coefficient is exactly
  `reliability_omega(denominator = "observed")`, and the “h” (for
  “hierarchical”) described no hierarchy in the single-factor model the
  function fits while inviting confusion with the bifactor
  omega-hierarchical. The
  [`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
  wrapper drops `type = "omega_h"` and instead forwards a `denominator`
  argument to
  [`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md).
  The documentation refers to the coefficient as robust omega, records
  the original motivation (model misfit as minor common factors, with
  the coefficient isolating the general factor’s variance against the
  observed composite variance) and the authors’ retrospective preference
  for a name stating the behavior, states the qualifications the word
  robust requires (robust to misspecification of the total variance
  only; distinct from outlier-robust estimation and from robust standard
  errors), and notes the design principle shared with categorical omega:
  in both, the total variance in the denominator is not taken from the
  fitted factor model.

- **[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  is the categorical omega function’s full name**, and its only name (an
  earlier `reliability_omega_c()` alias was removed). The
  [`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md)
  wrapper’s canonical type is `"omega_categorical"` with `"omega_c"`
  accepted as a shorthand. The coefficient attribute is now
  `"omega_categorical"`.

- **The reliability vignette was rewritten** around the framing of
  Kelley and Pornprasertmanit (2016): choosing the coefficient (a claim
  about the measurement model and the composite being scored) and
  choosing the interval (an empirical performance question their Monte
  Carlo studies answered, which is what the family’s defaults encode).
  One running example walks alpha versus omega, the two omega
  denominators under a minor-factor contamination, categorical omega
  under same versus differing threshold patterns, and coefficient H as a
  different composite rather than a different assumption.

- **[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)**
  now accepts a published k x k frequency table (`table =`) in place of
  raw rater vectors, and custom weight matrices in either scaling:
  agreement weights (diagonal 1) or Cohen’s (1968) ratio-scaled
  disagreement weights (`weight_scaling = "disagreement"`, zero
  diagonal, invariant to positive rescaling, converted internally via w
  = 1 - v / max(v)). Asymmetric weight matrices are supported for
  validity designs where the two directions of a confusion carry
  different costs. The help page replicates Cohen’s (1968) Table 1
  analyses in full: unweighted kappa .492, weighted kappa .348 under his
  disagreement weights, .574 with the 6 and 1 weights interchanged, his
  asymmetric computer-diagnosis validity example (.353), and his Formula
  10 and 13 standard errors (.0901, .0916, z = 3.80), which the examples
  reproduce for the historical record while the function reports the
  Fleiss, Cohen, and Everitt (1969) standard error that superseded them.
  When both raters are supplied as factors with the same level set,
  their level order is now respected instead of alphabetical sorting,
  which previously could silently misalign ordinal categories with
  linear, quadratic, or custom weights. Every result also carries a
  `cells` attribute holding the per-cell detail in the form of Cohen’s
  Table 1: observed proportion, chance-expected proportion, and the
  weights (both scalings when disagreement weights were supplied), one
  row per cell of the confusion matrix.

- **`diagnosis_agreement`** ships Cohen’s (1968) Table 1 as a data set
  in its original layout (Judge B in rows, Judge A in columns): one row
  per cell with the frequency, Cohen’s ratio-scaled disagreement weight,
  the observed proportion, and the chance-expected proportion. The
  reconstruction is verified against every quantity computed from the
  table in the paper.

- **A weighted kappa vignette** works Cohen’s (1968) illustration in
  full on the `diagnosis_agreement` data: unweighted and weighted kappa,
  both weight scalings, his Formula 10 and 13 standard errors beside the
  Fleiss-Cohen-Everitt interval, linear and quadratic weighting with the
  weighted-kappa-equals-r identity under equal marginals, and the
  asymmetric validity example, including a neutral working of how the
  printed weight display’s orientation relates to the published values.

- **[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)**
  fits a confirmatory factor analysis model with one or more factors,
  each specified by naming its indicators. The measurement structure is
  specified by describing what is constrained (`equal_loading`,
  `equal_intercept`, `equal_error`, each a single value or per-factor),
  and the function names the classical structure the description implies
  (congeneric, essentially tau-equivalent, tau-equivalent, essentially
  parallel, or parallel; Graham, 2006) in the printed header and the
  `"model"` attribute. The table reports every estimate with a
  confidence interval, and per factor coefficient omega, the average
  variance extracted, and coefficient H as `lavaan` defined parameters
  with delta method standard errors and intervals (no semTools
  involvement). `output = "measurement"` gathers the measurement
  properties, the latent correlations, and
  [`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md) per
  factor pair; `output = "fit"` hands back the `lavaan` object so two
  descriptor fits feed
  [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
  directly. The RMSEA interval level is reported explicitly as the
  `rmsea_ci_level` row.

- **[`plot_cfa_k()`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md)**
  displays the item-level loadings, error variances, or intercepts of a
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) fit
  with confidence intervals, one panel per factor, with a dashed
  reference line that shows the equated value (or, for free estimates,
  the informal “one common value” anchor), so the equality questions
  behind the classical structures can be seen before they are tested.

- **[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)**
  computes the bifactor dimensionality and reliability indices (ECV,
  omega, omega hierarchical and hierarchical subscale, PUC, and
  coefficient H) from a fitted bifactor `lavaan` model, with a guard
  that flags improper (Heywood) solutions (Rodriguez, Reise, & Haviland,
  2016).

- **[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)**
  quantifies Thurstonian simple structure in a loading matrix: Hofmann
  item complexity, the hyperplane proportion, and pure/complex item
  counts.

- **[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md)**
  gives the Browne and Cudeck (1989) expected cross-validation index for
  a covariance-structure model, with a confidence interval derived from
  the noncentral chi square
  ([`conf_limits_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md));
  accepts a `lavaan` fit or a published fit table.

- **[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md)**
  and
  **[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md)**
  implement the single-common-factor (Harman) screen and the Lindell and
  Whitney (2001) marker-variable adjustment for common method variance.
  The Harman screen is implemented factor analytically: a one-factor
  model is fit by maximum likelihood
  ([`stats::factanal()`](https://rdrr.io/r/stats/factanal.html)) and the
  reported proportion of variance is the common factor’s, not a
  principal component’s.

### ANCOVA multiple comparisons

- **Bryant–Paulson multiple comparisons for ANCOVA**: a new family for
  simultaneous inference on covariate-adjusted means when the covariates
  are random.
  [`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
  gives the simultaneous critical value (the ANCOVA member of the `cv_*`
  family, reducing to `sqrt(2) * cv_tukey_hsd()` when there are no
  covariates), and
  [`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
  places simultaneous (familywise) confidence intervals on contrasts of
  adjusted means, the familywise counterpart of the per-comparison
  [`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md).
  Both are computed exactly from the Bryant–Paulson generalized
  studentized range distribution
  ([`qbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
  /
  [`pbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
  /
  [`dbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)),
  evaluated as a Beta mixture of
  [`ptukey()`](https://rdrr.io/r/stats/Tukey.html) rather than read from
  a table. Three vignettes cover the critical values, an end-to-end
  ANCOVA workflow, and a simulation confirming exact familywise error
  control. Implements Bryant and Paulson (1976) and Bryant and Bruvold
  (1980).
- **[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md)**:
  the values of the covariate at which two groups differ significantly
  when the within-group regression slopes are *not* equal (heterogeneity
  of regression), where the group difference is a function of the
  covariate rather than a single number. The boundaries solve the
  quadratic that sets the squared group difference against its sampling
  variance at the critical value, the Johnson and Neyman (1936)
  procedure; with more than two groups the calculation is carried out
  for every pair.
  [`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md)
  draws the estimated difference across the covariate with the
  confidence band the region is read from, so the plot *is* the decision
  rule.

### The remaining textbook critical-value tables

- **[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
  [`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
  and
  [`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md)**
  complete the `cv_*` family’s coverage of the Maxwell, Delaney, and
  Kelley (2027) Appendix. The family had
  [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md) and
  [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md) but not
  the *F* or chi square counterparts, and the Bonferroni *F* table had
  no function at all.
  [`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md) covers
  Appendix Table A.2,
  [`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md)
  Table A.3, and
  [`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md)
  Table A.9; the tests assert each against the printed values.
- [`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md) and
  [`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md)
  default to `alternative = "greater"` rather than the `"not_equal"`
  that [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md)
  and [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)
  use. Neither distribution is symmetric, and both are used one-sided in
  the upper tail for the tests they serve: a restricted model fits worse
  than a full one, so evidence against a restriction is a large *F*,
  never a small one. Both tails remain available, through
  `alternative = "not_equal"` or through `alpha_lower` and
  `alpha_upper`, for an interval on a variance or on a ratio of
  variances. Both accept a noncentral parameter, as
  [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md) does.
- [`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md)
  reports the per-comparison rate alpha / C back in its `area_greater`
  column, which is the whole of what the adjustment does. Its help page
  separates it from the rank-sum procedure of
  [`dunn_test()`](https://yelleknek.github.io/DMAR/reference/dunn_test.md):
  Dunn (1961) is the Bonferroni procedure, Dunn (1964) is the
  nonparametric one.

### Pairwise comparisons without the usual assumptions

- **[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)**:
  simultaneous confidence intervals for all pairwise comparisons when
  homogeneity of variance is not assumed. Every other all-pairs
  procedure in the package pools the within-group variances into MS_W,
  so none of them is robust when that assumption fails. Games-Howell
  uses a separate error term and a Welch-Satterthwaite degrees of
  freedom for each pair, then takes its critical value from the
  studentized range, following Maxwell, Delaney, and Kelley (2027,
  Chapter 5, Equations 5.13 and 5.14). It is the heterogeneity-robust
  counterpart of
  [`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
  and handles unequal n as a matter of course. With two groups it is
  exactly Welch’s *t* test, which the tests assert against
  `t.test(var.equal = FALSE)`. Implements Games and Howell (1976).
- **[`dunn_test()`](https://yelleknek.github.io/DMAR/reference/dunn_test.md)**:
  Dunn’s rank-sum test of all pairwise differences, the follow-up to a
  significant Kruskal-Wallis test. It ranks all observations together,
  as the omnibus test does, and uses the variance of the ranks implied
  by the Kruskal-Wallis null (with the tie correction), so it stays
  coherent with the omnibus result in a way that running a Mann-Whitney
  test on each pair does not. The `method` argument passes the pairwise
  *p*-values to [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).
  Implements Dunn (1964). Note this is the *nonparametric* Dunn
  procedure, not the Bonferroni procedure of Dunn (1961) that Maxwell,
  Delaney, and Kelley
  2027. call Dunn’s procedure; the help page disambiguates the two.
- **[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
  and
  [`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md)**:
  randomization (permutation) tests for two independent groups and for
  paired observations. The two-group test refers the observed statistic
  to its distribution over reassignments of the observed scores to the
  groups, so the *p*-value needs no assumption about the population’s
  shape, and the result travels with the effect sizes the test only
  screens for: the mean difference with a randomization-based interval,
  the standardized mean difference with a noncentral *t* interval, the
  common language effect size, and Cliff’s delta. The paired test treats
  the within-pair sign of each difference as the randomization
  mechanism, enumerating all 2^n sign patterns exactly for small samples
  and sampling them otherwise. Implements the logic of Fisher
  1935. as developed by Edgington and Onghena (2007).
        [`plot_randomization_test()`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md)
        displays the randomization distribution the *p*-value is read
        from, with the observed statistic marked.

### API and Documentation Consistency

- **The Type I error rate is `alpha_level` everywhere.** Twenty-two
  functions took a bare `alpha` while twenty-eight already took
  `alpha_level`, and the bare name was carrying three incompatible
  meanings inside one package: the Type I error rate in the `cv_*`,
  equivalence, TOST, sequential, and Fisher-exact families; coefficient
  alpha in
  [`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md);
  and ggplot2 transparency in
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md).
  A reader who learned one meaning would misread the others. The Type I
  error uses are now `alpha_level` throughout, so fifty functions share
  one unambiguous name, and the argument that echoes it in a result
  table is named to match. The other two keep `alpha`, where it is
  unmistakable: `var_alpha(alpha = )` is the coefficient the function is
  named for and parallels `var_omega(omega = )`, and transparency is the
  universal ggplot2 convention in a plotting call where no error rate
  appears. Since the package is unreleased, the old name fails loudly
  rather than being aliased. `alpha_lower` and `alpha_upper`, which name
  the two tails of an asymmetric critical value, are unchanged.

### Multiple-comparison critical values computed exactly

- **[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md)
  and
  [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
  are now deterministic.** Both formerly obtained their multivariate-*t*
  quantile from
  [`mvtnorm::qmvt()`](https://rdrr.io/pkg/mvtnorm/man/qmvt.html), a
  Monte Carlo integrator whose sampling error reached a few hundredths
  at `alpha_level = .01`, enough to move the second decimal of a tabled
  critical value. They now evaluate the quantile by exact numerical
  quadrature and root finding: the studentized maximum modulus
  factorizes into a single one-dimensional integral over the shared
  scale (the *m* statistics are independent given it), and the
  balanced-design Dunnett statistic, with its constant correlation
  `1/2`, factorizes through a one-factor representation into two nested
  one-dimensional integrals. The returned values are reproducible to the
  solver tolerance, reproduce the published Dunnett and studentized
  maximum modulus tables more closely than the Monte Carlo path did, and
  no longer depend on `mvtnorm`. Both functions now accept `df = Inf`
  (the known-variance normal limit).
- The `seed` argument of
  [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) and
  [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
  is removed: the computation is no longer random, so there is nothing
  to seed. Call sites that passed `seed` should drop it.
- **[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
  adjusted *p*-values are now exact.** They were computed from
  [`mvtnorm::pmvt()`](https://rdrr.io/pkg/mvtnorm/man/pmvt.html), a
  Monte Carlo integrator (and, when was absent, from a conservative
  Sidak-Bonferroni fallback). They now use the same deterministic
  one-factor integral as the critical value, so the reported
  `p_adjusted` is reproducible and no longer depends on . The shared
  numeric engine lives in `R/dunnett_internals.R`.
- A new vignette, “Reproducing the Textbook Critical-Value Tables,”
  walks the `cv_*` family reproducing the Appendix critical-value tables
  and validates them by simulation.
- **[`qbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
  /
  [`pbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
  /
  [`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
  are now accurate at small error df.** They obtained the
  studentized-range part of the distribution from
  [`stats::ptukey()`](https://rdrr.io/r/stats/Tukey.html), whose
  algorithm loses accuracy at small df (at nu = 3, k = 20 by ~3e-4 in
  probability, enough to move the critical value by ~0.2, and more at nu
  = 2). For nu \< 7 the studentized-range distribution is now evaluated
  directly, without `ptukey`, by integrating the probability integral of
  the range against the chi squared error density. The range CDF is
  splined and cached per group count; the knots are placed at a fixed
  spacing of about 0.003 in the range argument, fine enough that the
  monotone (Fritsch and Carlson) interpolation error stays below 1e-9
  and the returned critical value matches the direct integral to about
  seven figures, while the small-df path still costs a fraction of a
  second per group count. The functions reproduce Bryant and
  Paulson’s (1976) Table 1 exactly, to the two decimal places tabled,
  over the whole of its range: both tail areas, all three covariate
  counts, every tabled group count, and every tabled error degrees of
  freedom from nu = 2 to nu = 120, which is 1188 critical values in all.
  Two entries, q\_.01;2,8,3 = 23.165013 and q\_.01;2,20,4 = 19.745008,
  sit about 1e-5 above the point where the second decimal turns over, so
  they round to 23.17 and 19.75 while the 1976 table rounds them down;
  both were confirmed to fourteen significant figures by two independent
  high-order quadrature engines that share no code with the package. A
  large-scale simulation of the statistic confirms the computed values
  independently. Values for nu \>= 7 are unchanged.

### Bug fixes

- **[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md)
  extracts the full adjusted (least-squares) means table.** From an
  `lm`/`aov` factorial or ANCOVA fit it returns every cell’s
  covariate-adjusted mean with SE and CI, or marginal means over chosen
  factors via `by`, with `weights = "equal"` (the population marginal
  means of Searle, Speed, and Milliken, 1980) or `"proportional"`
  (observed frequencies). The reference grid follows the model’s own
  terms, so transformed covariates and interactions are handled exactly,
  and nonestimable cells in rank-deficient designs error plainly by
  name. Validated against emmeans 2.0.3 at 1e-10 across fourteen fits
  (unbalanced factorials, multiple and transformed covariates, character
  and ordered predictors, weighted fits, an empty-cell design) for both
  weightings, with the comparisons pinned in tests and re-run live in
  tools/oracle_checks.R; also pinned to
  [`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md)’s
  adjusted means and to
  [`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)
  at 1e-10.

- **[`loa()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md)
  reports the Carkeet (2015) pair intervals by default.** The exact CIs
  for the limits of agreement considered as a pair, the construction
  Carkeet recommends for most uses, replace the per-limit intervals as
  the default; `method = "individual"` keeps the one-sided
  tolerance-factor form. Both constructions are anchored in tests to the
  paper’s printed coefficients and worked example, and a misattributed
  comment (the Bland and Altman 1999 approximate SE, previously credited
  to Carkeet) is corrected.

- **Sixteen packages leave Suggests.** DMAR’s test oracles (MBESS,
  emmeans, semTools, metafor, mirt, sirt, multcomp, performance, r2mlm,
  irr, irrCAC, psych, BayesFactor, gsl) are no longer dependencies:
  every live comparison was replaced by its pinned value, with
  provenance comments naming package, version, and date, and the live
  comparisons themselves moved to tools/oracle_checks.R, which re-runs
  all of them against the installed oracles at release time. BiasedUrn
  is replaced by a self-contained Fisher noncentral hypergeometric
  density (agreement 1.2e-13), kableExtra by knitr-only HTML and LaTeX
  table construction, and AMCP by shipping the depression_bdi data
  directly. Suggests drops from thirty packages to fourteen, none of
  them needing compiled system libraries beyond what the remaining
  features genuinely use.

- **The Bayes factor functions take informed priors, stated either
  way.** The prior on the standardized effect can now be an informed
  Cauchy, moving `prior_location` off zero (Gronau, Ly, & Wagenmakers,
  2020), or a normal with `prior_mean` and `prior_sd` for beliefs stated
  as moments; a Cauchy has no mean and no variance, so moment beliefs
  could not previously be expressed at all. The two families are
  exclusive, and the help pages explain what the Cauchy scale fixes (the
  quartiles: half the prior mass within one scale of the location) and
  the exact identity linking the families: a Cauchy is a normal prior
  whose variance is itself uncertain, and the test suite uses that
  identity to cross-validate the two code paths against each other with
  no external oracle. Every result now carries the full posterior of the
  effect in a `"posterior"` attribute (a data frame of `delta` and
  `density`), so any posterior probability can be computed, not only the
  reported ones; the prior is echoed in `prior_location` and
  `prior_scale` rows with the family recorded as an attribute. Default
  calls reproduce the previous JZS results exactly.

- **The Bayes factor functions accept summary statistics.** The JZS
  Bayes factor depends on the data only through the *t* statistic and
  the sample sizes, so
  [`bayes_one_sample_t()`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
  [`bayes_paired_t()`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md),
  and
  [`bayes_independent_t()`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md)
  now take the summary statistics a paper reports (`mean`, `sd`, and
  `n`; `mean_diff`, `sd_diff`, and `n`; or `mean_1`, `sd_1`, `n_1`,
  `mean_2`, `sd_2`, `n_2`) as readily as raw data, following the same
  exactly-one-path rule as the rest of the package. The summary form is
  exact, not an approximation, and the help pages show how a
  standardized effect size enters (a *d* of 0.5 is
  `mean_1 = 0.5, mean_2 = 0` with unit standard deviations). Tests pin
  the raw and summary forms to each other at machine precision.

- **The directional test family speaks snake_case.** The canonical
  `alternative` values of
  [`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
  [`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
  [`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
  [`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
  [`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
  and
  [`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)
  are now `"two_sided"`, `"less"`, and `"greater"`, with the base-R
  spelling `"two.sided"` accepted as an alias so existing calls keep
  working; anything stored on a returned object carries the underscore
  form. The `cv_*` critical value family keeps its deliberately wider
  synonym vocabulary with `"not_equal"` canonical, unchanged.

- **[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md)
  returned a variance that did not shrink with the sample size.** The
  worker multiplied the wrong combination of terms, and the consequence
  was not a small bias: the returned variance converged to a positive
  constant instead of to zero, so its error grew with N. Against a
  200,000-replication Monte Carlo of the sampling distribution it was
  too large by a factor of 1.4 at the smallest cell tested and by 115 at
  the largest, and N times the variance climbed from 0.8 at N = 30 to 51
  at N = 5100 where the true value settles near 0.31. Any standard error
  or Wald interval built on it was badly inflated, and worse in larger
  samples. The variance is now the delta-method transfer of Fleishman’s
  (1980, Eq. 22) exact variance of the unbiased estimator of the
  signal-to-noise ratio, carried to the omega squared scale by his Eq. 8
  with Jacobian (1 - omega²⁾2, which tracks the Monte Carlo to within a
  few percent across every cell tested and to 1.00 by N = 300. Two
  regression tests now guard it: agreement with a fixed-seed Monte
  Carlo, and the consistency check that N times the variance stays
  bounded, which the old formula would have failed. The help page also
  no longer attributes an omega squared variance to Fleishman, who gives
  one only for the signal-to-noise ratio and says explicitly that the
  correlation ratio has an interval and a median but not a variance.

- **Generalized eta squared now classifies interactions the way its
  sources do.** Olejnik and Algina (2003, Eq. 5) and Bakeman (2005)
  define an effect as a measured source of variance when *any* factor in
  its term is measured, so a measured-by-manipulated interaction belongs
  in the denominator automatically.
  [`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)
  and
  [`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)
  had instead included only the effects explicitly listed in `observed`,
  a rule the help page misattributed to Bakeman (2005), who recommends
  the opposite. Listing a factor in `observed` now also places every
  interaction containing it in the denominator, so `observed = "c"`
  reproduces the papers’ worked examples as written; an explicit
  interaction label is still honored as given. Values change only for
  model-interface calls on designs with measured-by-manipulated
  interactions, where the corrected values are smaller. Tests now anchor
  the Keppel and Kirk worked values printed in Olejnik and Algina
  (2003).

- [`signal_to_noise_R2()`](https://yelleknek.github.io/DMAR/reference/signal_to_noise_R2.md)’s
  nonlinear estimator now follows Muirhead’s

  1985. Equation 10 exactly: the c/Y correction is added to the
        untruncated linear estimate before truncation at zero, rather
        than after the linear estimate was itself truncated. The two
        orders agree whenever the sample R-squared is not tiny; for
        sample R-squared below p/(N-1) the reported nonlinear value is
        now smaller, as the paper defines it.

- In
  [`ss_aipe_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md)’s
  assurance search, the critical *t* now tracks each candidate sample
  size’s degrees of freedom instead of reusing the value frozen from the
  expected-width search. The frozen *t* was very slightly conservative;
  the correction moves nothing by as much as one subject, and every
  documented example returns the same sample size as before.

- The
  [`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)
  example carried the wrong sum of squares for its covariate. The
  Maxwell, Delaney, and Kelley Chapter 9 depression example needs the
  within-groups sum of squares of the pretest, which the chapter data
  give as 752.5 (and which the sibling
  [`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
  example already used); the example instead passed 313.37, the
  within-groups sum of squares the covariate explains in the posttest, a
  different row of the same ANCOVA decomposition. The interval moves
  only in the third decimal, but the two pages now agree with each other
  and with the chapter data, and a test recomputes 752.5 from the data
  (via the AMCP package) so the constant cannot drift again.

- [`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) now
  works when the noncentrality parameter is supplied directly. The `ncp`
  branch set the noncentrality but never the standardized contrast, so
  the documented call errored on the closing table construction (three
  terms against two values). The branch now inverts the same
  relationship the `means` and `psi` branches use (the noncentrality is
  the standardized contrast divided by `sqrt(sum(c_weights^2 / n))`), so
  the three parameterizations of one effect return identical tables;
  that equivalence and the ncp path’s coverage are now tested.

- [`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md) and
  [`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) now
  recycle a scalar `n` by the number of contrast weights. The scalar was
  recycled by `length(means)`, which is zero on the `psi`-only path
  (and, for
  [`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md), the
  `ncp`-only path), so those documented calls stopped on the `n` /
  `c_weights` length check. No result changes on the `means` path, where
  the two lengths agree; a `means` vector whose length differs from
  `c_weights` is now rejected with a clear error instead of being
  silently recycled.

- [`ss_aipe_cv()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  now accepts the documented `mu` / `sigma` parameterization, deriving
  `C_of_V = sigma / mu`; previously neither formal was read and the call
  stopped on “argument is of length zero”. Supplying `C_of_V` together
  with `mu` or `sigma` is rejected as conflicting. The `mu` / `sigma`
  path returns exactly the sample size of the equivalent `C_of_V` call
  (for `mu = 10`, `sigma = 1`, `width = .1`, `conf_level = .99`, both
  give `N = 20`). The example prose on the help page said the population
  coefficient of variation was .25 while the calls passed .1 (which
  plans `N = 20`, where .25 plans `N = 100`); the prose now matches the
  calls, and a `mu` / `sigma` example was added.

- [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) with
  `missing = "ml"` now estimates the mean structure with free item
  intercepts. The bare
  [`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html)
  interface it calls turns the mean structure on for FIML but leaves
  `int.ov.free` at `FALSE`, so every item intercept was fixed to zero
  and the loadings absorbed the item means; on any data not centered at
  zero the FIML estimates were wrong (the function’s own example,
  simulated with mean zero, could not show it). The fit now matches
  `lavaan::cfa(..., missing = "ml")` exactly. Listwise fits are
  unaffected.

- [`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)
  now uses the correct between-subjects noncentrality. The
  between-subjects test used the per-group sample size `n` where its
  noncentrality should use the total sample size `N = n a`, so the
  noncentrality was too small by the factor `a` (the number of
  between-subjects groups), halved for two groups. This understated
  between-subjects power and, planning in reverse, overstated the
  necessary per-group sample size; the within-subjects and interaction
  tests were unaffected. The corrected between-subjects test now agrees
  with a split-plot [`aov()`](https://rdrr.io/r/stats/aov.html)
  simulation and, for two groups, reproduces the two-level treatment
  test of
  [`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md)
  exactly (the between-subjects *F*(1, .) is the two-level treatment *t*
  squared), the identity that ties the two planners together. Both a
  value anchor and the cross-planner identity are now tested.

- [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) and
  [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
  now work at a large error df. Both integrated the chi squared error
  variate over `(0, Inf)`, and once df is large that density is a narrow
  bump far from the origin, so the adaptive rule sampled its way past
  the mass and returned zero: at df = 200 the integral evaluated to
  exactly 0 for every candidate quantile, which left the root finder
  with no bracket and the functions stopped with “f() values at end
  points not of opposite sign.” The limits are now the extreme quantiles
  of that variate, which puts the quadrature on the mass at any df. Both
  functions gained a large-df test; the values they already returned are
  unchanged.

- [`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md) no
  longer returns `NaN` when the numerator df is infinite. Supplying
  `ncp = 0` explicitly sends
  [`qf()`](https://rdrr.io/r/stats/Fdist.html) and
  [`pf()`](https://rdrr.io/r/stats/Fdist.html) down their noncentral
  algorithm, which does not admit an infinite numerator df, so the
  infinite-numerator column of Appendix Table A.2 came back `NaN`. Both
  functions now pass `ncp` only when it is nonzero, so the central
  algorithm serves the central case.
  [`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md)
  takes the same guard.

- [`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md) now
  reports the covariate-adjusted omnibus *F* for the treatment effect.
  It previously entered the treatment before the covariate and read the
  sequential (Type I) sum of squares, which is the *unadjusted*
  treatment *F*; the covariate is now entered first so the treatment’s
  sequential sum of squares is its adjusted sum of squares.

- [`ss_aipe_sm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md)
  corrects the noncentrality used in the assurance branch (it was the
  reciprocal `sm / sqrt(n)` instead of `sm * sqrt(n)`), which previously
  left the assurance target unreachable and pinned the search to a
  bracket endpoint.

- [`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md) now
  honors the documented `df_error` argument; supplying it previously
  left the error degrees of freedom undefined and errored.

- `tidy()` and `glance()` on
  [`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md)
  and
  [`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md)
  now return the user-supplied sample size on the realized-power path
  (previously `NA`).

- [`ss_aipe_reliability()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
  now computes confidence intervals when the default `interval = TRUE`
  is used (a string-only equality check had skipped the interval for the
  logical default).

- `ss_aipe_reliability(type = "Factor Analytic")` works again. The path
  read
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)’s
  legacy list layout (`$factor_loadings`, `$parameter_cov`), which the
  reworked
  [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) no
  longer returns, so it errored; coefficient omega and its delta method
  interval now route through the maintained `reliability` internals
  (`.omega_fit_cfa()`, `.ci_omega_delta()`).

### Highlight features

- **[`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)**
  are the new lm-like front end to full information maximum likelihood
  (FIML) regression. The univariate
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) mirrors
  the [`lm()`](https://rdrr.io/r/stats/lm.html) API (formula interface,
  `coef` / `vcov` / `confint` / `summary` / `anova` / `predict` /
  `update` S3 methods, profile / Wald / bootstrap CIs); the multivariate
  sibling
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)
  takes `cbind(y1, y2) ~ ...` and models the joint distribution of
  correlated outcomes, with the residual covariance among outcomes
  estimated as part of the fit. See
  [`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md)
  for when the FIML route actually buys you something over
  [`lm()`](https://rdrr.io/r/stats/lm.html) + listwise deletion.
- **Auxiliary variables in the FIML family.**
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)
  gained an `auxiliary` argument that brings variables related to the
  missingness or to the incomplete outcome into the model as saturated
  correlates (Graham, 2003): each auxiliary is correlated with the
  outcome residual, every predictor, and each other auxiliary, but never
  enters as a predictor, so the focal regression coefficients keep their
  meaning. This is the inclusive analysis strategy (Collins, Schafer, &
  Kam, 2001); it leaves complete-data estimates unchanged and, under
  MAR, recovers information that listwise deletion discards. The
  differential missingness scenario in
  [`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md)
  works through when it helps and when randomization already protects
  the estimand.
- **broom-style `tidy()` and `glance()`** integration via the `generics`
  package: a uniform interface to the broom ecosystem for `mlmr`,
  `mlmr_mv`, `cfa_1`, the reliability family, the long- format CI
  family, the ANOVA effect size CI family, and the power planner family.
  `purrr::map_dfr(fits, generics::tidy)` now works across DMAR outputs.
- **Performance**: inner-loop fast paths for the `convert_R2_*` family
  cut iterative AIPE / sensitivity planning calls by roughly 3x;
  representative
  [`ss_aipe_R2()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)
  calls go from ~2.7s to ~0.9s.
- **Numerical-correctness tests**: a dedicated test file compares
  marquee CI / variance / planner functions to MBESS and to closed-form
  derivations from the foundational papers, so silent numerical
  regressions are caught immediately.

### Data sets

- New data set `test_market`: a small balanced ANCOVA example (sales by
  promotion type, adjusted for a baseline covariate) used to illustrate
  the Bryant–Paulson simultaneous intervals in
  [`vignette("bryant_paulson_ancova", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/bryant_paulson_ancova.md)
  and
  [`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md).

- New data set `drinks_trial`: nine-month follow-up drinks-per-week
  outcomes for the *N* = 88 homeless alcohol-dependent participants in
  Smith, Meyers, and Delaney’s

  1998. randomized trial of the Community Reinforcement Approach at the
        Salvation Army Adult Rehabilitation Center in Albuquerque, New
        Mexico. Two consecutive cohorts: Cohort 1 compared Standard,
        CRA, and CRA + Disulfiram with cell sizes 17, 15, 19; Cohort 2
        dropped the disulfiram cell after Cohort 1 results and compared
        Standard against CRA with cell sizes 20,

  &nbsp;

  17. The outcome ships in raw form (heavily right-skewed, range 0 to
      624.6 drinks per week) and on the log10 scale used in the
      published analyses to recover approximate normality. Reproduced in
      Maxwell, Delaney, and Kelley (2027, *Designing Experiments and
      Analyzing Data: A Model Comparison Perspective*, 4th ed.,
      Routledge), Chapter 3, Section 3.10.4.

- New data set `bessel_errors`: Friedrich Wilhelm Bessel’s (1818) 9-bin
  grouped frequency distribution of the absolute errors of 300 stellar
  position observations made by British Astronomer Royal James Bradley
  at the Greenwich Observatory between 1750 and 1762. Both the observed
  and expected (normal-model) frequency columns sum to 300, matching
  Maxwell, Delaney, and Kelley (2027, 4th ed.), Table 1.4. Documented as
  a worked example for approximating moments from grouped frequency data
  (frequency-weighted means and variances using bin midpoints) and for
  plotting empirical-versus-theoretical frequency comparisons. Ships in
  the original grouped form Bessel reported; individual error values are
  not extant.

- New data set `prime_time_achievement` (also accessible via the short
  alias `Prime_Time`): the full Indiana Prime Time third grade
  achievement evaluation file (Lapsley, Daytner, Kelley, and Maxwell,
  2002, ERIC ED466679), built from the original Indiana Department of
  Education SPSS system file. 10,927 students nested in 586 classrooms
  in 163 schools in 61 school corporations (district x region
  combinations) in 9 educational service regions on 113 variables.
  Includes ISTEP+ NCE composites (the criterion in the published HLM
  analyses), Gates-MacGinitie and AANCE test scores, NPA cognitive
  ability scores, classroom enrollment, pupil to teacher ratio, Prime
  Time aide indicator and status, school and corporation demographics
  and finance, and three derived unique cluster identifiers (`corp_id`,
  `school_id`, `class_id`) that respect the nesting irregularities in
  the source file (one corporation ID spans two regions). Original
  Indiana DOE variable names and the original six-category race coding
  are preserved verbatim; the SPSS variable labels are retained as a
  `label` attribute on every column. The SPSS Select Cases artifact
  `FILTER_$` and the 888 “not applicable” codes have been dropped or
  recoded to `NA`. Includes documented examples that show level-1,
  level-2, and level-3 `lmer` fits mapped onto the multilevel framework
  used in Lapsley et al. (2002) and in Finch, Bolin, and Kelley (2019,
  *Multilevel Modeling Using R*, 2nd ed., CRC Press, chapters 3, 4, 6,
  9, 10).

- New data set `holzinger_swineford` (also accessible via the short
  alias `HS_Data`): the complete Holzinger and Swineford (1939) factor
  analysis data, 301 pupils on 26 ability tests from the Pasteur (n
  = 156) and Grant-White (n = 145) elementary schools in Chicago.
  Variable names follow the MBESS convention (e.g.,
  `t1_visual_perception` through `t26_flags`) so scripts written against
  `MBESS::HS` port over with a single rename of the object. The values
  are the corrected version of the data, identical to `MBESS::HS` as of
  MBESS 4.9.3 and to `psychTools::holzinger.raw`. The documentation
  discusses the bi- factor study design, the five ability blocks
  (spatial, verbal, mental speed, memory, reasoning), the
  Joreskog (1969) 9 test subset, and the silent post-4.6.0 MBESS
  correction that is the source of values still found in `sem::HS.data`
  and
  [`OpenMx::HS.ability.data`](https://rdrr.io/pkg/OpenMx/man/HS.ability.data.html).

- New data set `pygmalion`: the teacher-expectancy data from Rosenthal
  and Jacobson’s (1968) *Pygmalion in the Classroom*, 310 elementary
  school pupils in grades 1 to 6, of whom 64 were randomly designated to
  their teachers as likely intellectual “bloomers” and 246 served as
  controls, with pretest and follow-up IQ. This is the classic benchmark
  for analysis of covariance with *heterogeneity of regression*, and the
  running example for that topic in Maxwell, Delaney, and Kelley
  (*Designing Experiments and Analyzing Data*, Chapter 9). The
  within-group slopes (0.778 for controls, 0.969 for bloomers), the
  pooled residual variance (175.3251), and the covariate variance
  (348.91) reproduce the worked example for the variance of the
  estimated treatment effect at selected covariate values
  (cf. `MBESS::var.ete`). The same numbers ship with the book’s data
  companion **AMCP** as `chapter_9_exercise_15`; here the experimental
  condition is a labeled factor (`Control`, `Bloomer`) and the columns
  use DMAR’s descriptive names, with no measured value altered. Pairs
  with
  [`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md) and
  is documented in
  [`vignette("pygmalion", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/pygmalion.md).

### Migrating from MBESS

The largest visible change is a uniform `snake_case` interface for both
function names and argument names.

- Argument names that were `dot.case` in MBESS are `snake_case` in DMAR.
  The most common renames a user will hit when porting a script:

  | MBESS                | DMAR                                     |
  |----------------------|------------------------------------------|
  | `conf.level`         | `conf_level`                             |
  | `Random.Predictors`  | `random_predictors`                      |
  | `Specified.N`        | `specified_N`                            |
  | `alpha.lower`        | `alpha_lower`                            |
  | `alpha.upper`        | `alpha_upper`                            |
  | `degrees.of.freedom` | Use `df` or the explicit `df_1` / `df_2` |
  | `Group.1`, `Mean.1`  | `group_1`, `mean_1`                      |

  The rule of thumb: change every `.` between words to `_`, and
  lowercase the non-statistical prefix words. Capitals are preserved
  when the capital is statistically meaningful: `R2` (squared multiple
  correlation), `N` (sample size), `S` (a covariance matrix), `Lambda`
  (a factor-loadings matrix), `F_value` (an *F*-statistic).

- `aipe_smd()` is renamed to
  [`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
  to align it with the rest of the `ss_aipe_*` family.

- Every estimation, inference, and planning function returns a
  `data.frame` with `term` and `value` columns. (Plotting functions
  return a `ggplot` object, and a few utilities return their natural
  type.) Return objects in MBESS were sometimes named lists, sometimes
  data.frames, and sometimes vectors; the unified tidy return makes the
  package compose cleanly with the rest of the modern R ecosystem.

- Confidence intervals and sample size annotations travel with every
  effect size plot.
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  and
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)
  default to `show_ci = TRUE` and `show_n = TRUE`.

There are no deprecation shims. Old MBESS argument names will throw
“unused argument” errors at the call site, which is intentional; the fix
is mechanical, and a silent-forwarding shim would hide it.

### What is new in DMAR (beyond MBESS)

DMAR adds 40+ functions across families that were absent or
under-developed in MBESS. New families include:

- Mediation, both halves:
  [`mediate()`](https://yelleknek.github.io/DMAR/reference/mediate.md)
  analyzes the simple mediation model (optional covariates) with
  percentile bootstrap, BCa, Monte Carlo, and Sobel intervals for the
  indirect effect, and
  [`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)
  plans the study by joint significance with exact noncentral *t*
  component powers (Sobel power for comparison), validated against
  raw-data simulation.
  [`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
  extends the inference side to arbitrary mediation structures through
  the likelihood ratio model comparison framework (see its own section
  above).

- Measurement invariance:
  [`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  fits the configural, metric, scalar, and strict multi-group ladder for
  any measurement model, with a thresholds rung for ordered indicators,
  and reports the full comparison table (fit, the likelihood ratio test
  per step, and delta CFI / delta RMSEA), pinned to direct lavaan fits
  in the tests. See the entry at the top of this file for the arguments
  and the returned `"fits"` attribute.

- Construct validity:
  [`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md) (the
  Henseler heterotrait-monotrait ratio with an optional bootstrap upper
  bound) and
  [`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
  (Fornell-Larcker AVE from a lavaan fit or standardized loadings).
  Composite reliability needs no new function; it is coefficient omega
  ([`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).

- Clinical and behavioral endpoints:
  [`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md)
  (per-group responder proportions with Wilson intervals, the Newcombe
  risk difference interval, the number needed to treat, an omnibus chi
  square, and a threshold sweep) and
  [`ci_proportion()`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md),
  the package’s Wilson score interval for a single proportion.

- A meta-analysis family.
  [`meta_es()`](https://yelleknek.github.io/DMAR/reference/meta_es.md)
  pools any effect sizes given sampling variances (REML between-study
  variance by default, with Paule-Mandel, DerSimonian-Laird, and fixed
  effect options, the Hartung-Knapp adjustment on by default, Q-profile
  confidence intervals for tau-squared mapped to I-squared, and a
  prediction interval always reported);
  [`meta_smd()`](https://yelleknek.github.io/DMAR/reference/meta_smd.md)
  (exact Hedges correction by default) and
  [`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md)
  (Fisher’s Z pooling, with optional per-study attenuation corrections)
  are the metric front ends.
  [`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md)
  provides the four classical combined significance tests and
  [`meta_contrast()`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md)
  the Rosenthal and Rubin contrast among effect sizes;
  [`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)
  draws the studies, the pool, and the prediction interval. Validated
  against metafor, and against Raudenbush (1984) line by line in the
  teacher expectancy vignette; the 19 effect sizes of that synthesis
  ship as `teacher_expectancy`.

- Factorial ANCOVA power:
  [`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md)
  extends the factorial ANOVA planner to baseline covariates (error
  variance scaled by 1 - R2, one error df per covariate), with the
  complete 2 x 4 x 3 worked example, planning, simulation, Type III
  analysis, interaction plots, and focused complex comparisons, in the
  `ancova_2x4x3_power` vignette.

- A new vignette, “Power and Precision for the One-Way ANOVA: A Model
  Comparison Perspective,” works the one-way design from the comparison
  of a full model and a restricted model: the omnibus power and sample
  size through
  [`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
  a planned contrast through
  [`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
  effect size confidence intervals through
  [`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
  [`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  [`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
  and
  [`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
  the model comparison made literal with
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md), the
  Type S and Type M consequences of the design through
  [`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
  and accuracy in parameter estimation for a contrast through
  [`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
  following Maxwell, Delaney, and Kelley (2027).

- A Bayesian t family with the probability statement front and center:
  [`bayes_one_sample_t()`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
  [`bayes_paired_t()`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md),
  and
  [`bayes_independent_t()`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md)
  report the JZS posterior of the standardized effect (median, mean,
  credible interval, and P(delta \> 0 \| data)) with the default Bayes
  factor as a secondary row, computed by exact quadrature and validated
  against the BayesFactor package.

- [`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
  (renamed from its working name `correct_attenuation()`) documents and
  demonstrates the latent variable route it approximates: when
  item-level data exist, prefer the two-factor model’s latent
  correlation to the plug-in formula.

- Moments of the noncentral distributions:
  [`moments_nct()`](https://yelleknek.github.io/DMAR/reference/moments_nct.md),
  [`moments_ncf()`](https://yelleknek.github.io/DMAR/reference/moments_ncf.md),
  and
  [`moments_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/moments_nc_chisq.md)
  return the mean, variance, standard deviation, skewness, and excess
  kurtosis of the noncentral *t*, *F*, and chi square distributions,
  with `NA` for moments whose degrees of freedom conditions fail. The
  noncentral *t* mean is the quantity behind the upward bias of the
  standardized mean difference.

- A random-coefficients polynomial growth simulator:
  [`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md)
  generates longitudinal data of any polynomial order (order 0 is a flat
  line) for one or several groups, ties the level-one error to a target
  measurement reliability (reported per occasion), allows
  assessment-time jitter around the nominal schedule, and supports
  structured level-one error covariance (AR(1), compound symmetry,
  Toeplitz, heteroscedastic, or a full matrix). It is the Monte Carlo
  companion to
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md).

- [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
  now plans power for any polynomial change coefficient (intercept,
  linear, quadratic, cubic, …) through a `trend` argument, carrying the
  general Raudenbush and Liu (2001) sampling variance; the linear
  default reproduces the National Youth Survey benchmark.

- [`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)
  returns orthogonal polynomial trend contrast weights stored
  levels-by-trends (ready for
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html) and
  [`lm()`](https://rdrr.io/r/stats/lm.html)) and printed in the Maxwell,
  Delaney, and Kelley Table A.10 layout with a trailing
  sum-of-squared-weights column.

- [`unbiased_R2()`](https://yelleknek.github.io/DMAR/reference/unbiased_R2.md)
  gives the Olkin and Pratt (1958) exactly unbiased estimator of the
  population squared multiple correlation alongside the Ezekiel (1930)
  adjusted estimator (the `adj.r.squared` of `summary.lm`).

- Design consequences:
  [`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
  reports what a chosen design delivers under both of the package’s
  lenses, the significance lens (power, the Type S wrong-sign error
  rate, and the Type M exaggeration ratio, after Gelman and Carlin,
  2014, computed exactly rather than by their simulation) and the
  precision lens (the expected, median, and standard deviation of the
  realized confidence interval width, and the probability the realized
  interval beats a target width, the closed-form versions of the
  `ss_aipe_*_sensitivity()` Monte Carlo terms). Accepts a standard error
  directly or derives it from sd and per-group n. Every `ss_power_*` and
  `ss_aipe_*` help page points to it.

- The Spearman correction for attenuation:
  [`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
  disattenuates a correlation for measurement error in either or both
  variables, with a confidence interval when the sample size is
  supplied, using reliabilities from the `reliability_*` family.

- New parameterization conversions:
  [`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md)
  /
  [`convert_r_d()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md)
  (standardized mean difference and point-biserial correlation, with an
  unequal-group factor) and
  [`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md)
  /
  [`convert_or_d()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md)
  (the Hasselblad and Hedges logistic link to the odds ratio).

- **[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md)
  and
  [`convert_chisq_F()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md)**
  move a test statistic between the *F* and chi square metrics,
  returning the converted statistic itself (a value, not a *p*-value).
  Two conversions are offered through one argument, `df_denominator`.
  The default, `df_denominator = Inf`, treats the *F*’s error variance
  as known and returns the standard scaling, `df_numerator * F_value`
  (and back, `chi_square / df`); an *F* is a chi square whose error
  variance is estimated rather than known, and infinite denominator
  degrees of freedom is the case where it is known. A finite
  `df_denominator` instead returns the chi square value with the same
  upper-tail probability (the same *p*-value) as the *F*, which is exact
  at any denominator degrees of freedom; the two conversions agree as
  `df_denominator` grows, so the scaling default is the large-sample
  limit of the same family. The finite case is computed from ordinary
  upper-tail *p*-values (no logarithms): the upper tail is used because
  the lower-tail probability rounds to 1 in double precision by about
  *F* = 500 at small denominator degrees of freedom, which would send
  the result to infinity, whereas the upper-tail computation stays
  accurate past *F* = 1e20. Each help page states the exact computation.
  The map preserves the *p*-value but does not transport a noncentrality
  parameter, so noncentral work belongs in
  [`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
  and
  [`conf_limits_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md).

- Display helpers extending the package’s *p*-value convention to
  objects DMAR does not produce:
  [`format_p()`](https://yelleknek.github.io/DMAR/reference/format_p.md),
  [`print_anova()`](https://yelleknek.github.io/DMAR/reference/print_anova.md),
  and
  [`print_summary()`](https://yelleknek.github.io/DMAR/reference/print_summary.md).

- DMAR no longer requires the GSL system library: the Gauss
  hypergeometric function behind the exact squared multiple correlation
  moments is computed in base R (and remains accurate as the squared
  multiple correlation approaches 1, where the previous route lost
  precision). Installation now has no system prerequisites.

- Plot functions in the `plot_*` family:
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
  [`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).
  All built on `ggplot2` and colored by default with a neutral,
  colorblind-safe palette; each accepts a `palette` argument
  (`"okabe_ito"` or `"tableau"`) and, where applicable, a `colors`
  override for full manual control.

- Plot colors come from base R. The `plot_*` family colors itself with
  base R’s Okabe-Ito colorblind-safe palette by default, with base R’s
  Tableau 10 available through each plot’s `palette` argument. DMAR
  defines no palette of its own and adds no color dependency; a user who
  wants other colors adds an ordinary `ggplot2` scale to the plot.

- Ordinal and non-parametric effect sizes:
  [`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
  [`vargha_delaney_A()`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
  [`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
  [`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
  [`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)
  (sometimes called Cohen’s U3).

- Agreement and reliability extensions: `bland_altman_loa()`,
  [`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
  [`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
  [`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
  [`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
  [`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
  [`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
  [`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
  [`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  `reliability_omega_h()`,
  [`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
  [`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
  plus the `var_*` family of asymptotic variance utilities.

- Within-subjects, mixed, and multivariate ANOVA:
  [`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
  [`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
  [`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
  [`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
  [`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
  [`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
  [`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
  [`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md).

- Sample size planning has been broadened, particularly the `ss_power_*`
  family, which now covers between-subjects, within-subjects, mixed, and
  multi-level designs across the chapters of Maxwell, Delaney, and
  Kelley (2027).

- Parameterization conversions in the `convert_*` family
  (`convert_R2_f`, `convert_f_R2`, `convert_lambda_R2`,
  `convert_R2_lambda`, `convert_r_Z`, `convert_Z_r`,
  `convert_delta_lambda`, `convert_lambda_delta`, `convert_cor_cov`).
  Every conversion is exact-invertible and the inverse direction is
  shipped as a sibling function where it makes sense.

- A worked simulation study of the AIPE family reports expected and
  realized CI widths, realized coverage, and assurance-achievement rates
  across 10,000 Monte Carlo replications per cell. It establishes how
  the methods perform rather than how they are used, so it is maintained
  alongside the package rather than shipped in it.

### Internal changes

- [`ss_aipe_reliability()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
  no longer embeds a vendored copy of MBESS 3.2.0’s `ci.reliability()`
  for the Monte Carlo assurance search. The interval at each candidate
  sample size now comes from DMAR’s own reliability internals (the
  single-factor delta method interval for the factor analytic type; the
  van Zyl, Neudecker, and Nel closed form for the normal theory type). A
  before-and-after grid covering every model and interval type
  combination confirmed the computed interval widths, and therefore
  every planned sample size, are unchanged.

- [`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md)
  now reads its grouping factor from lme4’s grouping list, so
  interaction groupings such as `(1 | school:teacher)` work; and
  `mlmr_mv` gained the `tidy()` / `glance()` methods its documentation
  promised (one row per coefficient per outcome, with a `response`
  column).

- The test suite runs with `warnPartialMatchArgs`,
  `warnPartialMatchDollar`, and `warnPartialMatchAttr` enabled, so any
  internal reliance on partial matching fails loudly, and the display
  layer (the `dmar_tbl` print methods,
  [`format_p()`](https://yelleknek.github.io/DMAR/reference/format_p.md),
  [`print_anova()`](https://yelleknek.github.io/DMAR/reference/print_anova.md),
  [`print_summary()`](https://yelleknek.github.io/DMAR/reference/print_summary.md))
  is covered by testthat snapshots.

- The `seed` argument defaults to `NULL` everywhere it is exposed,
  meaning the function uses the caller’s current RNG state and does not
  seed; a user who wants reproducibility passes an explicit integer, and
  a supplied seed restores the caller’s RNG state on exit. The value
  `113` is used only in `@examples` and tests, never baked into a
  function default.

- The cluster-randomized design helpers under `ss_aipe_crd*` were
  consolidated into a single set of shared internals
  (`R/ss_aipe_crd_internals.R`); both the difference and effect size
  families now call the same `.find_*_crd_*()` back end.

- [`ss_aipe_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md)
  had a typo in the cubic change coefficient constant (`K_3` was coded
  as `1/1000800`, about one tenth of the correct `1/100800`). Because
  the constant sits in the denominator of the slope variance, this
  inflated the cubic variance, and hence the resolved sample size,
  roughly tenfold for `trend = "cubic"`. The constant now matches the
  closed form `K_p = (p!)^2 / [(2p)! (2p+1)!]` of Raudenbush and Liu
  (2001, p. 392), and a regression test recovers it directly from an OLS
  fit. The linear (`1/12`) and quadratic (`1/720`) constants were
  already correct.

- [`ss_aipe_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md)
  now honors its documented contract that a user may supply either
  `error_variance` or the converted variance
  `variance_true_minus_estimated_trend`. Previously `error_variance` was
  effectively required, and the consistency check between the two used
  `round(..., 3)`, which mishandled the small variances typical of
  slope-change designs. The cross-check now uses
  [`all.equal()`](https://rdrr.io/r/base/all.equal.html), and supplying
  neither argument raises an informative error.

- [`ss_aipe_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md)
  now scales the within-subject contribution to the slope variance by
  `frequency^(2p)`, matching
  [`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
  and Raudenbush and Liu (2001, p. 392). The change coefficient is a
  per-unit-time rate, so its sampling variance is
  `error_variance * frequency^(2p) / (sum of the squared polynomial weights)`;
  the factor was previously omitted. Every Kelley and Rausch (2011)
  benchmark uses `frequency = 1`, where the factor is 1, so all tabled
  sample sizes are unchanged; the correction only affects designs with
  `frequency != 1`, where it now agrees with the direct OLS slope
  variance on the actual time grid.

- [`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md)
  now simulates the same estimand its planner targets: the between-group
  difference in mean slopes (the group-by-time change parameter), using
  two independent groups, a pooled standard error, and `2n - 2` degrees
  of freedom. The previous simulator built a one-group mean-slope
  interval, whose width was systematically `1/sqrt(2)` of the planned
  target at every frequency. The realized mean CI width now tracks the
  planned width. The schema terms `mean_slope` / `median_slope` /
  `sd_slope` are renamed to `mean_slope_diff` / `median_slope_diff` /
  `sd_slope_diff` to reflect the difference estimand.

### Authorship

Ken Kelley (Department of Information Technology, Analytics, and
Operations; Mendoza College of Business; University of Notre Dame) is
the package author and maintainer. Bug reports and feature requests are
welcomed by email to <kkelley@nd.edu>; please put “DMAR” in the subject
line.
