# Live oracle comparison against semTools, run at release time after DMAR is
# installed (it calls DMAR::fa_bifactor), as the other snippets in this
# folder are.
## Written against semTools 0.5.9 and lavaan 0.7.3. The test pins the
## values semTools gives, since semTools is not a dependency of DMAR.
## from tests/testthat/test-fa_bifactor.R
local({
  # fa_bifactor() reports omega, omega_H, and omega_HS for the unweighted
  # sum of the items as scored, from the model implied variance of that
  # sum. semTools::compRelSEM() with obs.var = FALSE computes the same
  # coefficients from the same lavaan fit: with no `true` argument, the
  # reliability of each factor's composite with every factor counted as
  # true score (omega for the total score, the subscale omega for a group
  # factor); with true = list(<factor> = "<factor>"), the share due to that
  # factor alone (omega_H for the general factor, omega_HS for a group
  # factor). Its printed description of each composite is captured, and
  # the numbers must agree to the optimizer's precision.
  data(holzinger_swineford, package = "DMAR")
  groups <- list(
    spatial = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                "t4_lozenges"),
    verbal  = c("t5_general_information", "t6_paragraph_comprehension",
                "t7_sentence", "t8_word_classification", "t9_word_meaning"),
    memory  = c("t14_word_recognition", "t15_number_recognition",
                "t16_figure_recognition", "t17_object_number",
                "t18_number_figure", "t19_figure_word"))
  latent <- c("general", names(groups))
  tab <- DMAR::fa_bifactor(holzinger_swineford, groups,
                           output = "measurement")
  fit <- DMAR::fa_bifactor(holzinger_swineford, groups, output = "fit")
  invisible(utils::capture.output(
    total <- semTools::compRelSEM(fit, obs.var = FALSE, simplify = TRUE),
    hier <- semTools::compRelSEM(fit, obs.var = FALSE, simplify = TRUE,
                                 true = stats::setNames(as.list(latent),
                                                        latent))))
  est <- function(term) tab$estimate[match(term, tab$term)]
  stopifnot(isTRUE(all.equal(est(paste0("omega_", latent)),
                             as.numeric(total[latent]), tolerance = 1e-6)))
  stopifnot(isTRUE(all.equal(
    est(c("omega_H_general", paste0("omega_HS_", names(groups)))),
    as.numeric(hier[latent]), tolerance = 1e-6)))
  # The constants the test pins, so that a change in semTools that left
  # DMAR agreeing with it would still be seen.
  stopifnot(isTRUE(all.equal(
    as.numeric(c(total[latent], hier[latent])),
    c(0.8755302000, 0.6951105684, 0.9023000439, 0.7283031141,
      0.5523387266, 0.2660922789, 0.6124102558, 0.4816024383),
    tolerance = 1e-6)))
  # The value the help of fa_bifactor() quotes.
  stopifnot(isTRUE(all.equal(round(as.numeric(hier["general"]), 3), 0.552)))
})
