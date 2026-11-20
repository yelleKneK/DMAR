# Live oracle comparison against semTools, run at release time after DMAR is
# installed (it calls DMAR::fa_second_order), as the other snippets in this
# folder are.
## Written against semTools 0.5.9 and lavaan 0.7.3. The test pins the
## values semTools gives, since semTools is not a dependency of DMAR.
## from tests/testthat/test-fa_second_order.R
local({
  # fa_second_order() reports, for each first-order factor, omega for the
  # unweighted sum of its items and omega_HS, the share of that sum's
  # variance due to the factor's disturbance, so the share due to the
  # second-order factor is their difference. semTools::compRelSEM() with
  # obs.var = FALSE reads a second-order fit's first-order factors as
  # composites whose true score is the second-order factor, with the
  # disturbances counted as error, and computes that share from the same
  # lavaan fit. Its printed description of each composite is captured, and
  # the numbers must agree to the optimizer's precision.
  data(holzinger_swineford, package = "DMAR")
  factors <- list(
    spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
    verbal  = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
    speed   = c("t10_addition", "t12_counting_groups_of_dots",
                "t13_straight_and_curved_capitals"),
    memory  = c("t14_word_recognition", "t15_number_recognition",
                "t17_object_number"))
  second <- list(general = names(factors))
  tab <- DMAR::fa_second_order(holzinger_swineford, factors, second,
                               output = "measurement")
  fit <- DMAR::fa_second_order(holzinger_swineford, factors, second,
                               output = "fit")
  invisible(utils::capture.output(
    share <- semTools::compRelSEM(fit, obs.var = FALSE, simplify = TRUE)))
  est <- function(term) tab$estimate[match(term, tab$term)]
  f <- names(factors)
  stopifnot(isTRUE(all.equal(
    est(paste0("omega_", f)) - est(paste0("omega_HS_", f)),
    as.numeric(share[f]), tolerance = 1e-6)))
  # The constants the test pins.
  stopifnot(isTRUE(all.equal(
    as.numeric(share[f]),
    c(0.4263133166, 0.2504824422, 0.2305681256, 0.1234972521),
    tolerance = 1e-6)))
})
