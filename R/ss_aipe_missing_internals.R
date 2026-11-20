# Internal helpers for AIPE sample size planning under missing data.
#
# Shared by ss_aipe_reliability(), which plans for a fraction of missing
# information (FMI) in the coefficient, and ss_aipe_reliability_sensitivity(),
# which passes that planning value to the planner, simulates incomplete
# data, missing completely at random, analyzed by full information maximum
# likelihood, and records the FMI those data realized for the coefficient.
# The print line that names a plan's missing data assumptions is built here
# (.aipe_missing_phrase()) and called from .dmar_tbl_footer.default() in
# R/dmar_output_helpers.R, so the console, knit_print(), and as_kable() all
# show it.
#
# The arithmetic. The FMI of a parameter is one minus the ratio of its
# complete data sampling variance to its sampling variance with the data
# actually observed (Savalei & Rhemtulla, 2012, Eq. 4), so at a fixed N the
# missing values multiply the sampling variance by 1 / (1 - fmi). The
# reliability planner's sampling variance is V / (N - 1), which makes its
# information (N - 1) / V; missing information leaves (N - 1)(1 - fmi) / V.
# A plan whose complete data answer is N_c therefore needs
#
#   N = (N_c - 1) / (1 - fmi) + 1,
#
# and a sample of N under missing data carries the information of a complete
# sample of the effective size (N - 1)(1 - fmi) + 1. The factor scales
# N - 1, not N, because N - 1 is what the planner's information is
# proportional to. The target sample size of Savalei and Rhemtulla, N / (1 -
# fmi), is the same idea for an information proportional to N.


# .aipe_check_fmi(fmi)
#
# Refuses anything but a single finite number in [0, 1). An FMI of 1 would
# mean the data carry no information about the coefficient, for which no N
# is enough, and a negative FMI has no planning interpretation.
.aipe_check_fmi <- function(fmi) {
  if (!is.numeric(fmi) || length(fmi) != 1L || !is.finite(fmi) ||
      fmi < 0 || fmi >= 1) {
    .dmar_refuse_argument(c(
      "'fmi' must be a single number in [0, 1): the fraction of missing ",
      "information anticipated for the reliability coefficient, with 0 ",
      "(the default) for complete data. It is not the proportion of ",
      "missing values."), arg = "fmi")
  }
  invisible(fmi)
}


# .aipe_check_missing_prop(missing_prop)
#
# The probability with which each simulated item score is deleted. A
# probability of 1 deletes everything, so [0, 1) is the admissible range.
.aipe_check_missing_prop <- function(missing_prop) {
  if (!is.numeric(missing_prop) || length(missing_prop) != 1L ||
      !is.finite(missing_prop) || missing_prop < 0 || missing_prop >= 1) {
    .dmar_refuse_argument(c(
      "'missing_prop' must be a single number in [0, 1): the ",
      "probability that each simulated item score is missing, with 0 ",
      "(the default) for complete data."), arg = "missing_prop")
  }
  invisible(missing_prop)
}


# .aipe_fmi_inflate(N, fmi)
#
# The sample size under missing data that carries the information of a
# complete sample of size N: (N - 1) / (1 - fmi) + 1. Unrounded; the caller
# takes .aipe_ceiling() of it. Callers apply it only when fmi > 0, so a
# complete data plan never passes through the arithmetic.
.aipe_fmi_inflate <- function(N, fmi) {
  (N - 1) / (1 - fmi) + 1
}


# .aipe_fmi_effective(N, fmi)
#
# The inverse map: the complete data size whose information a sample of N
# carries under missing data, (N - 1)(1 - fmi) + 1. Unrounded.
.aipe_fmi_effective <- function(N, fmi) {
  (N - 1) * (1 - fmi) + 1
}


# .aipe_ceiling(x)
#
# ceiling() for the two maps above. Their value is often an integer in exact
# arithmetic (a complete data size of 50 and an FMI of .3 give exactly 71),
# but 1 - fmi is not exact in binary (1 - .3 is 0.7000000000000001 and
# 1 - .8 is 0.19999999999999996), so the computed value can land a rounding
# error above that integer, and a plain ceiling() would add a case: for
# complete data sizes 2 to 400 it does so 13 times at an FMI of .3 and every
# time at .8 and .9. A value within a relative sqrt(.Machine$double.eps) of
# the integer below it is taken to be that integer; a value further above
# it is a genuine fraction of a case and is rounded up as usual.
.aipe_ceiling <- function(x) {
  ceiling(x - sqrt(.Machine$double.eps) * abs(x))
}


# .mcar_delete(Y, missing_prop)
#
# Deletes each entry of the matrix Y independently with probability
# missing_prop, the missing completely at random (MCAR) mechanism. One
# uniform draw per entry, in column-major order, so the deletion pattern is
# reproducible under a seed. A row that loses every entry is kept here; the
# FIML analysis drops it, as lavaan does, because it carries no information.
# Called only when missing_prop > 0, so a complete data simulation draws
# nothing extra and its random stream is unchanged.
.mcar_delete <- function(Y, missing_prop) {
  Y[stats::runif(length(Y)) < missing_prop] <- NA
  Y
}


# .aipe_realized_fmi(data, estimator)
#
# The fraction of missing information that one simulated data set realized
# for the coefficient ss_aipe_reliability_sensitivity() estimates: the
# unstructured FMI of fraction_missing_information(), computed by .fmi_core()
# in R/missing_information_internals.R, for the defined parameter that is the
# coefficient in the FIML fit behind it. For alpha that fit is the saturated
# model of the items whose covariance matrix reliability_alpha(missing =
# "fiml") applies the classical coefficient to, with its `alpha :=`
# definition, the fit .fiml_cov() in R/reliability_internals.R returns (at
# the EM estimates); in a saturated model the three FMI estimators coincide
# (Chen & Savalei, 2021, p. 9). For omega it is the single-factor model of
# .omega_fit_cfa(), whose `omega_f1` is the model implied coefficient
# reliability_omega(denominator = "model_implied") reports. The FMI is a
# normal theory quantity at the maximum likelihood estimates, which MLR
# shares, so both fits are ML with no standard errors of their own:
# .fmi_core() makes the two refits it needs. The warnings either fit can
# raise about the simulated data repeat those of the reliability fit and are
# not wanted once per replication.
#
# The FMI is returned relative to all nrow(data) simulated cases, the N the
# planner's fmi refers to. lavaan sets aside a case with no observed item,
# and .fmi_core() gives the FMI relative to the n_used cases the fit kept,
# fmi_used; the information the complete sample of n cases would have
# carried is n / n_used times that of the n_used cases, so relative to n the
# FMI is 1 - (n_used / n)(1 - fmi_used). With few items and much missing
# data the difference is material (12.5 percent of the cases lose all three
# items at missing_prop = .5), and carrying fmi_used into a plan would
# understate the information the plan loses.
#
# Returns NA_real_ when the fit fails, does not converge, or the FMI cannot be
# computed (.fmi_core() stops when lavaan will not refit an improper solution
# at its estimates). Draws nothing from the random number generator, so the
# simulation's stream is the one it would be without this step.
.aipe_realized_fmi <- function(data, estimator) {
  data <- as.data.frame(data)
  fit <- tryCatch({
    if (estimator == "alpha") {
      suppressWarnings(.fiml_cov(data, items = names(data), fit = TRUE)$fit)
    } else {
      suppressWarnings(.omega_fit_cfa(data = data, estimator = "ML",
                                      se = "none", missing = "ml")$fit)
    }
  }, error = function(e) NULL)
  if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) {
    return(NA_real_)
  }
  pt <- tryCatch(suppressWarnings(.fmi_core(fit, "unstructured")),
                 error = function(e) NULL)
  if (is.null(pt)) return(NA_real_)
  name <- if (estimator == "alpha") "alpha" else "omega_f1"
  out <- pt$fmi[pt$op == ":=" & pt$lhs == name]
  if (length(out) != 1L || !is.finite(out)) return(NA_real_)
  n_used <- sum(lavaan::lavInspect(fit, "nobs"))
  1 - (n_used / nrow(data)) * (1 - out)
}


# .aipe_missing_phrase(x)
#
# The line written under a result that was planned or simulated under
# missing data, or NULL when the result carries no such record. Its caller
# is .dmar_tbl_footer.default() in R/dmar_output_helpers.R, through which
# print.dmar_tbl(), knit_print(), and as_kable() all reach it. Two shapes
# are recognized:
#
#   - ss_aipe_reliability() with fmi > 0 carries `fmi` and
#     `necessary_N_complete` (for a plan made for an assurance, the complete
#     data size at which the assurance was simulated);
#   - ss_aipe_reliability_sensitivity() with fmi > 0 or missing_prop > 0
#     carries `fmi` (NA when the sample size was specified rather than
#     planned) and `missing_prop`, and with missing_prop > 0 also
#     `mean_realized_fmi` and `n_realized_fmi`, the mean FMI the simulated
#     data realized for the coefficient and the number of replications it
#     averages, which add a second line.
#
# The numbers are formatted for reading only; the attributes keep full
# precision. The result is one string; the realized FMI, when there is one,
# starts a new line within it.
.aipe_missing_phrase <- function(x) {
  num <- function(v) trimws(formatC(v, digits = 3, format = "fg"))
  mp <- attr(x, "missing_prop", exact = TRUE)
  fmi <- attr(x, "fmi", exact = TRUE)
  if (is.numeric(mp) && length(mp) == 1L && is.numeric(fmi) &&
      length(fmi) == 1L) {
    planned <- if (is.na(fmi)) {
      "Evaluated at the specified N"
    } else if (fmi == 0) {
      "Planned for complete data"
    } else {
      paste0("Planned for a fraction of missing information of ", num(fmi))
    }
    simulated <- if (mp == 0) {
      "simulated with complete data"
    } else {
      paste0("simulated with each item score missing completely at random ",
             "with probability ", num(mp), " and analyzed by full ",
             "information maximum likelihood")
    }
    line <- paste0(planned, "; ", simulated, ".")
    realized <- attr(x, "mean_realized_fmi")
    n_realized <- attr(x, "n_realized_fmi")
    if (is.numeric(realized) && length(realized) == 1L &&
        is.numeric(n_realized) && length(n_realized) == 1L) {
      line <- paste0(line, "\n", if (n_realized == 0) {
        paste0("The coefficient's fraction of missing information could ",
               "not be computed in any replication.")
      } else {
        paste0("The fraction of missing information the simulated data ",
               "realized for the coefficient averaged ", num(realized),
               " over ", if (n_realized == 1) "the 1 replication" else
                 paste0("the ", format(n_realized, scientific = FALSE),
                        " replications"),
               " in which it was computed.")
      })
    }
    return(line)
  }
  n_c <- attr(x, "necessary_N_complete", exact = TRUE)
  if (is.numeric(n_c) && length(n_c) == 1L && is.numeric(fmi) &&
      length(fmi) == 1L && !is.na(fmi)) {
    line <- paste0("Planned for a fraction of missing information of ",
                   num(fmi), " for the coefficient; with complete data, N = ",
                   format(n_c, scientific = FALSE))
    if ("assurance" %in% x$term) {
      line <- paste0(line, ", the size at which the assurance was simulated")
    }
    return(paste0(line, "."))
  }
  NULL
}
