# Shared internals of the SEM fit family. The first part serves
# sem_fit_measures(), equivalence_rmsea(), and sem_residual_fit(); the
# information criterion helpers after it serve sem_information_criteria()
# and compare_sem(). None of them is exported.

# ---- validation ---------------------------------------------------------

.sem_scalar <- function(x, name, lower = -Inf, strict = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict) x <= lower else x < lower) ||
      (integer && x != floor(x))) {
    .dmar_refuse_argument(c(
      "'", name, "' must be a finite ", if (integer) "integer " else "number ",
      if (strict) "greater than " else "at least ", lower, "."), arg = NULL)
  }
  invisible(x)
}

.sem_probability <- function(x, name) {
  .sem_scalar(x, name, 0, strict = TRUE)
  if (x >= 1) .dmar_refuse_argument(c(
    "'", name, "' must be below 1."), arg = NULL)
}

# A validated scalar as a bare number. A value read from
# lavaan::fitMeasures() carries its name ("chisq", "df", ...), and a name
# carried into c() would rename the row it lands in.
.sem_number <- function(x) as.numeric(x)[1L]

# The checks every fitted model passes before its fit is evaluated: lavaan
# is installed, the object is a lavaan fit, the optimizer converged, and the
# solution is admissible (.sem_admissible(): lavaan's post-estimation check,
# which finds no negative variance estimate and no covariance matrix of the
# latent variables, of the residuals of the observed variables (theta), or
# of the residuals of the structural part of the model (psi) that is not
# positive definite; with auxiliary variables, the same check of the
# model's own parameters). Fit measures of an improper solution describe
# estimates outside the parameter space and are withheld. A fit that did
# not converge stops with an error of class dmar_not_converged, the class
# fa_confirmatory() gives such a fit, and an inadmissible one with a
# refusal of class dmar_inadmissible (.sem_refuse_inadmissible()) that
# names what failed. The two checks before them are argument checks and
# stop plainly. `label` names the object in the messages, so a caller that
# checks the elements of a list can say which element failed. Returns the
# fit's lavaan options.
.sem_fit_check <- function(fit, label = "'fit'") {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required for a fitted model. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  if (!inherits(fit, "lavaan")) .dmar_refuse_argument(c(
    label, " must be a lavaan model."), arg = NULL)
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    .dmar_refuse(
      paste0("The model in ", label, " did not converge."),
      "not_converged")
  }
  admissible <- .sem_admissible(fit)
  if (!admissible) {
    .sem_refuse_inadmissible(attr(admissible, "problems"),
                             attr(admissible, "auxiliary"), label)
  }
  invisible(lavaan::lavInspect(fit, "options"))
}

# ---- admissibility ------------------------------------------------------

# The auxiliary variables of a FIML fit that lavaan models as saturated
# correlates (its aux argument with missing = "ml"): observed variables of
# the fit's matrices, with residual variances and covariances of their own,
# that the model's parameters do not otherwise involve. lavaan 0.7-3 keeps
# their names among the fit's options (and in its call); the auxiliary
# names of its data slot (@Data@ov.names.aux) are those of a two-stage fit,
# whose auxiliary variables enter the first stage only and are not among
# the observed variables of the model. Only names among the fit's observed
# variables are returned, so a two-stage fit, and a fit without auxiliary
# variables, has none. A bootstrap refit made through lavaan's slots
# carries the options of the fit it copies, and so its auxiliary variables.
# Never stops: an element or slot that cannot be read counts as empty, and
# the call contributes a name only when it holds the names themselves.
.sem_aux_names <- function(fit) {
  read <- function(f) {
    x <- tryCatch(f(), error = function(e) NULL)
    x <- unlist(x, use.names = FALSE)
    if (is.character(x)) x else character(0)
  }
  aux <- unique(c(read(function() lavaan::lavInspect(fit, "options")$aux),
                  read(function() fit@Data@ov.names.aux),
                  read(function() lavaan::lavInspect(fit, "call")$aux)))
  intersect(aux, read(function() lavaan::lavNames(fit, "ov")))
}

# The SRMR of a fit with auxiliary variables, over the model's own observed
# variables. lavaan's srmr_bentler_nomean averages over every observed
# variable of the fit, the auxiliary variables among them, and an auxiliary
# variable is a saturated correlate whose variance and covariances are
# reproduced exactly or nearly so: each one adds residuals of about zero to
# the average and lowers the SRMR of the same model (by about 8 percent for
# one auxiliary variable and 15 percent for two in a nine-item model). The
# residuals here are lavaan's own (lavResiduals(), type "cor.bentler", the
# covariance residuals over the product of the two sample standard
# deviations, about the h1 moments under FIML), with the rows and columns of
# the auxiliary variables `aux` set aside; the groups are combined as lavaan
# combines them, weighted by their numbers of cases. With no auxiliary
# variable this is srmr_bentler_nomean.
.sem_srmr_own <- function(fit, aux) {
  r <- lavaan::lavResiduals(fit, type = "cor.bentler", se = FALSE,
                            zstat = FALSE, summary = FALSE)
  if (lavaan::lavInspect(fit, "ngroups") == 1L) r <- list(r)
  nobs <- unlist(lavaan::lavInspect(fit, "nobs"))
  srmr <- vapply(r, function(g) {
    m <- if (!is.null(g$cov)) g$cov else g$res.cov
    keep <- setdiff(rownames(m), aux)
    m <- m[keep, keep, drop = FALSE]
    sqrt(mean(m[lower.tri(m, diag = TRUE)]^2))
  }, numeric(1L))
  sum(nobs * srmr) / sum(nobs)
}

# Whether a matrix passes lavaan's test of positive (semi)definiteness in
# its post-estimation check (lav_object_post_check_psd(), lavaan 0.7-3): no
# eigenvalue below -eps^(3/4) times the larger of 1 and the largest
# absolute eigenvalue, eps being the machine epsilon.
.sem_psd <- function(m) {
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  !any(ev < -.Machine$double.eps^(3 / 4) * max(1, abs(ev[1L])))
}

# Whether the solution of a lavaan fit is admissible. Without auxiliary
# variables among the fit's matrices the verdict is lavaan's own
# post-estimation check (lavTech(fit, "post.check"), lav_object_post_check()),
# exactly as lavaan returns it, warnings included. With them (`aux`, read by
# .sem_aux_names()) it is the same check of the model's own parameters,
# .sem_inadmissible_conditions() with the auxiliary variables set aside:
# lavaan's check also tests the residual covariance matrix of the
# observed variables (theta) with the auxiliary variables in it, and in the
# saturated correlates parameterization an auxiliary variable's covariances
# with the items run entirely through their residuals, so that block is
# not positive definite once the auxiliary variable is related to a factor
# strongly enough, while the model's own parameters are proper. For one
# factor whose items have uncorrelated residuals, with the factor and the
# auxiliary variable in standard units, the auxiliary variable's row of
# theta is its correlation rho with the factor times the loadings, and the
# block is positive semidefinite only while rho^2 is at most (1 - H) / H,
# H being the maximal reliability of the items (coefficient H): a
# correlation above .33 is enough when H is .9, and a useful auxiliary
# variable is chosen for such a relation. On six of the
# holzinger_swineford tests with the cube and
# sentence scores missing at random given the general information score,
# that score as an auxiliary variable has residual covariances in theta
# with the three verbal tests that amount to correlations of 1.2 to 1.6,
# so lavaan's check fails a fit whose own estimates differ from those of
# the fit without the auxiliary variable by at most .003 in a
# standardized loading.
#
# Returns TRUE or FALSE with two attributes: `problems`, the data frame of
# .sem_inadmissible_conditions() (for a fit without auxiliary variables
# computed only when lavaan's check fails and `explain` is TRUE, and NULL
# otherwise), and `auxiliary`, the auxiliary variables set aside. Errors
# are not caught; the bootstrap counts a refit whose check fails to run as
# improper (.bs_refit_one()).
.sem_admissible <- function(fit, aux = .sem_aux_names(fit), explain = TRUE) {
  problems <- NULL
  if (length(aux)) {
    problems <- .sem_inadmissible_conditions(fit, aux)
    ok <- nrow(problems) == 0L
  } else {
    ok <- isTRUE(lavaan::lavTech(fit, "post.check"))
    if (!ok && explain) problems <- .sem_inadmissible_conditions(fit)
  }
  structure(ok, problems = problems, auxiliary = aux)
}

# What lavaan's post-estimation check (lav_object_post_check(), lavaan
# 0.7-3) finds wrong with a solution, with the auxiliary variables `aux` set
# aside: a data frame with one row per failed condition, in the order the
# check makes them, and columns `reason`, `group` (the label of the group
# where it failed; NA for a fit of one group), `level` (the level where it
# failed, 1 or 2, in a multilevel fit, whose blocks lavaan numbers level
# within group; NA for a fit of one level), `variable` (the variable whose
# variance estimate is negative; NA for a matrix), and `value` (that
# estimate, or the smallest eigenvalue of the matrix). The reasons, and the
# conditions under which lavaan makes each check, are lavaan's:
#   "negative_variance" a variance estimate below zero, of an observed and
#                       then of a latent variable, read from the parameter
#                       table; an NA among the variances of either kind
#                       leaves the matrices below unchecked (and those of
#                       that kind unjudged);
#   "cov.lv"            the covariance matrix of the latent variables is not
#                       positive definite, checked when no latent variance
#                       is negative and the model has latent variables;
#   "theta"             the residual covariance matrix of the continuous
#                       observed variables is not, checked when no observed
#                       variable's variance is negative;
#   "psi"               the residual covariance matrix of the structural
#                       part of the model is not, checked in the LISREL
#                       representation when no variance is negative and the
#                       block's cov.lv and theta passed.
# A matrix fails .sem_psd(). With `aux` empty the rows are exactly what
# lavaan's check flags, read with lavaan's own accessors (the parameter
# table, lavTech()'s cov.lv and theta, and the psi matrices of the model);
# with auxiliary variables their rows and columns are dropped from each
# matrix and their variances from the parameter table first. A fit whose
# matrices lavaan does not label cannot have its auxiliary variables set
# aside, and stops.
.sem_inadmissible_conditions <- function(fit, aux = character(0)) {
  pt <- fit@ParTable
  model <- fit@Model
  nblocks <- model@nblocks
  # The group and the level of each block: lavaan numbers the blocks level
  # within group, so block b is level (b - 1) %% L + 1 of group
  # (b - 1) %/% L + 1 for L levels. Should the blocks not be the groups
  # times the levels, each is named by its number.
  ngroups <- lavaan::lavInspect(fit, "ngroups")
  nlevels <- lavaan::lavInspect(fit, "nlevels")
  group_labels <- if (ngroups > 1L) {
    as.character(lavaan::lavInspect(fit, "group.label"))
  } else NA_character_
  if (nblocks == ngroups * nlevels && length(group_labels) == ngroups) {
    block_group <- rep(group_labels, each = nlevels)
    block_level <- if (nlevels > 1L) {
      rep(seq_len(nlevels), times = ngroups)
    } else rep(NA_integer_, nblocks)
  } else {
    block_group <- if (nblocks > 1L) as.character(seq_len(nblocks)) else
      NA_character_
    block_level <- rep(NA_integer_, nblocks)
  }
  drop_aux <- function(m) {
    if (!length(aux) || !nrow(m)) return(m)
    if (is.null(rownames(m))) {
      .dmar_refuse_argument(c(
        "The fit's matrices carry no variable names, so its auxiliary ",
        "variables cannot be set aside in the check of admissibility."),
        arg = NULL)
    }
    keep <- !rownames(m) %in% aux
    m[keep, keep, drop = FALSE]
  }
  rows <- list()
  add <- function(reason, block, variable = NA_character_, value) {
    data.frame(reason = reason, group = block_group[block],
               level = block_level[block], variable = variable,
               value = value, stringsAsFactors = FALSE)
  }
  # The variances, of the observed variables and then of the latent ones:
  # the rows of each that are negative, or NULL when one of them is NA.
  is_var <- pt$op == "~~" & pt$lhs == pt$rhs
  block_of <- if (!is.null(pt$block)) pt$block else
    if (!is.null(pt$group)) pt$group else rep(1L, length(pt$lhs))
  negative <- function(names) {
    idx <- which(is_var & pt$lhs %in% names)
    if (anyNA(pt$est[idx])) return(NULL)
    idx[pt$est[idx] < 0]
  }
  neg_ov <- negative(setdiff(lavaan::lavNames(fit, "ov"), aux))
  neg_lv <- negative(setdiff(lavaan::lavNames(fit, "lv"), aux))
  var_na <- is.null(neg_ov) || is.null(neg_lv)
  for (i in c(neg_ov, neg_lv)) {
    rows[[length(rows) + 1L]] <- add("negative_variance", block_of[i],
                                     pt$lhs[i], pt$est[i])
  }
  var_ov_ok <- !length(neg_ov)
  var_lv_ok <- !length(neg_lv)
  # The covariance matrix of the latent variables.
  veta_ok <- rep(TRUE, nblocks)
  if (!var_na && var_lv_ok && length(lavaan::lavNames(fit, "lv.regular"))) {
    eta <- lavaan::lavTech(fit, "cov.lv", add.labels = TRUE)
    for (b in seq_len(nblocks)) {
      m <- drop_aux(eta[[b]])
      if (!nrow(m)) next
      if (!.sem_psd(m)) {
        veta_ok[b] <- FALSE
        rows[[length(rows) + 1L]] <- add("cov.lv", b, value = min(
          eigen(m, symmetric = TRUE, only.values = TRUE)$values))
      }
    }
  }
  # The residual covariance matrix of the continuous observed variables.
  theta_ok <- rep(TRUE, nblocks)
  if (!var_na && var_ov_ok) {
    theta <- lavaan::lavTech(fit, "theta", add.labels = TRUE)
    for (b in seq_len(nblocks)) {
      num <- model@num.idx[[b]]
      if (!length(num)) next
      m <- drop_aux(theta[[b]][num, num, drop = FALSE])
      if (!nrow(m)) next
      if (!.sem_psd(m)) {
        theta_ok[b] <- FALSE
        rows[[length(rows) + 1L]] <- add("theta", b, value = min(
          eigen(m, symmetric = TRUE, only.values = TRUE)$values))
      }
    }
  }
  # The residual covariance matrix of the structural part.
  if (!var_na && var_ov_ok && var_lv_ok &&
      identical(model@representation, "LISREL")) {
    psi_idx <- which(names(model@GLIST) == "psi")
    for (b in seq_len(nblocks)) {
      if (!veta_ok[b] || !theta_ok[b] || b > length(psi_idx)) next
      m <- model@GLIST[[psi_idx[b]]]
      if (!nrow(m) || anyNA(m)) next
      dimnames(m) <- model@dimNames[[psi_idx[b]]]
      m <- drop_aux(m)
      if (!nrow(m)) next
      if (!.sem_psd(m)) {
        rows[[length(rows) + 1L]] <- add("psi", b, value = min(
          eigen(m, symmetric = TRUE, only.values = TRUE)$values))
      }
    }
  }
  if (!length(rows)) return(.sem_inadmissible_conditions_empty())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The data frame of .sem_inadmissible_conditions() when nothing fails.
.sem_inadmissible_conditions_empty <- function() {
  data.frame(reason = character(0), group = character(0),
             level = integer(0), variable = character(0),
             value = numeric(0), stringsAsFactors = FALSE)
}

# The refusal of an inadmissible solution: an error of class
# dmar_inadmissible (a dmar_refusal; .dmar_refuse(), R/dmar_conditions.R)
# whose message names what failed, from the rows of
# .sem_inadmissible_conditions(). The condition carries `reason` (the
# reasons that failed, each once, in the order of the check), `variables`
# (the variables with a negative variance estimate, empty when none),
# `problems` (the rows themselves), and `auxiliary` (the auxiliary
# variables set aside, empty when none). The message places each failure
# in its group and, in a multilevel fit, at its level ("at level 2"). A
# fit lavaan's check fails with nothing found here, which the check above
# rules out, is refused in lavaan's general terms.
.sem_refuse_inadmissible <- function(problems, aux, label) {
  if (is.null(problems)) problems <- .sem_inadmissible_conditions_empty()
  if (is.null(problems$level)) problems$level <- rep(NA_integer_, nrow(problems))
  where <- function(g, l) {
    paste0(ifelse(is.na(g), "", paste0(" in group ", g)),
           ifelse(is.na(l), "", paste0(" at level ", l)))
  }
  and <- function(x) {
    if (length(x) < 2L) return(x)
    if (length(x) == 2L) return(paste(x, collapse = " and "))
    paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  }
  parts <- character(0)
  neg <- problems[problems$reason == "negative_variance", , drop = FALSE]
  if (nrow(neg)) {
    parts <- c(parts, paste0(
      if (nrow(neg) == 1L) "a negative variance estimate for "
      else "negative variance estimates for ",
      and(paste0(neg$variable, where(neg$group, neg$level)))))
  }
  what <- c(cov.lv = "the latent variables (cov.lv)",
            theta = "the residuals of the observed variables (theta)",
            psi = "the residuals of the structural part of the model (psi)")
  for (r in names(what)) {
    hit <- problems[problems$reason == r, , drop = FALSE]
    if (!nrow(hit)) next
    at <- if (all(is.na(hit$level))) {
      groups <- hit$group[!is.na(hit$group)]
      if (length(groups)) {
        paste0(" in group", if (length(groups) > 1L) "s", " ", and(groups))
      } else ""
    } else {
      paste0(" ", and(sub("^ ", "", where(hit$group, hit$level))))
    }
    parts <- c(parts, paste0(
      "a covariance matrix of ", what[[r]], " that is not positive definite",
      at))
  }
  found <- if (length(parts)) {
    paste0(" finds ", paste(parts, collapse = "; "))
  } else " flags it"
  check <- if (length(aux)) {
    paste0("lavaan's post-estimation check of the model's own parameters, ",
           "with the auxiliary variable", if (length(aux) > 1L) "s", " ",
           and(aux), " set aside,")
  } else "lavaan's post-estimation check"
  .dmar_refuse(c("The model in ", label, " has an inadmissible solution (",
                 check, found, "); resolve it before evaluating fit."),
               "inadmissible", reason = unique(problems$reason),
               variables = unique(neg$variable), problems = problems,
               auxiliary = aux)
}

# lavResiduals() gained the unbiased SRMR interval and close fit test
# (usrmr_ci_level, usrmr_close_h0, h1_acov) in lavaan 0.7-2, and
# requireNamespace() does not check versions. The installed version is an
# argument so that the refusal can be tested.
.sem_require_lavaan_version <- function(fn, lavaan_version = utils::packageVersion("lavaan"),
                                        lavaan_floor = "0.7-2") {
  if (as.numeric_version(lavaan_version) < as.numeric_version(lavaan_floor)) {
    .dmar_refuse(c(
      fn, "() requires lavaan ", lavaan_floor, " or newer (installed: ",
      as.character(lavaan_version), "), whose lavResiduals() provides the ",
      "unbiased SRMR interval and close fit test. Update it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  invisible(TRUE)
}

# ---- the sample size multiplier -----------------------------------------

# The multiplier a in T = a F that turns the fit function into the test
# statistic and the population discrepancy into the noncentrality,
# lambda = a df RMSEA^2 / G for G groups (Steiger, 1998). Under "N-1" each
# group contributes its size less one, so a = N - G; under "N", a = N.
.sem_multiplier <- function(N, groups, convention) {
  if (identical(convention, "N")) N else N - groups
}

# The convention a lavaan fit used to form its test statistic, as lavaan
# itself reads it for its fit measures: T = N F for ML under the normal
# likelihood (the default), and T = (N - G) F under likelihood = "wishart"
# and for the least squares estimators (GLS, WLS, ULS, DWLS).
.sem_lavaan_convention <- function(fit) {
  o <- lavaan::lavInspect(fit, "options")
  N <- lavaan::lavInspect(fit, "ntotal")
  G <- lavaan::lavInspect(fit, "ngroups")
  normal <- o$estimator %in% c("ML", "PML", "FML", "catML") &&
    o$likelihood %in% c("default", "normal")
  convention <- if (normal) "N" else "N-1"
  list(convention = convention, multiplier = .sem_multiplier(N, G, convention),
       N = N, groups = G, likelihood = o$likelihood, estimator = o$estimator)
}

# ---- RMSEA noncentrality limits -----------------------------------------

# Limits for the noncentrality of a chi square statistic, inverted with
# ci_nc_chisq(), under the conventions of Browne and Cudeck (1993) and
# lavaan. When the statistic lies at or below the upper (1 - conf_level)/2
# quantile of the central chi square, no positive noncentrality reaches the
# lower limit and it is 0. When it lies at or below the lower
# (1 - conf_level)/2 quantile, no noncentrality is consistent with it at
# this level; the upper limit is then 0, the T-size convention of Yuan,
# Chan, Marcoulides, and Bentler (2016, Eq. 4), which lavaan shares, so that
# a test of not close fit and the upper limit always agree. ci_nc_chisq()
# warns in both cases (it clamps the first and returns NA for the second);
# the caller states the convention instead, so those two warnings are
# muffled here and any other warning passes.
.sem_ncp_bounds <- function(chisq, df, conf_level) {
  nc <- withCallingHandlers(
    ci_nc_chisq(chisq = chisq, df = df, conf_level = conf_level,
                verbose = FALSE),
    warning = function(cond) {
      if (grepl("clamped to 0|upper noncentrality limit is undefined",
                conditionMessage(cond)))
        invokeRestart("muffleWarning")
    })
  lower <- nc$value[nc$term == "ci_lower"]
  upper <- nc$value[nc$term == "ci_upper"]
  if (!isTRUE(lower > 0)) lower <- 0
  if (is.na(upper)) upper <- 0
  c(lower = lower, upper = upper)
}

# RMSEA from a noncentrality (Steiger, 1998, Eq. 23): sqrt(G lambda / (a df)).
.sem_rmsea_from_ncp <- function(lambda, df, multiplier, groups) {
  sqrt(groups * pmax(lambda, 0) / (multiplier * df))
}

# The RMSEA rows of one version, computed from the statistic T on df degrees
# of freedom that the version inverts, formed with multiplier a in G groups:
# the point estimate, the limits, and the p-values of the tests of close and
# not close fit. With c_hat = 1 these are the standard formulas, applied to
# the standard statistic for the standard measures and to the scaled one for
# the scaled measures (Brosseau-Liard, Savalei, & Li, 2012, Eqs. 11 and 13;
# Savalei, 2018, Table 1). With c_hat the scaling constant of a scaled
# statistic they are the robust interval of Brosseau-Liard et al. (2012,
# Eq. 10), which inverts the scaled statistic and multiplies its
# noncentrality limits by c_hat (Savalei, 2018, for a scaled and shifted
# statistic), and the tests lavaan pairs with it. RMSEA^2 = k lambda, so a
# bound epsilon has noncentrality epsilon^2 / k. The distribution function
# is .pchisq_nc(), accurate for a statistic of any size.
.sem_rmsea_rows <- function(chisq, df, multiplier, groups, conf_level,
                            close_fit, not_close_fit, c_hat = 1) {
  ncp <- .sem_ncp_bounds(chisq, df, conf_level)
  k <- c_hat * groups / (multiplier * df)
  c(rmsea = sqrt(k * max(chisq - df, 0)),
    rmsea_ci_lower = sqrt(k * ncp[["lower"]]),
    rmsea_ci_upper = sqrt(k * ncp[["upper"]]),
    p_close_fit = .pchisq_nc(chisq, df, close_fit^2 / k, lower.tail = FALSE),
    p_not_close_fit = .pchisq_nc(chisq, df, not_close_fit^2 / k))
}

# The noncentrality above which ci_nc_chisq() leaves R's pchisq(), which
# lavaan uses for its RMSEA limits and close fit tests (see
# R/ci_nc_chisq_internals.R). A fit whose statistic, or the noncentrality of
# a close fit bound, exceeds it has its limits and tests computed here.
.SEM_NCP_PCHISQ_LIMIT <- 1e5

# ---- restricted means ---------------------------------------------------

# Whether a lavaan fit restricts the means. D holds the derivatives of the
# model implied moments of every group, the means and the covariances,
# with respect to the parameters the constraints leave free
# (.sem_free_delta()). The parameter directions that leave every
# covariance unchanged move the means in rank(D) - rank(D_cov) independent
# directions, D_cov being the covariance rows of D, and the means are free
# exactly when that is the number of modeled means; the shortfall is the
# number of degrees of freedom the means contribute. (With fixed.x = TRUE
# the means of the exogenous observed variables are fixed at their sample
# values, so their rows of D are zero and they are not counted; a mean
# fixed by the model, such as an intercept fixed at a constant with the
# latent mean fixed too, is counted, and its zero row is a restriction.)
# A count of the free mean parameters, net of the equality constraints
# between two of them, missed a shared label under ceq.simple = TRUE (one
# parameter, no == row), an intercept constrained to a constant or in a
# linear constraint, an intercept tied to a loading, and an intercept held
# at the bound of an inequality; the rank sees each of them, and it does
# not depend on how a free mean structure is parameterized (an indicator's
# intercept fixed at zero with the latent mean free has as many directions
# as free intercepts). sem_fit_measures() reads it for gamma hat, which is
# defined for covariance structures, and for a ULS fit of continuous data,
# and measurement_invariance() for the rungs of a ULS ladder. Returns NA
# when lavaan's derivatives or constraints cannot be read, or when the
# fit is conditional on its exogenous covariates (conditional.x = TRUE),
# whose moments are those of the residuals.
.sem_means_restricted <- function(fit) {
  options <- lavaan::lavInspect(fit, "options")
  if (!isTRUE(options$meanstructure)) return(FALSE)
  if (isTRUE(options$conditional.x)) return(NA)
  D <- .sem_free_delta(fit)
  if (is.null(D)) return(NA)
  mean_rows <- endsWith(rownames(D), "~1")
  modeled <- sub("~1$", "", rownames(D)[mean_rows])
  if (isTRUE(options$fixed.x)) {
    modeled <- modeled[!modeled %in% lavaan::lavNames(fit, "ov.x")]
  }
  qr(D)$rank - qr(D[!mean_rows, , drop = FALSE])$rank < length(modeled)
}

# Whether the mean and the covariance parts of a lavaan fit share a
# parameter direction, read from the same D as .sem_means_restricted():
# the rows of D_cov and of D_mean (the covariance and the mean rows) span
# subspaces of the free parameter directions, and rank(D) = rank(D_cov) +
# rank(D_mean) exactly when the two subspaces meet only at zero, that is,
# when the free parameters split, locally, into directions that move only
# the covariances and directions that move only the means. A fit function
# that is a weighted sum of a covariance part and a mean part is then
# minimized by minimizing each part along its own directions, whatever the
# weights. Otherwise, as when the loadings enter the means through the
# latent means of a multiple group model (the scalar and strict rungs of an
# invariance ladder), or through a latent mean freed beside intercepts
# fixed at zero, the minimizer depends on the weights. Returns FALSE
# without a mean structure, and NA when the derivatives or constraints
# cannot be read or the fit is conditional on its exogenous covariates.
.sem_means_linked <- function(fit) {
  options <- lavaan::lavInspect(fit, "options")
  if (!isTRUE(options$meanstructure)) return(FALSE)
  if (isTRUE(options$conditional.x)) return(NA)
  D <- .sem_free_delta(fit)
  if (is.null(D)) return(NA)
  mean_rows <- endsWith(rownames(D), "~1")
  qr(D)$rank < qr(D[!mean_rows, , drop = FALSE])$rank +
    qr(D[mean_rows, , drop = FALSE])$rank
}

# Whether the estimates of a ULS fit of continuous data change with the
# units of the variables when every variable is multiplied by the same
# number k: TRUE, FALSE, or NA when that cannot be read. The ULS fit
# function with a mean structure is k^4 times its covariance part plus k^2
# times its mean part, in the rescaled parameters. When the model leaves
# the means free (.sem_means_restricted() FALSE) the mean part is zero at
# the minimum whatever the covariance parameters, and when the two parts
# share no parameter direction (.sem_means_linked() FALSE) each is
# minimized on its own; in both cases the estimates are equivariant (every
# loading multiplied by k, every variance by k^2). Otherwise the weights
# k^4 and k^2 trade the two parts off differently for every k, and the
# estimates are not equivariant. Browne's residual based statistic,
# N e'[G^-1 - G^-1 D (D'G^-1 D)^-1 D'G^-1] e (Browne, 1984), has no units
# only when e and D are evaluated at equivariant estimates, so this is also
# whether that statistic, and every index formed from it, changes with the
# units: on the scalar rung of the first example of
# measurement_invariance() fit by estimator = "ULS", a CFI of .8874 for the
# raw scores and .8861 with every score divided by 10.
.sem_uls_estimates_units <- function(fit) {
  restricted <- .sem_means_restricted(fit)
  if (isFALSE(restricted)) return(FALSE)
  linked <- .sem_means_linked(fit)
  if (isFALSE(linked)) return(FALSE)
  if (isTRUE(restricted) && isTRUE(linked)) TRUE else NA
}

# The derivatives of a fit's model implied moments, those of every group
# stacked, with respect to the parameters its constraints leave free:
# lavaan's Delta, whose columns are the unconstrained free parameters,
# times a basis of the directions in which every constraint that binds at
# the estimates holds. Under ceq.simple = TRUE, when every equality is a
# simple equality of two parameters, lavaan stores a basis K of those
# directions itself (its slot ceq.simple.K); otherwise, and for the other
# constraints, the basis spans the null space of the Jacobian of the
# constraints that bind (lavaan's con.jac, evaluated at the estimates): the
# == rows, and the inequality rows lavaan counts as active, the ones it
# keeps when it computes the standard errors (an inequality whose slack
# exceeds 1e-5 is inactive). An inequality held at its bound fixes its
# parameter there, as an equality would, although lavaan's degrees of
# freedom do not count it: on the three-factor model of nine
# holzinger_swineford tests with t1's intercept held above 100, the
# fitted mean misses the sample mean by 88 and the means are restricted.
# This is how .ts_constrained_delta() (R/two_stage_internals.R) reads the
# equality constraints for the residual based statistic, which refuses
# inequality constraints instead. Returns NULL when the derivatives or
# constraints cannot be read.
.sem_free_delta <- function(fit) {
  D <- tryCatch(lavaan::lavInspect(fit, "delta"), error = function(e) NULL)
  if (is.list(D) && !is.matrix(D)) {
    D <- tryCatch(do.call(rbind, lapply(D, unclass)), error = function(e) NULL)
  } else if (!is.null(D)) {
    D <- unclass(D)
  }
  if (!is.matrix(D) || is.null(rownames(D))) return(NULL)
  J <- tryCatch(fit@Model@con.jac, error = function(e) NULL)
  if (!is.matrix(J)) return(NULL)
  active <- c(attr(J, "ceq.idx"),
              setdiff(attr(J, "cin.idx"), attr(J, "inactive.idx")))
  J <- matrix(unclass(J)[active, , drop = FALSE], nrow = length(active))
  if (length(active) && ncol(J) != ncol(D)) return(NULL)
  if (isTRUE(tryCatch(fit@Model@ceq.simple.only, error = function(e) FALSE))) {
    K <- tryCatch(fit@Model@ceq.simple.K, error = function(e) NULL)
    if (!is.matrix(K) || nrow(K) != ncol(D)) return(NULL)
    D <- D %*% K
    J <- J %*% K
  }
  if (!nrow(J)) return(D)
  q <- qr(t(J))
  if (q$rank == 0L) return(D)
  D %*% qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
}

# Whether a fit with no degrees of freedom fits exactly, from its
# fitMeasures() `fm`: its df is zero and its statistic is zero to the
# accuracy of the optimizer (about 1e-12 on the fits checked: FIML, MLR,
# two-stage, MLMVS, ULSMV, two groups). Zero degrees of freedom alone does
# not make a model fit exactly. An inequality constraint that binds holds
# its parameter at the bound, as an equality would, though lavaan's count
# of the degrees of freedom does not include it, so the model can misfit
# with df = 0: one factor of three variables whose correlations are .8,
# .8, and .5 has a negative error variance, and with bounds = "pos.var"
# its statistic is 32.8 on 0 df (N = 200) and its CFI .925, not 1. So can a
# model that is not identified at its estimates. The tolerance, 1e-6 times
# the baseline statistic (or 1 when that is missing), is far above the
# optimizer's error and far below any misfit that would move the CFI:
# misclassifying a statistic that small moves the CFI by less than about
# 1e-6. A missing statistic with df = 0 counts as exact fit, as df = 0
# alone did.
.sem_exact_fit <- function(fm) {
  get <- function(nm) if (nm %in% names(fm)) as.numeric(fm[[nm]]) else NA_real_
  df <- get("df")
  if (!isTRUE(df == 0)) return(FALSE)
  chisq <- get("chisq")
  if (!is.finite(chisq)) return(TRUE)
  abs(chisq) <= 1e-6 * max(1, get("baseline.chisq"), na.rm = TRUE)
}

# Whether a fit has no degrees of freedom and does not fit exactly
# (.sem_exact_fit()), from its fitMeasures() `fm`. Its TLI (NNFI) and
# RMSEA, with the RMSEA's limits, divide by the degrees of freedom: with a
# statistic above zero they are undefined (T/df has no value, and they
# grow without bound as df falls to 0 with T fixed), not the values of
# exact fit, 1 and 0, that lavaan reports whenever df = 0. So is the exact
# fit p-value, which lavaan leaves NA. sem_fit_measures(), fa_confirmatory(), and
# measurement_invariance() report them as NA; the CFI, whose model term is
# T - c df = T, is formed from the statistic.
.sem_zero_df_misfit <- function(fm) {
  df <- if ("df" %in% names(fm)) as.numeric(fm[["df"]]) else NA_real_
  isTRUE(df == 0) && !.sem_exact_fit(fm)
}

# ---- robust measures lavaan does not compute ----------------------------

# The robust measures of a continuous ML fit of complete data whose
# corrected statistic is the older mean and variance adjusted one (test =
# "mean.var.adjusted", which estimator = "MLMVS" sets), for which lavaan
# 0.7-3 computes none (its lav_fit_rmsea_lavobject() and
# lav_fit_cfi_lavobject() make robust measures for the Satorra-Bentler,
# Yuan-Bentler, and scaled and shifted statistics only), under lavaan's
# names, so that sem_fit_measures(), fa_confirmatory(), and measurement_invariance()
# report the robust version for such a fit as for any other.
#
# The robust RMSEA and CFI of Brosseau-Liard, Savalei, and Li (2012) and
# Brosseau-Liard and Savalei (2014) subtract from the ML statistic T its
# asymptotic expectation under nonnormality, c df, where c = tr(U Gamma)/df,
# and so estimate the population values of the ML indices whichever
# corrected statistic accompanies them (Savalei, 2018, Eqs. 6 to 9, who
# shows this for the scaled and shifted statistic). The mean and variance
# adjusted statistic refers T/s to the chi square on df_s degrees of
# freedom, with s = tr((U Gamma)^2)/tr(U Gamma) and df_s = tr(U Gamma)^2 /
# tr((U Gamma)^2) (Satterthwaite's approximation), so s df_s = tr(U Gamma)
# and c = s df_s / df; the baseline model gives its own constant the same
# way. These are the constants of the Satorra-Bentler statistic of the
# same fit, so the values equal lavaan's robust measures for the model
# fitted with estimator = "MLM", whose information (expected) and U Gamma
# are the same; the tests check that. The interval and the tests of close
# and not close fit are those of Brosseau-Liard et al. (2012, Eq. 10),
# from the mean scaled statistic T/c on df with its noncentrality limits
# multiplied by c (.sem_rmsea_rows()); the TLI and RNI follow lavaan's
# robust formulas with the same constants.
#
# `conf_level`, `close_fit`, and `not_close_fit` are the level and the
# bounds of the interval and tests (lavaan's defaults, 0.90, 0.05, and
# 0.08, unless the caller passes its own; sem_fit_measures() and
# fa_confirmatory() pass theirs, and measurement_invariance(), which reads
# the point estimates alone, takes these). Returns a named numeric vector
# (cfi.robust, tli.robust, nnfi.robust, rni.robust, rmsea.robust,
# rmsea.ci.lower.robust, rmsea.ci.upper.robust, rmsea.pvalue.robust,
# rmsea.notclose.pvalue.robust), empty for any other fit. A fit with no
# degrees of freedom that fits exactly (.sem_exact_fit()), where c is 0/0,
# gets the values every version gives a
# model that fits exactly: its statistic is zero, and so is the model term
# T - c df of the robust CFI and RNI, which are 1; its TLI divides by df
# and is NA; and its RMSEA and limits take lavaan's 0 (with no test of
# close fit), as lavaan's robust RMSEA of the same model fitted by MLM
# does. Returning these rather than nothing keeps a saturated model or
# rung from changing the version its caller reports (a ladder of
# measurement_invariance() whose configural rung is saturated would
# otherwise report the standard indices at every rung).
.sem_robust_mva <- function(fit, conf_level = 0.90, close_fit = 0.05,
                            not_close_fit = 0.08) {
  o <- lavaan::lavInspect(fit, "options")
  corrected <- .lavaan_corrected_tests(o$test)
  if (!identical(o$estimator, "ML") || !length(corrected) ||
      corrected[[1L]] != "mean.var.adjusted" ||
      isTRUE(o$missing %in% c("ml", "ml.x", "two.stage", "robust.two.stage")) ||
      length(lavaan::lavNames(fit, "ov.ord")) > 0L ||
      lavaan::lavInspect(fit, "nlevels") > 1L) {
    return(numeric(0))
  }
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "chisq.scaling.factor",
                                   "df.scaled", "baseline.chisq",
                                   "baseline.df",
                                   "baseline.chisq.scaling.factor",
                                   "baseline.df.scaled", "rmsea.robust"))
  get <- function(nm) if (nm %in% names(fm)) as.numeric(fm[[nm]]) else NA_real_
  T <- get("chisq")
  df <- get("df")
  if ("rmsea.robust" %in% names(fm) || !is.finite(df) || df < 0 ||
      !is.finite(T)) {
    return(numeric(0))
  }
  if (df == 0) {
    # A model with no degrees of freedom that does not fit exactly (an
    # inequality bound that binds) has no constant c to correct by, and
    # gets no robust measures, as any fit whose constant is undefined.
    if (!.sem_exact_fit(fm)) return(numeric(0))
    return(c(cfi.robust = 1, tli.robust = NA_real_, nnfi.robust = NA_real_,
             rni.robust = 1, rmsea.robust = 0, rmsea.ci.lower.robust = 0,
             rmsea.ci.upper.robust = 0, rmsea.pvalue.robust = NA_real_,
             rmsea.notclose.pvalue.robust = NA_real_))
  }
  c_hat <- get("chisq.scaling.factor") * get("df.scaled") / df
  if (!is.finite(c_hat) || c_hat <= 0) return(numeric(0))
  conv <- .sem_lavaan_convention(fit)
  rows <- .sem_rmsea_rows(T / c_hat, df, conv$multiplier, conv$groups,
                          conf_level, close_fit, not_close_fit, c_hat = c_hat)
  T_b <- get("baseline.chisq")
  df_b <- get("baseline.df")
  c_b <- get("baseline.chisq.scaling.factor") * get("baseline.df.scaled") / df_b
  num <- T - c_hat * df
  den <- T_b - c_b * df_b
  incremental <- if (all(is.finite(c(num, den)))) {
    t2 <- max(num, den, 0)
    c(cfi = if (isTRUE(all.equal(t2, 0))) 1 else 1 - max(num, 0) / t2,
      tli = if (den != 0) 1 - (num * df_b) / (den * df) else NA_real_,
      rni = if (den != 0) 1 - num / den else NA_real_)
  } else c(cfi = NA_real_, tli = NA_real_, rni = NA_real_)
  out <- c(cfi.robust = incremental[["cfi"]], tli.robust = incremental[["tli"]],
           nnfi.robust = incremental[["tli"]], rni.robust = incremental[["rni"]],
           rmsea.robust = rows[["rmsea"]],
           rmsea.ci.lower.robust = rows[["rmsea_ci_lower"]],
           rmsea.ci.upper.robust = rows[["rmsea_ci_upper"]],
           rmsea.pvalue.robust = rows[["p_close_fit"]],
           rmsea.notclose.pvalue.robust = rows[["p_not_close_fit"]])
  attr(out, "scaling") <- c(model = c_hat, baseline = c_b)
  out
}

# Whether a lavaan fit is a ULS fit of continuous data (estimator =
# "ULS", "ULSM", or "ULSMV" without ordered variables). The ULS fit
# function, half the sum of the squared residual covariances (with the
# squared mean residuals added under a mean structure), weights every
# moment alike, so its value, the statistic T = N F, and the scaling
# constant c of a corrected statistic are in the units of the squared
# covariances, and ULS is not scale invariant: rescaling some variables
# and not others changes its estimates and every measure of its fit.
# sem_fit_measures() and measurement_invariance() report a measure formed
# from T only where it does not change when every variable is multiplied
# by the same number k, as a measure of fit must not:
#   - lavaan's robust RMSEA (lav_fit_rmsea_lavobject() forms it for any
#     continuous estimator with a Satorra-Bentler or scaled and shifted
#     test as sqrt(max(T - c df, 0)/(N df)), the correction of
#     Brosseau-Liard, Savalei, and Li, 2012, which rests on the ML fit
#     function and which Savalei, 2018, cautions does not carry over to
#     other estimators) is multiplied by k^2: 16.1 for the three-factor
#     model of nine holzinger_swineford tests fitted by ULSM to the scores
#     and 0.161 to the scores divided by 10. The standard RMSEA, formed from
#     T itself, moves the same way (17.1 and 0.161). Both are withheld.
#   - the standard CFI, TLI, RNI, and IFI set the degrees of freedom, a
#     number, against T (T - df), and change with k too (a CFI of .9767
#     and .9793); the NFI, 1 - T/T_b, does not. The first four are
#     withheld.
#   - the robust CFI, TLI, and RNI, ratios of T - c df to its baseline
#     analog, and the scaled measures, formed from T/c, do not change with
#     k while the model leaves the means free. ML-derived robust indices
#     are nevertheless withheld: this invariance does not establish their
#     validity outside ML. Eligible scaled measures remain qualified
#     diagnostics (the scaled RMSEA is 0.160 both times).
#   - when the model restricts the means (.sem_means_restricted()), the
#     fit function adds the squared mean residuals, in the squared units of
#     the variables, to the squared covariance residuals, in the squared
#     units of the covariances, so T is a sum of terms in different units
#     and no measure formed from it keeps its value (nor do the estimates,
#     when a parameter enters both the means and the covariances, as the
#     loadings do through the latent means of a multiple group model): with
#     t6's and t9's intercepts constrained equal in the same model, the
#     robust CFI is .979 for the scores and .972 for the scores divided by
#     10, and the scaled RMSEA .160 and .181. Every RMSEA and incremental
#     index is withheld; the test statistic is kept, whose reference
#     distribution under a corrected test does not depend on the units,
#     and so is the SRMR, which standardizes the covariance residuals of
#     the fit at hand.
# A fit by estimator = "ULS" alone is tested in lavaan 0.7-3 by Browne's
# normal theory residual based statistic, which weights the residuals by
# the inverse of their normal theory covariance matrix and has no units
# when it is evaluated at estimates that are equivariant, and its standard
# measures, formed from that statistic, are reported where they are: when
# the model leaves the means free, or restricts them through parameters
# the covariances do not share (.sem_uls_estimates_units()), as with t6's
# and t9's intercepts constrained equal above. When it restricts them
# through parameters the covariances share, the estimates change with the
# units, and so do that statistic and every RMSEA row and incremental
# index formed from it, which are withheld; the chi_square row, whose chi
# square reference distribution holds at any consistent estimates
# (Browne, 1984, Proposition 4), and the SRMR are kept.
# With ordered data the robust measures are the categorical correction of
# Savalei (2021), which lavaan computes from polychoric correlations,
# which have no units, and none of this applies.
.sem_uls_continuous <- function(fit) {
  identical(lavaan::lavInspect(fit, "options")$estimator, "ULS") &&
    !length(lavaan::lavNames(fit, "ov.ord"))
}

# The sentence that says why the robust RMSEA of a ULS fit of continuous
# data is not reported (.sem_uls_continuous()).
.sem_uls_robust_text <- function() {
  paste0("lavaan's robust RMSEA of a ULS fit of continuous data puts the ",
         "ULS statistic, which is in the units of the squared covariances, ",
         "into a correction derived for the ML fit function (Brosseau-Liard, ",
         "Savalei, & Li, 2012; Savalei, 2018), so it changes with the units ",
         "of the variables (multiplying every variable by k multiplies it ",
         "by k^2)")
}

# The sentence that says why no index of a ULS fit of continuous data
# that restricts the means is reported (.sem_uls_continuous()).
.sem_uls_means_text <- function() {
  paste0("the ULS fit function of a model that restricts the means adds the ",
         "squared residuals of the means, in the squared units of the ",
         "variables, to those of the covariances, in the squared units of ",
         "the covariances, so the statistic and every index formed from it ",
         "change with the units of the variables, and so do the estimates ",
         "when a parameter enters both, as the loadings do through the latent ",
         "means of a multiple group model")
}

# The sentence that says why no index formed from Browne's residual based
# statistic is reported for a ULS fit of continuous data whose estimates
# change with the units (.sem_uls_estimates_units()). `does` is that
# value for the fit at hand: TRUE, or NA when it could not be read.
.sem_uls_browne_text <- function(does = TRUE) {
  paste0("Browne's residual based statistic has no units, but it is ",
         "evaluated at the ULS estimates, which change with the units of ",
         "the variables when the model restricts the means through ",
         "parameters the covariances share, as the loadings enter the means ",
         "through the latent means of a multiple group model, ",
         if (isTRUE(does)) "and this model does" else {
           "which DMAR could not read for this model and treats as if it does"
         },
         ": the ULS fit function then trades the squared mean residuals, in ",
         "the squared units of the variables, against the squared ",
         "covariance residuals, in the squared units of the covariances, so ",
         "the statistic and every index formed from it change with the units ",
         "too")
}

# The line, under a table of fa_confirmatory() or measurement_invariance(), that
# names robust indices DMAR computed (.sem_robust_mva()).
.sem_robust_mva_line <- function(rows) {
  paste0(rows, ": the robust indices of Brosseau-Liard, Savalei, and Li ",
         "(2012) and Brosseau-Liard and Savalei (2014), computed by DMAR ",
         "from the scaling constant and degrees of freedom of the mean and ",
         "variance adjusted statistic, for which lavaan computes none; they ",
         "equal the robust indices under estimator = \"MLM\" (Savalei, 2018)")
}

# ---- printing -----------------------------------------------------------

# The three fit functions return a dmar_tbl with the leading class
# dmar_sem_fit. Every display of it (the console, knit_print(), as_kable())
# writes, after the lines any table carries, the lines in the footer
# element of its sem_fit attribute, which say which version of each
# measure was reported, the sample size convention, and the verdict of a
# test, so that the table cannot be read without them. The lines are
# built where the table is (sem_fit_measures(), equivalence_rmsea(),
# sem_residual_fit()); see .dmar_tbl_footer() in R/dmar_output_helpers.R.
#' @exportS3Method .dmar_tbl_footer dmar_sem_fit
#' @noRd
.dmar_tbl_footer.dmar_sem_fit <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "sem_fit")$footer, "sem_fit")
}

# ---- information criteria (sem_information_criteria(), compare_sem()) ----

# The criteria in the order of every table that reports them. The names are
# the row names (the case cfa_1() and fa_confirmatory() use for AIC and BIC); the
# arguments that select criteria accept them in either case.
.SEM_IC_NAMES <- c("AIC", "BIC", "SABIC", "HBIC", "IBIC", "SPBIC")

# Matches criterion names given in either case to .SEM_IC_NAMES. `several`
# allows a vector of distinct names (sem_information_criteria()'s
# `measures`); otherwise exactly one name (compare_sem()'s `criterion`).
.sem_ic_match <- function(x, arg, several = TRUE) {
  ok <- is.character(x) && length(x) >= 1L && !anyNA(x) &&
    (several || length(x) == 1L)
  m <- if (ok) match(toupper(x), .SEM_IC_NAMES) else NA_integer_
  if (!ok || anyNA(m) || anyDuplicated(m)) {
    .dmar_refuse_argument(c(
      "'", arg, "' must be ", if (several) "distinct names" else "one name",
      " from ", paste(.SEM_IC_NAMES, collapse = ", "),
      " (in either case)."), arg = NULL)
  }
  .SEM_IC_NAMES[m]
}

# The checks that make a fitted model's log likelihood the maximized log
# likelihood the criteria are defined on, with a message that names the
# model (`label`) and says what to do. `measures` decides whether the
# sample size enters (every criterion but AIC). Returns a list that
# describes the likelihood, for the attribute and the printed footer.
.sem_likelihood_check <- function(fit, options, label, measures) {
  if (options$missing %in% c("two.stage", "robust.two.stage")) {
    diagnostic <- .fiml_fit_h1_diagnostic(fit)
    if (!is.null(diagnostic)) {
      .dmar_refuse(c("The model in ", label, " uses two-stage estimation, ",
                     "whose unrestricted first-stage moments are unavailable. ",
                     diagnostic),
                   "h1_unavailable", diagnostic = diagnostic)
    }
  }
  ord <- lavaan::lavNames(fit, "ov.ord")
  if (length(ord)) {
    .dmar_refuse_argument(c(
      "The model in ", label, " has ordered categorical indicators (",
      paste(ord, collapse = ", "), "). Their estimators (WLSMV and the ",
      "other least squares estimators, or a pairwise likelihood) give no ",
      "full log likelihood, so information criteria are not defined for it."),
      arg = NULL)
  }
  if (!identical(options$estimator, "ML")) {
    .dmar_refuse_argument(c(
      "The model in ", label, " was fitted with estimator = \"",
      options$estimator.orig, "\", whose fit function is not a log ",
      "likelihood. Information criteria need a maximum likelihood fit ",
      "(estimator = \"ML\", or its robust variants MLM, MLMV, and MLR, ",
      "which maximize the same likelihood)."), arg = NULL)
  }
  if (!identical(options$likelihood, "normal")) {
    .dmar_refuse_argument(c(
      "The model in ", label, " was fitted with likelihood = \"",
      options$likelihood, "\", whose log likelihood is that of the sample ",
      "covariance matrix under the Wishart distribution rather than the ",
      "normal log likelihood of the data the criteria are defined on. ",
      "Refit it with likelihood = \"normal\" (lavaan's default)."), arg = NULL)
  }
  if (isTRUE(options$.sampling.weights)) {
    .dmar_refuse_argument(c(
      "The model in ", label, " was fitted with sampling weights. lavaan ",
      "then maximizes a weighted pseudo log likelihood (Asparouhov, 2005), ",
      "for which neither the penalties of AIC and the BIC family nor the ",
      "information matrices of IBIC and SPBIC have their likelihood ",
      "justification, so information criteria are not computed for it."),
      arg = NULL)
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L && any(measures != "AIC")) {
    .dmar_refuse_argument(c(
      "The model in ", label, " is a two-level model, for which the ",
      "sample size in the penalty of the BIC family is ambiguous: the ",
      "number of level 1 units, the number of clusters, or something ",
      "between (Bollen, Ray, Zavisca, & Harden, 2012, p. 295). AIC, which ",
      "does not involve the sample size, is available with measures = ",
      "\"AIC\"; for the others, supply log_likelihood, npar, and the N you ",
      "choose as summary input."), arg = NULL)
  }
  conditioning <- if (isTRUE(options$fixed.x) || isTRUE(options$conditional.x)) {
    sort(lavaan::lavNames(fit, "ov.x"))
  } else character(0)
  missing <- options$missing
  dropped <- sum(lavaan::lavInspect(fit, "norig")) - lavaan::lavInspect(fit, "ntotal")
  likelihood <- if (missing %in% c("ml", "ml.x")) {
    "the observed data (FIML) normal log likelihood"
  } else if (missing %in% c("two.stage", "robust.two.stage")) {
    paste0("two-stage: the complete data normal log likelihood at the stage ",
           "one EM estimates of the means and covariances (Lai, 2021)")
  } else if (dropped > 0) {
    paste0("the normal log likelihood of the complete cases (", dropped,
           " case", if (dropped == 1) "" else "s", " removed by listwise deletion)")
  } else "the normal log likelihood of the complete data"
  list(likelihood = likelihood, missing = missing, conditioning = conditioning,
       conditional_x = isTRUE(options$conditional.x),
       meanstructure = isTRUE(options$meanstructure),
       groups = lavaan::lavInspect(fit, "ngroups"),
       estimator = options$estimator.orig)
}

# How each latent variable's scale is set, read from the parameter table of
# a single group fit: "variance fixed" when its variance (or disturbance
# variance) is fixed and every loading on it is free (std.lv = TRUE),
# "marker indicator" when a loading is fixed to a nonzero value and its
# variance is free (lavaan's default), and "other" for anything else
# (effects coding, a latent variable without indicators, both fixed). IBIC,
# and SPBIC under its default prior, depend on this choice.
.sem_identification <- function(fit) {
  pt <- lavaan::parTable(fit)
  lv <- lavaan::lavNames(fit, "lv")
  if (!length(lv)) return(character(0))
  pt <- pt[pt$block == 1L, , drop = FALSE]
  vapply(lv, function(f) {
    load <- pt[pt$op == "=~" & pt$lhs == f, , drop = FALSE]
    v <- pt[pt$op == "~~" & pt$lhs == f & pt$rhs == f, , drop = FALSE]
    var_fixed <- nrow(v) == 1L && v$free == 0L && isTRUE(v$est > 0)
    load_fixed <- load$free == 0L & load$est != 0
    if (!nrow(load)) return("other")
    if (var_fixed && !any(load_fixed)) return("variance fixed")
    if (!var_fixed && sum(load_fixed) == 1L) return("marker indicator")
    "other"
  }, character(1))
}

# The conditions under which the information matrices lavaan reports are the
# ones IBIC and SPBIC are defined on (Bollen et al., 2014, Eqs. 9, 16, and
# 18): a single group, single level normal theory ML fit of complete data
# with unconstrained free parameters, the expected information evaluated at
# the fitted moments for IBIC, and the observed information as the negative
# Hessian for SPBIC. Returns the parameter table; refuses with every reason
# that applies.
.sem_information_fit_check <- function(fit, options, measures, label = "'fit'") {
  pt <- lavaan::parTable(fit)
  # A fit to summary statistics has no case indices and no raw data to
  # inspect; its moments are complete by construction.
  raw <- !is.null(lavaan::lavInspect(fit, "case.idx"))
  dat <- if (raw) lavaan::lavInspect(fit, "data") else NULL
  ll_df <- as.numeric(attr(lavaan::logLik(fit), "df"))
  n_coef <- length(lavaan::coef(fit))
  why <- character(0)
  if (options$estimator.orig != "ML") {
    why <- c(why, paste0("it was fitted with estimator = \"", options$estimator.orig,
                         "\", and IBIC and SPBIC rest on the normal theory ",
                         "information matrices a robust fit was chosen not to trust ",
                         "(refit with estimator = \"ML\")"))
  }
  if (lavaan::lavInspect(fit, "ngroups") != 1L) why <- c(why, "it has more than one group")
  if (lavaan::lavInspect(fit, "nlevels") > 1L) why <- c(why, "it is a two-level model")
  if ((isTRUE(options$fixed.x) || isTRUE(options$conditional.x)) &&
      length(lavaan::lavNames(fit, "ov.x"))) {
    why <- c(why, paste0("its log likelihood is conditional on the exogenous ",
                         "observed variables (refit with fixed.x = FALSE)"))
  }
  if (any(pt$op %in% c("==", "<", ">")) || anyDuplicated(pt$free[pt$free > 0L])) {
    why <- c(why, "it has equality or inequality constraints, so its free coefficients are not independent")
  }
  if (!is.null(pt$efa) && any(nzchar(pt$efa))) {
    why <- c(why, paste0("it has exploratory (rotated) factor blocks, whose ",
                         "information matrices are singular, of rank the number ",
                         "of coefficients less the rotational constraints"))
  } else if (!isTRUE(ll_df == n_coef) ||
             !isTRUE(lavaan::lavInspect(fit, "npar") == n_coef)) {
    why <- c(why, "its number of parameters differs from its number of free coefficients")
  }
  incomplete <- raw && anyNA(unlist(dat))
  if (incomplete) {
    why <- c(why, paste0("its analysis data have missing values (missing = \"",
                         options$missing, "\")"))
  }
  if ("SPBIC" %in% measures && options$observed.information[1L] != "hessian") {
    why <- c(why, paste0("its observed information is lavaan's observed.information = \"",
                         options$observed.information[1L], "\" approximation rather ",
                         "than the negative Hessian SPBIC is defined on (refit with ",
                         "observed.information = \"hessian\", lavaan's default)"))
  }
  # lavaan sets h1.information = "unstructured" itself for the two-stage
  # method, whose refusal the missing values already explain.
  if ("IBIC" %in% measures && !incomplete &&
      options$h1.information[1L] != "structured") {
    why <- c(why, paste0("its expected information is evaluated at the sample ",
                         "moments (h1.information = \"", options$h1.information[1L],
                         "\") rather than at the fitted moments IBIC is defined on ",
                         "(refit with h1.information = \"structured\", lavaan's default)"))
  }
  if (length(why)) {
    which <- paste(intersect(c("IBIC", "SPBIC"), measures), collapse = " and ")
    .dmar_refuse_argument(c(
      "Automatic ", which, " extraction is not available for the model in ",
      label, ": ", paste(why, collapse = "; "), ". It needs a single group, ",
      "single level, normal theory ML fit of complete continuous data whose ",
      "free parameters are unconstrained. In another parameterization, ",
      "supply log_likelihood, npar, N, and the information matrices as ",
      "summary input."), arg = NULL)
  }
  invisible(pt)
}

.sem_information_matrix <- function(x, d, name, scale, N, estimate_names = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(d, d)) ||
      any(!is.finite(x)) || !isSymmetric(unname(x), tol = 1e-8)) {
    .dmar_refuse_argument(c(
      "'", name, "' must be a finite symmetric ", d, " by ", d,
      " matrix in the same parameter order as 'parameter_estimate'."),
      arg = "parameter_estimate")
  }
  if (!is.null(estimate_names) && !is.null(rownames(x)) &&
      !identical(rownames(x), estimate_names)) {
    .dmar_refuse_argument(c(
      "Parameter names/order in '", name, "' and 'parameter_estimate' must ",
      "agree."), arg = "parameter_estimate")
  }
  if (scale == "per_observation") x <- x * N
  ch <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(ch) || any(!is.finite(ch)) || any(diag(ch) <= 0)) {
    .dmar_refuse_argument(c(
      "'", name, "' must be positive definite in independent free parameters."),
      arg = NULL)
  }
  list(matrix = x, chol = ch, logdet = 2 * sum(log(diag(ch))))
}

.sem_prior <- function(parameter_estimate, prior_mean, zero_parameters,
                       pt = NULL) {
  d <- length(parameter_estimate)
  if (!is.null(prior_mean) && !is.null(zero_parameters)) {
    .dmar_refuse_argument(
      "Supply 'prior_mean' or 'zero_parameters', not both.",
      arg = c("prior_mean", "zero_parameters"))
  }
  if (!is.null(prior_mean)) {
    if (!is.numeric(prior_mean) || length(prior_mean) != d ||
        any(!is.finite(prior_mean))) {
      .dmar_refuse_argument(
        "'prior_mean' must have one finite value per free parameter.",
        arg = "prior_mean")
    }
    if (!is.null(names(prior_mean)) &&
        !identical(names(prior_mean), names(parameter_estimate))) {
      .dmar_refuse_argument(c(
        "Named 'prior_mean' must match the names and order of ",
        "'parameter_estimate'."),
        arg = c("prior_mean", "parameter_estimate"))
    }
    return(list(mean = prior_mean, rule = "explicit prior mean"))
  }
  if (is.null(zero_parameters)) {
    if (is.null(pt)) {
      .dmar_refuse_argument(
        "SPBIC summary input requires 'prior_mean' or explicit 'zero_parameters'.",
        arg = c("prior_mean", "zero_parameters"))
    }
    zero_parameters <- unique(pt$free[pt$free > 0L & pt$op %in% c("=~", "~", "<~")])
    rule <- paste0("Bollen et al. (2014, p. 6): free loadings and regression ",
                   "coefficients centered at zero, other free parameters at ",
                   "their estimates")
  } else rule <- "selected parameters centered at zero, other free parameters at their estimates"
  if (is.character(zero_parameters)) {
    if (is.null(names(parameter_estimate)) ||
        anyDuplicated(names(parameter_estimate)) ||
        any(!zero_parameters %in% names(parameter_estimate))) {
      .dmar_refuse_argument(c(
        "'zero_parameters' names must uniquely identify entries of ",
        "'parameter_estimate'."),
        arg = c("zero_parameters", "parameter_estimate"))
    }
    zero_parameters <- match(zero_parameters, names(parameter_estimate))
  }
  if (!is.numeric(zero_parameters) || !length(zero_parameters) ||
      any(!is.finite(zero_parameters)) || any(zero_parameters != floor(zero_parameters)) ||
      any(zero_parameters < 1 | zero_parameters > d)) {
    .dmar_refuse_argument(c(
      "'zero_parameters' must select at least one free parameter. For a saturated ",
      "covariance model, explicitly select the off-diagonal covariances."),
      arg = "zero_parameters")
  }
  mean <- parameter_estimate
  mean[zero_parameters] <- 0
  list(mean = mean, rule = rule)
}

# What compare_sem() requires of every model, recorded per fit: whether it
# was fitted to raw data or to summary statistics, its groups, its observed
# variables, its cases and data values (columns put in name order, since
# the likelihood does not depend on the order of the variables), the
# conventions that define its likelihood, and for a two-stage fit the
# auxiliary variables of its first stage, which move the EM estimates its
# log likelihood is evaluated at. `lik` is the result of
# .sem_likelihood_check().
.sem_comparison_signature <- function(fit, lik) {
  ov <- sort(lavaan::lavNames(fit, "ov"))
  cases <- lavaan::lavInspect(fit, "case.idx")
  raw <- !is.null(cases)
  by_name <- function(m) m[, ov, drop = FALSE]
  values <- if (raw) {
    d <- lavaan::lavInspect(fit, "data")
    if (is.list(d)) lapply(d, by_name) else by_name(d)
  } else {
    s <- lavaan::lavInspect(fit, "sampstat")
    sort_moments <- function(m) {
      out <- list(cov = unclass(m$cov)[ov, ov])
      if (!is.null(m$mean)) out$mean <- unclass(m$mean)[ov]
      out
    }
    if (!is.null(s$cov)) sort_moments(s) else lapply(s, sort_moments)
  }
  aux <- if (lik$missing %in% c("two.stage", "robust.two.stage")) {
    sort(.ts_fit_aux(fit)$names)
  } else character(0)
  # The cases a two-stage first stage read, which a fit lavaan made counts
  # short by those observed on an auxiliary variable alone and a fit
  # fa_confirmatory() or measurement_invariance() returned does not
  # (.ts_rows_read()), so that two such fits of the same data compare.
  rows <- .ts_rows_read(fit)
  list(raw = raw, groups = lavaan::lavInspect(fit, "ngroups"),
       group_labels = lavaan::lavInspect(fit, "group.label"),
       N = if (is.null(rows)) lavaan::lavInspect(fit, "ntotal")
           else sum(rows$read),
       observed = ov,
       cases = cases, values = values, missing = lik$missing,
       aux = aux,
       meanstructure = lik$meanstructure, conditioning = lik$conditioning,
       conditional_x = lik$conditional_x && length(lik$conditioning) > 0L)
}

# The missing data method a `missing` setting stands for, as the likelihood
# sees it: lavaan's "two.stage" and "robust.two.stage" fit the same EM
# moments by the same complete data likelihood (Lai, 2021) and differ only
# in their standard errors and test statistic.
.sem_missing_method <- function(missing) {
  if (missing %in% c("two.stage", "robust.two.stage")) "two.stage" else missing
}

# Refuses a set of models whose likelihoods are not comparable, naming the
# first part of the signature that differs, the two models, and what to do.
.sem_comparison_check <- function(sig) {
  nm <- names(sig)
  first <- sig[[1L]]
  show <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
  for (i in seq_along(sig)[-1L]) {
    s <- sig[[i]]
    pair <- paste0("fits$", nm[1L], " and fits$", nm[i])
    fail <- function(...) .dmar_refuse_argument(c(pair, " ", ...), arg = NULL)
    if (!identical(first$raw, s$raw)) {
      fail("were fitted differently, one to raw data and one to summary ",
           "statistics (sample.cov); fit every model the same way.")
    }
    if (!identical(first$groups, s$groups) ||
        !identical(first$group_labels, s$group_labels)) {
      fail("have different groups (", first$groups, " and ", s$groups,
           "); the models must be fitted to the same groups.")
    }
    if (!identical(first$observed, s$observed)) {
      fail("use different observed variables (", show(setdiff(first$observed, s$observed)),
           " only in the first; ", show(setdiff(s$observed, first$observed)),
           " only in the second). A log likelihood is a density of the variables ",
           "it models, so models of different variables cannot be compared by it.")
    }
    # The two two-stage settings share their EM moments, estimates, and log
    # likelihood and differ only in their standard errors and test, which
    # the criteria do not use, so they count as one method.
    if (!identical(.sem_missing_method(first$missing), .sem_missing_method(s$missing))) {
      fail("handle missing data differently (missing = \"", first$missing,
           "\" and \"", s$missing, "\"), which gives different likelihoods; fit ",
           "every model with the same missing data method.")
    }
    if (!identical(first$aux, s$aux)) {
      fail("use different auxiliary variables in their first stage (",
           show(first$aux), " and ", show(s$aux), "). A two-stage log ",
           "likelihood is evaluated at the EM estimates of the means and ",
           "covariances, which the auxiliary variables change; fit every ",
           "model with the same aux.")
    }
    if (!identical(first$N, s$N)) {
      fail("were fitted to different numbers of cases (", first$N, " and ", s$N,
           "); compare models fitted to the same cases.")
    }
    if (!identical(first$cases, s$cases)) {
      fail("were fitted to different cases; compare models fitted to the same cases.")
    }
    if (!identical(first$values, s$values)) {
      fail("were fitted to different data values for the same ",
           if (first$raw) "cases" else "summary statistics", "; compare models ",
           "fitted to the same data.")
    }
    if (!identical(first$meanstructure, s$meanstructure)) {
      fail("differ in their mean structure (meanstructure = ", first$meanstructure,
           " and ", s$meanstructure, "). Without a mean structure lavaan leaves ",
           "the saturated means out of the parameter count, which moves every ",
           "criterion by an amount that depends on the setting; fit every model ",
           "with the same meanstructure.")
    }
    if (!identical(first$conditioning, s$conditioning)) {
      fail("condition their log likelihoods on different exogenous observed ",
           "variables ({", paste(first$conditioning, collapse = ", "), "} and {",
           paste(s$conditioning, collapse = ", "), "}). With fixed.x = TRUE (the ",
           "default of lavaan::sem()) or conditional.x = TRUE, lavaan's log ",
           "likelihood is the conditional density given the exogenous observed ",
           "variables, and conditional densities given different variables are ",
           "not comparable. Refit every model with fixed.x = FALSE so that all ",
           "the log likelihoods are joint densities of the same variables.")
    }
    if (!identical(first$conditional_x, s$conditional_x)) {
      fail("differ in conditional.x (", first$conditional_x, " and ",
           s$conditional_x, "), which lavaan's parameter count follows; fit ",
           "every model with the same conditional.x, or every model with ",
           "fixed.x = FALSE.")
    }
  }
  invisible(TRUE)
}

# A fit convention is separate from uncertainty in estimated parameters.
.sem_saturated_text <- function() {
  paste0("The model is saturated (zero degrees of freedom and zero discrepancy). ",
    "RMSEA = 0 with limits [0, 0] and TLI = 1 follow lavaan's saturated-model conventions; ",
    "they describe reproduction of the sample moments, not evidence of structural fit. ",
    "Their conventional limits do not estimate sampling precision, and no fit-index standard error is supplied. ",
    "Parameter standard errors and confidence intervals still describe uncertainty in the estimated parameters. ",
    "The exact-fit p-value and close-fit tests are unavailable.")
}
