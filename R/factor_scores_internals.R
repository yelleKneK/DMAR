# =============================================================================
#  Internal helpers for factor_scores()
# =============================================================================
#
# factor_scores() returns, for every row given to a fitted lavaan model and
# every factor, the regression score and the Bartlett score (for ordered
# items, lavaan's empirical Bayes modal and maximum likelihood scores), each
# with a per-case standard error and interval, the number of the factor's
# items the case answered, its missing data pattern, and the determinacy of
# each score type for that pattern. The helpers below do the parts that
# lavaan::lavPredict() does not: the per-pattern determinacy, the scores
# from the indicator block when auxiliary variables were carried as
# saturated correlates, the scores of a fit of continuous items that keeps
# incomplete cases without FIML (two-stage ML, pairwise deletion, which
# lavPredict() scores as complete data, leaving every incomplete case
# unscored), the factors whose Bartlett or ML score a case's observed items
# do not determine, the alignment of lavaan's per-group output with the
# rows given to the fit, the check that `data` holds those rows, and the
# reasons attached to every value that is missing or borrowed.
#
# Notation, per group: Lambda (p x m) is the loading matrix with lavaan's
# dummy latent variables included, Theta (p x p) the residual covariance
# matrix, B (m x m) the regressions among latent variables, Psi (m x m)
# their residual covariance matrix, and
#   V(eta) = (I - B)^{-1} Psi (I - B)^{-T},
#   Sigma  = Lambda V(eta) Lambda' + Theta.
# lavaan carries an observed variable that takes part in a regression, or
# covaries with a latent variable or with another such variable, as a dummy
# latent variable measured by itself without error: an exogenous covariate
# (an ov.x dummy), and any other (an ov.y dummy), such as an outcome the
# factor predicts, an item with a direct effect of a covariate
# (differential item functioning in a MIMIC model), or a criterion that
# covaries with a factor (f ~~ z). An ov.y dummy's loadings on the latent
# variables are then in its row of B, and its residual (co)variances in
# Psi. In a model with a B matrix, lavPredict() works with the loading
# matrix in which each ov.y dummy's row is replaced by its row of B
# (lav_lisrel_lambda()); without B (a fit conditional on the covariates
# with no regression among the variables lavaan holds as latent) every
# dummy keeps its unit row, one regressed only on exogenous covariates
# (through Gamma) included. lavaan's check of the solution and its EBM
# and ML scores use the residual covariance matrix into which the
# dummies' block of Psi is folded (lav_lisrel_theta()). For an ov.y dummy
# regressed on nothing, the row of B is zero, and the variable is better
# kept as lavaan holds it, measured by its own dummy latent variable: then
# the scores condition on it as on an exogenous covariate. .fs_dummies()
# builds both matrices that way, and every computation below uses them.
# Sigma and E(y) are lavaan's model implied moments, always formed from the
# raw matrices (.fs_moments()). Lambda_eff and Theta_eff are not in general
# a decomposition of Sigma, and nothing below forms Lambda_eff V(eta)
# Lambda_eff' + Theta_eff: Theta_eff folds in the Psi block of every dummy,
# as lav_lisrel_theta() does, so a dummy that keeps its unit row (an
# exogenous covariate, or an observed variable regressed on nothing or,
# without B, only on exogenous covariates) would have its variance, and
# its covariances with other such dummies, counted twice (f ~~ z doubles
# the variance of z); the covariance of z1 and z2 in a chain z1 ~ f,
# z2 ~ z1 would not be reproduced; and in a fit .fs_check_dummies()
# refuses, nor would the covariances of a factor's items with an observed
# variable regressed on others that predicts the factor. The decomposition
# the scores rest on is the one given at .fs_check_dummies(). For the
# observed items o of a missing data pattern (Lambda_o, Sigma_oo), with
# Lambda the loading matrix of .fs_dummies(),
#   regression: A_R = V(eta) Lambda_o' Sigma_oo^{-1},
#               posterior covariance V(eta) - A_R Lambda_o V(eta);
#   Bartlett:   A_B = (Lambda_o' Sigma_oo^{-1} Lambda_o)^{-1}
#                     Lambda_o' Sigma_oo^{-1},
#               error covariance (Lambda_o' Sigma_oo^{-1} Lambda_o)^{-1}
#                     - V(eta) = (Lambda_o' Theta_oo^{-1} Lambda_o)^{-1}.
# These are the matrices lav_predict_eta_normal() in lavaan 0.7-3 uses
# pattern by pattern under missing = "ml", except where its loading matrix
# gives an ov.y dummy regressed on nothing its zero row of B (whenever the
# model has a B matrix), which drops the variable's covariances with the
# factors (such fits are scored here; see factor_scores()); written with
# Sigma^{-1}, the Bartlett matrix also covers a dummy latent variable whose
# residual variance is zero (an exogenous covariate, or an observed
# variable regressed on nothing, which the Bartlett score so conditions
# on). The determinacy of a score is its correlation with the
# factor, [A Lambda_o V(eta)]_kk / sqrt([A Sigma_oo A']_kk V(eta)_kk), the
# formula lavaan uses for fs.determinacy, here evaluated for each pattern.
#
# Layout:
#   (A) Checks on the fit          .fs_check_fit, .fs_inadmissible,
#                                  .fs_check_dummies, .fs_reaches_lv,
#                                  .fs_residual_pairs
#   (B) Model structure            .fs_as_list, .fs_items, .fs_held_cells,
#                                  .fs_holds, .fs_zero_held,
#                                  .fs_zero_held_fit,
#                                  .fs_auxiliary, .fs_dummies, .fs_moments
#   (C) Per-pattern algebra        .fs_patterns, .fs_identified,
#                                  .fs_pattern_info
#   (D) lavaan's scores            .fs_lavpredict, .fs_ml_finite,
#                                  .fs_in_cone
#   (E) Rows given to the fit      .fs_data_check
#   (F) Reasons                    .fs_add_reason
#
# Nothing in this file is exported.
# =============================================================================


# ---------------------------------------------------------------------------
# (A) Checks on the fit
# ---------------------------------------------------------------------------

# Refuse what cannot be scored, each with the remedy: an object that is not
# a lavaan fit, a fit to summary statistics (no rows), a fit that did not
# converge, and a two-level fit (lavPredict() scores each level separately,
# which is a different table).
.fs_check_fit <- function(fit) {
  if (!inherits(fit, "lavaan") || inherits(fit, "blavaan")) {
    .dmar_refuse_argument(c(
      "'fit' must be a lavaan fit object, for example from ",
      "lavaan::cfa(), lavaan::sem(), or fa_confirmatory(output = \"fit\")."),
      arg = "fit")
  }
  if (!length(unlist(lavaan::lavInspect(fit, "case.idx")))) {
    .dmar_refuse_argument(c(
      "The model was fitted to summary statistics (a covariance ",
      "matrix), so it holds no rows to score. Fit it to the raw data ",
      "to obtain factor scores."), arg = NULL)
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    .dmar_refuse(c(
      "The model did not converge, so its estimates are not a ",
      "solution to score with. Resolve the convergence problem ",
      "(starting values, a simpler model, more iterations) and refit."),
      "not_converged")
  }
  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(c(
      "factor_scores() scores single-level models. For a two-level ",
      "model, lavaan::lavPredict(fit, level = 1) and level = 2 give ",
      "the within and between scores."), arg = NULL)
  }
  if (!length(lavaan::lavNames(fit, "lv"))) {
    .dmar_refuse_argument(
      "The model has no latent variables to score.", arg = NULL)
  }
  invisible(TRUE)
}

# Ordered items are scored by lavaan's empirical Bayes modal (EBM) and
# maximum likelihood (ML) scores, which lav_predict_eta_ebm_ml() finds case
# by case by minimizing, over the factor value eta, the negative log
# likelihood of the ordered responses (plus, for EBM, the normal prior
# centered at the factor's model implied mean mu). A response's probability
# compares each threshold tau with the mean of its latent response,
# nu + Lambda eta. In lavaan 0.7-2 the objective took its thresholds from
# lav_model_th(), the model implied thresholds tau - (nu + Lambda mu), and
# compared them with nu + Lambda (eta - mu) (it passed `x2 <- x - mu_i` to
# lav_predict_fy_eta_i()): the factor mean cancels, but the intercepts nu
# of the ordered indicators enter twice. The scores were therefore wrong in
# any model whose ordered indicators have intercepts other than zero, a
# multiple group model with equal thresholds for one. lavaan 0.7-3
# compares the thresholds of the tau matrix with nu + Lambda eta (it passes
# `x2 <- x`), and it also absorbs the scale factors of the delta
# parameterization first (lav_model_delta_absorb()); the functions that
# evaluate each case's response probabilities and latent response means
# (lav_predict_fy_eta_i(), lav_lisrel_eyetax()) are the same in both
# versions. SIEVER's test of its categorical scores against values
# computed directly from the fitted ordered probit model, for a two-group
# model with equal loadings and thresholds, found the 0.7-2 scores off by up
# to 0.47 (EBM) and 0.55 (ML) and the 0.7-3 scores in agreement. lavaan sits
# in DMAR's Suggests, whose version floor R does not enforce when DMAR
# loads, so the floor for ordered items is enforced here. It applies to
# every ordered fit rather than only those with nonzero intercepts: 0.7-3
# changed the treatment of the scale factors as well, and which other
# models 0.7-2 scored correctly was not established. The installed version
# is an argument so that the refusal can be tested.
.fs_require_lavaan_ordered <- function(lavaan_version =
                                         utils::packageVersion("lavaan"),
                                       lavaan_floor = "0.7-3") {
  if (as.numeric_version(lavaan_version) < as.numeric_version(lavaan_floor)) {
    .dmar_refuse(c(
      "Factor scores for ordered items require lavaan ", lavaan_floor,
      " or newer (installed: ", as.character(lavaan_version), "). ",
      "lavaan 0.7-2 counted the intercepts of the ordered indicators ",
      "twice in its EBM and ML scores, which are then wrong in any model ",
      "where those intercepts are not zero, such as a multiple group ",
      "model with equal thresholds. Update it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  invisible(TRUE)
}

# Describe an inadmissible solution, or return character(0): a negative
# residual variance of an observed variable, a negative variance of a
# latent variable, a covariance matrix of the factors that is not positive
# definite (a correlation beyond 1, for one), a residual covariance matrix
# of the continuous observed variables that is not positive definite, or,
# when none of these holds, a residual covariance matrix of the structural
# part of the model (Psi) that is not. These are the conditions lavaan's
# post-estimation check (lavInspect(fit, "post.check"),
# lav_object_post_check()) flags, which sem_fit_measures() refuses, and a
# matrix is judged with that check's tolerance: an eigenvalue below
# -eps^(3/4) times the larger of 1 and the largest eigenvalue. As there, a
# matrix is judged only when no variance on its diagonal is negative, and
# the residual covariance matrix of the observed variables is lavaan's,
# with the residual (co)variances of an observed variable regressed on
# others (a dummy latent variable, whose own Theta entries are zero) taken
# from Psi (`theta_eff`, from .fs_dummies()). Read from Theta alone, a
# negative residual variance of an item with a direct effect of a
# covariate (y1 ~ x) went unseen and was scored without a warning. The
# regression weights use the factors' covariance matrix and both weights
# the loadings and the residual covariance matrix, so the scores inherit
# the improper solution. Auxiliary variables (`aux`) are left out of the
# residual covariance matrices judged: the scores set them aside, and a
# saturated correlate strongly related to the factor routinely makes the
# whole matrix indefinite while that of the indicators is not.
.fs_inadmissible <- function(est, glab, lv, ordered = character(0),
                             aux = character(0)) {
  out <- character(0)
  indefinite <- function(s) {
    if (!length(s)) return(FALSE)
    ev <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
    any(ev < -.Machine$double.eps^(3 / 4) * max(1, abs(ev[1L])))
  }
  for (g in seq_along(est)) {
    where <- if (length(est) > 1L) sprintf(" in group %s", glab[g]) else ""
    theta <- est[[g]]$theta_eff %||% est[[g]]$theta
    th <- diag(theta)
    names(th) <- rownames(theta)
    neg_th <- names(th)[th < 0]
    if (length(neg_th)) {
      out <- c(out, sprintf("a negative residual variance for %s%s",
                            paste(neg_th, collapse = ", "), where))
    }
    ps <- diag(est[[g]]$psi)
    names(ps) <- rownames(est[[g]]$psi)
    neg_lv <- intersect(names(ps)[ps < 0], lv)
    if (length(neg_lv)) {
      out <- c(out, sprintf("a negative variance for the factor %s%s",
                            paste(neg_lv, collapse = ", "), where))
    }
    lv_g <- intersect(lv, colnames(est[[g]]$psi))
    veta_bad <- FALSE
    if (!length(neg_lv) && length(lv_g) > 1L) {
      m <- ncol(est[[g]]$psi)
      beta <- if (is.null(est[[g]]$beta)) matrix(0, m, m) else est[[g]]$beta
      ib_inv <- solve(diag(m) - beta)
      veta <- ib_inv %*% est[[g]]$psi %*% t(ib_inv)
      dimnames(veta) <- dimnames(est[[g]]$psi)
      veta <- veta[lv_g, lv_g, drop = FALSE]
      if (indefinite(veta)) {
        veta_bad <- TRUE
        far <- NULL
        if (all(diag(veta) > 0)) {
          r <- stats::cov2cor(veta)
          far <- which(abs(r) > 1 & upper.tri(r), arr.ind = TRUE)
        }
        out <- c(out, if (length(far) && nrow(far)) {
          sprintf("a correlation of %s between the factors %s and %s%s",
                  format(r[far], digits = 3), lv_g[far[, 1L]],
                  lv_g[far[, 2L]], where)
        } else {
          sprintf(paste0("a covariance matrix of the factors that is not ",
                         "positive definite%s"), where)
        })
      }
    }
    num <- setdiff(rownames(theta), c(ordered, aux))
    theta_bad <- !length(neg_th) &&
      indefinite(theta[num, num, drop = FALSE])
    if (theta_bad) {
      out <- c(out, sprintf(paste0("a residual covariance matrix of the ",
                                   "continuous observed variables that is ",
                                   "not positive definite%s"), where))
    }
    # lavaan judges Psi as a whole only when the checks above all pass.
    psi_v <- setdiff(rownames(est[[g]]$psi), aux)
    if (!length(neg_th) && !length(neg_lv) && !veta_bad && !theta_bad &&
        indefinite(est[[g]]$psi[psi_v, psi_v, drop = FALSE])) {
      out <- c(out, sprintf(paste0("a residual covariance matrix of the ",
                                   "structural part of the model (Psi) ",
                                   "that is not positive definite%s"), where))
    }
  }
  out
}

# Refuse a fit whose observed variables in regressions lavaan scores
# wrongly. In a model with a B matrix, lavaan's scores treat an observed
# variable regressed on others (a variable of `dummy_b`, see .fs_dummies())
# as a measure of the latent variables in its row of B, with its
# disturbance zeta as the measurement error: y = nu + Lambda_eff eta +
# (zeta on those rows) + epsilon, where the variables regressed on nothing
# and the exogenous covariates keep their unit rows and so are conditioned
# on. Both scores rest on that decomposition. Without B (a fit conditional
# on the covariates with no regression among the variables lavaan holds as
# latent), lavaan replaces no row: an observed variable regressed only on
# exogenous covariates (z ~ x, through Gamma) keeps its unit row, as one
# regressed on nothing does, and is conditioned on, so `dummy_b` is empty
# and there is nothing to check. lavPredict()'s regression score is then
# E(eta | y, z, x) and its Bartlett score conditionally unbiased given the
# factors, z, and x, also when z covaries with a factor (f ~ x, z ~ x,
# f ~~ z, where lavaan frees the covariance by default), which was
# refused. The regression weights V(eta) Lambda_eff' Sigma^{-1} give
# a factor's conditional mean only when Lambda_eff V(eta) gives the
# covariances of the observed variables with the factor, that is, when the
# factor is uncorrelated with every such disturbance; and the Bartlett score
# is conditionally unbiased given the factors and the variables it
# conditions on only when E(y | those) = Lambda_eff eta, that is, when
# every such disturbance is uncorrelated with all of them as well. So the
# condition is that Cov(eta, zeta) = (I - B)^{-1} Psi vanish in the rows of
# the regular latent variables, the exogenous covariates, and the observed
# variables regressed on nothing (auxiliary variables, which the scores set
# aside, excepted), in the columns of the variables of `dummy_b`.
# It does not ask the rows of those variables themselves: in a chain such
# as z1 ~ f, z2 ~ z1, Cov(eta, zeta) is not zero in the row of z2 and the
# column of z1, and lavaan's regression score is still E(f | y) and its
# Bartlett score conditionally unbiased. The condition fails when such a
# variable predicts a factor, directly or through others (an observed
# mediator, f ~ a with a ~ x), or its residual covaries with a factor or
# with a variable the scores condition on; in such a fit lavaan's
# regression score of the factor misses its conditional mean (by 1.08 in
# the Holzinger and Swineford example of the tests) when a factor is
# involved, and its Bartlett score is biased given the variables it
# conditions on, so the fit is refused.
#
# Ordered items are scored by lavaan's EBM and ML scores, which take each
# case's likelihood from the latent response means of lav_lisrel_eyetax()
# and the factors' prior from their model implied means and covariance
# matrix given the exogenous covariates (lav_model_eetax(),
# lav_model_vetax()). Three parts of a model are left out there, and a fit
# in which leaving one out changes the scores is refused. (1) The direct
# effect of an exogenous covariate on an observed variable (the gamma row
# of an ov.y dummy): lav_lisrel_eyetax() never reads Gamma. For a variable
# that depends on a factor directly (its row of B reaches one: y1 ~ x for
# an ordered item, the standard MIMIC model of differential item
# functioning), the thresholds are estimated given the covariate, so the
# scores are those of a case with the covariate at zero (the ML scores off
# by 0.11 to 0.44 in the example of the tests), and the fit is refused.
# For a variable that depends on no factor directly (w ~ x), the term
# left out shifts the mean of a density that does not involve the factors,
# which changes each case's likelihood by a factor constant in them, so
# the scores are the model's and the fit is scored. (2) A regression on
# an observed variable that is not an exogenous covariate (a nonzero entry
# B[r, d] in the column of an ov.y dummy d): lav_lisrel_eyetax() sets the
# ov.y dummies to zero before it multiplies by B, which leaves B[r, d] d
# out of the mean of r. That changes the scores, and the fit is refused,
# when r depends on a factor directly (B[r, lv] nonzero: y4 ~ y3 for local
# dependence between two items, the ML scores off by up to 0.18 in the
# example of the tests, or a partially mediated w ~ m + f), when d is an
# ordered item (its latent response is not observed, so the model's
# distribution of r given d's category involves the factors), or when r
# is a factor (whose prior in the EBM score then has the mean that d's
# mean implies rather than the mean its value implies; the check of
# Cov(eta, zeta) above refuses such a fit first). Otherwise, as in a fully
# mediated chain of continuous variables (m ~ f, w ~ m), the density of r
# given its parents involves no factor (the local Markov property of a
# recursive model with uncorrelated disturbances), so the term left out
# is constant in the factors wherever d is observed, and the scores are
# the model's (to 1.5e-4 for the EBM and 9e-4 for the ML score, against a
# direct maximization of the model's likelihood, in the example of the
# tests); such a fit is scored, unless a row the fit keeps has d missing
# and r observed while a path of B leads from d to a factor through no
# variable the row has observed (pairwise deletion): the model's
# likelihood of r then runs through d's dependence on the factors, which
# lavaan leaves out, and the fit is refused. A variable the row has
# observed on each such path (m in m ~ f, w ~ m, s ~ w, with w missing)
# cuts that dependence, and the row is scored. Only r's own row matters,
# since lav_lisrel_eyetax() applies B once. (3) An
# observed variable whose row lavaan keeps and whose disturbance covaries
# with a factor (regressed on nothing, f ~~ x; or, in a fit without B,
# regressed only on exogenous covariates, w ~ x with f ~~ w, which
# lavaan::sem() frees by default when x also predicts f): its own
# likelihood given the factors is constant, and the prior leaves its value
# out, so the EBM score is the posterior mode given the items (and the
# covariates) alone (for f ~~ x, the same for every value of x within a
# response pattern). An ordered fit with exogenous covariates
# that is not conditional on them (conditional.x = FALSE) is refused too,
# since lavPredict() stops on it with "replacement has length zero". The
# case likelihood is also the product of the observed variables'
# univariate probabilities given the factors, which leaves their residual
# covariances out; factor_scores() warns of those (.fs_residual_pairs()).
.fs_check_dummies <- function(est, lv, glab, categorical, conditional_x,
                              ov_x, aux = character(0),
                              ordered = character(0), data = NULL,
                              group_equal = character(0)) {
  if (categorical && !conditional_x && length(ov_x)) {
    .dmar_refuse_argument(c(
      "lavaan::lavPredict() cannot score a fit with ordered items and ",
      "exogenous covariates (", paste(ov_x, collapse = ", "), ") that ",
      "is not conditional on them (conditional.x = FALSE). Refit with ",
      "conditional.x = TRUE, lavaan's default for ordered items, to ",
      "score."), arg = NULL)
  }
  listed <- function(x) {
    if (length(x) < 3L) return(paste(x, collapse = " and "))
    paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  }
  bartlett_issues <- rep(NA_character_, length(est))
  for (g in seq_along(est)) {
    where <- if (length(est) > 1L) sprintf(" in group %s", glab[g]) else ""
    e <- est[[g]]
    if (!length(e$dummy_y)) next
    m <- ncol(e$psi)
    beta <- if (is.null(e$beta)) matrix(0, m, m) else e$beta
    cz <- solve(diag(m) - beta) %*% e$psi
    dimnames(cz) <- dimnames(e$psi)
    tol <- 1e-8 * max(1, abs(cz))
    d_b <- e$dummy_b
    if (length(d_b)) {
      cond <- setdiff(colnames(e$psi), c(d_b, aux))
      bad <- abs(cz[cond, d_b, drop = FALSE]) > tol
      # With continuous items the regression score remains the factor's
      # conditional mean when Cov(factor, disturbance) is zero, even if
      # that disturbance covaries with an observed conditioning variable.
      # The stronger Bartlett condition then fails. Keep the valid method
      # and explicitly withhold the other; ordered scores remain refused.
      if (!categorical && !conditional_x && any(bad) &&
          !any(bad[cond %in% lv, , drop = FALSE])) {
        bartlett_issues[g] <- paste0(
          "Bartlett scores withheld", where, ": the disturbance of ",
          listed(d_b[colSums(bad) > 0L]), " covaries with ",
          listed(cond[rowSums(bad) > 0L]),
          ", a variable the scores condition on. Conditional unbiasedness ",
          "does not hold. Regression scores retain the fitted normal ",
          "conditional mean because the factor-disturbance covariances are zero.")
        bad[,] <- FALSE
      }
      if (any(bad)) {
        vars <- d_b[colSums(bad) > 0L]
        with <- cond[rowSums(bad) > 0L]
        one <- length(vars) == 1L
        # In a multiple group fit, group.equal holds parameters equal by
        # kind, and the remedy moves the variable's intercept, residual
        # variance, and covariances to a latent variable: with equal
        # intercepts (and free means) lavaan frees the new latent
        # variable's mean after the first group, and the equalities of
        # "residuals" and "residual.covariances" no longer reach the
        # variable while those of "lv.variances" and "lv.covariances" now
        # do, so the refit as first worded had another log likelihood and
        # number of degrees of freedom (-4246.416 on 15 against -4246.508
        # on 16 for two schools with equal loadings and intercepts). Fixing
        # the latent variable's intercept at zero restores the model in
        # the first case; in the others the equalities must be carried
        # over by hand, which depends on the identification (under
        # std.lv = TRUE with "lv.variances", the latent variable's
        # variance is fixed at 1 in every group, and freeing it takes more
        # than group.partial), so the message asks for them and for the
        # check.
        multi <- length(est) > 1L
        fix_int <- multi && "intercepts" %in% group_equal &&
          !("means" %in% group_equal)
        fix_lab <- multi && any(c("residuals", "residual.covariances",
                                  "lv.variances", "lv.covariances") %in%
                                  group_equal)
        across <- if (fix_int || fix_lab) {
          paste0(
            ", once it carries this fit's constraints across the groups: ",
            "group.equal holds parameters equal by kind, and the variable's ",
            "intercept, residual variance, and covariances become a latent ",
            "variable's, so ",
            paste(c(if (fix_int) {
              sprintf(paste0(
                "fix the latent variable's intercept at zero in every group ",
                "(%s_lv ~ 0*1), which lavaan otherwise frees after the ",
                "first group when intercepts are held equal"), vars[1L])
            }, if (fix_lab) {
              paste0(
                "give the variable's residual variance and covariances the ",
                "equalities group.equal placed on them here and no others ",
                "(with labels, or group.partial), which the refitted model's ",
                "log likelihood and degrees of freedom confirm")
            }), collapse = ", and "))
        } else ""
        .dmar_refuse_argument(
          sprintf(paste0(
       "The disturbance of %s, %s regressed on other variables, ",
       "covaries with %s%s%s. lavaan's scores treat an observed variable ",
       "regressed on others as a measure of the latent variables it ",
       "depends on, with its disturbance as the measurement error, which ",
       "then gives %s, so the fit is not scored. Carried as a single ",
       "indicator latent variable (%s_lv =~ %s, with %s ~~ 0*%s), such a ",
       "variable is scored correctly: written so, lavaan fixes one scale ",
       "constraint (the loading, or with std.lv = TRUE the latent ",
       "variable's variance), and the refitted model has this fit's log ",
       "likelihood and degrees of freedom%s."),
       listed(vars), if (one) "an observed variable" else
         "observed variables",
       listed(with), where,
       if (any(with %in% lv)) {
         paste0(" (the variable predicts a factor, directly or through ",
                "others, or its residual covaries with one)")
       } else {
         paste0(", ", if (sum(!(with %in% lv)) == 1L) "a variable" else
           "variables", " the scores condition on (an exogenous ",
           "covariate, or an observed variable regressed on nothing)")
       },
       if (any(with %in% lv)) {
         paste0("neither the factors' conditional means nor ",
                "conditionally unbiased Bartlett scores")
       } else {
         paste0("Bartlett scores that are not conditionally unbiased ",
                "given the variables they condition on")
       }, vars[1L], vars[1L], vars[1L], vars[1L], across), arg = NULL)
        # The remedy names no fixed loading: a single indicator latent
        # variable with zero error variance reproduces the observed
        # variable's model only with exactly one scale constraint (Bollen,
        # 1989, ch. 7), and with std.lv = TRUE lavaan fixes the latent
        # variable's (residual) variance at 1 even when a loading is fixed,
        # so =~ 1*z there fixed both and gave a different, misspecified model
        # (its log likelihood -4502.0 against -4262.4 for a verbal factor
        # predicting the counting test, with a residual covariance with age).
      }
    }
    if (!categorical) next
    lv_g <- intersect(lv, colnames(e$psi))
    # The dummies whose row of B reaches a factor directly: lavaan's term of
    # such a variable in a case's likelihood varies with the factors.
    on_lv <- function(v) {
      if (is.null(e$beta)) return(rep(FALSE, length(v)))
      vapply(v, function(r) {
        r %in% rownames(e$beta) && any(e$beta[r, lv_g] != 0)
      }, logical(1))
    }
    if (!is.null(e$gamma)) {
      gy <- e$gamma[intersect(e$dummy_y, rownames(e$gamma)), , drop = FALSE]
      gy <- gy[rowSums(gy != 0) > 0L & on_lv(rownames(gy)), , drop = FALSE]
      if (nrow(gy)) {
        vars <- rownames(gy)
        .dmar_refuse_argument(
          sprintf(paste0(
       "%s %s regressed on an exogenous covariate (%s)%s and also %s ",
       "on a factor directly. lavaan's EBM and ML scores of ordered items ",
       "leave ",
       "such a direct effect out of each case's likelihood, so they are ",
       "those of a case with the covariate at zero rather than this ",
       "model's scores, and the fit is not scored. A model without the ",
       "direct effect, or one that treats the items as continuous, can ",
       "be scored."),
       listed(vars), if (length(vars) == 1L) "is" else "are",
       paste(colnames(gy)[colSums(gy != 0) > 0L], collapse = ", "),
       where, if (length(vars) == 1L) "depends" else "depend"), arg = NULL)
      }
    }
    if (!is.null(e$beta)) {
      bc <- e$beta[, e$dummy_y, drop = FALSE] != 0
      at <- which(bc, arr.ind = TRUE)
      if (nrow(at)) {
        r_v <- rownames(bc)[at[, 1L]]
        d_v <- colnames(bc)[at[, 2L]]
        why <- character(length(r_v))
        why[on_lv(r_v)] <- "the variable predicted depends on a factor"
        why[d_v %in% ordered] <- "the predictor is an ordered item"
        why[r_v %in% lv_g] <- "a factor is predicted"
        bad <- nzchar(why)
        if (any(bad)) {
          reg <- paste(r_v[bad], "~", d_v[bad])
          one <- length(unique(d_v[bad])) == 1L
          .dmar_refuse_argument(
            sprintf(paste0(
         "The regression%s %s%s %s on %s that lavaan holds among the ",
         "latent variables rather than as %s. lavaan's EBM and ML scores ",
         "of ordered items take such a variable as zero where it predicts ",
         "another, so they leave the regression out of each case's ",
         "likelihood (and, where it predicts a factor, out of the ",
         "factor's prior). That leaves the scores unchanged only when the ",
         "variable predicted is an observed variable that depends on no ",
         "factor directly and the predictor is a continuous variable (where ",
         "a case lacks it, one that reaches no factor except through a ",
         "variable the case has observed); here %s, so they are not this ",
         "model's scores, and the fit is not scored. ",
         "A model without the regression, or one that treats the items as ",
         "continuous, can be scored."),
         if (length(reg) == 1L) "" else "s", listed(reg), where,
         if (length(reg) == 1L) "is" else "are",
         if (one) "an observed variable" else "observed variables",
         if (one) "an exogenous covariate" else "exogenous covariates",
         listed(unique(why[bad]))), arg = NULL)
        }
        # The remaining regressions leave out a term that is constant in the
        # factors where the predictor is observed. Where it is missing and
        # the variable predicted observed, the model's likelihood of the
        # variable predicted integrates the predictor out given the
        # predictor's own predictors, so it depends on the factors when a
        # path of B from the predictor reaches a factor through variables
        # the row also lacks (or directly), which lavaan leaves out: such
        # a row refuses the fit. A variable the row has observed on every
        # such path cuts the dependence (by the local Markov property cited
        # above, given it the variables beyond it are independent of the
        # factors), and the row is scored: in m ~ f, w ~ m, s ~ w with w
        # missing and m observed, the likelihood of s given m does not
        # involve f, and lavaan's EBM scores of such rows agreed with a
        # direct maximization of the model's posterior to 5.5e-5, as over
        # complete rows. Each row is judged on its own variables
        # (.fs_reaches_lv()); a variable absent from the data is taken as
        # missing.
        x_g <- if (length(data) >= g) data[[g]] else NULL
        if (!is.null(x_g)) {
          nb <- e$beta != 0
          dum <- intersect(rownames(nb), c(e$dummy_y, e$dummy_x))
          for (k in seq_along(r_v)) {
            if (!all(c(r_v[k], d_v[k]) %in% colnames(x_g))) next
            rows <- which(is.na(x_g[, d_v[k]]) & !is.na(x_g[, r_v[k]]))
            if (!length(rows)) next
            seen <- intersect(dum, colnames(x_g))
            miss <- is.na(x_g[rows, seen, drop = FALSE])
            key <- apply(miss, 1L, function(z) paste(as.integer(z),
                                                     collapse = ""))
            first <- !duplicated(key)
            gap_key <- key[first][vapply(which(first), function(i) {
              lacks <- c(setdiff(dum, seen), seen[miss[i, ]])
              .fs_reaches_lv(nb, d_v[k], lv_g, lacks)
            }, logical(1))]
            n_gap <- sum(key %in% gap_key)
            if (!n_gap) next
            .dmar_refuse_argument(
              sprintf(paste0(
           "The regression %s ~ %s%s is on an observed variable that ",
           "lavaan holds among the latent variables, and %s is missing in ",
           "%d %s the fit keeps where %s is observed and where %s depends ",
           "on a factor through no variable the row has observed ",
           "(directly, or only through variables the row also lacks). ",
           "lavaan's EBM and ML scores of ordered items take %s as zero ",
           "where it predicts %s, which leaves the scores unchanged where ",
           "%s is observed, or where every path from it to the factors ",
           "runs through a variable the row has observed; in those rows ",
           "the model's likelihood of %s depends on the factors through ",
           "%s, and lavaan leaves that out, so the fit is not scored. ",
           "Refit without those rows (missing = \"listwise\", for one), or ",
           "without the regression, to score."),
           r_v[k], d_v[k], where, d_v[k], n_gap,
           if (n_gap == 1L) "row" else "rows", r_v[k], d_v[k], d_v[k],
           r_v[k], d_v[k], r_v[k], d_v[k]), arg = NULL)
          }
        }
      }
    }
    # Every ov.y dummy whose row lavaan keeps (regressed on nothing, or,
    # without B, only on exogenous covariates): for one regressed on
    # nothing, Cov(eta, zeta) is its covariance with the latent variables,
    # V(eta)[, d]; for one regressed only on covariates, the covariance of
    # its disturbance with theirs, Psi[, d].
    d_0 <- setdiff(e$dummy_y, e$dummy_b)
    if (length(d_0) && length(lv_g)) {
      bad <- abs(cz[lv_g, d_0, drop = FALSE]) > tol
      # Those regressed on nothing are named first.
      on_x <- d_0 %in% e$dummy_r
      if (any(bad[, !on_x])) bad[, on_x] <- FALSE
      if (any(bad)) {
        vars <- d_0[colSums(bad) > 0L]
        one <- length(vars) == 1L
        with <- listed(lv_g[rowSums(bad) > 0L])
        fac <- if (sum(rowSums(bad) > 0L) == 1L) "factor" else "factors"
        if (!any(vars %in% e$dummy_r)) {
          .dmar_refuse_argument(
            sprintf(paste0(
         "The observed variable%s %s, regressed on nothing, %s with the ",
         "%s %s%s, and lavaan holds %s among the latent variables rather ",
         "than as %s. lavaan's EBM score of ordered items takes the ",
         "factors' prior from their model implied means, which leave out ",
         "the value of such a variable, so it is the posterior mode given ",
         "the items alone, and the fit is not scored. Written as a ",
         "predictor of the factor (f ~ %s) that covaries with nothing ",
         "else, the variable is an exogenous covariate, which the EBM ",
         "score conditions on."),
         if (one) "" else "s", listed(vars),
         if (one) "covaries" else "covary", fac, with, where,
         if (one) "it" else "them",
         if (one) "an exogenous covariate" else "exogenous covariates",
         vars[1L]), arg = NULL)
        }
        .dmar_refuse_argument(
          sprintf(paste0(
       "The %s %s, regressed on ",
       "exogenous covariates only, %s with the %s %s%s. lavaan's EBM ",
       "score of ordered items takes the factors' prior from their model ",
       "implied means given the covariates, which leave out the value of ",
       "such a variable, and the variable's own likelihood given the ",
       "factors leaves out that covariance, so the score is the posterior ",
       "mode given the items and the covariates alone, and the fit is not ",
       "scored. A model in which that disturbance is uncorrelated with the ",
       "factors, or one that treats the items as continuous, can be ",
       "scored."),
       if (one) "disturbance of the observed variable" else
         "disturbances of the observed variables", listed(vars),
       if (one) "covaries" else "covary", fac, with, where), arg = NULL)
      }
    }
  }
  bartlett_issues
}

# Whether the variable `from` reaches one of the latent variables `lv`
# through the regressions of B (`nb`, B != 0, each row regressed on its
# columns), following a predictor only when it is one of `lv` or one of
# `lacks`, the variables a row has not observed: a path through a variable
# the row has observed is cut there.
.fs_reaches_lv <- function(nb, from, lv, lacks) {
  todo <- from
  done <- character(0)
  while (length(todo)) {
    v <- todo[1L]
    todo <- todo[-1L]
    done <- c(done, v)
    pa <- colnames(nb)[nb[v, ]]
    if (any(pa %in% lv)) return(TRUE)
    todo <- union(todo, setdiff(intersect(pa, lacks), done))
  }
  FALSE
}


# The residual covariances that lavaan's EBM and ML scores of ordered items
# leave out, as "a ~~ b" (with the group when there are several). Each
# case's likelihood there is the product of its observed variables'
# univariate probabilities (or densities) given the factors
# (lav_predict_fy_eta_i()), so a nonzero off-diagonal entry of the
# residual covariance matrix (`theta_eff`, the dummies' block of Psi folded
# in) is left out whenever either variable depends on the factors: the
# scores are those of the model without it. In a model of four Holzinger
# and Swineford tests cut at their quartiles with a residual correlation of
# -0.12 between two of them, the ML scores of the 160 response patterns
# differed from the scores that maximize the model's own likelihood by up
# to 0.11, 15% of their standard error (0.03 at the median), and agreed
# with those that maximize the product to 2e-4. Groups in `skip` (whose cases lavaan does
# not score) are passed over.
.fs_residual_pairs <- function(est, lv, glab, skip = logical(length(est))) {
  out <- character(0)
  for (g in seq_along(est)) {
    if (skip[g]) next
    th <- est[[g]]$theta_eff %||% est[[g]]$theta
    lam <- est[[g]]$lambda_eff %||% est[[g]]$lambda
    v <- rownames(th)
    dep <- rowSums(lam[v, intersect(lv, colnames(lam)), drop = FALSE] !=
                     0) > 0L
    at <- which(th != 0 & upper.tri(th) & outer(dep, dep, `|`),
                arr.ind = TRUE)
    if (!nrow(at)) next
    out <- c(out, paste0(v[at[, 1L]], " ~~ ", v[at[, 2L]],
                         if (length(est) > 1L) {
                           sprintf(" (group %s)", glab[g])
                         } else ""))
  }
  out
}


# ---------------------------------------------------------------------------
# (B) Model structure
# ---------------------------------------------------------------------------

# lavaan returns a bare matrix (or a list of model matrices) for a single
# group and a list of them for several; this wraps the single-group case so
# that every per-group loop indexes a list. `inner` names an element that a
# single group's object has directly (lavInspect(fit, "est") gives a list of
# matrices with $lambda for one group, a list of such lists for several).
.fs_as_list <- function(x, inner = NULL) {
  if (!is.null(inner)) {
    if (!is.null(x[[inner]])) list(x) else x
  } else if (is.list(x) && !is.data.frame(x)) {
    x
  } else {
    list(x)
  }
}

# The items of each factor, per group: the observed variables reached from
# the factor through =~ relations, so a higher order factor's items are
# those of the factors it measures. A loading the model holds at zero
# (.pt_zero_rows() in R/reliability_internals.R: 0*x, or a*x with a == 0,
# directly or through a chain of equalities) does not make an item.
# Returns a list over groups of named lists over factors.
.fs_items <- function(fit, lv, ngroups) {
  pt <- lavaan::parTable(fit)
  ov <- lavaan::lavNames(fit, "ov")
  zero <- .pt_zero_rows(pt)
  lapply(seq_len(ngroups), function(g) {
    m <- pt[pt$op == "=~" & pt$group == g & !zero, c("lhs", "rhs")]
    reach <- function(f, seen = character(0)) {
      if (f %in% seen) return(character(0))
      kids <- m$rhs[m$lhs == f]
      c(intersect(kids, ov),
        unlist(lapply(setdiff(kids, ov), reach, seen = c(seen, f))))
    }
    out <- lapply(lv, function(f) unique(reach(f)))
    names(out) <- lv
    out
  })
}

# The entries of the model matrices that the model holds at zero
# (.pt_zero_rows()), and the matrices with those entries set to exactly
# zero. A row fixed at zero is zero already, but lavaan returns a row
# constrained to zero through a chain of equalities (a == b with b == 0)
# at its optimizer's rounding, -1.4e-15 for one loading, and every
# structural test of factor_scores() (which factors a pattern's observed
# items load on, which factors have items of their own, which dummies are
# regressed, which residual covariances are left out) reads the matrices
# by value, as does lavPredict()'s construction of a higher order factor's
# Bartlett score (which looks for a column of zeros among the loadings).
# The constraint, not the estimate, decides: such a row describes the
# same model as its absence, with the same chi square and degrees of
# freedom, so setting it to the value the constraint fixes makes every
# spelling of the model (0*x, a*x with a == 0, a chain, or no row) reach
# the same scores and reasons. No tolerance on the estimate is used, since
# a genuinely small free loading must stay a loading.
#
# .fs_held_cells() places each such row in the cell of the model matrices
# that holds it in this fit. Only a free row can differ from zero (a row
# fixed at zero is zero in the matrices already), and lavaan's own map of
# the fit's free parameters places it: fit@Model@x.free.idx gives, for
# each model matrix, the number of the free parameter (the `free` column of
# the parameter table) that each of its free cells holds, and
# fit@Model@m.free.idx the cells themselves, in the order of the fit's own
# matrices (fit@Model@GLIST, which lavInspect(fit, "est") returns with
# names, block by block). Placing a row by the parameter table alone
# (lavMatrixRepresentation()) orders the observed variables as the model's
# syntax does, which is not the fit's order when it was made with
# ov_order = "data": the cell found then belonged to another variable,
# and a real loading (24.1, of the counting test on the speed factor) was
# set to zero, which moved the scores by up to 3.8. The cells are returned
# as a data frame (block, mat, mm, the index of the matrix in GLIST, row,
# col, and est, the row's estimate), or NULL when there are none.
# .fs_zero_held() sets them to zero in each group's model matrices (`est`,
# lavInspect(fit, "est") through .fs_as_list()), a symmetric matrix in
# both of its cells. .fs_zero_held_fit() does the same in the model
# matrices a copy of the fit carries for lavPredict(), and returns the fit
# unchanged when no held entry differs from zero. Both set a cell to zero
# only when it holds the row's own estimate, so a cell placed wrongly is
# left as it is rather than erasing another parameter. The copy keeps the
# model implied moments the fit cached (fit@implied), which lavPredict()
# reads for the model implied covariance matrix in its regression and
# Bartlett weights rather than recomputing it from the matrices. They
# differ from the moments of the zeroed matrices by terms linear in such
# residues, since a residue multiplies the other loadings and the factors'
# covariances in Sigma (2.1e-14 for the nine tests with a cross-loading
# held at zero through a chain), far below the precision of the scores.
.fs_held_cells <- function(fit) {
  pt <- lavaan::parTable(fit)
  held <- .pt_zero_rows(pt) & pt$free > 0L
  if (!any(held)) return(NULL)
  mod <- fit@Model
  glist <- mod@GLIST
  block <- rep(seq_along(mod@nmat), mod@nmat)
  cells <- lapply(seq_along(glist), function(mm) {
    x_idx <- mod@x.free.idx[[mm]]
    hit <- x_idx %in% pt$free[held]
    if (!any(hit)) return(NULL)
    at <- arrayInd(mod@m.free.idx[[mm]][hit], dim(glist[[mm]]))
    data.frame(block = block[mm], mat = names(glist)[mm], mm = mm,
               row = at[, 1L], col = at[, 2L],
               est = pt$est[match(x_idx[hit], pt$free)],
               stringsAsFactors = FALSE)
  })
  cells <- do.call(rbind, cells)
  if (!NROW(cells)) NULL else cells
}

# Whether a cell holds a held row's estimate (to the rounding of a copy).
.fs_holds <- function(value, est) {
  isTRUE(abs(value - est) <= sqrt(.Machine$double.eps) * max(1, abs(est)))
}

.fs_zero_cell <- function(m, mat, row, col) {
  m[row, col] <- 0
  if (mat %in% c("theta", "psi")) m[col, row] <- 0
  m
}

.fs_zero_held <- function(est, cells) {
  for (i in seq_len(NROW(cells))) {
    b <- cells$block[i]
    m <- cells$mat[i]
    if (b > length(est) || is.null(est[[b]][[m]])) next
    if (!.fs_holds(est[[b]][[m]][cells$row[i], cells$col[i]],
                   cells$est[i])) next
    est[[b]][[m]] <- .fs_zero_cell(est[[b]][[m]], m, cells$row[i],
                                   cells$col[i])
  }
  est
}

.fs_zero_held_fit <- function(fit, cells) {
  if (!NROW(cells)) return(fit)
  glist <- fit@Model@GLIST
  changed <- FALSE
  for (i in seq_len(nrow(cells))) {
    mm <- cells$mm[i]
    if (mm > length(glist) || names(glist)[mm] != cells$mat[i]) next
    value <- glist[[mm]][cells$row[i], cells$col[i]]
    if (value == 0 || !.fs_holds(value, cells$est[i])) next
    glist[[mm]] <- .fs_zero_cell(glist[[mm]], cells$mat[i], cells$row[i],
                                 cells$col[i])
    changed <- TRUE
  }
  if (changed) fit@Model@GLIST <- glist
  fit
}

# Auxiliary variables carried as saturated correlates (Graham, 2003): an
# observed variable that indicates no factor, takes part in no regression,
# and covaries (~~) with at least one indicator. Under that
# parameterization the auxiliary is uncorrelated with every factor and
# correlated with the item residuals, so a score that conditions on it
# subtracts from each item the part the auxiliary predicts, which includes
# whatever the auxiliary shares with the factor. The scores are then
# computed from the indicators alone. A covariate that predicts a factor
# (op ~) is told apart by its operator and left to lavPredict(). A row the
# model holds at zero (.pt_zero_rows(): 0*x, or a*x with a == 0, directly
# or through a chain of equalities) is read as absent, as everywhere in
# factor_scores(): it describes the same model as no row, so an auxiliary
# written with a loading held at zero (f =~ 0*a) is set aside as one
# written without it is. Read from lavNames(fit, "ov.ind"), such a
# variable was an indicator, the fit was scored with it conditioned on,
# and its scores differed from those of the same model written without
# the row by up to 1.9.
.fs_auxiliary <- function(fit) {
  pt <- lavaan::parTable(fit)
  pt <- pt[!.pt_zero_rows(pt), , drop = FALSE]
  ov <- lavaan::lavNames(fit, "ov")
  ind <- intersect(pt$rhs[pt$op == "=~"], ov)
  in_model <- unique(c(pt$lhs[pt$op %in% c("~", "=~")],
                       pt$rhs[pt$op %in% c("~", "=~")]))
  candidates <- setdiff(ov, c(ind, in_model))
  cov_rows <- pt[pt$op == "~~" & pt$lhs != pt$rhs, c("lhs", "rhs")]
  keep <- vapply(candidates, function(v) {
    any((cov_rows$lhs == v & cov_rows$rhs %in% ind) |
          (cov_rows$rhs == v & cov_rows$lhs %in% ind))
  }, logical(1))
  candidates[keep]
}

# The dummy latent variables of one group's model matrices (lavInspect(fit,
# "est")), and the two matrices built from them (see the notation above):
# `est_g` is returned with `dummy_x` (the exogenous covariates carried as
# latent variables), `dummy_y` (lavaan's other dummies, its ov.y dummies),
# `dummy_r` (those of them regressed on something, with a nonzero row of B
# or, in a fit conditional on the covariates, of Gamma), `dummy_b` (those
# of `dummy_r` whose row lav_lisrel_lambda() replaces by its row of B for
# lavPredict(): all of them in a model with a B matrix, none without one),
# `lambda_eff` (the loading matrix with the row of each variable in
# `dummy_b` replaced by its row of B, as lav_lisrel_lambda() builds it), and
# `theta_eff` (the residual covariance matrix with the dummies' block of
# Psi folded in, as lav_lisrel_theta() builds it for the post-estimation
# check and the EBM and ML scores). A dummy latent variable is a column of
# Lambda that is not one of the model's latent variables (`lv`); lavaan
# names it after its observed variable, the row it measures. Exogenous
# covariates are in the model matrices only when the fit is not
# conditional on them.
#
# An ov.y dummy regressed on nothing (a criterion that covaries with a
# factor, f ~~ z; a covariate that predicts a factor and also covaries with
# another variable, so that lavaan does not count it among the exogenous
# covariates; an auxiliary variable) keeps its unit row, as an exogenous
# covariate's: its value is its latent variable's, Cov(eta, y_d) is
# V(eta)[, d] exactly, and the scores condition on it. lav_lisrel_lambda()
# replaces its row by its row of B as well whenever B exists, a row of
# zeros, which drops the variable's covariances with the factors from the
# weights: lavPredict()'s regression score then missed the factor's
# conditional mean by up to 0.03 for f ~ x with x ~~ a and 1.9 for f ~ x
# with f ~~ z, in the Holzinger and Swineford examples of the tests.
# Without B lavaan keeps the unit row of every dummy, of one regressed on
# nothing and of one regressed only on exogenous covariates (its row of
# Gamma) alike, and the scores condition on each.
.fs_dummies <- function(est_g, lv, ov_x = character(0)) {
  lambda <- est_g$lambda
  dummy <- intersect(setdiff(colnames(lambda), lv), rownames(lambda))
  est_g$dummy_x <- intersect(dummy, ov_x)
  est_g$dummy_y <- setdiff(dummy, ov_x)
  regressed <- function(m) {
    if (is.null(m) || !length(est_g$dummy_y)) {
      return(rep(FALSE, length(est_g$dummy_y)))
    }
    rows <- match(est_g$dummy_y, rownames(m))
    out <- rep(FALSE, length(est_g$dummy_y))
    out[!is.na(rows)] <- rowSums(m[rows[!is.na(rows)], , drop = FALSE] !=
                                   0) > 0L
    out
  }
  est_g$dummy_r <- est_g$dummy_y[regressed(est_g$beta) |
                                   regressed(est_g$gamma)]
  est_g$dummy_b <- if (is.null(est_g$beta)) character(0) else est_g$dummy_r
  lambda_eff <- lambda
  if (length(est_g$dummy_b)) {
    lambda_eff[est_g$dummy_b, ] <- est_g$beta[est_g$dummy_b,
                                              colnames(lambda), drop = FALSE]
  }
  theta_eff <- est_g$theta
  if (length(dummy)) {
    theta_eff[dummy, dummy] <- est_g$psi[dummy, dummy, drop = FALSE]
  }
  est_g$lambda_eff <- lambda_eff
  est_g$theta_eff <- theta_eff
  est_g
}

# The model matrices of one group, completed: the latent covariance and
# mean, V(eta) and E(eta), and the model implied covariance and mean of the
# observed variables. Sigma and E(y) come from lavaan's matrices as they
# are; `lambda` is returned as lavPredict() uses it (`lambda_eff`, from
# .fs_dummies(), when present), for the weights of the scores.
#
# Without a mean structure the means of the observed variables are not
# modeled, and they are the sample means of the rows the fit used, as in
# lavaan. The latent variables' means then follow lavaan's rule
# (lav_lisrel_alpha0()): the intercepts of the latent variables are zero,
# and those of the dummy latent variables are the values that make the
# model implied means of their observed variables the sample means. A
# factor predicted by a covariate (f ~ x) then has the mean beta times the
# covariate's mean, which is lavInspect(fit, "mean.lv"), the mean of the
# same fit with a mean structure, and the mean FIML gives. Taking every
# intercept as zero instead put the factor's mean at zero while the
# covariate was centered at its mean, which shifted every score of such a
# factor, in a fit by pairwise deletion or of the available cases and with
# auxiliary variables set aside, by -beta times the covariate's mean
# (-0.26 in the example of the tests).
.fs_moments <- function(est_g, data_g) {
  lambda <- est_g$lambda
  m <- ncol(lambda)
  beta <- if (is.null(est_g$beta)) matrix(0, m, m) else est_g$beta
  ib_inv <- solve(diag(m) - beta)
  veta <- ib_inv %*% est_g$psi %*% t(ib_inv)
  sigma <- lambda %*% veta %*% t(lambda) + est_g$theta
  if (is.null(est_g$nu)) {
    ey <- colMeans(data_g[, rownames(lambda), drop = FALSE], na.rm = TRUE)
    alpha <- numeric(m)
    dummy <- c(est_g$dummy_y, est_g$dummy_x)
    if (length(dummy)) {
      d_lv <- match(dummy, colnames(lambda))
      d_ov <- match(dummy, rownames(lambda))
      alpha[d_lv] <- solve((lambda %*% ib_inv)[d_ov, d_lv, drop = FALSE],
                           ey[d_ov])
    }
  } else {
    alpha <- if (is.null(est_g$alpha)) numeric(m) else as.numeric(est_g$alpha)
  }
  eeta <- as.numeric(ib_inv %*% alpha)
  if (!is.null(est_g$nu)) {
    ey <- as.numeric(est_g$nu) + as.numeric(lambda %*% eeta)
  }
  dimnames(veta) <- list(colnames(lambda), colnames(lambda))
  dimnames(sigma) <- list(rownames(lambda), rownames(lambda))
  names(eeta) <- colnames(lambda)
  names(ey) <- rownames(lambda)
  list(lambda = est_g$lambda_eff %||% lambda, theta = est_g$theta,
       veta = veta, eeta = eeta, sigma = sigma, ey = ey)
}


# ---------------------------------------------------------------------------
# (C) Per-pattern algebra
# ---------------------------------------------------------------------------

# Missing data patterns of a logical matrix of observed entries (rows are
# cases). Patterns are numbered from the most complete down (pattern 1 is
# the complete one when any case is complete), ties in the order the
# patterns first appear. Returns the pattern of every row and the patterns
# themselves as a logical matrix.
.fs_patterns <- function(obs) {
  key <- apply(obs, 1L, function(r) paste(as.integer(r), collapse = ""))
  first <- !duplicated(key)
  keys <- key[first]
  n_obs <- rowSums(obs[first, , drop = FALSE])
  ord <- order(-n_obs, seq_along(keys))
  keys <- keys[ord]
  observed <- obs[first, , drop = FALSE][ord, , drop = FALSE]
  rownames(observed) <- NULL
  list(id = match(key, keys), observed = observed)
}

# Which factors the observed items of a pattern identify: given their rows
# `lo` of the loading matrix (items by factors), a factor is identified when
# it has a nonzero loading among them and its unit vector e_k lies in the
# row space of `lo`, that is, when the projector P onto that row space has
# P_kk = 1. Only then is the factor's value a function of the items'
# expected values, Lambda_o eta, so only then does an estimate of it that is
# unbiased given the factor values exist (an estimable function, in linear
# model terms). With simple structure every factor with an observed item is
# identified; a pattern whose observed items load on more factors than they
# can separate (one item loading on two factors, say) is not. The rank and
# the projector come from the singular value decomposition with the
# tolerance of MASS::ginv(), the generalized inverse lavaan falls back on.
.fs_identified <- function(lo, tol = sqrt(.Machine$double.eps)) {
  has <- colSums(lo != 0) > 0
  ok <- has
  if (any(has)) {
    s <- svd(lo[, has, drop = FALSE], nu = 0L)
    r <- sum(s$d > tol * max(s$d))
    if (r < sum(has)) {
      v <- s$v[, seq_len(r), drop = FALSE]
      ok[has] <- 1 - rowSums(v^2) < tol
    }
  }
  ok
}

# Everything the scores need for one pattern of observed items `o`
# (logical over the rows of mom$lambda): the regression and Bartlett weight
# matrices, the variances behind the standard errors, and the two
# determinacies, for every latent variable in mom$veta. The Bartlett score
# exists only for a latent variable that the pattern's observed items
# identify (`bart_ok`, from .fs_identified()): one with an observed item
# (`has`) whose value the items separate from the other factors'. When
# they all are, the information matrix of the factors with an observed item
# is inverted, and lavaan's generalized inverse over every factor gives the
# same. When some are not (an item shared by two factors and nothing else
# observed of either), that matrix is singular and lavaan's generalized
# inverse returns the minimum norm solution, which is not an estimate of
# those factors' values and whose "standard error" is not the error of
# anything; their Bartlett score is NA here. For an identified factor in
# such a pattern the generalized inverse gives the unique estimate and its
# error variance, (Lambda' Sigma^{-1} Lambda)^- - V(eta) evaluated for that
# factor, since both are invariant to the choice of generalized inverse for
# an estimable function, so it is kept and agrees with lavaan. A variance
# below `floor` gives an NA standard error, as lavPredict() does.
.fs_pattern_info <- function(mom, o, floor = 1e-5) {
  m <- ncol(mom$veta)
  lv <- colnames(mom$veta)
  na_m <- stats::setNames(rep(NA_real_, m), lv)
  no_m <- stats::setNames(rep(FALSE, m), lv)
  out <- list(reg_w = NULL, reg_var = na_m, reg_se = na_m, det_reg = na_m,
              bart_w = NULL, bart_var = na_m, bart_se = na_m,
              det_bart = na_m, has = no_m, bart_ok = no_m)
  if (!any(o)) return(out)
  veta <- mom$veta
  lo <- mom$lambda[o, , drop = FALSE]
  s_oo <- mom$sigma[o, o, drop = FALSE]
  s_inv <- tryCatch(solve(s_oo), error = function(e) NULL)
  if (is.null(s_inv)) return(out)

  a_r <- veta %*% t(lo) %*% s_inv
  post <- veta - a_r %*% lo %*% veta
  cov_r <- diag(a_r %*% lo %*% veta)
  var_r <- diag(a_r %*% s_oo %*% t(a_r))
  det_r <- cov_r / sqrt(var_r * diag(veta))
  det_r[!is.finite(det_r)] <- NA_real_
  out$reg_w <- a_r
  out$reg_var <- stats::setNames(diag(post), lv)
  out$det_reg <- stats::setNames(det_r, lv)

  has <- colSums(lo != 0) > 0
  ok <- .fs_identified(lo)
  out$has <- stats::setNames(has, lv)
  out$bart_ok <- stats::setNames(ok, lv)
  if (any(has)) {
    lh <- lo[, has, drop = FALSE]
    info <- t(lh) %*% s_inv %*% lh
    info_inv <- if (all(ok[has])) {
      tryCatch(solve(info), error = function(e) NULL)
    } else {
      MASS::ginv(info)
    }
    if (!is.null(info_inv)) {
      a_b <- matrix(NA_real_, m, sum(o), dimnames = list(lv, NULL))
      a_b[has, ] <- info_inv %*% t(lh) %*% s_inv
      a_b[!ok, ] <- NA_real_
      bart_var <- diag(info_inv) - diag(veta)[has]
      bart_var[!ok[has]] <- NA_real_
      out$bart_w <- a_b
      out$bart_var[has] <- bart_var
      if (any(ok)) {
        ab_h <- a_b[ok, , drop = FALSE]
        cov_b <- diag(ab_h %*% lo %*% veta[, ok, drop = FALSE])
        var_b <- diag(ab_h %*% s_oo %*% t(ab_h))
        det_b <- cov_b / sqrt(var_b * diag(veta)[ok])
        det_b[!is.finite(det_b)] <- NA_real_
        out$det_bart[ok] <- det_b
      }
    }
  }
  out$reg_se <- sqrt(ifelse(out$reg_var < floor, NA_real_, out$reg_var))
  out$bart_se <- sqrt(ifelse(out$bart_var < floor, NA_real_, out$bart_var))
  out
}


# ---------------------------------------------------------------------------
# (D) lavaan's scores
# ---------------------------------------------------------------------------

# Run lavPredict() for one method with its standard errors and align the
# result to the rows of each group's data: a list over groups of
# list(score, se), each an n_g x K matrix over the reported factors. Under
# complete data or listwise deletion lavaan returns one row of standard
# errors per group (they do not vary across cases), which is recycled; a
# matrix that omits the rows with no observed item is padded. A row without
# a score gets no standard error, so a recycled row never stands beside a
# score that does not exist. `optim_method` (lavaan 0.7-3, which ordered
# items require) is passed on only when given, for the ML scores of ordered
# items, which lavaan finds case by case by optimization.
.fs_lavpredict <- function(fit, method, lv, data_list, optim_method = NULL) {
  args <- list(fit, type = "lv", method = method, se = "standard",
               assemble = FALSE)
  if (!is.null(optim_method)) args$optim_method <- optim_method
  res <- do.call(lavaan::lavPredict, args)
  se <- attr(res, "se", exact = TRUE)
  res <- .fs_as_list(res)
  if (!is.null(se)) se <- .fs_as_list(se)
  lapply(seq_along(res), function(g) {
    n_g <- nrow(data_list[[g]])
    sc <- as.matrix(res[[g]])
    if (nrow(sc) != n_g || !all(lv %in% colnames(sc))) {
      .dmar_refuse_argument(c(
        "lavaan::lavPredict() returned scores that cannot be aligned ",
        "with the rows of the fit."), arg = NULL)
    }
    sc <- sc[, lv, drop = FALSE]
    s <- matrix(NA_real_, n_g, length(lv), dimnames = list(NULL, lv))
    if (length(se) >= g && !is.null(se[[g]])) {
      e <- as.matrix(se[[g]])
      if (is.null(colnames(e))) colnames(e) <- colnames(res[[g]])
      nonempty <- rowSums(!is.na(data_list[[g]])) > 0L
      if (nrow(e) == 1L) {
        e <- e[rep(1L, n_g), , drop = FALSE]
      } else if (nrow(e) == sum(nonempty) && nrow(e) != n_g) {
        pad <- matrix(NA_real_, n_g, ncol(e),
                      dimnames = list(NULL, colnames(e)))
        pad[nonempty, ] <- e
        e <- pad
      }
      if (nrow(e) == n_g && all(lv %in% colnames(e))) {
        s <- e[, lv, drop = FALSE]
      }
    }
    s[!is.finite(sc)] <- NA_real_
    list(score = sc, se = s)
  })
}

# Ordered items: which factors of each case have a finite, unique maximum
# likelihood score. lavaan's ML score maximizes over the factor values eta
# the log likelihood of the case's responses, the sum over its observed
# items of log P(y_i | eta), less a ridge of 1e-5 (eta - mu)'(eta - mu) / 2
# that keeps its optimizer finite (lav_predict_eta_ebm_ml()); the EBM score
# adds the factors' normal prior and is always finite. An item depends on
# eta only through its latent response mean nu_i + lambda_i'eta, with
# lambda_i its row of the loading matrix lavaan's likelihood uses
# (lav_lisrel_eyetax()): for an observed variable regressed on others (an
# item with a direct effect, y1 ~ f2, or an outcome the factor predicts,
# z ~ f), lavaan holds the variable as a dummy latent variable whose
# loadings are its row of B, which `lambda_eff` (.fs_dummies()) puts in its
# row. Read from Lambda alone, such an item had no loading at all, and a
# case whose other items sat at an end was told it had no finite ML score
# although the item, in a middle category, bounds the likelihood. Exogenous
# covariates shift the mean by a constant, and the regressions among the
# latent variables do not enter. The probability of the lowest category
# rises as lambda_i'eta falls, that of the highest as it rises, and that of
# any other category, like the density of a continuous item, falls toward
# zero as lambda_i'eta moves away in either direction. So moving eta along
# a direction d raises the likelihood without end, or leaves it unchanged,
# exactly when d lies in the cone
#   C = {d : -lambda_i'd >= 0 for an item in its lowest category,
#             lambda_i'd >= 0 for an item in its highest category,
#             lambda_i'd  = 0 for every other observed item},
# and, the log likelihood being concave in eta, factor k has a finite and
# unique ML score exactly when every d in C has d_k = 0. By Farkas' lemma,
# some d in C has d_k < 0 if and only if e_k is not a nonnegative
# combination of the rows that define C, and some has d_k > 0 if and only
# if -e_k is not (.fs_in_cone()). When each observed item loads on one
# factor, C is a product of one-factor cones and the rule reduces to the
# familiar one: a factor whose observed items all sit at the end its
# loadings point to (all in the lowest category, or all in the highest, for
# positive loadings) has no finite ML score, and an observed continuous
# item of the factor makes it finite. lavaan returns where its ridge stops
# the optimizer for such a factor (or NA for every factor of the case, when
# the optimizer stops on its iteration limit), so these cases are flagged.
# Returns an n_g x K matrix of codes: "extreme" for the one-factor rule;
# "flat" when the observed items do not separate the factor from another
# (e_k outside their row space, .fs_identified(): the likelihood is
# constant along a direction that moves it); "joint" when the likelihood
# rises without bound along a direction that moves it with other factors
# whose items it shares; and "" otherwise, a finite and unique ML score.
# That includes a factor that shares observed items with a factor coded
# "extreme" or "joint": as the other factor moves out, the shared items'
# probabilities tend to 1 whatever this factor's value, so the likelihood's
# maximum in this factor is that of the case's other items, and lavaan's
# optimizer, whose ridge stops the other factor where those probabilities
# differ from 1 by about 1e-5, reaches it within its tolerance (within
# 1.2e-4, and its standard error within 6e-5, of a brute force maximization
# in the cross-loading example of the tests, closer than on rows without a
# divergent factor). Only factors with items of their own (`own`) and an
# observed item in the row are judged. Category codes run from 1 to the
# number of thresholds plus 1 in lavInspect(fit, "data").
.fs_ml_finite <- function(x_g, est_g, lv, ordered, own) {
  K <- length(lv)
  out <- matrix("", nrow(x_g), K, dimnames = list(NULL, lv))
  lam <- (est_g$lambda_eff %||% est_g$lambda)[, lv, drop = FALSE]
  lam <- lam[rownames(lam) %in% colnames(x_g) & rowSums(lam != 0) > 0, ,
             drop = FALSE]
  if (!nrow(lam)) return(out)
  items <- rownames(lam)
  x <- x_g[, items, drop = FALSE]
  tau_names <- rownames(est_g$tau)
  # The state of each item in each row: 0 not observed, 1 lowest category,
  # 3 highest category, 2 any other observed value (a middle category, or a
  # continuous item).
  st <- matrix(2L, nrow(x), ncol(x))
  for (j in which(items %in% ordered)) {
    top <- sum(startsWith(tau_names, paste0(items[j], "|"))) + 1
    st[!is.na(x[, j]) & x[, j] == 1, j] <- 1L
    st[!is.na(x[, j]) & x[, j] == top, j] <- 3L
  }
  st[is.na(x)] <- 0L
  key <- apply(st, 1L, paste, collapse = ".")
  first <- which(!duplicated(key))
  codes <- matrix("", length(first), K)
  for (u in seq_along(first)) {
    s <- st[first[u], ]
    o <- s > 0L
    if (!any(o)) next
    lo <- lam[o, , drop = FALSE]
    so <- s[o]
    nz <- lo != 0
    judge <- colSums(nz) > 0 & own
    if (!any(judge)) next
    code <- rep("", K)
    # The one-factor rule: every observed item of the factor at the end its
    # loading on the factor points to (`so` recycles down the columns).
    up <- (so == 3L & lo > 0) | (so == 1L & lo < 0)
    down <- (so == 1L & lo > 0) | (so == 3L & lo < 0)
    ext <- judge & (colSums(nz & !up) == 0L | colSums(nz & !down) == 0L)
    code[ext] <- "extreme"
    if (any(rowSums(nz) > 1L)) {
      ok <- .fs_identified(lo)
      code[judge & !ext & !ok] <- "flat"
      b <- rbind(-lo[so == 1L, , drop = FALSE], lo[so == 3L, , drop = FALSE],
                 lo[so == 2L, , drop = FALSE], -lo[so == 2L, , drop = FALSE])
      for (k in which(judge & !ext & ok)) {
        e_k <- replace(numeric(K), k, 1)
        if (!.fs_in_cone(b, e_k) || !.fs_in_cone(b, -e_k)) code[k] <- "joint"
      }
    }
    codes[u, ] <- code
  }
  out[] <- codes[match(key, key[first]), , drop = FALSE]
  out
}

# Is the vector `v` a nonnegative combination of the rows of `b`? By
# nonnegative least squares, min over y >= 0 of ||b'y - v||, computed by the
# active set method of Lawson and Hanson (1974, Chapter 23), which is finite
# and exact up to rounding; `v` lies in the cone when the residual is zero.
# A row added to the passive set whose coefficient is not positive (which
# exact arithmetic rules out) is set aside until the solution moves, so the
# method cannot cycle; the iteration count is capped as well.
.fs_in_cone <- function(b, v, tol = 1e-8) {
  a <- t(b)
  p <- ncol(a)
  if (!p) return(sqrt(sum(v^2)) <= tol)
  eps <- 1e-12 * max(1, sqrt(max(colSums(a^2))))
  y <- numeric(p)
  passive <- logical(p)
  barred <- logical(p)
  for (iter in seq_len(5L * p + 10L)) {
    w <- drop(crossprod(a, v - a %*% y))
    cand <- which(!passive & !barred & w > eps)
    if (!length(cand)) break
    j <- cand[which.max(w[cand])]
    passive[j] <- TRUE
    repeat {
      z <- numeric(p)
      if (!any(passive)) break
      z[passive] <- qr.coef(qr(a[, passive, drop = FALSE]), v)
      z[is.na(z)] <- 0
      if (all(z[passive] > 0)) break
      if (z[j] <= 0 && y[j] == 0) {
        passive[j] <- FALSE
        barred[j] <- TRUE
        z <- y
        break
      }
      neg <- passive & z <= 0
      alpha <- min(y[neg] / (y[neg] - z[neg]))
      y <- y + alpha * (z - y)
      passive <- passive & y > 0
    }
    if (!identical(z, y)) barred[] <- FALSE
    y <- z
  }
  sqrt(sum((v - a %*% y)^2)) <= tol
}


# ---------------------------------------------------------------------------
# (E) Rows given to the fit
# ---------------------------------------------------------------------------

# Does `data` hold the rows the model was fitted to, in their order? The
# scores come from the fit and the counts of items answered, the group of
# an unscored row, and the reasons from `data`, so the two must be the
# same rows. `data` must hold every observed variable of the model; its
# row count must be the fit's original count (lavaan's norig, which
# includes rows deleted listwise) plus the rows whose group is missing
# (lavaan drops them before counting) plus, for a fit restricted to some
# groups with group.label, the rows of the other groups; and in the rows
# the fit used, the values must be those it used. A continuous variable
# agrees to 1e-8 (or, under std.ov = TRUE, correlates perfectly); an
# ordered variable is held by lavaan as category codes, so its supplied
# values must map one to one and in order onto the codes. Missing values
# must fall in the same places. Returns the row count, the rows without a
# group, the rows of groups the fit left out, and the group variable.
.fs_data_check <- function(fit, data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    .dmar_refuse_argument(c(
      "'data' must be the data frame of rows given to the model, in ",
      "their order."), arg = "data")
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  ov <- lavaan::lavNames(fit, "ov")
  absent <- setdiff(ov, names(data))
  if (length(absent)) {
    .dmar_refuse_argument(
      sprintf(paste0("'data' does not hold the model's observed ",
                     "variable%s %s, so it is not the data the model ",
                     "was fitted to."),
              if (length(absent) == 1L) "" else "s",
              paste(absent, collapse = ", ")), arg = "data")
  }
  index <- .fs_as_list(lavaan::lavInspect(fit, "case.idx"))
  norig <- as.integer(unlist(lavaan::lavInspect(fit, "norig")))
  group <- as.character(lavaan::lavInspect(fit, "group"))
  has_group <- length(group) == 1L && nzchar(group)
  unassigned <- integer(0)
  unselected <- integer(0)
  if (has_group && !(group %in% names(data))) {
    .dmar_refuse_argument(
      sprintf(paste0("'data' does not hold the group variable %s; ",
                     "supply the rows given to the model with it, so ",
                     "that rows the model excluded for a missing group ",
                     "can be identified."), group), arg = "data")
  }
  if (has_group) {
    labels <- as.character(lavaan::lavInspect(fit, "group.label"))
    held <- as.character(data[[group]])
    selected <- length(lavaan::lavInspect(fit, "options")$group.label) > 0L
    outside <- !is.na(held) & !(held %in% labels)
    foreign <- if (selected) character(0) else unique(held[outside])
    if (length(foreign)) {
      .dmar_refuse_argument(
        sprintf(paste0("'data' holds %s %s of the group variable %s, ",
                       "which the model was not fitted to; supply the ",
                       "rows given to the model, in their order."),
                if (length(foreign) == 1L) "the value" else "the values",
                paste(foreign, collapse = ", "), group), arg = "data")
    }
    sizes <- vapply(labels, function(l) sum(held == l, na.rm = TRUE),
                    integer(1))
    placed <- all(vapply(seq_along(index), function(g) {
      all(held[as.integer(index[[g]])] %in% labels[g])
    }, logical(1)))
    if (!identical(unname(sizes), norig) || !placed) {
      .dmar_refuse_argument(
        sprintf(paste0("'data' does not hold the rows the model was ",
                       "fitted to: the rows of its groups (%s) differ ",
                       "from those the model used. Supply the rows ",
                       "given to the model, in their order."), group),
        arg = "data")
    }
    unassigned <- which(is.na(held))
    if (selected) unselected <- which(outside)
    expected <- sum(norig) + length(unassigned) + length(unselected)
  } else {
    expected <- max(c(sum(norig), unlist(index)))
  }
  if (nrow(data) != expected) {
    .dmar_refuse_argument(
      sprintf(paste0("'data' has %d rows, but the model was fitted to ",
                     "%d; supply the rows given to the model, in their ",
                     "order."), nrow(data), expected), arg = "data")
  }
  model_data <- .fs_as_list(lavaan::lavInspect(fit, "data"))
  ordered <- lavaan::lavNames(fit, "ov.ord")
  standardized <- isTRUE(lavaan::lavInspect(fit, "options")$std.ov)
  for (g in seq_along(index)) {
    ids <- as.integer(index[[g]])
    x <- as.matrix(model_data[[g]])
    for (v in colnames(x)) {
      fitted <- as.numeric(x[, v])
      supplied <- data[[v]][ids]
      if (v %in% ordered) {
        s <- if (is.factor(supplied)) as.integer(supplied) else supplied
        same <- all(is.na(s) == is.na(fitted))
        if (same && any(!is.na(fitted))) {
          seen <- !is.na(fitted)
          pairs <- unique(data.frame(code = fitted[seen], value = s[seen],
                                     stringsAsFactors = FALSE))
          same <- !anyDuplicated(pairs$code) &&
            !anyDuplicated(pairs$value) &&
            (!is.numeric(pairs$value) ||
               !is.unsorted(pairs$value[order(pairs$code)], strictly = TRUE))
        }
      } else {
        s <- suppressWarnings(as.numeric(supplied))
        same <- all(is.na(s) == is.na(fitted))
        seen <- !is.na(fitted)
        if (same && sum(seen) > 0L) {
          same <- if (standardized) {
            sum(seen) < 3L ||
              isTRUE(stats::cor(s[seen], fitted[seen]) > 1 - 1e-8)
          } else {
            isTRUE(all(abs(s[seen] - fitted[seen]) <=
                         1e-8 * pmax(1, abs(fitted[seen]))))
          }
        }
      }
      if (!same) {
        .dmar_refuse_argument(
          sprintf(paste0("'data' does not hold the rows the model was ",
                         "fitted to: the values of %s differ from those ",
                         "the model used. Supply the rows given to the ",
                         "model, in their order."), v), arg = "data")
      }
    }
  }
  list(n = nrow(data), unassigned = unassigned, unselected = unselected,
       group = if (has_group) group else NULL, data = data)
}


# ---------------------------------------------------------------------------
# (F) Reasons
# ---------------------------------------------------------------------------

# Append a reason to the cells [rows, cols] of a character matrix of
# reasons, separated by "; ", without repeating one already there.
.fs_add_reason <- function(reason, rows, cols, text) {
  if (!length(rows) || !length(cols)) return(reason)
  cur <- reason[rows, cols, drop = FALSE]
  has <- vapply(cur, function(r) text %in% strsplit(r, "; ", fixed = TRUE)[[1]],
                logical(1))
  new <- ifelse(nzchar(cur), paste(cur, text, sep = "; "), text)
  new[has] <- cur[has]
  reason[rows, cols] <- new
  reason
}
