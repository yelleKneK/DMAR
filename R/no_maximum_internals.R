# A likelihood with no maximum, and the refusal it leads to. DMAR fits its
# factor analysis models with lavaan: the CFA engine (.fa_cfa_engine(),
# behind fa_confirmatory(), cfa_k(), cfa_1(), cfa_2(), fa_bifactor(),
# fa_second_order(), and the single-factor fits of the reliability family),
# fa_exploratory(), measurement_invariance(), and the configural fit of
# measurement_alignment(). lavaan lets an error
# variance fall below zero, and for some models and data the likelihood
# then has no maximum anywhere: it keeps rising, by less and less, as one
# variable's error variance falls without bound and its loading grows to
# match. lavaan's optimizer stops somewhere on that ridge, where its steps
# stop gaining, and reports convergence; the point depends on the starting
# values and the tolerance, not on the data, and so does every estimate,
# index, and fit row read from it. On holzinger_swineford, a bifactor
# model with the group factors spatial, verbal, and speed and a fourth on
# t24_woody_mccall, t25_paper_form_board_r, and t26_flags stops with the
# error variance of t25 at -1163 under FIML: the residual covariances of
# those three items are all negative, which no group factor can reproduce,
# and the likelihood rises toward the limit in which t25's loading is
# infinite. A finite negative estimate is a different case: the nine tests
# in three groups of three that fa_bifactor()'s help describes give visual
# perception an error variance of -41, at a maximum of the likelihood,
# which falls on both sides of it.
#
# The helpers below decide which case a fit is (.lavaan_no_maximum()) and
# stop with a classed refusal that says so and shows the lavaan code of the
# fit (.no_max_refuse()), or, inside a resampling loop, with a short
# condition the loop counts (.no_max_replication()). None is exported.
# The last helper of the file, .lavaan_check_converged(), comes before all
# of this: it refuses a fit given by the user whose optimizer did not
# converge at all.
#
# The structural equation models DMAR fits with lavaan are checked as
# well. The SEM planners read the population values of their parameters
# from a fit of the analysis model to the population moments
# (ss_aipe_composite_sem(), ss_power_composite_sem(), ss_aipe_sem_path(),
# and ss_aipe_sem_path_sensitivity()), and mediation_mbco() starts its
# OpenMx fits, which maximize the same likelihood, from a lavaan fit; each
# of these fits is refused when it has no maximum (.no_max_check()), with
# the user's model syntax printed whole. A replication of those planners'
# Monte Carlo loops with no maximum is set aside as one that did not
# converge (to a maximum, which lavaan's report of convergence on a ridge
# does not mean).
#
# A function that reads only the chi square of a fit is not checked:
# ss_aipe_rmsea_sensitivity(), whose RMSEA is formed from it, and the
# bootstrap of sem_exact_fit_bootstrap(). Along a ridge the likelihood
# levels off, rising by less and less toward a limit or flat, and lavaan's
# optimizer stops only where its steps no longer gain, so the chi square,
# unlike the estimates, does not depend on where on the ridge it stops: it
# moves by twice the gain the check's refit finds further out. Of 60
# bootstrap samples of a factor with four weak indicators (weak(2) of the
# no-maximum tests, resampled under seed 113), 7 ran off, and their
# statistics, 0.37 to 3.96 on 2 df, would move by 0.00005 to 0.0011;
# setting them aside moved the p-value from .125 (of 56) to .122 (of 49).
# Setting such a fit aside changes the reference distribution of the
# statistic, or refuses a population RMSEA that is determined, and so
# those fits are used as they are.


# .lavaan_no_maximum(fit, step = 10)
#
# Whether the lavaan fit `fit` sits on a ridge of its likelihood, with no
# maximum. The question arises only when some observed variable's error
# (residual) variance is estimated below zero, and only for a likelihood:
# the check applies to a maximum likelihood fit (lavaan's estimator "ML",
# which covers "MLR" and the other robust corrections, with complete data,
# FIML, or two-stage estimation, whose second stage maximizes the normal
# likelihood of the first stage's moments). It does not apply to least
# squares (ULS, GLS, WLS) or to the diagonally weighted fits behind WLSMV:
# their discrepancy is not a likelihood, lavaan computes no log likelihood
# for them, and a ridge of a least squares discrepancy is a question this
# check does not answer. A multilevel fit, and a rotated exploratory fit
# (whose parameter table holds the rotated solution, not the one
# estimated; the caller checks the unrotated fit), are not checked either.
# A structural model is checked as a factor analysis model is, with one
# limit. lavaan carries an observed variable that enters a regression as a
# latent variable of its own (.no_max_matrices()): an outcome or a
# covariate of the structural model, and an indicator that also enters a
# regression (a direct effect of a covariate on it, as in a MIMIC model
# with a direct effect, or its effect on another variable), whose loading
# then sits among the regressions. The variance of such a variable is not
# held, and its loading is not moved (.no_max_start()). For a covariate,
# and an outcome of observed variables alone, that variance is not an
# error variance; but an indicator so carried, or an outcome regressed on
# a latent variable (the same model as one more indicator), can carry a
# ridge as any indicator can, and the check does not look for one there:
# the fit is reported as it was before the check existed. On the ridges
# of that kind tried, lavaan reported no convergence, and a fit it reports
# as not converged is not checked: the three items of the tests whose
# correlations have a negative product, whose first runs off (N = 300),
# with a direct effect of a covariate on that item, with that item an
# outcome of the factor, or with that item a predictor of the covariate,
# under either scaling (the one factor model alone reports convergence on
# its ridge with the factor's variance fixed at one, and no convergence
# with a marker).
#
# The negative error variances are fixed further out, each at ten times
# its estimate (`step`) or one implied variance below it, whichever is
# further, and the model is fit again with every other parameter free,
# from a start that keeps every variable's implied variance
# (.no_max_start()): first all of them together, and then, when that does
# not show a ridge (or there is one), each alone, the most negative
# relative to its variable's model implied variance first. Together comes
# first because two variances can run off at once (listwise deletion on
# the bifactor example above sends t2_cubes and t25 off together), and
# then a refit that holds one of them leaves the other running and does
# not converge. At a finite maximum the log likelihood falls (by 0.012 for
# visual perception in the nine-test bifactor model above); on a ridge it
# does not (it rises by 0.0008 for t25 at ten times -1163). The fit is
# taken to have no maximum when the log likelihood of a fit further out is
# no more than 1e-6 per case below the fit's (`tol`, 1e-6 times the number
# of cases, the tolerance .fa_second_order_maximum() compares two fits
# with: the log likelihood is a sum over the cases, and so are the gaps
# lavaan's optimizer accepts). The floor of one implied variance matters
# for an estimate just below zero, which ten times itself barely moves: a
# finite maximum there would not lose 1e-6 per case over so short a step.
# What the comparison establishes is that the likelihood is flat to
# within the tolerance, or still rising, from the fit to the point further
# out, and the refusal says so in those words. A finite maximum that
# shallow is taken as none, and it can occur: of 23 resamples of the
# bifactor models of fa_bifactor()'s help that the check flagged, 22 rose
# at every multiple of the estimate tried (two, ten, a hundred, and a
# thousand times, where the refit converged), by 1.8e-5 at the least, and
# one (t4_lozenges at -1683 in the 15-test model, seed 113, the third
# replication) rose by 4e-6 to about -2300 and fell by 3.5e-5 from there
# toward its limit, against a tolerance of 3e-4. Its estimate is as much a
# matter of where the optimizer stopped as a ridge's is.
#
# Every row that shares a variance's parameter is fixed with it (equal
# error variances within a factor, or across the groups of a strict
# invariance model), so its equality constraints stay satisfied, and the
# constraints that tied those rows are dropped (.no_max_tied_rows()). The
# check is inconclusive, and the function returns NULL rather than declare
# a ridge, when a variance enters any other constraint (an inequality, a
# bound, or an equation with another parameter), where fixing it would
# change the model rather than move along it (that variance is not taken),
# or when every refit stops with an error or does not converge: a check
# that could not be made is no evidence of a ridge, and the fit is then
# reported as it was before the check existed. The refit reuses the fit's
# options, data, and sample statistics (lavaan's slotOptions, slotData,
# and slotSampleStats), with no standard errors, test, baseline, or
# unrestricted model, a single optimizer attempt, and no check of the
# starting values. lavaan reads the start from the `est` column of the
# parameter table. Along a ridge that start implies a positive definite
# covariance matrix; at a finite maximum the move toward a ridge can
# disturb the other covariances so that it does not, and lavaan was seen
# to set the start aside for values of its own then (on a two-group fit of
# a saturated triad), the variances still held. The comparison holds
# either way, since it reads only where the refit ends. One refit is made
# when a single variance is negative, and with several, one with all of
# them held and then one for each, until one shows a ridge: 0.03 s on the
# nine-test bifactor example above, and 0.08 s on the ridge, whose own fit
# takes 5 s; a proper solution costs nothing beyond reading its parameter
# table. Inside a bootstrap a replication with no maximum costs more than
# that: its stop is counted among the errors, whose messages the engine
# recovers by fitting each such replication once more, on the same
# processes as the replications (.boot_refit_errors()), and a fit that
# runs off is among the slowest lavaan makes (100 replications of the
# 15-test bifactor model of fa_bifactor()'s help took 143 s, against
# 99 s without the check).
#
# Returns NULL when the fit has a maximum, the check does not apply, or it
# is inconclusive; otherwise a list with `variable` (the observed variable
# or variables whose error variance was held), `group` (the label of each
# one's group, NULL with one group), `estimate` (their error variances at
# the fit), `fixed_at` (their values in the refit), `loglik` (the log
# likelihoods of the fit and of the refit, named "fit" and "further"),
# `tol`, and `N`.
.lavaan_no_maximum <- function(fit, step = 10) {
  if (!inherits(fit, "lavaan")) return(NULL)
  # Every fit of the families above passes through here, a bootstrap's
  # replications included, and nearly all have no negative error
  # variance, so that is asked first, of the parameter table lavaan keeps
  # in the fit (a list, read without building the data frame parTable()
  # returns): a free variance of a variable that is not a factor (the
  # factors are the left sides of "=~").
  raw <- fit@ParTable
  lv <- unique(raw$lhs[raw$op == "=~"])
  if (!any(raw$op == "~~" & raw$lhs == raw$rhs & !raw$lhs %in% lv &
           raw$free > 0L & is.finite(raw$est) & raw$est < 0)) {
    return(NULL)
  }
  opts <- lavaan::lavInspect(fit, "options")
  if (!identical(opts$estimator, "ML")) return(NULL)
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) return(NULL)
  if (isTRUE(lavaan::lavInspect(fit, "nlevels") > 1L)) return(NULL)
  pt <- lavaan::parTable(fit)
  if (!is.null(pt$efa) && any(nzchar(pt$efa)) &&
      !identical(opts$rotation, "none")) {
    return(NULL)
  }
  ov <- lavaan::lavNames(fit, "ov")
  neg <- which(pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% ov &
                 pt$free > 0L & is.finite(pt$est) & pt$est < 0)
  if (!length(neg)) return(NULL)
  ll0 <- .no_max_loglik(fit)
  if (!is.finite(ll0)) return(NULL)
  N <- lavaan::lavInspect(fit, "ntotal")
  tol <- 1e-6 * N
  mats <- .no_max_matrices(fit)
  # The variance of an observed variable that lavaan carries as a latent
  # variable of its own (an outcome or a covariate of a structural model,
  # or an indicator that also enters a regression; .no_max_matrices()) is
  # not a candidate: the start of the refit does not move such a variable
  # (.no_max_start()), the limit described above.
  neg <- neg[!mapply(function(v, b) v %in% mats[[b]]$dummy, pt$lhs[neg],
                     pt$block[neg])]
  if (!length(neg)) return(NULL)
  implied_var <- vapply(neg, function(r) {
    S <- mats[[pt$block[r]]]$implied
    if (pt$lhs[r] %in% rownames(S)) S[pt$lhs[r], pt$lhs[r]] else NA_real_
  }, numeric(1L))
  scale <- ifelse(is.finite(implied_var) & implied_var > 0, implied_var, 1)
  # The candidates, most negative first: each variance with the rows tied
  # to it and the value it is held at. One already tied to an earlier
  # candidate is held with it.
  cand <- list()
  for (k in order(pt$est[neg] / scale)) {
    r <- neg[k]
    if (any(vapply(cand, function(cc) r %in% cc$rows, logical(1L)))) next
    tied <- .no_max_tied_rows(pt, r)
    if (is.null(tied)) next
    to <- step * pt$est[r]
    if (is.finite(implied_var[k]) && implied_var[k] > 0) {
      to <- min(to, pt$est[r] - implied_var[k])
    }
    cand <- c(cand, list(list(r = r, rows = tied, to = to)))
  }
  if (!length(cand)) return(NULL)
  sets <- if (length(cand) > 1L) c(list(cand), lapply(cand, list)) else
    list(cand)
  groups <- lavaan::lavInspect(fit, "group.label")
  for (cs in sets) {
    rows <- unlist(lapply(cs, `[[`, "rows"))
    to <- unlist(lapply(cs, function(cc) rep(cc$to, length(cc$rows))))
    refit <- .no_max_refit(fit, pt, rows, to, mats, opts)
    if (is.null(refit)) next
    ll1 <- .no_max_loglik(refit)
    if (!is.finite(ll1) || ll1 < ll0 - tol) next
    r <- vapply(cs, `[[`, integer(1L), "r")
    return(list(variable = pt$lhs[r],
                group = if (length(groups) > 1L) unname(groups[pt$group[r]]),
                estimate = pt$est[r],
                fixed_at = vapply(cs, `[[`, numeric(1L), "to"),
                loglik = c(fit = ll0, further = ll1), tol = tol, N = N))
  }
  NULL
}

# .no_max_loglik(fit)
#
# The log likelihood of a lavaan fit, NA when lavaan gives none. The one
# warning lavaan's logLik() gives, "logLik only available if estimator is
# ML" for a fit without a likelihood, comes with the NA that already says
# so, and is muffled at this one call (the check reads only fits under
# lavaan's estimator "ML", but a caller's fit is read as it comes).
.no_max_loglik <- function(fit) {
  tryCatch(suppressWarnings(as.numeric(lavaan::logLik(fit))),
           error = function(e) NA_real_)
}

# .no_max_matrices(fit)
#
# For each group of a lavaan fit, the matrices the start of the refit is
# built from: the loadings (`lambda`, observed variables by latent
# variables), the regressions among the latent variables (`beta`, NULL
# when there are none), the model implied covariance matrix of the latent
# variables (`cov_lv`), that of the observed variables (`implied`), and
# the observed variables lavaan carries as latent variables of their own
# (`dummy`). In a structural model lavaan represents an observed variable
# that enters a regression (an outcome, a covariate, or an indicator that
# is also regressed on a covariate or predicts another variable) by a
# latent variable of the same name with a loading of one and no error
# variance; such an indicator's loading on its factor then sits in the
# regressions (`beta`), and its row of the loadings holds a zero there.
# The loadings and the regressions carry these, and lavaan's covariance
# matrix of the latent variables (lavInspect(fit, "cov.lv")) leaves them
# out, so with them that matrix is formed from the regressions and the
# disturbances, (I - B)^-1 Psi (I - B)^-T. The variance of such a
# variable, which lavaan holds among the disturbances, is neither held nor
# moved by the check, and an indicator's loading held among the
# regressions is not moved (.lavaan_no_maximum(), .no_max_start()).
.no_max_matrices <- function(fit) {
  one <- lavaan::lavInspect(fit, "ngroups") == 1L
  wrap <- function(x) if (one) list(x) else x
  est <- wrap(lavaan::lavInspect(fit, "est"))
  lv <- wrap(lavaan::lavInspect(fit, "cov.lv"))
  implied <- wrap(lavaan::lavInspect(fit, "implied"))
  lapply(seq_along(est), function(b) {
    lambda <- est[[b]]$lambda
    beta <- est[[b]]$beta
    cov_lv <- lv[[b]]
    dummy <- intersect(colnames(lambda), rownames(lambda))
    if (!is.null(lambda) && !identical(colnames(lambda), colnames(cov_lv))) {
      psi <- unclass(est[[b]]$psi)
      IB <- if (is.null(beta)) {
        diag(nrow(psi))
      } else {
        solve(diag(nrow(psi)) - unclass(beta))
      }
      cov_lv <- IB %*% psi %*% t(IB)
      dimnames(cov_lv) <- dimnames(psi)
    }
    list(lambda = lambda,
         beta = if (!is.null(beta) && any(beta != 0)) beta,
         cov_lv = cov_lv, implied = implied[[b]]$cov, dummy = dummy)
  })
}

# .no_max_tied_rows(pt, r)
#
# The rows of the parameter table `pt` that share the parameter of row
# `r`: rows with its label, rows with its free parameter number (lavaan's
# ceq.simple = TRUE gives equated rows one number), and rows tied to any
# of these by a simple equality constraint ("a == b", both sides naming
# parameters), taken to closure. Returns NULL when one of them enters any
# other constraint (an inequality, an equation with an expression or a
# constant) or carries a finite bound: fixing the variance would then
# change the model rather than move along it, and the check is
# inconclusive.
.no_max_tied_rows <- function(pt, r) {
  con <- pt$op %in% c("==", "<", ">")
  param <- !con & pt$op != ":="
  label <- if (is.null(pt$label)) rep("", length(pt$lhs)) else pt$label
  plabel <- if (is.null(pt$plabel)) rep("", length(pt$lhs)) else pt$plabel
  keys_of <- function(rows) {
    k <- c(plabel[rows], label[rows])
    unique(k[nzchar(k)])
  }
  tied <- r
  repeat {
    keys <- keys_of(tied)
    free <- pt$free[tied]
    free <- free[free > 0L]
    add <- which(param & ((nzchar(label) & label %in% keys) |
                            (nzchar(plabel) & plabel %in% keys) |
                            (pt$free > 0L & pt$free %in% free)))
    simple <- which(pt$op == "==" & (pt$lhs %in% keys | pt$rhs %in% keys))
    for (s in simple) {
      sides <- c(pt$lhs[s], pt$rhs[s])
      add <- c(add, which(param & ((nzchar(label) & label %in% sides) |
                                     (nzchar(plabel) & plabel %in% sides))))
    }
    new <- setdiff(add, tied)
    if (!length(new)) break
    tied <- c(tied, new)
  }
  keys <- keys_of(tied)
  for (s in which(con)) {
    words <- unlist(regmatches(c(pt$lhs[s], pt$rhs[s]),
                               gregexpr("[A-Za-z0-9_.]+", c(pt$lhs[s],
                                                            pt$rhs[s]))))
    if (!any(words %in% keys)) next
    if (!(pt$op[s] == "==" && pt$lhs[s] %in% keys && pt$rhs[s] %in% keys)) {
      return(NULL)
    }
  }
  if (!is.null(pt$lower) && any(is.finite(pt$lower[tied]))) return(NULL)
  if (!is.null(pt$upper) && any(is.finite(pt$upper[tied]))) return(NULL)
  sort(tied)
}

# .no_max_refit(fit, pt, tied, to, mats, opts)
#
# The fit of the model of `fit` with the rows `tied` of its parameter table
# `pt` fixed at `to` (one value per row), from the start .no_max_start()
# gives, reusing the fit's options (`opts`), data, and sample statistics;
# NULL when lavaan stops with an error or the refit does not converge.
# lavaan's warnings about the refit (a covariance matrix of the estimates
# it is not asked for, a variance that is negative by construction, a
# solution that has run away) concern a fit no one sees, and are muffled
# at this one call.
.no_max_refit <- function(fit, pt, tied, to, mats, opts) {
  start <- .no_max_start(pt, tied, to, mats)
  keys <- c(pt$plabel[tied], pt$label[tied])
  keys <- keys[nzchar(keys)]
  m <- pt
  m$est <- start
  m$start <- start
  m$free[tied] <- 0L
  m$ustart[tied] <- to
  m$se <- NULL
  # The equalities among the fixed rows hold by construction. The defined
  # parameters (the indices DMAR's syntax defines) do not enter the
  # likelihood, and a definition that reads a fixed row is not needed, so
  # they are left out unless a constraint reads one of them.
  drop <- which(m$op == "==" & m$lhs %in% keys & m$rhs %in% keys)
  defined <- which(m$op == ":=")
  con <- m$op %in% c("==", "<", ">")
  con_words <- unlist(regmatches(c(m$lhs[con], m$rhs[con]),
                                 gregexpr("[A-Za-z0-9_.]+",
                                          c(m$lhs[con], m$rhs[con]))))
  if (!any(m$lhs[defined] %in% con_words)) drop <- c(drop, defined)
  if (length(drop)) m <- m[-drop, , drop = FALSE]
  free <- m$free > 0L
  m$free[free] <- match(m$free[free], sort(unique(m$free[free])))
  m$id <- seq_len(nrow(m))
  o <- opts
  o$se <- "none"
  o$test <- "none"
  o$baseline <- FALSE
  o$h1 <- FALSE
  o$check.start <- FALSE
  o$check.post <- FALSE
  o$optim.attempts <- 1L
  if (identical(o$optim.method, "none")) o$optim.method <- "nlminb"
  o$start <- "default"
  o$verbose <- FALSE
  refit <- tryCatch(suppressWarnings(lavaan::lavaan(
    model = m, slotOptions = o, slotSampleStats = fit@SampleStats,
    slotData = fit@Data)), error = function(e) NULL)
  if (is.null(refit) || !isTRUE(lavaan::lavInspect(refit, "converged"))) {
    return(NULL)
  }
  refit
}

# .no_max_start(pt, tied, to, mats)
#
# The starting values of the refit, one per row of `pt`, with the rows
# `tied` held at `to` (one value per row): the fit's estimates, moved
# toward the ridge the check follows. Fixing an error
# variance far below zero alone would leave its variable a negative
# implied variance, where the likelihood is not defined, so each fixed
# variable keeps its implied variance. In the metric of the factors'
# standard deviations, its largest loading grows by the factor s that
# makes up the difference, the other loadings on that factor shrink by
# 1/s, and, in a model without regressions among the factors, so do the
# factor's correlations with the other factors, which keeps most of the
# implied covariances where they were. The new loadings are then put back
# on the factor's own scale: its standard deviation as it was, or, when a
# loading fixes the scale (a marker variable), the one that keeps that
# loading at its fixed value, the factor's free variance taking the
# difference. Each moved row passes its value to the rows it is equated
# with (.no_max_tied_rows()), and every other free error variance then
# takes the value that keeps its own variable's implied variance. On a
# ridge this is the direction the fit runs off in (the loading of the
# variable grows without bound and those of the factor's other items
# shrink); at a finite maximum it is a start from which nlminb finds the
# best fit with the variance held, and that is all the check needs.
#
# An indicator that lavaan carries as a latent variable of its own
# (.no_max_matrices()) has its loading among the regressions and a zero in
# its row of the loadings, so its loading row is left at its estimate, not
# moved and not read as a marker, and its variance, held among the
# disturbances, is left as it is. With regressions in the model the
# factor's variance is not changed either, and so such an indicator keeps
# its implied variance as every other variable does; a factor whose scale
# is fixed by such an indicator's loading keeps that scale.
.no_max_start <- function(pt, tied, to, mats) {
  st <- pt$est
  st[tied] <- to
  for (r in tied) {
    b <- pt$block[r]
    m <- mats[[b]]
    item <- pt$lhs[r]
    L <- m$lambda
    V <- m$cov_lv
    if (is.null(L) || !item %in% rownames(L) || !ncol(L)) next
    sdv <- sqrt(pmax(diag(V), 0))
    lam <- L[item, ] * sdv
    k <- which.max(abs(lam))
    if (!length(k) || lam[k] == 0 || sdv[k] <= 0) next
    lv <- colnames(L)[k]
    gap <- pt$est[r] - to[match(r, tied)]
    if (is.null(m$beta)) {
      # The cross terms of the variable's variance through factor k keep
      # their size when the loading grows by s and the correlations shrink
      # by 1/s, so s has a closed form.
      s <- sqrt(1 + gap / lam[k]^2)
    } else {
      R <- V / outer(sdv, sdv)
      bq <- sum(R[k, -k] * lam[-k])
      cq <- if (length(lam) > 1L) {
        drop(lam[-k] %*% R[-k, -k, drop = FALSE] %*% lam[-k])
      } else {
        0
      }
      disc <- bq^2 - (cq - drop(lam %*% R %*% lam) - gap)
      if (!is.finite(disc) || disc < 0) next
      side <- if (lam[k] + bq != 0) sign(lam[k] + bq) else sign(lam[k])
      s <- (-bq + side * sqrt(disc)) / lam[k]
    }
    if (!is.finite(s) || s <= 0) next
    in_b <- pt$block == b
    load_rows <- which(in_b & pt$op == "=~" & pt$lhs == lv &
                         pt$rhs %in% rownames(L) & !pt$rhs %in% m$dummy)
    col <- L[, k] * sdv[k]
    col_new <- col / s
    col_new[item] <- col[item] * s
    marker <- load_rows[pt$free[load_rows] == 0L & pt$est[load_rows] != 0]
    d <- sdv[k]
    if (length(marker)) {
      if (!is.null(m$beta)) next
      d <- col_new[pt$rhs[marker[1L]]] / pt$est[marker[1L]]
    }
    if (!is.finite(d) || d <= 0) next
    free_load <- load_rows[pt$free[load_rows] > 0L]
    st[free_load] <- col_new[pt$rhs[free_load]] / d
    moved <- free_load
    if (is.null(m$beta)) {
      var_row <- which(in_b & pt$op == "~~" & pt$lhs == lv & pt$rhs == lv &
                         pt$free > 0L)
      st[var_row] <- d^2
      cov_rows <- which(in_b & pt$op == "~~" & pt$lhs != pt$rhs &
                          (pt$lhs == lv | pt$rhs == lv) &
                          pt$lhs %in% colnames(L) & pt$rhs %in% colnames(L) &
                          pt$free > 0L)
      other <- ifelse(pt$lhs[cov_rows] == lv, pt$rhs[cov_rows],
                      pt$lhs[cov_rows])
      st[cov_rows] <- V[lv, other] / (sdv[k] * sdv[other]) / s * d *
        sdv[other]
      moved <- c(moved, var_row, cov_rows)
    }
    for (q in moved) {
      eq <- .no_max_tied_rows(pt, q)
      if (!is.null(eq)) st[eq] <- st[q]
    }
  }
  # The free error variances, each keeping its variable's implied
  # variance under the moved loadings (and, without regressions among the
  # factors, the moved factor variances and covariances).
  for (b in seq_along(mats)) {
    m <- mats[[b]]
    L <- m$lambda
    if (is.null(L) || !ncol(L)) next
    in_b <- pt$block == b
    ld <- which(in_b & pt$op == "=~" & pt$rhs %in% rownames(L) &
                  pt$lhs %in% colnames(L))
    L_new <- L
    L_new[cbind(match(pt$rhs[ld], rownames(L)),
                match(pt$lhs[ld], colnames(L)))] <- st[ld]
    V_new <- m$cov_lv
    if (is.null(m$beta)) {
      cv <- which(in_b & pt$op == "~~" & pt$lhs %in% colnames(L) &
                    pt$rhs %in% colnames(L))
      V_new[cbind(match(pt$lhs[cv], colnames(L)),
                  match(pt$rhs[cv], colnames(L)))] <- st[cv]
      V_new[cbind(match(pt$rhs[cv], colnames(L)),
                  match(pt$lhs[cv], colnames(L)))] <- st[cv]
    }
    err <- which(in_b & pt$op == "~~" & pt$lhs == pt$rhs &
                   pt$lhs %in% rownames(L) & !pt$lhs %in% m$dummy &
                   pt$free > 0L)
    err <- setdiff(err, tied)
    if (!length(err)) next
    common <- rowSums((L_new %*% V_new) * L_new)
    st[err] <- diag(m$implied)[pt$lhs[err]] - common[pt$lhs[err]]
  }
  st
}

# .no_max_call_text(fun, args, data = "my_data", select = list())
#
# The lavaan call that reproduces a fit, as text: `fun` (for example
# "lavaan::lavaan") called with the model syntax as `model` (written
# beside the call, so named and not repeated) and the arguments `args` as
# DMAR passed them, the data replaced by a placeholder name (`data`, or
# my_S and my_M for a covariance matrix and its means, beside the number
# of cases as given, in lavaan's spelling or the one lavaan::efa() also
# takes), which comes first. The arguments are joined into lines of about
# 72 characters, broken only between arguments, the lines after the first
# indented.
#
# The placeholder stands for what the user passed, and DMAR may have
# passed lavaan only part of it. lavaan::lavaan() and lavaan::cfa() read
# the model's variables from a data frame, and from a covariance matrix by
# its names, whatever else either holds; but lavaan::efa() takes every
# column of a data frame and every row of a covariance matrix (it reads
# the variables from its row names, setting its own ov_names aside), and
# lavaan pairs means with the rows of the covariance matrix by position,
# not by name. `select` names, for such an argument ("data", "sample_cov",
# "sample.cov", "sample.mean", ...), the variables DMAR passed, and its
# placeholder then selects them: my_data[c(...)], my_S[c(...), c(...)],
# or my_M[c(...)]. An argument not named in `select` is printed as the
# bare placeholder.
.no_max_call_text <- function(fun, args, data = "my_data",
                              select = list()) {
  first <- if (!is.null(args$model)) "model"
  args$model <- NULL
  input <- c("data", "sample.cov", "sample_cov", "sample.mean",
             "sample_mean", "sample.nobs", "sample_nobs")
  holder <- function(arg, name, matrix = FALSE) {
    x <- as.name(name)
    sel <- select[[arg]]
    if (is.null(sel)) return(x)
    if (matrix) call("[", x, sel, sel) else call("[", x, sel)
  }
  if (!is.null(args$data)) args$data <- holder("data", data)
  for (nm in intersect(c("sample.cov", "sample_cov"), names(args))) {
    args[[nm]] <- holder(nm, "my_S", matrix = TRUE)
  }
  for (nm in intersect(c("sample.mean", "sample_mean"), names(args))) {
    args[[nm]] <- holder(nm, "my_M")
  }
  args <- args[!vapply(args, is.null, logical(1L))]
  args <- args[order(!names(args) %in% input)]
  # deparse() breaks a long value into lines, indented after the first;
  # they are joined again with one space.
  pieces <- c(first, vapply(names(args), function(nm) {
    paste0(nm, " = ", paste(trimws(deparse(args[[nm]], width.cutoff = 500L)),
                            collapse = " "))
  }, character(1L), USE.NAMES = FALSE))
  pieces <- paste0(pieces, c(rep(",", length(pieces) - 1L), ")"))
  lines <- paste0(fun, "(", pieces[1L])
  for (piece in pieces[-1L]) {
    last <- lines[length(lines)]
    if (nchar(last) + 1L + nchar(piece) <= 72L) {
      lines[length(lines)] <- paste(last, piece)
    } else {
      lines <- c(lines, paste0("  ", piece))
    }
  }
  paste(lines, collapse = "\n")
}

# .no_max_refuse(fn, nm, syntax = nm$syntax, lavaan_call = nm$lavaan_call,
#                where = nm$where, ..., code_intro = NULL,
#                placeholder = NULL, engine = NULL, defined = "drop",
#                extra_class = NULL)
#
# Stops with the refusal of a fit whose likelihood has no maximum
# (.lavaan_no_maximum()'s list `nm`): an error of class c("dmar_no_maximum",
# "dmar_refusal", "error", "condition") with code "no_maximum", whose
# message names the function (`fn`, without parentheses), says that the fit
# failed because the likelihood has no maximum, names the variable whose
# error variance ran off below zero with its value and what the check
# found there, says that DMAR fits the model with lavaan, the estimation
# engine, and to see lavaan for more details, with the options there
# (bounds = "pos.var", or a different model), and prints the lavaan code
# that reproduces the fit (the model syntax DMAR built, `syntax`, when
# there is one, and the call, `lavaan_call`, from .no_max_call_text()).
# `where` names the part of the analysis the fit belongs to (a rung of an
# invariance ladder, a planner's fit to the population moments), NULL for
# none; `code_intro`, the words that introduce the code when they are not
# "This lavaan code reproduces the fit"; `placeholder`, the words that
# say what the placeholders of the code stand for when they are not
# "your data in place of my_data" (or its covariance matrix counterpart);
# `engine`, the sentences that say how DMAR fits the model and what the
# choices are when DMAR does not fit it with lavaan alone; and `defined`,
# whether the syntax is printed without its ":=" lines ("drop", for the
# syntax DMAR builds, whose ":=" lines define its indices) or whole
# ("keep", for a model the user wrote, whose ":=" lines define the user's
# own quantities) (.no_max_message()). The condition carries `syntax` (in
# full), `lavaan_call`, `variables`, `group`, `estimate`, `fixed_at`,
# `loglik` (the log likelihoods of the fit and of the fit further out,
# named "fit" and "further"), `tolerance` (the check's, 1e-6 per case) and
# `N` (the cases), `fn`, and `where`, with any further fields in `...`;
# `extra_class` is passed to .dmar_refuse().
#
# R prints at most getOption("warning.length") bytes of an error's
# message, 1000 by default, and the lavaan code runs past that, so the
# limit is raised to its largest value, 8170, while the error is signaled
# and restored on the way out (after R has printed the message, when no
# caller catches it). A message longer still leaves the model syntax out
# (.no_max_message()), so that nothing R would cut is lost.
.no_max_refuse <- function(fn, nm, syntax = nm$syntax,
                           lavaan_call = nm$lavaan_call, where = nm$where,
                           ..., code_intro = NULL, placeholder = NULL,
                           engine = NULL, defined = "drop",
                           extra_class = NULL) {
  old <- options(warning.length = 8170L)
  on.exit(options(old), add = TRUE)
  .dmar_refuse(.no_max_message(fn, nm, syntax, lavaan_call, where,
                               code_intro = code_intro,
                               placeholder = placeholder, engine = engine,
                               defined = defined),
               "no_maximum", variables = nm$variable, group = nm$group,
               estimate = nm$estimate, fixed_at = nm$fixed_at,
               loglik = nm$loglik, tolerance = nm$tol, N = nm$N,
               syntax = syntax, lavaan_call = lavaan_call, fn = fn,
               where = where, ..., extra_class = extra_class)
}

# .no_max_message(fn, nm, syntax, lavaan_call, where = NULL,
#                 code_intro = NULL, placeholder = NULL, engine = NULL,
#                 defined = "drop", limit = 8000L)
#
# The message of .no_max_refuse(). It says what the check established:
# with the variance held further out, the log likelihood is higher (the
# likelihood is still rising as the variance falls) or lower by no more
# than the tolerance of 1e-6 per case (flat). The options in lavaan come
# before the code, and the code last: by default, that DMAR fits the model
# with lavaan and that bounds = "pos.var" or a different model are the
# choices there, which `engine` replaces where lavaan is not the engine
# (mediation_mbco() fits its model with OpenMx from a lavaan fit's
# estimates). `placeholder` replaces the words that say what my_data, or
# my_S and my_M, stand for, where they are not the user's data or
# covariance matrix and its means as given (a planner's population
# moments, which may come from cov_sem()). The model syntax DMAR builds is
# printed without the lines that define its indices (lavaan's ":=" lines,
# `defined` = "drop"), which do not change the fit and can run to
# thousands of characters, and the message says so; a model the user wrote
# (`defined` = "keep", the SEM planners and mediation_mbco()) is printed
# whole, since its ":=" lines define the user's own quantities, such as an
# indirect effect. The condition keeps the syntax whole either way. A
# message that would still pass `limit` bytes (R cuts a printed error at
# 8170, less the "Error: " before it) leaves the syntax out and says where
# it is, and keeps the call; the 12-test bifactor model's message runs to
# 3,400 characters, about 200 for each item with a mean structure. The
# numbers are formatted for the message only.
.no_max_message <- function(fn, nm, syntax, lavaan_call, where = NULL,
                            code_intro = NULL, placeholder = NULL,
                            engine = NULL, defined = "drop", limit = 8000L) {
  defined <- match.arg(defined, c("drop", "keep"))
  num <- function(v) {
    .no_max_and(trimws(formatC(v, format = "fg", digits = 4,
                               big.mark = "")))
  }
  ll <- function(v) formatC(v, format = "f", digits = 4)
  small <- function(v) format(signif(v, 2), scientific = FALSE)
  one <- length(nm$variable) == 1L
  it <- if (one) "it" else "them"
  falls <- if (one) "the variance falls" else "the variances fall"
  gain <- nm$loglik[["further"]] - nm$loglik[["fit"]]
  compare <- if (gain >= 0) {
    paste0("higher, ", ll(nm$loglik[["further"]]), " against ",
           ll(nm$loglik[["fit"]]), ": the likelihood is still rising as ",
           falls)
  } else {
    paste0("lower by only ", small(-gain),
           ", within the tolerance of 1e-6 per case",
           if (!is.null(nm$tol) && !is.null(nm$N)) {
             paste0(" (", small(nm$tol), " for these ",
                    format(nm$N, big.mark = ",", scientific = FALSE),
                    " cases)")
           },
           ": the likelihood is flat as ", falls)
  }
  head <- paste0(
    fn, "() stopped: the fit failed because the likelihood has no ",
    "maximum", if (!is.null(where)) paste0(" (", where, ")"), ". The ",
    "error ", if (one) "variance" else "variances", " of ",
    .no_max_variables(nm), if (one) " runs" else " run",
    " off below zero: lavaan's optimizer stopped with ", it, " at ",
    num(nm$estimate), ", and with ", it, " held at ", num(nm$fixed_at),
    " the log likelihood is ", compare, ". The estimates are where the ",
    "optimizer stopped, not at a maximum, so no result is returned. ",
    if (is.null(engine)) {
      paste0("DMAR fits the model with lavaan, the estimation engine; see ",
             "lavaan for more details. There, bounds = \"pos.var\" keeps ",
             "every error variance at zero or above, or a different model ",
             "can be fit.")
    } else {
      engine
    })
  if (is.null(syntax) && is.null(lavaan_call)) return(head)
  if (is.null(placeholder)) {
    placeholder <- if (!is.null(lavaan_call) &&
                       grepl("my_S", lavaan_call, fixed = TRUE)) {
      paste0("your covariance matrix in place of my_S",
             if (grepl("my_M", lavaan_call, fixed = TRUE)) {
               " and its means in place of my_M"
             })
    } else {
      "your data in place of my_data"
    }
  }
  intro <- if (is.null(code_intro)) {
    "This lavaan code reproduces the fit"
  } else {
    code_intro
  }
  call_line <- if (!is.null(lavaan_call)) paste0("fit <- ", lavaan_call)
  fields <- c(if (!is.null(syntax)) "syntax", if (!is.null(lavaan_call)) {
    "lavaan_call"
  })
  carried <- paste0("; the error carries it as its ",
                    paste(fields, collapse = " and "),
                    if (length(fields) > 1L) " elements" else " element",
                    " (see ?dmar_metadata)")
  drop <- defined == "drop"
  model <- if (!is.null(syntax)) {
    lines <- strsplit(syntax, "\n")[[1L]]
    kept <- lines[nzchar(trimws(lines)) &
                    !(drop & grepl(":=", lines, fixed = TRUE))]
    c("model <- '", kept, "'")
  }
  left_out <- drop && !is.null(syntax) && grepl(":=", syntax, fixed = TRUE)
  full <- paste0(
    head, " ", intro, ", with ", placeholder,
    if (left_out) {
      paste0(" (the lines of DMAR's syntax that define its indices are ",
             "left out; they do not change the fit)")
    }, carried, ":\n\n", paste(c(model, call_line), collapse = "\n"))
  if (nchar(full, type = "bytes") <= limit || is.null(model) ||
      is.null(call_line)) {
    return(full)
  }
  paste0(
    head, " The model syntax", if (drop) " DMAR built", ", ",
    length(model) - 2L, " lines, is too long to print here; the error ",
    "carries it as its syntax element, and the call below as its ",
    "lavaan_call element (see ?dmar_metadata). ", intro, ", with the ",
    "syntax as model and ", placeholder, ":\n\n", call_line)
}

# .no_max_variables(nm)
#
# The variables of .lavaan_no_maximum()'s list, each with its group when
# there are groups, joined into a phrase ("a", "a and b", "a, b, and c").
.no_max_variables <- function(nm) {
  v <- nm$variable
  if (!is.null(nm$group)) v <- paste0(v, " in group ", nm$group)
  .no_max_and(v)
}

# .no_max_and(x)
#
# The strings `x` joined into a phrase: "a", "a and b", "a, b, and c".
.no_max_and <- function(x) {
  n <- length(x)
  if (n <= 1L) return(paste(x, collapse = ""))
  if (n == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-n], collapse = ", "), ", and ", x[n])
}

# .no_max_short(nm)
#
# The short statement of .lavaan_no_maximum()'s list `nm`, with no values:
# "the likelihood has no maximum: the error variance of x runs off below
# zero". The message of a replication's stop (.no_max_replication()), and
# the reason a caller records where the whole message, with its lavaan
# code, would not fit (factor_count()'s reasons, the reliability family's
# fit_error).
.no_max_short <- function(nm) {
  one <- length(nm$variable) == 1L
  paste0("the likelihood has no maximum: the error ",
         if (one) "variance" else "variances", " of ", .no_max_variables(nm),
         if (one) " runs" else " run", " off below zero")
}

# .no_max_replication(nm)
#
# Inside a resampling loop (a bootstrap of the reliability family, of
# ci_bifactor_indices(), or of reliability_H() on a fit), the stop of a
# replication whose likelihood has no maximum: the same class as
# .no_max_refuse()'s, with the fields of `nm`, and a short message that
# names the variable but no value, so that the bootstrap engine, which
# records each replication's message and counts the replication among
# those that raised an error (.boot_refit_errors()), gathers the
# replications that ran off on the same variable under one message rather
# than one message apiece, each with the lavaan code of a resample.
.no_max_replication <- function(nm) {
  .dmar_refuse(.no_max_short(nm),
               "no_maximum", variables = nm$variable, group = nm$group,
               estimate = nm$estimate, fixed_at = nm$fixed_at,
               loglik = nm$loglik, tolerance = nm$tol, N = nm$N,
               syntax = NULL, lavaan_call = NULL)
}

# .no_max_fields(cond)
#
# The list .lavaan_no_maximum() returns, with `syntax`, `lavaan_call`, and
# `where`, read back from a dmar_no_maximum condition, so that a caller
# that catches the refusal of a fit it made for its own purposes can raise
# it again under its own name (.no_max_rethrow()) or keep it for later
# (the reliability family's single-factor fits, .omega_fit_cfa()).
.no_max_fields <- function(cond) {
  list(variable = cond$variables, group = cond$group,
       estimate = cond$estimate, fixed_at = cond$fixed_at,
       loglik = cond$loglik, tol = cond$tolerance, N = cond$N,
       syntax = cond$syntax, lavaan_call = cond$lavaan_call,
       where = cond$where)
}

# .no_max_renamed(nm, lavaan_names, names)
#
# The list `nm` (.no_max_fields()) with the variables given the caller's
# names, in the variables, the syntax, and the call: the reliability
# family enters items into lavaan as y1, ..., yJ (and auxiliary variables
# as a1, ..., aK), and the user knows them by their own names
# (.relia_lavaan_names_back()).
.no_max_renamed <- function(nm, lavaan_names, names) {
  back <- function(x) {
    if (is.null(x)) x else .relia_lavaan_names_back(x, lavaan_names, names)
  }
  nm$variable <- back(nm$variable)
  nm$syntax <- back(nm$syntax)
  nm$lavaan_call <- back(nm$lavaan_call)
  nm
}

# .no_max_rethrow(cond, fn)
#
# The refusal `cond` raised again under the name of the function the user
# called, `fn`: cfa_1() and cfa_2() reach the fit through
# fa_confirmatory(), and the message should name the one in the user's
# call.
.no_max_rethrow <- function(cond, fn) {
  .no_max_refuse(fn, .no_max_fields(cond))
}

# .no_max_as(expr, fn)
#
# Evaluates `expr` and raises a dmar_no_maximum refusal from it again
# under `fn` (.no_max_rethrow()); anything else passes as it is.
.no_max_as <- function(expr, fn) {
  tryCatch(expr, dmar_no_maximum = function(e) .no_max_rethrow(e, fn))
}

# .no_max_check(fit, fn, syntax = NULL, fun = "lavaan::sem", args = NULL,
#               select = list(), where = NULL, placeholder = NULL,
#               engine = NULL, defined = "keep")
#
# The check and the refusal together, for a function that fits one model
# and reads its result from that fit: when .lavaan_no_maximum() finds that
# `fit` has no maximum, stops with .no_max_refuse() under `fn`, with the
# model syntax `syntax` and the call to `fun` with the arguments `args` as
# DMAR passed them (.no_max_call_text(), with `select`), or with no code
# when `args` is NULL (a fit the user made, which has no DMAR syntax);
# `where`, `placeholder`, `engine`, and `defined` go to the message, the
# syntax printed whole by default, since it is the model the user wrote.
# Returns NULL, invisibly, when the fit has a maximum or the check does
# not apply. The SEM planners check their fits to the population moments
# with it; mediation_mbco() calls .no_max_refuse() itself, with the code
# it builds for its own fit (.mbco_no_max_code()).
.no_max_check <- function(fit, fn, syntax = NULL, fun = "lavaan::sem",
                          args = NULL, select = list(), where = NULL,
                          placeholder = NULL, engine = NULL,
                          defined = "keep") {
  nm <- .lavaan_no_maximum(fit)
  if (is.null(nm)) return(invisible(NULL))
  .no_max_refuse(fn, nm, syntax = syntax,
                 lavaan_call = if (!is.null(args)) {
                   .no_max_call_text(fun, args, select = select)
                 },
                 where = where, placeholder = placeholder, engine = engine,
                 defined = defined)
}


# .lavaan_check_converged(fit, fn, label = "'fit'")
#
# The check a function given a fitted lavaan model makes before it reads
# an estimate from it: lavaan's optimizer met its convergence criterion.
# A fit that did not converge holds the point where the optimizer stopped
# (at its iteration limit, for instance), not a solution, and every
# loading, variance, and index read from it would describe that point.
# Through DMAR 1.0.1 and the first builds of 1.1.0,
# average_variance_extracted() reported an AVE of .412 for a one-factor
# fit stopped after one iteration (with control = list(iter.max = 1)),
# bifactor_indices() a full table for a nine-item bifactor fit stopped
# the same way, and dmacs() the dMACS values of a stopped multiple-group
# fit. Such a fit is refused with class dmar_not_converged, the class the
# same refusal has elsewhere in the package (fa_confirmatory(),
# factor_scores(), sem_fit_measures()); the condition carries `fn` and
# the number of iterations lavaan ran as `iterations`. A fit that
# converged to an improper solution passes: each function treats a
# negative variance or a standardized loading beyond one under its own
# documented rule, and a fit whose likelihood has no maximum is the
# business of .lavaan_no_maximum(). `fn` names the function in the
# message, `label` the object. Returns `fit`, invisibly.
.lavaan_check_converged <- function(fit, fn, label = "'fit'") {
  if (isTRUE(lavaan::lavInspect(fit, "converged"))) return(invisible(fit))
  iterations <- tryCatch(as.integer(lavaan::lavInspect(fit, "iterations")),
                         error = function(e) NA_integer_)
  if (length(iterations) != 1L) iterations <- NA_integer_
  ran <- if (is.na(iterations)) {
    ""
  } else if (iterations == 0L) {
    " (lavaan's optimizer was not run, as with do.fit = FALSE)"
  } else {
    paste0(" (lavaan stopped its optimizer after ", iterations,
           if (iterations == 1L) " iteration" else " iterations", ")")
  }
  .dmar_refuse(c(
    "The model in ", label, " did not converge", ran, ", so its ",
    "estimates mark where the optimizer stopped and are not a solution; ",
    fn, " reports nothing from such a fit. Resolve the convergence ",
    "problem (better starting values, a simpler model, or more ",
    "iterations) and refit."),
    "not_converged", fn = fn, iterations = iterations)
}
