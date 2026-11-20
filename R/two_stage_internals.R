# Shared helpers for two-stage maximum likelihood with incomplete data
# (Yuan & Bentler, 2000; Savalei & Bentler, 2009), the missing data option
# "two_stage" (and "robust_two_stage") of fa_confirmatory(), measurement_invariance(),
# reliability_omega(), reliability_alpha(), and reliability(), and the
# two-stage branches of sem_fit_measures(), sem_information_criteria(), and
# compare_sem(). None of them is exported.
#
# Two-stage estimation first estimates the means and covariances of the
# variables by maximum likelihood under the saturated normal model (the EM
# estimates, stage one, with any auxiliary variables), and then fits the
# structured model to those moments with the complete data maximum
# likelihood discrepancy function, as if they were the moments of complete
# data (stage two). lavaan does both (missing = "two.stage" or
# "robust.two.stage") and supplies the two-stage standard errors and the
# scaled statistic. What it does not supply is the residual based statistic
# of Savalei and Bentler (2009, Eqs. 8 and 14), the test their simulations
# found best calibrated after normal theory two-stage estimation; it is
# computed here from lavaan's exported accessors and the observed
# information of the saturated model under missing data, which this file
# computes itself. Three traps in lavaan 0.7-3 shape the code: through the
# bare lavaan() interface the intercepts must be written into the syntax
# (fa_confirmatory() does, as for FIML); with auxiliary variables lavInspect(fit,
# "gamma") is not the matrix lavaan's own standard errors use, so the
# information is never read from it; and a case observed on an auxiliary
# variable alone is read by the first stage but left out of lavaan's N,
# so a fit with such cases is computed again with every case it read
# (section C2).


# ---------------------------------------------------------------------------
# (A) Spellings
# ---------------------------------------------------------------------------

# The values of `missing` that request two-stage estimation, in DMAR's
# snake_case spelling and in lavaan's dotted one (accepted by fa_confirmatory() and
# measurement_invariance(), which pass `missing` to lavaan, as they accept
# lavaan's "ml" beside "fiml").
.TWO_STAGE_VALUES <- c("two_stage", "robust_two_stage",
                       "two.stage", "robust.two.stage")

# .ts_lavaan_missing(missing)
#
# lavaan's value for a two-stage request: "two.stage" for normal theory
# standard errors (Savalei & Bentler, 2009) and "robust.two.stage" for
# standard errors robust to nonnormality (Yuan & Bentler, 2000; Savalei &
# Falk, 2014). NULL when `missing` does not request two-stage estimation.
.ts_lavaan_missing <- function(missing) {
  switch(missing,
         two_stage = , two.stage = "two.stage",
         robust_two_stage = , robust.two.stage = "robust.two.stage",
         NULL)
}

# .ts_is_two_stage_fit(fit)
#
# Whether a lavaan fit was estimated in two stages.
.ts_is_two_stage_fit <- function(fit) {
  isTRUE(lavaan::lavInspect(fit, "options")$missing %in%
           c("two.stage", "robust.two.stage"))
}


# ---------------------------------------------------------------------------
# (B) The saturated model under missing data
# ---------------------------------------------------------------------------

# .ts_dup_matrix(p)
#
# The duplication matrix D_p, with vec(S) = D_p vech(S) for a symmetric
# p x p matrix S and vech(S) its lower triangle taken column by column,
# lavaan's order for the covariance moments.
.ts_dup_matrix <- function(p) {
  idx <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  D <- matrix(0, p * p, nrow(idx))
  D[cbind((idx[, 2L] - 1L) * p + idx[, 1L], seq_len(nrow(idx)))] <- 1
  D[cbind((idx[, 1L] - 1L) * p + idx[, 2L], seq_len(nrow(idx)))] <- 1
  D
}

# .ts_sat_information(X, mu, S)
#
# The observed information matrix of the saturated normal model under
# missing data, summed over the cases (Savalei & Bentler, 2009, Eq. 5, times
# N), evaluated at the means `mu` and covariance matrix `S`, in lavaan's
# order of the moments: the p means, then vech(S). A case contributes
# through the variables it observes; for a pattern with observed set o, at
# Sigma_o^-1 = Si and the centered observations r_i,
#   d2/dmu dmu'   : Si
#   d2/dmu dvec'  : (Si r_i)' (x) Si
#   d2/dvec dvec' : Si (x) (Si r_i r_i' Si - Si / 2)
# each padded to the p variables and summed over the cases of the pattern,
# with the covariance block reduced to vech(S) by D_p once at the end. This
# is lavaan's h1.information.observed of a FIML fit with h1.information =
# "unstructured" and observed.information = "h1", times N (the tests check
# the agreement); on a two-stage fit lavInspect() returns the complete data
# information there instead, which is the stage two view.
#
# Inputs:
#   X  : an N x p numeric matrix, NA where a value is missing. Cases that
#        observe nothing contribute nothing.
#   mu : the p means at which the information is evaluated.
#   S  : the p x p covariance matrix at which it is evaluated.
#
# Returns:
#   The (p + p(p + 1)/2) square information matrix, summed over the cases.
.ts_sat_information <- function(X, mu, S) {
  X <- as.matrix(X)
  p <- ncol(X)
  Dp <- .ts_dup_matrix(p)
  Amm <- matrix(0, p, p)
  Ams <- matrix(0, p, p * p)
  Ass <- matrix(0, p * p, p * p)
  obs <- !is.na(X)
  key <- apply(obs, 1L, function(r) paste(which(r), collapse = ","))
  for (k in unique(key)) {
    rows <- which(key == k)
    idx <- which(obs[rows[1L], ])
    if (!length(idx)) next
    nk <- length(rows)
    Si <- solve(S[idx, idx, drop = FALSE])
    SR <- sweep(X[rows, idx, drop = FALSE], 2L, mu[idx]) %*% Si
    Sip <- matrix(0, p, p)
    Sip[idx, idx] <- Si
    M <- matrix(0, p, p)
    M[idx, idx] <- crossprod(SR) - 0.5 * nk * Si
    Amm <- Amm + nk * Sip
    cs <- numeric(p)
    cs[idx] <- colSums(SR)
    Ams <- Ams + kronecker(t(cs), Sip)
    Ass <- Ass + kronecker(Sip, M)
  }
  AmsD <- Ams %*% Dp
  rbind(cbind(Amm, AmsD), cbind(t(AmsD), crossprod(Dp, Ass %*% Dp)))
}

# .ts_schur_items(A, p_all, p_items)
#
# The information about the moments of the first p_items variables once the
# moments of the remaining (auxiliary) variables are estimated with them:
# the Schur complement A_bb - A_bg A_gg^-1 A_gb of the information over all
# p_all variables, with b the item means and the vech elements whose two
# variables are both items (Savalei & Bentler, 2009, Eqs. 11--14). Its
# inverse is the item block of the inverse of the full information, which
# is how lavaan forms the two-stage standard errors with auxiliary
# variables. The item elements keep lavaan's order.
.ts_schur_items <- function(A, p_all, p_items) {
  vi <- which(lower.tri(diag(p_all), diag = TRUE), arr.ind = TRUE)
  keep <- c(seq_len(p_items),
            p_all + which(vi[, 1L] <= p_items & vi[, 2L] <= p_items))
  other <- setdiff(seq_len(nrow(A)), keep)
  A[keep, keep, drop = FALSE] -
    A[keep, other, drop = FALSE] %*%
    solve(A[other, other, drop = FALSE], A[other, keep, drop = FALSE])
}

# .ts_em_moments(Y)
#
# The stage one estimates of the means and covariances of the columns of Y,
# by lavaan's own EM algorithm with its own tolerance, so that they are the
# moments a lavaan two-stage fit uses: a saturated model of the columns with
# missing = "two.stage" and do.fit = FALSE, which computes the sample
# statistics (the EM step) and fits nothing. With auxiliary variables this
# reproduces the augmented moments of lavaan's stage one to about 1e-11,
# where a FIML fit of the saturated model stops at its own tolerance (about
# 1e-5 apart on unscaled data). The standard errors and test are named as
# lavaan sets them for a two-stage fit, which keeps lavaan from saying it
# replaced them. The warnings such a call can raise concern the data (for
# example variances of very different sizes, or an EM that stopped at its
# iteration limit) and were raised already by the fit being evaluated, so
# they are muffled here by default; `quiet = FALSE` lets them through, for
# a caller whose estimates rest on these moments alone (the bootstrap
# replications of .omega_fit_cfa_ts_estimates()).
#
# Returns:
#   A list with `mean` (a numeric vector) and `cov` (a matrix), in the
#   column order of Y.
.ts_em_moments <- function(Y, quiet = TRUE) {
  Y <- as.data.frame(Y)
  v <- paste0("v", seq_len(ncol(Y)))
  names(Y) <- v
  sat_call <- function() {
    lavaan::lavaan(
      paste(v, "~~", v, collapse = "\n"), data = Y,
      missing = "two.stage", se = "two.stage", test = "satorra.bentler",
      meanstructure = TRUE, int.ov.free = TRUE, do.fit = FALSE)
  }
  sat <- if (quiet) suppressWarnings(sat_call()) else sat_call()
  ss <- lavaan::lavInspect(sat, "sampstat")
  list(mean = as.numeric(ss$mean[v]), cov = unclass(ss$cov)[v, v])
}


# ---------------------------------------------------------------------------
# (C) Auxiliary variables of a fit
# ---------------------------------------------------------------------------

# .ts_fit_aux(fit)
#
# The auxiliary variables of a lavaan fit (lavaan's aux = argument), which
# enter stage one only, with their data, group by group in the rows of the
# fit's data. lavaan 0.7-3 exposes neither through lavInspect(), so they are
# read from the fit's data slots, guarded: when the call carried aux = and
# the slots are not where this version of lavaan keeps them, the function
# stops rather than compute a statistic that ignores the auxiliaries.
#
# Returns:
#   NULL when the fit has no auxiliary variables; otherwise a list with
#   `names` (character) and `data` (one numeric matrix per group).
.ts_fit_aux <- function(fit) {
  call_aux <- tryCatch(lavaan::lavInspect(fit, "call")$aux,
                       error = function(e) NULL)
  nm <- tryCatch(fit@Data@ov.names.aux, error = function(e) NULL)
  dat <- tryCatch(fit@Data@aux, error = function(e) NULL)
  ng <- lavaan::lavInspect(fit, "ngroups")
  ok <- is.list(nm) && length(nm) == ng && is.list(dat) &&
    length(dat) == ng
  if (!ok) {
    if (!is.null(call_aux)) {
      .dmar_refuse_argument(c(
        "This two-stage fit has auxiliary variables (aux =), which the ",
        "installed version of lavaan does not keep where DMAR reads them, ",
        "so the residual based statistic cannot account for them. Refit ",
        "with missing = \"robust.two.stage\", whose scaled statistic ",
        "lavaan computes itself."), arg = NULL)
    }
    return(NULL)
  }
  names_aux <- unique(unlist(nm, use.names = FALSE))
  if (!length(names_aux)) return(NULL)
  data <- lapply(seq_along(dat), function(g) {
    x <- as.matrix(dat[[g]])
    if (length(nm[[g]]) != ncol(x)) {
      .dmar_refuse_argument(c(
        "The two-stage auxiliary names do not match their data columns. ",
        "Refit from the original raw data before reporting corrected ",
        "statistics."), arg = NULL)
    }
    colnames(x) <- nm[[g]]
    x
  })
  list(names = names_aux, data = data)
}

# .ts_dots_aux(dots)
#
# The auxiliary variables a call passes to lavaan through `...` (its aux
# argument), as fa_confirmatory() and measurement_invariance() receive them: the
# element of `dots` whose name is aux or a prefix of it, which lavaan's
# partial matching would take for aux. NULL when there is none.
.ts_dots_aux <- function(dots) {
  nm <- names(dots)
  if (is.null(nm)) return(NULL)
  hit <- which(nzchar(nm) & startsWith("aux", nm))
  if (!length(hit)) return(NULL)
  aux <- dots[[hit[1L]]]
  if (is.character(aux) && length(aux)) aux else NULL
}

# ---------------------------------------------------------------------------
# (C2) The cases the first stage reads
# ---------------------------------------------------------------------------
#
# A case observed on an auxiliary variable but on none of the model's
# variables is read by the first stage: the EM estimates of the means and
# covariances of the items and the auxiliaries together use its auxiliary
# values, and through them it informs the item moments the second stage
# fits. This is what makes two-stage estimation with the auxiliary
# consistent when the item block goes missing at random given the
# auxiliaries; the first stage likelihood is summed over every case
# (Savalei & Bentler, 2009, Eq. 10). lavaan 0.7-3 reads such a case in its
# first stage but counts it as empty and leaves it out of its N
# (lavInspect(fit, "nobs")), and N is what turns the per case quantities
# of the two stages into the totals behind every inferential number:
#
#   - the covariance matrix of the estimates is NVarCov / N, with NVarCov
#     the sandwich of the per case information of the second stage and
#     the asymptotic covariance matrix of the stage one moments, a per case
#     quantity that lavaan averages over every case the first stage read
#     (lav_mvnorm_missing_information_observed_samplestats() divides by the
#     number of those cases), so every variance is N_read / N times too
#     large and every standard error sqrt(N_read / N) times;
#   - the naive statistic is 2 N F, so it and the scaled statistic built on
#     it, the corrected indices, the baseline statistic, and the RMSEA
#     multiplier rest on N;
#   - the log likelihood is N times the complete data log likelihood per
#     case at the EM estimates, and BIC takes log N;
#   - in a multiple group fit the groups are weighted by N_g / N in the
#     discrepancy function the second stage minimizes, F = sum_g (N_g / N)
#     F_g, in the information, and in the trace behind the scaling
#     constant. When a constraint holds parameters of several groups
#     together, the estimates depend on those weights: lavaan's minimize
#     the discrepancy with the groups weighted by their cases with an
#     observed variable of the model, and the two-stage estimates with
#     every case read weight them by the cases the first stage read. When
#     the groups' counts grow in different ratios the two differ, and with
#     the estimates every value computed from them.
#
# On six Holzinger and Swineford tests with the whole block missing for 77
# of 301 pupils, at random given general information (the auxiliary),
# every standard error, normal theory and robust, was 1.1592 =
# sqrt(301/224) times the two-stage sandwich computed from all 301 rows
# (Savalei & Bentler, 2009, Eqs. 6 and 13). Setting such cases aside, which
# the first versions of 1.1.0 did, made the counts agree by conditioning on
# an observed item, which makes the estimates inconsistent in exactly the
# case an auxiliary is named for: on that example it moved the intercept
# of paragraph comprehension from 9.21 (rows kept; 9.18 with complete
# data) to 8.75, about two standard errors, and in 60 replications of the
# mechanism its mean error was -0.49 against -0.01 with the rows kept.
# The cases are therefore kept, and every value is computed with the
# number of cases the first stage read (.ts_all_rows()). In several groups
# with constraints across them the estimates are computed again too: on
# the same data by school, 130 and 94 cases with an observed item grow to
# 156 and 145 read, and the scalar rung's naive statistic is 64.139 at the
# estimates with every case read, the value semTools' twostage() gives,
# against 64.501 with the new counts at lavaan's estimates.

# .ts_rows_read(fit)
#
# The cases the first stage of a two-stage fit read and the cases lavaan's
# N counts, group by group. A case is read when any of the model's
# variables or of the auxiliary variables is observed on it; a case with
# nothing observed is dropped by lavaan from both stages. The auxiliary
# data are read with .ts_fit_aux(), whose rows are those of lavInspect(fit,
# "data").
#
# Returns:
#   NULL when the fit is not a two-stage fit or has no auxiliary variables,
#   and when they cannot be read or their rows do not line up with the
#   data's; otherwise a list with `read` and `nobs`, integer vectors with
#   one element per group.
.ts_rows_read <- function(fit) {
  if (!.ts_is_two_stage_fit(fit)) return(NULL)
  aux <- tryCatch(.ts_fit_aux(fit), error = function(e) NULL)
  if (is.null(aux)) return(NULL)
  X <- lavaan::lavInspect(fit, "data")
  ng <- lavaan::lavInspect(fit, "ngroups")
  if (ng == 1L) X <- list(X)
  if (length(X) != length(aux$data)) return(NULL)
  read <- vapply(seq_len(ng), function(g) {
    items <- as.matrix(X[[g]])
    extra <- as.matrix(aux$data[[g]])
    if (nrow(items) != nrow(extra)) return(NA_integer_)
    sum(rowSums(!is.na(items)) > 0L | rowSums(!is.na(extra)) > 0L)
  }, integer(1L))
  if (anyNA(read)) return(NULL)
  list(read = read,
       nobs = as.integer(unlist(lavaan::lavInspect(fit, "nobs"))))
}

# .ts_aux_only_cases(fit)
#
# The cases a two-stage fit's first stage read that its N leaves out: in a
# fit lavaan made, those observed on an auxiliary variable but on none of
# the model's variables; in a fit .ts_all_rows() computed again, none. The
# count decides only a warning or a line, so a fit whose auxiliary data
# cannot be read is not checked rather than refused.
#
# Returns:
#   The number of such cases, summed over the groups (0 when there are
#   none or they cannot be counted).
.ts_aux_only_cases <- function(fit) {
  rows <- .ts_rows_read(fit)
  if (is.null(rows)) return(0L)
  as.integer(sum(pmax(rows$read - rows$nobs, 0L)))
}

# .ts_all_rows(fit)
#
# A two-stage fit computed again with its N set to the cases its first
# stage read (.ts_rows_read()), group by group, so that every estimate,
# standard error, test statistic, fit index, and information criterion
# lavaan derives from the fit rests on the cases the moments describe.
# lavaan is called with the fit's own data, model, parameter table, and
# options, its sample statistics carrying the counts of the cases read,
# and computes the standard errors, the tests, the baseline model, and the
# log likelihoods as it does for any fit. lavaan's own two-stage formulas
# then use the right N; nothing is rescaled here.
#
# The estimates are those of the second stage with every case read, the
# minimum of sum_g (read_g / N_read) F_g. In one group, and in several
# whose parameters no constraint holds together (.ts_groups_linked()),
# they do not depend on the counts, so lavaan's are kept: its optimizer is
# switched off (optim.method = "none"), and lavaan starts at the estimates
# and stays there (they come back exactly, or under equality constraints
# to its rounding error of about 1e-14). In several groups held together
# by a constraint they move with the counts, and lavaan's optimizer starts
# at its estimates and goes on to the minimum with the new weights; these
# are the estimates of semTools' twostage(), whose first stage counts
# every case. The tests check the result against the two-stage sandwich
# and the trace behind the scaling constant built from the data over every
# case read, in one group and in several with constraints across groups,
# against an ordinary maximum likelihood fit of the stage one moments with
# the counts of the cases read, and against the naive statistic of
# semTools' twostage().
#
# The fit keeps lavaan's record of the optimizer at its estimates (the
# function value, which fitMeasures() reports as fmin, uses the new
# weights, so that the naive statistic is 2 N fmin in several groups too),
# with the iterations of the fit's estimation, added to those of the new
# one when the estimates moved; its optimizer and call are put back. The
# data slot counts the cases read too, so lavaan's summary() says that
# they were used, and lavTestLRT() and anova(), which compare those
# counts, warn that a difference test between the fit and one lavaan made
# of the same data on its smaller N is invalid rather than report it as if
# the two shared their cases; the rows of the data, and with them
# lavPredict(), are unchanged.
#
# The new fit is accepted only when lavaan computed it, converged, with the
# counts in place, and, when the estimates do not depend on the counts,
# with them unchanged; otherwise the fit lavaan made is returned, and
# `corrected` says so, for the caller's warning and footer. A fit without
# such cases is returned as it is.
#
# Returns:
#   A list with `fit` (the fit to report), `read` and `nobs` (the counts
#   per group, NULL when the fit has no auxiliary variables), `extra` (the
#   cases lavaan's N left out, summed over the groups), `corrected` (TRUE
#   when the fit was computed again, FALSE when it could not be, NA when
#   there was nothing to correct), and `reestimated` (TRUE when the
#   estimates depend on the counts and the fit had such cases, so that
#   the estimates reported are not lavaan's when `corrected` is TRUE, or
#   are lavaan's in place of the right ones when it is FALSE).
.ts_all_rows <- function(fit) {
  rows <- .ts_rows_read(fit)
  extra <- if (is.null(rows)) 0L else sum(pmax(rows$read - rows$nobs, 0L))
  out <- list(fit = fit, read = rows$read, nobs = rows$nobs,
              extra = as.integer(extra), corrected = NA, reestimated = FALSE)
  if (extra == 0L) return(out)
  linked <- .ts_groups_linked(fit)
  out$reestimated <- linked
  refit <- tryCatch(.ts_refit_n(fit, rows$read, optimize = linked),
                    error = function(e) NULL)
  est <- lavaan::coef(fit)
  ok <- inherits(refit, "lavaan") &&
    isTRUE(lavaan::lavInspect(refit, "converged")) &&
    identical(as.integer(unlist(lavaan::lavInspect(refit, "nobs"))),
              rows$read) &&
    (linked || isTRUE(all(abs(lavaan::coef(refit) - est) <=
                            1e-10 * pmax(1, abs(est)))))
  if (!ok) {
    out$corrected <- FALSE
    return(out)
  }
  refit@optim$iterations <- fit@optim$iterations +
    if (linked) refit@optim$iterations else 0L
  refit@Options$optim.method <- fit@Options$optim.method
  refit@call <- fit@call
  refit@Data@nobs <- as.list(rows$read)
  out$fit <- refit
  out$corrected <- TRUE
  out
}

# .ts_groups_linked(fit)
#
# Whether a constraint holds parameters of several groups of a fit
# together, so that the estimates depend on the weights N_g / N of the
# groups in the discrepancy function: a free parameter shared by several
# groups (lavaan's ceq.simple = TRUE), or an equality or inequality
# constraint whose sides name parameters of several groups (group.equal,
# a label shared across groups). A constraint that names a parameter this
# cannot place in a group, a defined parameter among them, is taken to
# link the groups, which costs only a search that stays where it starts.
#
# Returns:
#   TRUE or FALSE (FALSE in one group).
.ts_groups_linked <- function(fit) {
  if (lavaan::lavInspect(fit, "ngroups") < 2L) return(FALSE)
  pt <- lavaan::parTable(fit)
  free <- pt$free > 0L & pt$group > 0L
  if (any(tapply(pt$group[free], pt$free[free],
                 function(g) length(unique(g))) > 1L)) {
    return(TRUE)
  }
  con <- which(pt$op %in% c("==", "<", ">"))
  if (!length(con)) return(FALSE)
  nm <- c(pt$plabel, pt$label)
  gr <- c(pt$group, pt$group)
  named <- nzchar(nm)
  nm <- nm[named]
  gr <- gr[named]
  for (i in con) {
    vars <- tryCatch(unique(c(all.vars(str2lang(pt$lhs[i])),
                              all.vars(str2lang(pt$rhs[i])))),
                     error = function(e) NA_character_)
    if (anyNA(vars) || !all(vars %in% nm)) return(TRUE)
    g <- unique(gr[nm %in% vars])
    if (length(g) > 1L || any(g < 1L)) return(TRUE)
  }
  FALSE
}

# .ts_refit_n(fit, n, optimize = FALSE)
#
# The lavaan call behind .ts_all_rows(): the fit computed again from its
# own slots, with `n` cases in its groups, and its optimizer off, or on
# (`optimize`), started at the fit's estimates, when the estimates depend
# on the counts. Its data and model are those of the fit, so the warnings
# it could give (lavaan's post-estimation check of the estimates, notes
# about the data) repeat ones the fit gave already and are muffled, as are
# its messages; .ts_all_rows() checks that it converged.
.ts_refit_n <- function(fit, n, optimize = FALSE) {
  ss <- fit@SampleStats
  ss@nobs <- as.list(as.integer(n))
  ss@ntotal <- as.integer(sum(n))
  opt <- fit@Options
  if (!optimize) opt$optim.method <- "none"
  suppressMessages(suppressWarnings(lavaan::lavaan(
    slotOptions = opt, slotParTable = fit@ParTable, slotSampleStats = ss,
    slotData = fit@Data, slotModel = fit@Model, slotCache = fit@Cache)))
}

# .ts_aux_only_note(fun, extra, read, nobs, what)
#
# The message fa_confirmatory(), measurement_invariance(), and the reliability
# family give when the first stage of a two-stage fit reads cases observed
# on an auxiliary variable alone: how many, that they are kept, and what
# rests on every case read. fa_confirmatory() and measurement_invariance() give it
# once the fit has been computed again with every case
# (.ts_all_rows()), and only then; when it could not be, the warning of
# .ts_all_rows_warning() says what rests on lavaan's N instead. The
# reliability family gives it once per call, before any fit, and says
# only what holds either way. It has the class
# dmar_two_stage_aux_only_message, so that the family can muffle the one
# fa_confirmatory() would give for each of its fits. Its text ends with a newline,
# as the text base R's message() builds from its arguments does, so that
# the console's next output starts on a line of its own.
#
# Inputs:
#   fun   : the calling function, named in the message.
#   extra : the number of such cases.
#   read  : the cases the first stage read (total).
#   nobs  : the cases lavaan's N counts (total).
#   what  : what rests on every case read, a clause that follows "and".
.ts_aux_only_note <- function(fun, extra, read, nobs,
                              what = paste0("the standard errors, test ",
                                            "statistics, fit indices, and ",
                                            "information criteria are ",
                                            "computed on all ", read,
                                            " cases it reads")) {
  one <- extra == 1L
  message(messageCondition(paste0(
    fun, "(): ", extra, if (one) " case is" else " cases are",
    " observed on an auxiliary variable but on none of the model's ",
    "variables. Two-stage estimation keeps ", if (one) "it" else "them",
    ": the first stage reads ", if (one) "its" else "their",
    " auxiliary values, which inform the estimates of the item means and ",
    "covariances, and ", what, " (lavaan's two-stage estimator counts ",
    "only the ", nobs, " with an observed variable of the model).\n"),
    class = "dmar_two_stage_aux_only_message"))
}

# .ts_rows_line(rows, kind)
#
# The line under a fa_confirmatory() or measurement_invariance() table that names
# the N of a two-stage fit whose first stage read cases observed on an
# auxiliary variable alone (`rows`, from .ts_all_rows(), whose counts are
# summed over the groups here), or NULL when there were none. `kind` is
# that of .ts_footer(), which decides what the line says uses the N. When
# the estimates of a rung that constrains across groups depend on the
# counts (`rows$reestimated`), the line says that they were computed
# again too, or, when the fit could not be, that they weight the groups by
# lavaan's counts. `lavaan_missing` is the fit's ("two.stage" or
# "robust.two.stage"): the chi_square row of a normal theory two-stage fit
# is the residual based statistic, which sums the information of the
# saturated model over every case the first stage read and so does not
# rest on lavaan's N when the fit could not be computed again (its value
# is the same either way, at a rung whose estimates do not move); only the
# scaled statistic of a robust two-stage fit does.
.ts_rows_line <- function(rows, kind = c("cfa", "measurement",
                                         "invariance"),
                          lavaan_missing = "robust.two.stage") {
  kind <- match.arg(kind)
  residual <- identical(lavaan_missing, "two.stage")
  if (is.null(rows) || !isTRUE(rows$extra > 0L)) return(NULL)
  read <- sum(rows$read)
  nobs <- sum(rows$nobs)
  moved <- isTRUE(rows$reestimated)
  if (isTRUE(rows$corrected)) {
    paste0("N: the ", read, " cases the first stage read, ", rows$extra,
           " of them observed on an auxiliary variable alone, which ",
           "lavaan's own N (", nobs, ") leaves out; ",
           switch(kind,
                  cfa = "the standard errors, tests, indices, and criteria",
                  measurement = "the standard errors",
                  invariance = paste0(
                    if (moved) paste0("the estimates of the rungs that ",
                                      "constrain across groups and "),
                    "the statistics, indices, and difference tests of ",
                    "every rung")),
           " use all ", read)
  } else {
    paste0("N: lavaan's ", nobs, ", which leaves out ", rows$extra,
           " cases the first stage read (observed on an auxiliary variable ",
           "alone); DMAR could not compute the fit again with all ", read,
           ", so ",
           switch(kind,
                  cfa = paste0("the standard errors are about ",
                               format(signif(sqrt(read / nobs), 4)),
                               " times too large and the ",
                               if (residual) "indices and criteria" else
                                 "scaled statistic, indices, and criteria",
                               " rest on ", nobs, " cases",
                               if (residual) {
                                 paste0("; the residual based statistic in ",
                                        "chi_square uses all ", read)
                               }),
                  measurement = paste0("the standard errors are about ",
                                       format(signif(sqrt(read / nobs), 4)),
                                       " times too large"),
                  invariance = paste0(
                    "the ",
                    if (residual) "indices and difference tests" else
                      "scaled statistics, indices, and difference tests",
                    " rest on ", nobs, " cases",
                    if (residual) {
                      paste0(" (the residual based statistics in chi_square ",
                             "use every case read)")
                    },
                    if (moved) paste0(", and the estimates of the rungs ",
                                      "that constrain across groups weight ",
                                      "the groups by lavaan's counts",
                                      if (residual) {
                                        paste0(", which moves the residual ",
                                               "based statistics of those ",
                                               "rungs as well")
                                      }))))
  }
}

# .ts_all_rows_warning(fun, rows)
#
# The warning of a call whose two-stage fit could not be computed again
# with every case its first stage read (`rows$corrected` FALSE): which
# values then rest on lavaan's smaller N. The counts travel on the
# condition as `read` and `nobs`, so that the reliability family can say
# the same in its own terms. Class dmar_two_stage_n. `lavaan_missing` is
# the fit's, as in .ts_rows_line(): after normal theory two-stage
# estimation the residual based statistic uses every case read, and only
# the indices and the criteria (the difference tests of a ladder) rest on
# lavaan's N.
.ts_all_rows_warning <- function(fun, rows,
                                 lavaan_missing = "robust.two.stage") {
  if (!identical(rows$corrected, FALSE)) return(invisible(NULL))
  read <- sum(rows$read)
  nobs <- sum(rows$nobs)
  residual <- identical(lavaan_missing, "two.stage")
  ladder <- identical(fun, "measurement_invariance")
  .dmar_warn(
    paste0(
      fun, "(): the first stage of the two-stage fit read ", read, " cases, ",
      rows$extra, " of them observed on an auxiliary variable alone, which ",
      "lavaan's N of ", nobs, " leaves out, and DMAR could not compute the ",
      "fit again with all ", read, ". ",
      if (isTRUE(rows$reestimated)) {
        paste0("The estimates of a model that constrains parameters across ",
               "groups then weight the groups by lavaan's counts rather than ",
               "by the cases read, the ")
      } else {
        "The estimates are right, but the "
      },
      "standard errors are about ", format(signif(sqrt(read / nobs), 4)),
      " times too large, and ",
      if (!residual) {
        if (ladder) "the scaled statistics, " else "the scaled statistic, "
      },
      "the fit indices",
      if (ladder) ", and the difference tests" else
        if (residual) " and the information criteria" else
          ", and the information criteria",
      " rest on ", nobs, " cases",
      if (residual) {
        paste0("; the residual based statistic", if (ladder) "s",
               " in chi_square use", if (!ladder) "s", " every case read")
      },
      ". missing = \"fiml\" with the same auxiliary variables keeps every ",
      "case."),
    "two_stage_n", read = read, nobs = nobs, extra = rows$extra)
}

# .ts_coefficient_n_warning(read, nobs, extra)
#
# The warning the reliability family gives in place of fa_confirmatory()'s
# (.ts_all_rows_warning()) when the two-stage fit of a coefficient could
# not be computed again with every case its first stage read: the
# coefficient, a function of the estimates, is right, and its standard
# error is too large. A bootstrap interval refits both stages on each
# resample of every row and does not use lavaan's N. Class
# dmar_two_stage_n.
.ts_coefficient_n_warning <- function(read, nobs, extra) {
  .dmar_warn(
    paste0(
      "The first stage of the two-stage fit read ", read, " cases, ", extra,
      " of them observed on an auxiliary variable alone, which lavaan's N of ",
      nobs, " leaves out, and DMAR could not compute the fit again with all ",
      read, ". The coefficient is right, but its analytic standard error, ",
      "where one is reported, is about ", format(signif(sqrt(read / nobs), 4)),
      " times too large, which widens an interval built on it. A bootstrap ",
      "interval counts every case, as does missing = \"fiml\" with the same ",
      "auxiliary variables."),
    "two_stage_n", read = read, nobs = nobs, extra = extra)
}

# .ts_aux_only_warning(fit, rest, consequence, aside = "",
#                      where = "this two-stage fit")
#
# The warning sem_fit_measures() gives for a two-stage fit the user made
# with lavaan whose first stage read cases its N leaves out
# (.ts_aux_only_cases()): the rows it reads from lavaan that rest on the
# smaller N, which DMAR does not change in a fit it did not make.
# sem_information_criteria() and compare_sem() compute their values
# themselves and correct them instead (.ts_all_rows()), and give it only
# for a fit that could not be computed again; a fit that fa_confirmatory() or
# measurement_invariance() returns has no such cases left to warn about.
# The warning has the class dmar_two_stage_aux_only, so that compare_sem()
# can muffle the one each of its sem_information_criteria() calls would
# give and give its own, which names the fit. Nothing is warned for a fit
# that is not a two-stage fit.
#
# Inputs:
#   fit         : a lavaan fit.
#   rest        : what rests on the smaller N, with its verb ("the N row
#                 rests").
#   consequence : the clause after "so", what that means for the result.
#   aside       : a clause set off after "the smaller N", or "".
#   where       : how the fit is named.
#
# Returns:
#   The number of such cases, invisibly (0 when nothing was warned).
.ts_aux_only_warning <- function(fit, rest, consequence, aside = "",
                                 where = "this two-stage fit") {
  if (!.ts_is_two_stage_fit(fit)) return(invisible(0L))
  n <- .ts_aux_only_cases(fit)
  if (n == 0L) return(invisible(0L))
  one <- n == 1L
  .dmar_warn(
    sprintf(paste0(
      "The first stage of %s read %d %s observed on an auxiliary variable ",
      "but on none of the model's variables, which lavaan leaves out of its ",
      "N of %s. The EM estimates of the means and covariances include %s, ",
      "while %s on the smaller N%s, so %s. fa_confirmatory(missing = \"two_stage\") ",
      "and measurement_invariance() keep such cases and compute every value ",
      "with all of them; for another model, FIML (missing = \"ml\") with the ",
      "same auxiliary variables keeps them as well."),
      where, n, if (one) "case" else "cases",
      format(lavaan::lavInspect(fit, "ntotal")), if (one) "it" else "them",
      rest, aside, consequence),
    "two_stage_aux_only")
  invisible(n)
}


# ---------------------------------------------------------------------------
# (D) The residual based statistic
# ---------------------------------------------------------------------------

# .ts_residual_test(fit)
#
# The residual based statistic of Savalei and Bentler (2009, Eq. 8; Eq. 14
# with auxiliary variables) for a two-stage fit,
#
#   T = e' [W - W Delta (Delta' W Delta)^-1 Delta' W] e,
#
# with e the stage one moments less the moments the fitted model implies
# (lavaan's wls.obs - wls.est, means then vech), Delta their derivatives
# with respect to the free parameters, and W the observed information of
# the saturated model under missing data, summed over the cases (with
# auxiliary variables, its Schur complement on the item moments). In a
# multiple group fit e and Delta are stacked over the groups and W is block
# diagonal. Because the bracket annihilates Delta, T is referred to chi
# square on the model's degrees of freedom whatever weight stage two used,
# under normal data missing at random. The information is evaluated at the
# stage one estimates: lavaan's own (sampstat) without auxiliary variables,
# and with them the augmented EM estimates from .ts_em_moments().
#
# The statistic is the normal theory one. Its asymptotically distribution
# free form (W replaced by the inverse of A^-1 B A^-1) is not reported: no
# study evaluates it after two-stage estimation, and the nonnormal versions
# of residual based statistics "have been found to perform very poorly in
# small sample sizes" (Savalei & Bentler, 2009, citing Yuan & Bentler,
# 1998).
#
# A fit that conditions on observed exogenous variables (fixed.x = TRUE,
# lavaan::sem()'s default when the model regresses on an observed
# variable, or conditional.x = TRUE) is refused, with the advice to refit
# it with fixed.x = FALSE (and conditional.x = FALSE, since lavaan 0.7-3
# stops at fixed.x = FALSE beside conditional.x = TRUE).
# The statistic is derived for cases drawn from one multivariate normal
# distribution of every variable (Savalei & Bentler, 2009), and semTools'
# twostage() sets fixed.x = FALSE before it computes it. When every
# exogenous variable is complete, the estimates and the residuals are the
# same under both settings. When one has missing values they are not:
# under fixed.x = TRUE lavaan deletes every case with a missing exogenous
# variable before either stage (with its warning at the fit), robust
# two-stage estimation included, and under fixed.x = FALSE the first stage
# keeps those cases and estimates the exogenous moments with the others,
# so the refit has more cases and other estimates (on three verbal tests
# regressed on age with 8 ages and 90 sentence scores missing, 293 cases
# against 301, and a slope of -0.746 against -0.737). The message says so,
# and says how to keep the fit's estimates: the refit with fixed.x = FALSE
# on the rows the fit used, which are the rows on which every exogenous
# variable is observed (lavInspect(fit, "case.idx"), group by group with
# several groups), has the fit's 293 cases and its estimates to within the
# optimizer's tolerance (2e-6 there). Under conditional.x = TRUE lavaan
# deletes the same cases and also leaves out of its sample size a case
# that observes the exogenous variables alone, which case.idx still lists:
# with two such cases added to the data above, the refit with
# conditional.x = FALSE and fixed.x = FALSE on the rows case.idx lists
# (those on which age is observed) has the fit's estimates to within 1e-5,
# though its sample size counts the two, and the refit on the rows without
# them misses by 0.003. The message names the rows in words, since
# fa_confirmatory() and measurement_invariance() give it for a fit the
# caller does not hold. The refusal is an error of class
# dmar_two_stage_fixed_x (a dmar_refusal) that carries the exogenous
# variables as `variables`. Either way the projection
# needs the derivatives of the moments with respect to the means and
# covariances of the exogenous variables, which move the moments of the
# variables regressed on them too. Unit columns for the exogenous moments
# in their place reproduce the statistic only in special cases (the
# exogenous variables all predicting one factor with free slopes), and
# otherwise give another statistic that passes the rank check below (with
# two correlated covariates predicting different factors, or a covariate
# predicting a factor and one of its indicators, or a fixed slope).
#
# Returns:
#   A list with `statistic`, `df`, and `p_value`.
.ts_residual_test <- function(fit) {
  o <- lavaan::lavInspect(fit, "options")
  if (!isTRUE(o$missing %in% c("two.stage", "robust.two.stage"))) {
    .dmar_refuse_argument(
      "The residual based statistic is defined for a two-stage fit.",
      arg = NULL)
  }
  if (isTRUE(o$.sampling.weights) ||
      length(lavaan::lavInspect(fit, "cluster")) ||
      lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(c(
      "The residual based statistic of Savalei and Bentler (2009) ",
      "assumes independent, equally weighted cases; this fit has sampling ",
      "weights or clustered cases. Use missing = \"robust.two.stage\" ",
      "(DMAR's \"robust_two_stage\"), whose scaled statistic lavaan ",
      "computes."), arg = NULL)
  }
  ov_x <- lavaan::lavNames(fit, "ov.x")
  if ((isTRUE(o$fixed.x) || isTRUE(o$conditional.x)) && length(ov_x)) {
    conditional <- isTRUE(o$conditional.x)
    setting <- if (conditional) "conditional.x = TRUE" else "fixed.x = TRUE"
    .dmar_refuse(c(
      "The residual based statistic treats every observed variable as ",
      "random, and this fit conditions on its exogenous observed variables ",
      "(", paste(ov_x, collapse = ", "), "). Refit with ",
      if (conditional) {
        c("conditional.x = FALSE and fixed.x = FALSE (lavaan does not take ",
          "fixed.x = FALSE beside conditional.x = TRUE)")
      } else "fixed.x = FALSE",
      ". When those variables are complete the refit leaves the estimates ",
      "as they are; when one has missing values, lavaan has deleted the ",
      "cases missing it under ", setting, ", and the refit keeps them, so ",
      "its estimates differ. To keep the fit's estimates, refit ",
      if (conditional) "that way" else "with fixed.x = FALSE",
      " on the rows the fit used, those on which every one of those ",
      "variables is observed (for a lavaan fit, the rows ",
      "lavaan::lavInspect(fit, \"case.idx\") lists, one set per group). ",
      "missing = \"robust.two.stage\", whose scaled statistic lavaan ",
      "computes, is the alternative, but under ", setting, " it deletes ",
      "the same cases."), "two_stage_fixed_x", variables = ov_x)
  }
  ng <- lavaan::lavInspect(fit, "ngroups")
  as_list <- function(x) if (ng == 1L) list(x) else x
  X  <- as_list(lavaan::lavInspect(fit, "data"))
  ss <- as_list(lavaan::lavInspect(fit, "sampstat"))
  wo <- as_list(lavaan::lavInspect(fit, "wls.obs"))
  we <- as_list(lavaan::lavInspect(fit, "wls.est"))
  D  <- do.call(rbind, lapply(as_list(lavaan::lavInspect(fit, "delta")),
                              unclass))
  aux <- .ts_fit_aux(fit)
  W <- vector("list", ng)
  e <- numeric(0)
  for (g in seq_len(ng)) {
    Xg <- as.matrix(X[[g]])
    ov <- colnames(Xg)
    p <- ncol(Xg)
    if (length(wo[[g]]) != p + p * (p + 1L) / 2L) {
      .dmar_refuse_argument(c(
        "The residual based statistic needs the means and covariances of ",
        "every observed variable among the moments of the fit, which this ",
        "fit does not model. Use missing = \"robust.two.stage\"."), arg = NULL)
    }
    W[[g]] <- if (is.null(aux)) {
      .ts_sat_information(Xg, as.numeric(ss[[g]]$mean[ov]),
                          unclass(ss[[g]]$cov)[ov, ov, drop = FALSE])
    } else {
      Y <- cbind(Xg, aux$data[[g]])
      m <- .ts_em_moments(Y)
      .ts_schur_items(.ts_sat_information(Y, m$mean, m$cov),
                      p_all = ncol(Y), p_items = p)
    }
    e <- c(e, as.numeric(wo[[g]]) - as.numeric(we[[g]]))
  }
  Wb <- W[[1L]]
  if (ng > 1L) {
    n_all <- sum(vapply(W, nrow, integer(1L)))
    Wb <- matrix(0, n_all, n_all)
    at <- 0L
    for (g in seq_len(ng)) {
      k <- nrow(W[[g]])
      Wb[at + seq_len(k), at + seq_len(k)] <- W[[g]]
      at <- at + k
    }
  }
  df <- as.numeric(lavaan::fitMeasures(fit, "df"))
  D <- .ts_constrained_delta(fit, D)
  # The projection must be onto the moments the model can move, of
  # dimension the number of moments less the degrees of freedom; a
  # derivative matrix of another rank would give a statistic on other
  # degrees of freedom than it is referred to.
  if (qr(D)$rank != length(e) - df) {
    .dmar_refuse_argument(c(
      "The residual based statistic is not available for this model: ",
      "its derivative matrix does not have the rank its degrees of ",
      "freedom imply. Use missing = \"robust.two.stage\"."), arg = NULL)
  }
  WD <- Wb %*% D
  P <- Wb - WD %*% solve(crossprod(D, WD), t(WD))
  # The bracket is nonnegative definite, so the statistic is nonnegative,
  # and with zero degrees of freedom it is the zero matrix (Delta spans
  # every moment), so the statistic is zero; what remains of either in
  # floating point is a rounding error of the projection.
  stat <- if (df > 0) max(drop(crossprod(e, P %*% e)), 0) else 0
  list(statistic = stat, df = df,
       p_value = if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE)
                 else NA_real_)
}


# .ts_constrained_delta(fit, D)
#
# The derivatives of the moments with respect to the parameters the model
# leaves free once its equality constraints hold. lavInspect(fit, "delta")
# differentiates with respect to the unconstrained free parameters, so
# that a label shared by several parameters, or the equality constraints
# of group.equal, contribute one column per parameter. A projection on that
# matrix removes the directions the constraints forbid as well, so the
# statistic falls on the degrees of freedom of the unconstrained model: on
# the two-group, metric model of the tests (6 constraints) it gave 81.55,
# barely above the configural 81.41, against 89.81 with the constraints
# respected, a difference of 8.39 on 6 df beside the scaled difference of
# 7.92. (semTools 0.5-9's twostage() reports the 81.55.) lavaan 0.7-3
# keeps the constraints in one of two forms, and each is undone here.
#
#   - By default (ceq.simple = FALSE) it writes them as == rows over the
#     unconstrained parameters, and D is multiplied by a basis K of the
#     null space of the constraints' Jacobian, the directions in which the
#     parameters can move with every equality still holding.
#   - Under ceq.simple = TRUE, when every constraint is a simple equality
#     of two parameters, it writes no == row and instead stores the matrix
#     K that maps the parameters it estimates to the unconstrained ones
#     (fit@Model@ceq.simple.K, unconstrained by estimated); D is
#     multiplied by that K. Read from the == rows alone, such a fit kept
#     a column of D per unconstrained parameter, and the rank check
#     refused it (the metric rung of measurement_invariance() with
#     ceq.simple = TRUE passed through ...). With any general constraint
#     (b == 1.5*a) lavaan ignores ceq.simple and writes every constraint
#     as an == row, the first case. The slots are read guarded, as in
#     .ts_fit_aux(): a fit whose lavaan does not have them leaves D
#     unconstrained, and the rank check refuses it.
#
# On complete data the statistic then equals lavaan's own
# browne.residual.nt for the constrained model, under either setting,
# which the tests check. Inequality constraints, which hold with equality
# or not depending on the solution, are refused.
.ts_constrained_delta <- function(fit, D) {
  pt <- lavaan::parTable(fit)
  # lavaan writes an inequality on one parameter and a number (a > 10) as a
  # bound on that parameter rather than as a constraint; a parameter that
  # sits on its bound is as good as fixed, and the projection would treat
  # it as free.
  free <- pt$free > 0L
  on_bound <- function(b) {
    if (is.null(b)) return(FALSE)
    b <- b[free]
    is.finite(b) & abs(pt$est[free] - b) <= 1e-6 * pmax(1, abs(b))
  }
  if (any(pt$op %in% c("<", ">")) || any(on_bound(pt$lower)) ||
      any(on_bound(pt$upper))) {
    .dmar_refuse_argument(c(
      "The residual based statistic is not available for a model with ",
      "inequality constraints or a parameter on a bound. Use missing = ",
      "\"robust.two.stage\"."), arg = NULL)
  }
  if (isTRUE(tryCatch(fit@Model@ceq.simple.only, error = function(e) FALSE))) {
    K <- tryCatch(fit@Model@ceq.simple.K, error = function(e) NULL)
    if (!is.matrix(K) || nrow(K) != ncol(D)) {
      .dmar_refuse_argument(c(
        "The residual based statistic could not read the equality ",
        "constraints of this model (lavaan's ceq.simple = TRUE). Refit ",
        "it with ceq.simple = FALSE, which gives the same estimates, or ",
        "use missing = \"robust.two.stage\"."), arg = NULL)
    }
    return(D %*% K)
  }
  if (!any(pt$op == "==")) return(D)
  J <- unclass(lavaan::lavInspect(fit, "con.jac"))
  if (!is.matrix(J) || ncol(J) != ncol(D)) {
    .dmar_refuse_argument(c(
      "The residual based statistic could not read the constraints of ",
      "this model. Use missing = \"robust.two.stage\"."), arg = NULL)
  }
  q <- qr(t(J))
  K <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
  D %*% K
}


# ---------------------------------------------------------------------------
# (E) What a table says about two-stage estimation
# ---------------------------------------------------------------------------

# .ts_test_label(lavaan_missing)
#
# The name of the single-model test a two-stage fit reports: the residual
# based statistic after normal theory two-stage estimation (Savalei &
# Bentler, 2009) and lavaan's scaled statistic after robust two-stage
# estimation (Yuan & Bentler, 2000; Savalei & Falk, 2014).
.ts_test_label <- function(lavaan_missing) {
  if (identical(lavaan_missing, "robust.two.stage")) {
    "two-stage scaled statistic (Satorra-Bentler type; Savalei & Falk, 2014)"
  } else {
    "two-stage residual based statistic (Savalei & Bentler, 2009)"
  }
}

# .ts_indices_label(lavaan_missing)
#
# What lavaan's robust RMSEA and CFI are after two-stage estimation, in the
# words a footer prints. Zhang and Savalei (2023) correct the complete data
# indices of the second stage for the missing data by subtracting from the
# naive statistic T an estimate c of its expectation under the model,
# c = tr(U~ Gamma~), with U~ the residual weight matrix of the second stage
# and Gamma~ the asymptotic covariance matrix of the stage one moments.
# Their TS V2 (Table 2, p. 268) takes Gamma~ = W~^-1 V~ W~^-1, the sandwich
# of the FIML observed (W~) and first order (V~) information of the
# saturated model, which does not assume normality, and they recommend it
# after two-stage estimation, for incomplete nonnormal data as well
# (p. 280). lavaan uses that sandwich after robust two-stage estimation.
# After normal theory two-stage estimation it uses Gamma~ = W~^-1, the form
# the sandwich takes under normality (their Appendix A, p. 282, and the
# note to their Table 1, where the same substitution turns their FIML-C
# versions 1--3 into 4--6), so the indices are the normal theory analog of
# TS V2, which Zhang and Savalei did not study. The two estimate the same
# population values and agree closely with normal data; on nonnormal data
# they can differ (the tests build both from the information matrices and
# compare them with lavaan's). Zhang and Savalei define only the RMSEA and
# the CFI; lavaan applies the same correction to the TLI (and NNFI), and
# forms the RMSEA interval by inverting T/c and rescaling its limits by c.
.ts_indices_label <- function(lavaan_missing) {
  if (identical(lavaan_missing, "robust.two.stage")) {
    "the TS V2 indices of Zhang and Savalei (2023)"
  } else {
    paste0("the normal theory analog of the TS V2 indices of Zhang and ",
           "Savalei (2023)")
  }
}

# .ts_footer(lavaan_missing, kind)
#
# The lines printed under a table estimated in two stages, stored in the
# table's `two_stage` attribute at construction and read by
# .two_stage_footer_lines(). `kind` is "cfa" (fa_confirmatory(), with its fit rows
# and information criteria), "measurement" (fa_confirmatory(output = "measurement"),
# which has no fit rows, so only the first line), or "invariance"
# (measurement_invariance(), with its difference tests). `rows` is the
# .ts_all_rows() record of the fit (summed over the rungs by
# measurement_invariance()); when the first stage read cases observed on an
# auxiliary variable alone, a last line names the N (.ts_rows_line()).
.ts_footer <- function(lavaan_missing,
                       kind = c("cfa", "measurement", "invariance"),
                       rows = NULL) {
  kind <- match.arg(kind)
  n_line <- .ts_rows_line(rows, kind, lavaan_missing)
  robust <- identical(lavaan_missing, "robust.two.stage")
  head <- paste0(
    "Two-stage maximum likelihood: the model fit to the EM estimates of the ",
    "means and covariances, with ",
    if (robust) "standard errors robust to nonnormality (Yuan & Bentler, 2000)"
    else "normal theory standard errors (Savalei & Bentler, 2009)")
  chi <- if (robust) {
    "chi_square: the two-stage scaled statistic (Savalei & Falk, 2014)"
  } else {
    "chi_square: the residual based statistic (Savalei & Bentler, 2009)"
  }
  ind <- paste0("cfi, rmsea: ", .ts_indices_label(lavaan_missing),
                if (kind == "cfa") {
                  paste0("; tli, nnfi, and the RMSEA interval: lavaan's ",
                         "extension of the same correction")
                })
  if (kind == "measurement") {
    return(c(two_stage = head,
             if (!is.null(n_line)) c(two_stage_n = n_line)))
  }
  lines <- c(two_stage = head, two_stage_test = chi,
             two_stage_indices = ind)
  if (kind == "cfa") {
    lines["two_stage_ic"] <- paste0(
      "AIC, BIC, H0, H1: the complete data normal log likelihood at the EM ",
      "estimates; AIC and BIC are the two-stage criteria of Lai (2021), ",
      "comparable only with two-stage fits of the same variables, and ",
      "-2(H0 - H1) is not a chi square statistic")
  } else {
    lines["two_stage_diff"] <- paste0(
      "delta_chi_square: the scaled difference test (Satorra & Bentler, ",
      "2001), not the difference of the chi_square entries")
  }
  if (!is.null(n_line)) lines["two_stage_n"] <- n_line
  lines
}

# .two_stage_footer_lines(x)
#
# The two-stage lines of a table's printed footer (fa_confirmatory() and
# measurement_invariance() under missing = "two_stage" or
# "robust_two_stage"), called by .dmar_tbl_footer(). NULL for any other
# table.
.two_stage_footer_lines <- function(x) {
  ts <- attr(x, "two_stage", exact = TRUE)
  if (!is.list(ts) || !length(ts$footer)) return(NULL)
  ts$footer
}
