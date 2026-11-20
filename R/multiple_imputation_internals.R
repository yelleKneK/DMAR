# =============================================================================
#  Internal helpers for multiple imputation
# =============================================================================
#
# Multiple imputation of incomplete ordered-categorical items by chained
# equations (van Buuren & Groothuis-Oudshoorn, 2011, the mice package),
# and the bootstrap then impute interval built on it (Schomaker &
# Heumann, 2018, their Method 4, "Boot MI"): B bootstrap samples are
# drawn from the incomplete rows, each is imputed M_b times, a
# replication's estimate is the mean of its M_b completed data
# estimates, and the percentile interval of those means is the
# confidence interval. reliability_omega_categorical(missing = "mi") is
# the first caller. The imputation and the bootstrap take the analysis
# as a function of one completed data set (the point function contract
# of .bootstrap_ci() in R/reliability_internals.R), so another analysis
# of ordered items can use them unchanged. Nothing here is exported, and
# mice is in Suggests, so every path that imputes checks for it first
# (.mi_require_mice()).
#
# Layout:
#   (A) Imputation               .mi_have_mice, .mi_require_mice,
#                                .mi_methods, .mi_plain_columns,
#                                .mi_redundant_aux, .mi_complete,
#                                .mi_logged_counts
#   (B) Number of imputations    .MI_MAXIT, .MI_POINT_FLOOR,
#                                .mi_imputations_floor,
#                                .mi_imputations_rule,
#                                .mi_widening_imputations
#   (C) Bootstrap then impute    .boot_mi, .mi_replicate_fn,
#                                .mi_replicate, .boot_mi_tally,
#                                .mi_variance_components,
#                                .boot_mi_interval
#   (D) Categorical omega        .omega_c_mi, .omega_c_point_fn,
#                                .omega_c_mi_point, .omega_c_mi_combine,
#                                .omega_c_mi_check
#   (E) Printed lines            .missing_treatment_footer_lines
# =============================================================================


# -----------------------------------------------------------------------------
# (A) Imputation
# -----------------------------------------------------------------------------

# .mi_have_mice() / .mi_require_mice()
#
# mice is in Suggests. .mi_require_mice() stops with an installation
# message when it is absent; the check itself is the one-line
# .mi_have_mice(), kept apart so a test can stand in for a machine
# without mice.
.mi_have_mice <- function() {
  requireNamespace("mice", quietly = TRUE)
}

.mi_require_mice <- function() {
  if (!.mi_have_mice()) {
    .dmar_refuse(c(
      "missing = \"mi\" imputes the items with the 'mice' package, ",
      "which is not installed. Install it with ",
      "install.packages(\"mice\"), or use missing = \"pairwise\", ",
      "which assumes the values are missing completely at random."),
      "package_required", package = "mice")
  }
  invisible(TRUE)
}


# .mi_methods(d, items, impute)
#
# The imputation method of every column of `d` for mice::mice(): "" for
# a column with no missing value (mice imputes nothing there and uses it
# only as a predictor), the chosen model for an incomplete item, and
# predictive mean matching for an incomplete auxiliary variable, which
# is numeric. Under impute = "polr" the items arrive as ordered factors;
# an item with two categories gets "logreg", since the proportional odds
# model of two ordered categories is the logistic regression (mice's
# own choice for a binary factor, mice::make.method()).
#
# Inputs:
#   d      : the data frame mice will receive (items and auxiliaries).
#   items  : the item column names.
#   impute : "pmm" or "polr".
# Returns:
#   A character vector named by the columns of d.
.mi_methods <- function(d, items, impute) {
  incomplete <- vapply(d, anyNA, logical(1L))
  meth <- stats::setNames(ifelse(incomplete, "pmm", ""), names(d))
  if (identical(impute, "polr")) {
    for (j in intersect(items, names(d)[incomplete])) {
      meth[[j]] <- if (nlevels(d[[j]]) > 2L) "polr" else "logreg"
    }
  }
  meth
}


# .mi_plain_columns(d, items)
#
# The columns of `d` as plain vectors, before mice sees them: each item
# as integer category codes and every other column (an auxiliary
# variable) as a double, with every class and attribute removed. A
# column read with haven (class haven_labelled, a double carrying value
# labels) passes the checks of the callers, which ask only for numbers,
# but mice's predictive mean matching cannot take it: its regression
# transposes the residuals, which inherit the class, and the vctrs
# method behind haven refuses ("`t.haven_labelled()` not supported").
# A factor would reach mice as a factor, so the other items would enter
# its regressions through K - 1 contrasts and the item itself would be
# optimally scaled (mice.impute.pmm(quantify = TRUE)), a different
# imputation model from the integer codes the help page describes; its
# labels are read as the codes, as the caller reads them before any
# treatment (.omega_c_item_codes(), which also refuses labels that are
# not whole numbers), so this is a backstop. impute = "polr" makes its
# own ordered factors from these codes (.mi_complete()).
#
# Inputs:
#   d     : data frame of the items and any auxiliary variables.
#   items : the item column names.
# Returns:
#   `d` with plain integer items and plain double auxiliaries.
.mi_plain_columns <- function(d, items) {
  for (j in names(d)) {
    x <- d[[j]]
    if (is.factor(x)) x <- as.character(x)
    x <- as.vector(unclass(x))
    d[[j]] <- if (j %in% items) as.integer(x) else as.double(x)
  }
  d
}


# .mi_redundant_aux(d, items, aux)
#
# The auxiliary variables that add nothing to the imputation model: an
# auxiliary variable whose multiple correlation with the columns
# observed on every row that observes it (the items and the auxiliary
# variables kept before it, in the order given) reaches .999, mice's own
# threshold for a collinear pair (mice's find.collinear()), extended from
# one other column to several. On the rows where it is observed it is
# then a linear function of predictors already in the model, and where
# it is missing mice would impute it from them, so it carries no
# information about the missing item values, and it is left out.
#
# Why this is needed. mice's setup check (remove.collinear = TRUE) sets
# aside a variable whose pairwise correlation with another reaches .999;
# it is turned off (.mi_complete()) because items observed together on
# two or three rows reach it by arithmetic. Without it, two auxiliary
# variables that duplicate each other (x and 2x + 1, age in years and in
# months) both reach the regression of every item, and when the
# predictors carry a second exact dependency (a second duplicated pair,
# or an auxiliary variable that is a linear function of complete items)
# mice's remove.lindep() (mice 3.19.0) removes only one: after its first
# removal it indexes the correlation matrix with keep[keep], which is
# all TRUE and recycles over the full matrix, so it reads the next
# eigenvector of the unreduced matrix. mice then falls back to its ridge
# penalty and the imputations become unstable: on the test fixture with
# x = the standardized mean of y1 to y3 and x2 = 2x + 1, categorical
# omega was 0.768 against 0.800 with x alone, and the between-imputation
# standard deviation 0.062 against 0.008. A rank two dependency with no
# duplicated pair (x and the standardized mean of y1 and y2) did the
# same before remove.collinear was turned off, and the screen covers it
# too. A screen of the auxiliary variables alone leaves the items as
# they are, so a planned missing design is still imputed in full.
#
# The check needs more rows than columns; an auxiliary variable observed
# on too few rows to tell is kept. It draws no random number, so a call
# with no redundant auxiliary variable imputes exactly as before.
#
# Inputs:
#   d     : data frame of the items and the auxiliary variables.
#   items : the item column names.
#   aux   : the auxiliary column names, in the order given.
# Returns:
#   A list with `aux`, the auxiliary variables left out (character(0)
#   when none), and `why`, for each of them the columns it is a linear
#   function of.
.mi_redundant_aux <- function(d, items, aux) {
  out <- list(aux = character(), why = list())
  if (!length(aux)) return(out)
  # Plain numbers, as mice will see them (a labeled column refuses
  # as.numeric() when vctrs is loaded).
  d <- .mi_plain_columns(as.data.frame(d)[c(items, aux)], items)
  kept <- character()
  for (a in aux) {
    rows <- which(!is.na(d[[a]]))
    cand <- c(items, kept)
    cols <- cand[vapply(cand, function(v) !anyNA(d[[v]][rows]),
                        logical(1L))]
    y <- as.numeric(d[[a]][rows])
    if (!length(cols) || length(rows) <= length(cols) + 1L) {
      kept <- c(kept, a)
      next
    }
    X <- cbind(1, vapply(cols, function(v) as.numeric(d[[v]][rows]),
                         numeric(length(rows))))
    fit <- stats::lm.fit(X, y)
    tss <- sum((y - mean(y))^2)
    r2 <- 1 - sum(fit$residuals^2) / tss
    if (is.finite(r2) && r2 >= 0.999^2) {
      beta <- fit$coefficients[-1L]
      sd_x <- apply(X[, -1L, drop = FALSE], 2L, stats::sd)
      std <- abs(beta) * sd_x / stats::sd(y)
      out$aux <- c(out$aux, a)
      out$why[[a]] <- cols[!is.na(std) & std > 1e-6]
    } else {
      kept <- c(kept, a)
    }
  }
  out
}


# .mi_redundant_aux_message(screen)
#
# The message naming the auxiliary variables .mi_redundant_aux() left
# out, each with the columns it is a linear function of. Every such
# variable carries no information beyond the predictors already in the
# model, which is the reason given. Only with two or more left out does
# the message add that mice would have met more than one exact linear
# dependency, which it does not resolve: with one, mice's remove.lindep()
# removes it and the imputations are as without it (0.8001 with a
# between-imputation standard deviation of 0.008 on the test fixture,
# with or without the variable).
.mi_redundant_aux_message <- function(screen) {
  each <- vapply(screen$aux, function(a) {
    by <- screen$why[[a]]
    paste0(a, " is a linear function of ",
           if (length(by) > 2L) {
             paste0(paste(by[-length(by)], collapse = ", "), ", and ",
                    by[length(by)])
           } else {
             paste(by, collapse = " and ")
           })
  }, character(1L))
  one <- length(each) == 1L
  paste0(if (one) "The auxiliary variable " else "Among the auxiliary variables, ",
         paste(each, collapse = "; "),
         " (a multiple correlation of .999 or more on the rows that ",
         "observe ", if (one) "it" else "each", "), so ",
         if (one) "it carries" else "these carry",
         " no information beyond the predictors already in the imputation ",
         "model and ", if (one) "is" else "are", " left out of it.",
         if (!one) {
           paste0(" Kept, they would leave mice more than one exact linear ",
                  "dependency among the predictors, which it does not ",
                  "resolve: it falls back to a ridge penalty, and its ",
                  "imputations become unstable.")
         })
}


# .mi_logged_counts(events, m)
#
# Two kinds of mice's logged events, counted by imputation (mice logs the
# imputation each happened in): a fallback of the proportional odds
# model to a multinomial one ("polr falls back"), and the ridge penalty
# mice applies when the predictors of a variable are (nearly) linearly
# dependent ("multi-collinear"), after which, in mice's own words, the
# results can be unstable.
#
# Inputs:
#   events : mice's loggedEvents data frame, or NULL.
#   m      : the number of imputations.
# Returns:
#   A list with `polr` and `ridge`, integer vectors of length m.
.mi_logged_counts <- function(events, m) {
  count <- function(pattern) {
    if (is.null(events) || is.null(events$out)) return(integer(m))
    hit <- grepl(pattern, events$out, fixed = TRUE)
    if (!any(hit)) return(integer(m))
    at <- if (is.null(events$im)) rep(1L, sum(hit)) else events$im[hit]
    tabulate(at[at >= 1L & at <= m], nbins = m)
  }
  list(polr = count("polr falls back"), ridge = count("multi-collinear"))
}


# .mi_complete(d, items, impute = "pmm", m = 5L, maxit = .MI_MAXIT)
#
# m completed copies of the items of `d`, imputed by chained equations
# (mice::mice()) with every other column of `d` as a predictor (mice's
# default predictor matrix), auxiliary variables included. The columns
# are made plain first (.mi_plain_columns()).
#
#   impute = "pmm"  : predictive mean matching. The items enter as
#                     integer codes, each imputed value is an observed
#                     value of that item drawn from the donors whose
#                     predicted means are nearest, and the regression
#                     parameters behind the predictions are drawn, so the
#                     imputations are proper and never contain a category
#                     the item does not show.
#   impute = "polr" : the proportional odds model (MASS::polr() inside
#                     mice), with the items as ordered factors, so every
#                     other item enters as a factor (all K - 1 polynomial
#                     contrasts). mice draws categories at the polr point
#                     estimates, with no draw of its parameters. When
#                     polr fails mice falls back to a multinomial model;
#                     polr.to.loggedEvents = TRUE logs each fallback, and
#                     they are counted here.
#
# maxit is the number of chained equations iterations (mice's default is
# 5; 10 is the value of the Monte Carlo behind the default model, see
# ?reliability_omega_categorical). mice's own seed argument is left at
# NA, so the draws come from the R session's stream, which the caller
# has set (.dmar_local_seed()) or which a bootstrap replication has set
# for itself (.mi_replicate_fn()).
#
# mice sets aside, before it imputes, a variable that is constant or
# whose pairwise correlation with another reaches .999 in absolute value
# (mice's find.collinear(), with use = "pairwise.complete.obs"): it stops
# using it as a predictor and, when it is incomplete, stops imputing it,
# so its missing values survive into every completed data set. For items
# observed together on only two or three rows, as in a planned missing
# design with a small overlap, that correlation is one by arithmetic,
# not by redundancy; left to mice's default, such an item went
# unimputed, and lavaan then fitted each completed data set to the rows
# that observe it, by listwise deletion and without a word, while the
# result reported every row. So remove.collinear = FALSE is passed and
# every incomplete item is imputed. The check also set aside an
# auxiliary variable that duplicates another, which mice's own
# remove.lindep() cannot always do in its place (it removes one exact
# dependency of several); the callers screen the auxiliary variables
# first (.mi_redundant_aux()). When the predictors of an item are still
# (nearly) linearly dependent, mice applies its ridge penalty (1e-5)
# and logs it; those events are counted (.mi_logged_counts()) and
# reported, since the imputations can then be unstable. A constant item
# is still set aside by mice; the callers
# refuse an item with fewer than two observed values before imputing,
# and a completed data set that still has a missing item value stops
# this function, naming the item and the reason mice logged, rather
# than reaching a fit that would silently drop the incomplete rows.
#
# Two kinds of warning are muffled, and only these: mice's notice that
# it logged events (the polr fallbacks and ridge penalties among them are
# counted and reported by the caller, an item left unimputed stops the
# call as just described, and the rest, a predictor set aside as
# constant or as linearly dependent on the others, are mice's routine
# bookkeeping), and R's partial matching diagnostics ("partial
# match of 'coef' to 'coefficients'"), which mice's own code raises once
# per imputed variable and iteration when a session turns on the
# warnPartialMatch options, as this package's tests do; they concern
# mice's source, not the imputation.
#
# Inputs:
#   d      : data frame of the items (integer category codes) and any
#            auxiliary variables, at least one item value missing.
#   items  : the item column names.
#   impute : "pmm" or "polr".
#   m      : number of imputations.
#   maxit  : number of chained equations iterations.
# Returns:
#   A list with `data`, a list of m data frames of the items as integer
#   category codes (the codes of `d`); `polr_fallbacks`, the number of
#   fallbacks from polr to a multinomial model mice logged;
#   `polr_fallbacks_each`, the same count for each of the m imputations
#   (mice logs the imputation a fallback happened in); and `ridge` and
#   `ridge_each`, the same for mice's ridge penalty.
.mi_complete <- function(d, items, impute = "pmm", m = 5L,
                         maxit = .MI_MAXIT) {
  d <- .mi_plain_columns(as.data.frame(d), items)
  if (identical(impute, "polr")) {
    for (j in items) d[[j]] <- factor(d[[j]], ordered = TRUE)
  }
  meth <- .mi_methods(d, items, impute)
  muffle_logged <- function(w) {
    msg <- conditionMessage(w)
    if (startsWith(msg, "Number of logged events") ||
        startsWith(msg, "partial match of") ||
        startsWith(msg, "partial argument match of")) {
      invokeRestart("muffleWarning")
    }
  }
  imp <- withCallingHandlers(
    if (identical(impute, "polr")) {
      mice::mice(d, m = m, method = meth, maxit = maxit, printFlag = FALSE,
                 remove.collinear = FALSE, polr.to.loggedEvents = TRUE)
    } else {
      mice::mice(d, m = m, method = meth, maxit = maxit, printFlag = FALSE,
                 remove.collinear = FALSE)
    },
    warning = muffle_logged)
  events <- imp$loggedEvents
  logged <- .mi_logged_counts(events, m)
  fallbacks_each <- logged$polr
  completed <- withCallingHandlers(lapply(seq_len(m), function(k) {
    dk <- mice::complete(imp, k)[items]
    dk[] <- lapply(dk, function(x) as.integer(as.character(x)))
    dk
  }), warning = muffle_logged)
  left <- items[vapply(items, function(j) {
    any(vapply(completed, function(dk) anyNA(dk[[j]]), logical(1L)))
  }, logical(1L))]
  if (length(left)) {
    why <- vapply(left, function(j) {
      hit <- if (is.null(events) || is.null(events$out)) {
        character()
      } else {
        as.character(events$meth[events$out == j])
      }
      if (length(hit)) paste0(j, " (set aside as ", hit[1L], ")") else j
    }, character(1L))
    .dmar_refuse_argument(c(
      "mice did not impute ", paste(why, collapse = " or "), ", so ",
      "the completed data sets would still have missing values."), arg = NULL)
  }
  list(data = completed, polr_fallbacks = as.integer(sum(fallbacks_each)),
       polr_fallbacks_each = as.integer(fallbacks_each),
       ridge = as.integer(sum(logged$ridge)),
       ridge_each = as.integer(logged$ridge))
}


# -----------------------------------------------------------------------------
# (B) Number of imputations
# -----------------------------------------------------------------------------

# The chained equations iterations behind every imputation (documented on
# the help page, not exposed), and the fewest imputations of the observed
# data behind a point estimate. Twenty keeps Rubin's relative efficiency,
# (1 + FMI / M)^-1, at or above .976 for a fraction of missing
# information up to .5.
.MI_MAXIT <- 10L
.MI_POINT_FLOOR <- 20L


# .mi_imputations_floor(imputations)
#
# The number of imputations of the observed data made before anything
# else: max(imputations, 20). Without a bootstrap this is the point
# estimate's number of imputations; with one, von Hippel's rule may add
# more after the bootstrap (.omega_c_mi()).
.mi_imputations_floor <- function(imputations) {
  max(as.integer(imputations), .MI_POINT_FLOOR)
}


# .mi_imputations_rule(fmi_upper, cv = 0.05)
#
# von Hippel's (2020) quadratic rule, his Eq. 1, M = 1 + (1/2)(FMI/CV)^2,
# evaluated at the upper limit of a 95% confidence interval for the
# fraction of missing information, as his two stage procedure does (with
# the point estimate there is about an even chance of too few
# imputations). His target is the replicability of a standard error, CV
# its coefficient of variation over sets of imputations; CV = .05 is the
# criterion his Figure 1 fits to Bodner's (2008) results, under which the
# rule is close to M = 200 FMI^2. Under bootstrap then impute there is no
# Rubin's rules standard error, so the rule is borrowed for the number
# of imputations behind the point estimate, with the bootstrap as the
# pilot; the help page says so.
#
# Returns:
#   An integer, or NA when fmi_upper is not finite.
.mi_imputations_rule <- function(fmi_upper, cv = 0.05) {
  if (length(fmi_upper) != 1L || !is.finite(fmi_upper)) return(NA_integer_)
  as.integer(ceiling(1 + 0.5 * (fmi_upper / cv)^2))
}


# .mi_widening_imputations(fmi, limit = 1.05, kept = 1, even = FALSE)
#
# The fewest imputations per bootstrap replication that keep the
# widening of the percentile interval by imputation noise,
# sqrt(1 + FMI / M_b), at or below `limit` (see .omega_c_mi()). When a
# replication averages only the completed data sets that return a value,
# a share `kept` of them on the harmonic average (M_b mean(1 / n_b) is
# 1 / kept), the widening is sqrt(1 + FMI / (kept M_b)), and the number
# is raised in proportion; the share is taken to stay as it was. With
# `even` (some completed data sets returned no value) an odd number is
# raised by one: a replication is set aside when fewer than half of its
# sets return a value, and an odd number of imputations sets aside more
# replications than the even number above it, which also widens the
# interval less (.mi_replicate()).
.mi_widening_imputations <- function(fmi, limit = 1.05, kept = 1,
                                     even = FALSE) {
  n <- max(2L, as.integer(ceiling(fmi / kept / (limit^2 - 1) - 1e-9)))
  if (even && n %% 2L == 1L) n <- n + 1L
  n
}


# -----------------------------------------------------------------------------
# (C) Bootstrap then impute
# -----------------------------------------------------------------------------

# .boot_mi(data, items, estimate_fn, impute, m, maxit, B, ncpus,
#          progress = NULL)
#
# The B replications of the bootstrap then impute interval. Each draws n
# rows of `data` with replacement (the incomplete rows, auxiliaries
# included), imputes them m times (.mi_complete()), and applies
# estimate_fn to each completed data set.
#
# Random numbers. The replications draw random numbers of their own (the
# imputations), unlike the replications of .bootstrap_ci(), whose point
# functions are deterministic given the resample. So that a seed gives
# the same interval whatever ncpus is, every draw that decides a
# replication is made here, in this process, before any replication
# runs: the B resamples (B x n row indices) and one integer seed per
# replication. A replication then imputes under its own seed, set and
# restored around it by withr::with_seed(), with the generator kinds of
# this session, so forked workers, socket workers, and this process
# impute the same values. (boot::boot() could not be used: its parallel
# replications draw from whatever stream each worker inherits, which a
# fork reseeds unless the generator is L'Ecuyer's, and then by worker,
# not by replication.) The caller owns the seed (.dmar_local_seed()).
#
# Parallel replications follow the engine's rule (.boot_parallel_type()):
# forked workers where the platform has fork, a socket cluster on
# Windows, whose workers load the installed DMAR.
#
# A progress callback (see .bootstrap_ci()) is called after each
# replication, or each block of them on several processes
# (.boot_lapply()), and a FALSE stops the run with an error of class
# dmar_bootstrap_cancelled. Every draw is made before the first
# replication and each replication imputes under its own seed, so the
# callback cannot change a result. On a socket cluster the function of b,
# which holds the data, the resamples, and the seeds, is sent to each
# worker once, and each block sends only its replication numbers
# (.boot_worker_stage()).
#
# Inputs:
#   data        : data frame of the items and any auxiliaries, the rows
#                 with at least one observed item.
#   items       : the item column names.
#   estimate_fn : function of one completed data frame of the items,
#                 returning a number under the point function contract
#                 of .bootstrap_ci().
#   impute, m, maxit : as for .mi_complete().
#   B           : number of bootstrap replications.
#   ncpus       : number of processes.
#   progress    : NULL, or a function(done, total), as above.
# Returns:
#   A list with `t`, a B x (3 + 2m) matrix whose first three columns are
#   those of the engine's replicate matrix (the value the interval uses,
#   the status code, and the estimate; .boot_classify()), whose next m
#   columns are the replication's m completed data estimates (NA where
#   one returned no value), and whose last m are the status codes of
#   those completed data sets (.mi_replicate()); `messages`, the error
#   message of each replication that has one (NA otherwise);
#   `polr_fallbacks` and `ridge`, the totals of the events logged
#   (.mi_complete()); and `ncpus`, the number of processes used. Status
#   code 7, "imputation", marks a replication whose imputation itself
#   failed (codes 1 to 6 are the engine's, .BOOT_STATUS in
#   R/reliability_internals.R, whose sixth, "not_converged", this loop
#   never assigns).
.boot_mi <- function(data, items, estimate_fn, impute, m, maxit, B, ncpus,
                     progress = NULL) {
  n <- nrow(data)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- .mi_replicate_fn(data, items, estimate_fn, impute, m, maxit,
                          idx, seeds, RNGkind())
  par_type <- .boot_parallel_type(ncpus)
  cl <- NULL
  if (par_type == "snow") {
    cl <- parallel::makePSOCKcluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  res <- .boot_lapply(seq_len(B), one, par_type, ncpus, cl,
                      progress = progress)
  width <- 3L + 2L * m
  # A worker that died returns no list; count it as an error.
  res <- lapply(res, function(r) {
    if (is.list(r) && is.numeric(r$row) && length(r$row) == width) {
      r
    } else {
      list(row = c(NA_real_, 3, NA_real_, rep(NA_real_, m), rep(3, m)),
           message = if (inherits(r, "try-error")) {
             trimws(as.character(r))
           } else {
             "the worker process returned no result"
           },
           fallbacks = 0L, ridge = 0L)
    }
  })
  list(t = matrix(unlist(lapply(res, `[[`, "row")), ncol = width,
                  byrow = TRUE),
       messages = vapply(res, function(r) {
         if (is.character(r$message) && length(r$message) == 1L) {
           r$message
         } else {
           NA_character_
         }
       }, character(1L)),
       polr_fallbacks = sum(vapply(res, function(r) {
         as.integer(r$fallbacks)
       }, integer(1L))),
       ridge = sum(vapply(res, function(r) {
         as.integer(r$ridge %||% 0L)
       }, integer(1L))),
       ncpus = if (par_type == "no") 1L else ncpus)
}


# .mi_replicate_fn(data, items, estimate_fn, impute, m, maxit, idx,
#                  seeds, rng_kind)
#
# The function of b that .boot_mi() hands its workers: replication b on
# the rows idx[b, ], imputed under seeds[b] with the generator kinds
# rng_kind (RNGkind() of the calling session). Built in its own small
# closure with its arguments forced, like .boot_statistic(), so a socket
# worker receives values rather than promises it cannot evaluate.
.mi_replicate_fn <- function(data, items, estimate_fn, impute, m, maxit,
                             idx, seeds, rng_kind) {
  force(data)
  force(items)
  force(estimate_fn)
  force(impute)
  force(m)
  force(maxit)
  force(idx)
  force(seeds)
  force(rng_kind)
  function(b) {
    withr::with_seed(
      seeds[b],
      .mi_replicate(data[idx[b, ], , drop = FALSE], items, estimate_fn,
                    impute, m, maxit),
      .rng_kind = rng_kind[1L], .rng_normal_kind = rng_kind[2L],
      .rng_sample_kind = rng_kind[3L])
  }
}


# .mi_replicate(d, items, estimate_fn, impute, m, maxit)
#
# One replication: impute the resample m times and apply estimate_fn to
# each completed data set, classifying each result as the engine does
# (.boot_classify(), with improper estimates set aside). The replication
# computes the statistic the point estimate is: the mean of the
# completed data sets that return a value, defined when at least half of
# them do (.omega_c_mi_check(), Lim and Cheung's, 2022, criterion). A
# bootstrap replicates the estimator, so the replication follows the
# point estimate's rule, set by set. Through the first versions of 1.1.0
# a replication returned a value only when all m of its completed data
# sets did, a different statistic: with a share p of sets returning
# none, a replication was set aside with probability 1 - (1 - p)^m,
# which grows with m, so raising `imputations` set more replications
# aside, and their completed data sets that did return a value were
# discarded. On the help page's example (B = 200, the model implied
# denominator), that rule applied to the first 1 to 5 imputations of
# each replication sets aside 24, 35.5, 45.5, 48.5, and 52.5 percent of
# the replications.
#
# Under the point estimate's rule a replication is set aside with the
# probability that fewer than half of its m sets return a value,
# P(Bin(m, 1 - p) < m / 2) for independent sets, which goes to zero as
# m grows when p is below one half, but not steadily. An odd m needs
# (m + 1) / 2 sets to return a value, one more than the even m - 1 needs
# of one set fewer, so the share rises from each even m to the next odd
# one, and each even m sets aside fewer than either odd neighbor: at
# p = .2, 4.0, 10.4, 2.7, 5.8, and 1.7 percent at m = 2 to 6. On the help
# page's example (a quarter of the sets returned no value) the rule sets
# aside 24, 12.5, 21, 14.5, and 20 percent at m = 1 to 5, more than
# binomial because the sets of one resample fail together. So the
# number of imputations the warning names is even whenever sets were
# set aside (.mi_widening_imputations()), apart from replications whose
# imputation itself failed (.boot_mi_notes()). A replication set aside takes
# the status of its first completed data set that returned no value,
# or 7 ("imputation") when the imputation itself raised an error. The
# replicate means then average unequal numbers of imputations, which
# .mi_variance_components() allows for.
#
# Returns:
#   A list with `row` (c(value, status, estimate, the m estimates with NA
#   where a set returned no value, the m status codes of the sets)),
#   `message` (the error message behind status 3 or 7, NA otherwise),
#   `fallbacks` (polr fallbacks logged), and `ridge` (the times mice
#   applied its ridge penalty; .mi_complete()).
.mi_replicate <- function(d, items, estimate_fn, impute, m, maxit) {
  comp <- tryCatch(.mi_complete(d, items, impute, m, maxit),
                   error = function(e) e)
  if (inherits(comp, "error")) {
    return(list(row = c(NA_real_, 7, NA_real_, rep(NA_real_, m),
                        rep(7, m)),
                message = conditionMessage(comp), fallbacks = 0L,
                ridge = 0L))
  }
  vals <- lapply(comp$data, function(dk) {
    tryCatch(estimate_fn(dk), error = function(e) e)
  })
  cls <- vapply(vals, .boot_classify, numeric(3L), improper = "drop")
  code <- cls[2L, ]
  ok <- code == 1
  est <- ifelse(ok, cls[3L, ], NA_real_)
  msg <- NA_character_
  if (sum(ok) >= m / 2) {
    value <- mean(est[ok])
    row <- c(value, 1, value, est, code)
  } else {
    first <- which(!ok)[1L]
    row <- c(NA_real_, code[first], NA_real_, est, code)
    if (code[first] == 3) msg <- conditionMessage(vals[[first]])
  }
  list(row = row, message = msg, fallbacks = comp$polr_fallbacks,
       ridge = comp$ridge %||% 0L)
}


# .boot_mi_tally(boot, B)
#
# The counts of the bootstrap attribute for a bootstrap then impute run,
# in the form .bootstrap_ci() returns them (.boot_tally()), with a third
# reason, `imputation`, for the replications whose imputation failed
# (taken out of `non_finite`, where .boot_tally() would count them), the
# distinct messages behind those failures in `imputation_messages`
# (beside `error_messages`), and `status` with a seventh level,
# "imputation", after the engine's six (.BOOT_STATUS). The engine's
# .boot_tally() would count code 7 nowhere but in non_finite, so its
# not_converged count stays that of code 6, which is zero here.
.boot_mi_tally <- function(boot, B) {
  tt <- boot$t
  code <- as.integer(tt[, 2L])
  counts <- .boot_tally(tt[, 1:3, drop = FALSE], B)
  n_imp <- sum(code == 7L)
  counts$reasons <- c(error = counts$reasons[["error"]],
                      non_finite = counts$reasons[["non_finite"]] - n_imp,
                      imputation = n_imp)
  tab <- function(sel) {
    msgs <- boot$messages[sel & !is.na(boot$messages)]
    if (!length(msgs)) return(integer())
    tb <- table(msgs)
    stats::setNames(as.integer(tb), names(tb))
  }
  counts$error_messages <- tab(code == 3L)
  counts <- append(counts, list(imputation_messages = tab(code == 7L)),
                   after = match("error_messages", names(counts)))
  counts
}


# .mi_variance_components(Tm, level = 0.95)
#
# The one-way decomposition of the replicate by imputation matrix of a
# bootstrap then impute run (rows the replications that returned a
# value, columns their M_b imputations, NA where a completed data set
# returned no value and was set aside). Each replication's estimates
# share a bootstrap sample and differ only by imputation, so the one-way
# random effects model applies. When every replication averages all M_b
# of its imputations the design is balanced, with
#
#   MSB = M_b var(row means),  MSW = mean within-row variance,
#   E(MSW) = s2_imp,           E(MSB) = s2_imp + M_b s2_boot,
#
# where s2_imp is the between-imputation variance and s2_boot the
# sampling variance of the estimator with infinitely many imputations,
# the total variance of Rubin's rules (within plus between). Their ratio
# estimates the fraction of missing information,
#   FMI = s2_imp / s2_boot,
# bounded at 1. Its upper limit comes from the usual interval for the
# ratio of variance components, rho = s2_boot / s2_imp, whose lower
# limit at the 1 - (1 - level)/2 quantile of F on (B - 1, B (M_b - 1))
# degrees of freedom is (F_obs / F_crit - 1) / M_b, so FMI_upper =
# 1 / rho_lower (1 when rho_lower is not positive). The level is fixed at
# 95%, von Hippel's (2020), not the level of the reported interval. This
# is a construction from the variance components model, not a published
# procedure for omega; two cautions go on the help page: with polr's
# imputations, which ignore the uncertainty of the imputation model's
# parameters, s2_imp understates the proper between-imputation variance,
# so the FMI is too low, and the FMI is imprecise when B is small.
#
# The variance of a replicate mean is s2_boot + s2_imp / M_b, so imputation
# noise widens the percentile interval by sqrt(1 + FMI / M_b) relative to
# the interval with infinitely many imputations per replication.
#
# A replication that sets some of its completed data sets aside averages
# n_b < M_b of them (.mi_replicate()), and the design is unbalanced. The
# analysis of variance estimators then use the sums of squares about the
# replication means and about the grand mean weighted by n_b, with
# N = sum(n_b),
#
#   MSW = SSW / (N - B),       MSB = SSB / (B - 1),
#   E(MSW) = s2_imp,           E(MSB) = s2_imp + n0 s2_boot,
#   n0 = (N - sum(n_b^2) / N) / (B - 1)
#
# (Searle, Casella, & McCulloch, 1992, Sec. 3.6), and the interval for
# the ratio takes n0 in place of M_b with F on (B - 1, N - B) degrees of
# freedom, which is then approximate (MSB is no longer a multiple of a
# chi square variable; Thomas & Hultquist, 1978); with n_b = M_b for
# every replication, N = B M_b and n0 = M_b, and each formula is the
# balanced one above, which is computed as before. A replicate mean of
# n_b imputations has variance s2_boot + s2_imp / n_b, so the widening
# is sqrt(1 + FMI mean(1 / n_b)), sqrt(1 + FMI / M_b) when balanced;
# `inv_n` returns mean(1 / n_b).
#
# Inputs:
#   Tm    : numeric matrix, at least two rows and two columns, NA where
#           a replication's completed data set was set aside.
#   level : confidence level of the interval for FMI.
# Returns:
#   A list with s2_boot, s2_imp, fmi, fmi_upper, sd_replicates (the
#   standard deviation of the replicate means), and inv_n (the mean of
#   1 / n_b), each NA when fewer than two rows have a value or no row has
#   two.
.mi_variance_components <- function(Tm, level = 0.95) {
  Tm <- as.matrix(Tm)
  out <- list(s2_boot = NA_real_, s2_imp = NA_real_, fmi = NA_real_,
              fmi_upper = NA_real_, sd_replicates = NA_real_,
              inv_n = NA_real_)
  if (ncol(Tm) < 2L) return(out)
  n_b <- rowSums(is.finite(Tm))
  Tm <- Tm[n_b > 0L, , drop = FALSE]
  n_b <- n_b[n_b > 0L]
  B <- nrow(Tm)
  m <- ncol(Tm)
  N <- sum(n_b)
  if (B < 2L || N - B < 1L) return(out)
  if (all(n_b == m)) {
    means <- rowMeans(Tm)
    msb <- m * stats::var(means)
    msw <- sum((Tm - means)^2) / (B * (m - 1))
    n0 <- m
    df_within <- B * (m - 1)
  } else {
    means <- rowSums(Tm, na.rm = TRUE) / n_b
    grand <- sum(n_b * means) / N
    msb <- sum(n_b * (means - grand)^2) / (B - 1)
    msw <- sum((Tm - means)^2, na.rm = TRUE) / (N - B)
    n0 <- (N - sum(n_b^2) / N) / (B - 1)
    df_within <- N - B
  }
  s2_boot <- max((msb - msw) / n0, 0)
  if (msw <= 1e-20 * max(Tm^2, na.rm = TRUE)) {
    # No imputation noise beyond the rounding of the row means: every
    # imputation of a replication gave the same estimate, so nothing
    # missing mattered to it.
    fmi <- 0
    fmi_upper <- 0
  } else {
    fmi <- if (s2_boot > 0) min(msw / s2_boot, 1) else 1
    f_crit <- stats::qf(1 - (1 - level) / 2, B - 1, df_within)
    rho_lower <- (msb / msw / f_crit - 1) / n0
    fmi_upper <- if (rho_lower > 0) min(1 / rho_lower, 1) else 1
  }
  list(s2_boot = s2_boot, s2_imp = msw, fmi = fmi, fmi_upper = fmi_upper,
       sd_replicates = stats::sd(means), inv_n = mean(1 / n_b))
}


# .boot_mi_interval(t_means, estimate, se, conf_level, kind, B)
#
# The interval of a bootstrap then impute run. "percentile" is the Boot MI
# interval of Schomaker and Heumann (2018, Eq. 3.5), formed as the
# engine forms a percentile interval, by boot::boot.ci() (with t and t0
# supplied, it needs no boot object beyond the number of replications),
# so the order statistics and their interpolation are boot's. The Wald
# variants are centered at the point estimate, the mean of the
# imputations of the observed data, with its standard error `se`.
#
# Returns:
#   c(lower, upper); NA when fewer than two replications returned a
#   value or every replicate mean is the same.
.boot_mi_interval <- function(t_means, estimate, se, conf_level, kind, B) {
  if (sum(is.finite(t_means)) < 2L || !is.finite(estimate)) {
    return(c(NA_real_, NA_real_))
  }
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  switch(
    kind,
    percentile = {
      bci <- boot::boot.ci(list(R = B), conf = conf_level, type = "perc",
                           t0 = estimate, t = t_means)
      if (is.null(bci)) c(NA_real_, NA_real_) else bci$percent[4:5]
    },
    bootstrap_se = estimate + c(-1, 1) * crit * se,
    bootstrap_se_logistic = {
      lg <- .logistic_ci(estimate, se, crit)
      c(lg$lower, lg$upper)
    })
}


# -----------------------------------------------------------------------------
# (D) Categorical omega under multiple imputation
# -----------------------------------------------------------------------------

# .omega_c_point_fn(denominator)
#
# The point function a bootstrap replication applies to a completed data
# set, built in its own small closure so a socket worker receives the
# denominator and nothing of the caller's frame.
.omega_c_point_fn <- function(denominator) {
  force(denominator)
  # Completed resamples obey the same column checks as other
  # reliability resamples, before lavaan is called.
  .relia_checked_point(function(d) {
    .omega_c_from_data(d, denominator = denominator)
  })
}


# .omega_c_mi(data, items, impute, imputations, denominator, kind, B,
#             conf_level, seed, ncpus, progress = NULL)
#
# Categorical omega under missing = "mi": the point estimate from the
# imputations of the observed data and, when kind is not "none", the
# bootstrap then impute interval. This function owns the seed: every
# random draw (the imputations of the observed data, the resamples, and
# the replications' seeds) is made between its .dmar_local_seed() and
# its exit, in this order:
#
#   1. max(imputations, 20) imputations of the observed data
#      (.mi_imputations_floor()), so that without a bootstrap the point
#      estimate is the mean of those, and a call with a bootstrap and the
#      same seed starts from the same imputations;
#   2. with a bootstrap, the reserve: the further imputations of the
#      observed data that von Hippel's (2020) rule could ask for, up to
#      the number it asks for at a fraction of missing information of one
#      (201), from one further mice run (imputations from separate runs
#      are as valid as those of one run, each run's chains starting
#      afresh);
#   3. the B replications of the bootstrap, each with `imputations`
#      imputations (.boot_mi()).
#
# The rule, applied to the upper limit of the fraction of missing
# information the replications estimate, then decides how many of the
# imputations made in steps 1 and 2 the point estimate uses: the first
# max(floor, rule), in the order made. Which ones are used depends only
# on that number, never on their estimates, and the reserve beyond it is
# left unused. Every completed data set of the observed data is fitted
# and checked before the bootstrap starts, so nothing found in one can
# discard a bootstrap already run; the reserve costs at most 181
# imputations and fits, about 4% of the 5000 of the default bootstrap.
#
# A completed data set of the observed data that returns no categorical
# omega (an improper solution above all, see .omega_c_mi_check()) is set
# aside and counted, as a bootstrap replication is, and the point
# estimate is the mean of those that returned a value; the call stops,
# before the bootstrap, only when fewer than half of the completed data
# sets made return one. The standard error of the point estimate is
# sqrt(s2_boot + s2_imp / M) (.mi_variance_components()), with M the
# completed data sets it averages, so it includes the imputation noise of
# its own imputations; the Wald variants use it, and the percentile
# interval is that of the replicate means. The residual polychoric
# correlations are the mean, over the completed data sets the estimate
# averages, of each set's residuals (Rubin's point estimate), NA for a
# pair that one of them leaves undetermined.
#
# Returns:
#   A list with estimate, se, lower, upper, ci_method (the label the
#   result carries), bootstrap (the counts, or NULL), mi (the attribute
#   `imputation`; see ?reliability_omega_categorical), residual_cor
#   (.residual_cor_summary()), undetermined (the pairs a completed data
#   set the estimate averages leaves undetermined), and notes (the
#   sentences of the one warning the caller gives, empty when there is
#   nothing to say).
.omega_c_mi <- function(data, items, impute, imputations, denominator,
                        kind, B, conf_level, seed, ncpus, progress = NULL) {
  .mi_require_mice()
  if (kind != "none" && !requireNamespace("boot", quietly = TRUE)) {
    .dmar_refuse(c(
      "Bootstrap CIs require the 'boot' package. ",
      "Install with install.packages(\"boot\")."),
      "package_required", package = "boot")
  }
  # An auxiliary variable that is a linear function of columns already
  # among the predictors adds nothing and can leave mice with more exact
  # dependencies than it removes (.mi_redundant_aux()); it is left out of
  # the imputation of the observed data and of every replication.
  screen <- .mi_redundant_aux(data, items, setdiff(names(data), items))
  if (length(screen$aux)) {
    message(.mi_redundant_aux_message(screen))
    data <- data[setdiff(names(data), screen$aux)]
  }
  .dmar_local_seed(seed)
  maxit <- .MI_MAXIT
  m_first <- .mi_imputations_floor(imputations)
  made <- .omega_c_mi_point(data, items, impute, m_first, maxit,
                            denominator)
  .omega_c_mi_check(made)
  if (kind != "none") {
    m_most <- .mi_imputations_rule(1)
    if (!is.na(m_most) && m_most > m_first) {
      made <- .omega_c_mi_combine(
        made, .omega_c_mi_point(data, items, impute, m_most - m_first,
                                maxit, denominator))
      .omega_c_mi_check(made)
    }
  }

  vc <- .mi_variance_components(matrix(numeric(0), 0L, 0L))
  rule <- NA_integer_
  counts <- NULL
  boot <- NULL
  sets_aside <- NULL
  if (kind != "none") {
    boot <- .boot_mi(data, items, .omega_c_point_fn(denominator), impute,
                     imputations, maxit, B, ncpus, progress = progress)
    counts <- .boot_mi_tally(boot, B)
    used <- is.finite(boot$t[, 1L])
    # A replication averages the completed data sets that returned a
    # value (.mi_replicate()), so the matrix has NA where one did not,
    # and the variance components allow for the unequal numbers.
    vc <- .mi_variance_components(
      boot$t[used, 3L + seq_len(imputations), drop = FALSE])
    rule <- .mi_imputations_rule(vc$fmi_upper)
    sets_aside <- .boot_mi_sets_set_aside(
      boot$t[used, 3L + imputations + seq_len(imputations), drop = FALSE])
  }
  m_point <- m_first
  if (!is.na(rule) && rule > m_first) {
    m_point <- min(as.integer(rule), length(made$status))
  }
  point <- .omega_c_mi_take(made, seq_len(m_point))
  ok <- point$status == "proper"
  m_ok <- sum(ok)

  estimate <- mean(point$estimates[ok])
  se <- if (is.finite(vc$s2_boot)) {
    sqrt(vc$s2_boot + vc$s2_imp / m_ok)
  } else {
    NA_real_
  }
  limits <- c(NA_real_, NA_real_)
  label <- kind
  notes <- .omega_c_mi_point_notes(point)
  # sqrt(1 + FMI mean(1 / n_b)), which is sqrt(1 + FMI / M_b) when every
  # replication averages all M_b of its imputations.
  widening <- if (is.finite(vc$fmi)) {
    sqrt(1 + vc$fmi * vc$inv_n)
  } else {
    NA_real_
  }
  if (kind != "none") {
    # Fewer than half of the replications returning a value leaves no
    # interval and no standard error, as for the bootstrap of the other
    # treatments (.bootstrap_ci()): the rest are a
    # minority of the resamples, chosen by whether they returned a value.
    short_of_half <- counts$returned >= 2L && counts$returned < B / 2
    # The count rule of ?dmar_bootstrap, as in .bootstrap_ci(): it is
    # held to the number of replications requested, B
    # (.dmar_boot_count_met()), not to those that returned a value, so
    # it withholds only for a caller that did not hold B to it.
    boot_kind <- .boot_method_kind(kind)
    required <- .dmar_boot_required(boot_kind, tail = (1 - conf_level) / 2)
    short_of_count <- counts$returned >= 2L && !short_of_half &&
      !.dmar_boot_count_met(B, required)
    counts$required <- required
    counts$B_used <- as.integer(counts$returned)
    counts$bootstrap_status <- if (short_of_half) "insufficient_share" else
      if (short_of_count || counts$returned < 2L) "insufficient_count" else
        if (!is.na(required) && B < required) "exploratory" else
          if (counts$failed > 0L) "dropped" else "complete"
    if (short_of_half || short_of_count) {
      se <- NA_real_
    } else {
      limits <- .boot_mi_interval(boot$t[, 1L], estimate, se, conf_level,
                                  kind, B)
    }
    label <- paste0(kind, " (bootstrap then impute)")
    counts$improper_rule <- "drop"
    counts$caution <- .boot_caution(counts$failed, B, conf_level)
    counts$require_half <- TRUE
    counts$acceleration <- NA_real_
    counts$acceleration_method <- NA_character_
    counts$jackknife <- .boot_tally(matrix(numeric(0), 0L, 3L), 0L)
    counts$ncpus <- boot$ncpus
    counts$t <- as.numeric(boot$t[, 3L])
    levels_status <- c(.BOOT_STATUS, "imputation")
    counts$status <- factor(levels_status[as.integer(boot$t[, 2L])],
                            levels = levels_status)
    notes <- c(notes, .boot_mi_notes(
      counts, B, conf_level, imputations, vc$fmi, widening, denominator,
      sets_aside = sets_aside,
      kept = if (is.finite(vc$inv_n)) 1 / (imputations * vc$inv_n) else 1))
  }

  replicate_estimates <- if (is.null(boot)) {
    NULL
  } else {
    boot$t[, 3L + seq_len(imputations), drop = FALSE]
  }
  mi <- list(
    impute = impute,
    imputations_point = m_point,
    imputations_returned = m_ok,
    imputations_set_aside = .omega_c_mi_tally(point$status),
    imputations_replicate = as.integer(imputations),
    imputations_floor = m_first,
    imputations_rule = rule,
    maxit = maxit,
    estimates = point$estimates,
    status = point$status,
    between_sd = stats::sd(point$estimates[ok]),
    sd_replicates = vc$sd_replicates,
    s2_boot = vc$s2_boot,
    s2_imp = vc$s2_imp,
    fmi = vc$fmi,
    fmi_ci_upper = vc$fmi_upper,
    widening = widening,
    replicate_estimates = replicate_estimates,
    replicate_sets_set_aside = sets_aside,
    polr_fallbacks = c(
      point = sum(point$polr_fallbacks),
      replicates = if (is.null(boot)) 0L else boot$polr_fallbacks),
    ridge = c(point = sum(point$ridge),
              replicates = if (is.null(boot)) 0L else boot$ridge),
    aux_set_aside = screen$aux)
  ridge <- mi$ridge
  if (any(ridge > 0L)) {
    notes <- c(notes, sprintf(paste0(
      "mice applied its ridge penalty %s, which it does when the ",
      "predictors of an item are (nearly) linearly dependent, and the ",
      "imputations can then be unstable; items that duplicate each other ",
      "are the usual cause"),
      paste(c(if (ridge[["point"]] > 0L) {
        sprintf("%d %s in imputing the observed data", ridge[["point"]],
                if (ridge[["point"]] == 1L) "time" else "times")
      }, if (ridge[["replicates"]] > 0L) {
        sprintf("%d %s in imputing the bootstrap replications",
                ridge[["replicates"]],
                if (ridge[["replicates"]] == 1L) "time" else "times")
      }), collapse = " and ")))
  }

  resid <- Reduce(`+`, point$resid[ok]) / m_ok
  list(estimate = estimate, se = se, lower = limits[1L], upper = limits[2L],
       ci_method = label, bootstrap = counts, mi = mi,
       residual_cor = .residual_cor_summary(resid, df = point$df),
       undetermined = .omega_c_mi_pairs(point$undetermined[ok]),
       notes = notes)
}


# .boot_mi_notes(counts, B, conf_level, imputations, fmi, widening,
#                denominator = "observed", sets_aside = NULL, kept = 1)
#
# The sentences of the one warning a bootstrap then impute run gives:
# the replications set aside and why, the completed data sets set aside
# within the replications that returned a value (`sets_aside`, from
# .boot_mi_sets_set_aside()), what the undetermined and improper ones do
# to the interval (the undetermined sentence depends on the
# denominator), the caution when those set aside reach the tail area
# beyond a limit (.boot_caution()), and the widening of the interval by
# imputation noise when it exceeds 5% (with the number of imputations
# per replication that would bring it under, given the share `kept` of
# completed data sets that return a value, and even when any were set
# aside; .mi_widening_imputations()).
.boot_mi_notes <- function(counts, B, conf_level, imputations, fmi,
                           widening, denominator = "observed",
                           sets_aside = NULL, kept = 1) {
  notes <- character()
  # No interval when fewer than half of the replications returned a
  # value (.omega_c_mi()); the caution about the limits then has none
  # to speak of.
  short_of_half <- isTRUE(counts$require_half) && counts$returned >= 2L &&
    counts$returned < B / 2
  short_of_count <- identical(counts$bootstrap_status, "insufficient_count") &&
    counts$returned >= 2L
  if (counts$failed > 0L) {
    top <- ""
    if (length(counts$imputation_messages)) {
      msg <- names(counts$imputation_messages)[
        which.max(counts$imputation_messages)]
      top <- sprintf(" (the imputation failed most often with: %s)",
                     sub("[.]+$", "", trimws(gsub("\\s+", " ",
                         strsplit(msg, "\n")[[1L]][1L]))))
    }
    notes <- c(notes, sprintf(paste0(
      "%d of %d bootstrap replications %s set aside (%s)%s; a ",
      "replication averages those of its %d completed data sets that ",
      "return a value, as the point estimate does, and is set aside when ",
      "fewer than half do%s"),
      counts$failed, B, if (counts$failed == 1L) "was" else "were",
      .bootstrap_reason_phrase(counts$reasons, counts$error_messages,
                               counts$undetermined, counts$improper),
      top, imputations,
      if (counts$returned >= 2L && !short_of_half && !short_of_count) {
        sprintf(paste0("; the interval is computed from the %d that ",
                       "returned a value%s"),
                counts$returned,
                # The count is held to the request; the note says so and
                # says to request more where replications can fail.
                sub("\\.$", "", sub(
                  "^ The", "; the",
                  .dmar_boot_request_note(B, counts$returned,
                                          counts$required, "replications"))))
      } else if (short_of_count) {
        sprintf(paste0(
          "; the %d requested are fewer than the %d required at this ",
          "level, so no interval is reported (see ?dmar_bootstrap). The ",
          "point estimate is not affected"),
          B, counts$required)
      } else if (short_of_half) {
        sprintf(paste0(
          "; fewer than half returned a value, so no interval is reported: ",
          "the %d that did are a minority of the resamples, chosen by ",
          "whether they returned a value, and their spread describes those ",
          "resamples, not the sampling distribution of the coefficient"),
          counts$returned)
      } else {
        "; fewer than two returned a value, so no interval is reported"
      }))
  }
  n_sets <- sum(sets_aside)
  if (length(sets_aside) && n_sets > 0L) {
    notes <- c(notes, sprintf(paste0(
      "In the %d replications that returned a value, %d of their %d ",
      "completed data sets returned none (%s) and were set aside, and ",
      "each of those replications averages the rest, as the point ",
      "estimate does%s"),
      counts$returned, n_sets, counts$returned * as.integer(imputations),
      .omega_c_mi_why(rep(names(sets_aside), sets_aside)),
      if (isTRUE(sets_aside[["improper"]] > 0L)) {
        paste0("; an improper set has the largest loadings, so setting it ",
               "aside lowers its replication's mean, as it lowers the ",
               "point estimate")
      } else {
        ""
      }))
  }
  if (counts$undetermined > 0L || isTRUE(sets_aside[["undetermined"]] > 0L)) {
    notes <- c(notes, .omega_c_undetermined_note(denominator))
  }
  if (counts$improper > 0L) {
    notes <- c(notes, paste0(
      "An improper replication is one in which an error variance fell ",
      "below zero (a standardized loading past one), so setting these ",
      "aside can trim the upper tail of the bootstrap distribution"))
  }
  if (isTRUE(counts$caution) && counts$returned >= 2L && !short_of_half &&
      !short_of_count) {
    notes <- c(notes, sprintf(paste0(
      "Caution: %d of %d replications (%s) were set aside, at least the ",
      "%s of the bootstrap distribution that lies beyond each limit of a ",
      "%s%% interval, so the limits depend on how those replications ",
      "were treated; report the interval with these counts"),
      counts$failed, B, .pct(counts$failed / B),
      .pct((1 - conf_level) / 2), .format_conf_pct(conf_level)))
  }
  # The widening is that of the interval, so it is given only when an
  # interval was formed; with fewer than half of the replications
  # returning a value (or fewer than two) there is none, and the fraction
  # of missing information behind the widening rests on a minority of the
  # resamples, chosen by whether they returned a value. Through the first
  # versions of 1.1.0 the sentence followed "no interval is reported".
  if (!short_of_half && counts$returned >= 2L && is.finite(widening) &&
      widening > 1.05) {
    # When completed data sets returned no value, the number named is
    # even: an odd number sets aside more replications than the even
    # number above it (.mi_widening_imputations()). That holds for the
    # sets of the replications kept (n_sets) and for the replications set
    # aside because fewer than half of their sets returned a value, whose
    # reasons are those of a set, but not for those whose imputation
    # itself failed (reason "imputation"): no completed data set exists
    # there, the half rule never applies, and one call to mice over m
    # chains fails with a probability that rises with m with no parity to
    # it. Through the first versions of 1.1.0 those were counted too, and
    # the warning named an even number for a reason that did not apply.
    n_imp <- sum(counts$reasons[names(counts$reasons) == "imputation"])
    even <- n_sets > 0L || (counts$failed - n_imp) > 0L
    notes <- c(notes, sprintf(paste0(
      "Imputation noise widens the interval by about %s (the fraction of ",
      "missing information is about %s, with %d imputations in each ",
      "replication); imputations = %d would bring the widening under 5%%%s"),
      .pct(widening - 1), formatC(fmi, format = "f", digits = 2),
      as.integer(imputations),
      .mi_widening_imputations(fmi, kept = kept, even = even),
      if (even) {
        paste0(" (an even number, since a replication is set aside when ",
               "fewer than half of its completed data sets return a value, ",
               "and an odd number sets more aside than the even number ",
               "above it)")
      } else {
        ""
      }))
  }
  notes
}


# .boot_mi_sets_set_aside(codes)
#
# The completed data sets set aside within the replications that
# returned a value, counted by reason, from the status codes
# .mi_replicate() records for each set (the engine's codes, .BOOT_STATUS:
# 2 improper, 3 error, 4 non_finite, 5 undetermined, 6 not_converged),
# named by the reasons of the observed data's sets (.OMEGA_C_MI_REASONS).
#
# Inputs:
#   codes : the B_used x M_b matrix of status codes.
# Returns:
#   A named integer vector over .OMEGA_C_MI_REASONS.
.boot_mi_sets_set_aside <- function(codes) {
  reason <- c("2" = "improper", "3" = "error", "4" = "non_finite",
              "5" = "undetermined", "6" = "non_converged")
  codes <- as.integer(codes)
  .omega_c_mi_tally(c(unname(reason[as.character(codes[codes != 1L])]),
                      character()))
}


# .omega_c_mi_point(data, items, impute, m, maxit, denominator)
#
# m imputations of the observed data and the categorical omega fit on
# each completed data set. Every fit's outcome is recorded, set by set,
# rather than acted on, so the caller can decide on all of them at once
# (.omega_c_mi_check()) and later take the ones the point estimate uses
# (.omega_c_mi_take()). An imputation that raises an error stops the
# call: it concerns the data, not one completed data set.
#
# Returns:
#   A list whose first nine elements have one entry per completed data
#   set: `estimates` (NA where the fit returned none), `status`
#   ("proper", "improper", "undetermined", "error", "non_converged", or
#   "non_finite"; "error" when the fit, lavaan's included, stopped with
#   an error), `messages` (the error message behind "error", NA
#   otherwise), `improper` (the items found improper, a character vector
#   per set), `undetermined` (the item pairs whose polychoric the set
#   leaves undetermined, a two-column matrix per set), `resid` (the
#   residual polychoric correlations of a proper set, NA for an
#   undetermined pair; NULL for the others), `polr_fallbacks`,
#   `unidentified` (TRUE for a set set aside as "undetermined" because
#   its determined correlations do not identify the loadings, under
#   either denominator), and `ridge` (the times mice applied its ridge
#   penalty in imputing the set; .mi_logged_counts()); then `df`,
#   `items`, and `denominator`.
.omega_c_mi_point <- function(data, items, impute, m, maxit, denominator) {
  comp <- tryCatch(.mi_complete(data, items, impute, m, maxit),
                   error = function(e) e)
  if (inherits(comp, "error")) {
    .dmar_refuse(c(
      "The imputation of the items failed: ", conditionMessage(comp)),
      "computation_failed")
  }
  status <- character(m)
  estimates <- rep(NA_real_, m)
  messages <- rep(NA_character_, m)
  improper <- rep(list(character()), m)
  undetermined <- rep(list(matrix(character(0), ncol = 2L)), m)
  resid <- vector("list", m)
  unidentified <- logical(m)
  df <- NULL
  for (k in seq_len(m)) {
    parts <- tryCatch(.omega_c_parts(comp$data[[k]]), error = function(e) e)
    if (inherits(parts, "error")) {
      status[k] <- "error"
      messages[k] <- conditionMessage(parts)
      next
    }
    # A fit lavaan stopped with an error is counted, with its message, as
    # an error, apart from one that ran and did not converge.
    if (!isTRUE(parts$converged)) {
      fit_error <- parts[["fit_error"]]
      if (is.null(fit_error)) {
        status[k] <- "non_converged"
      } else {
        status[k] <- "error"
        messages[k] <- fit_error
      }
      next
    }
    df <- parts$df
    und <- parts$undetermined
    # A set whose determined correlations do not identify the loadings
    # (.omega_c_identification()) has no categorical omega under either
    # denominator, and its improper verdict would read arbitrary
    # loadings, so it is set aside first, among the undetermined ones.
    # Through the first versions of 1.1.0 the model implied denominator
    # kept it as proper.
    if (isFALSE(parts$identified)) {
      status[k] <- "undetermined"
      undetermined[[k]] <- und
      unidentified[k] <- TRUE
      next
    }
    if (length(parts$improper)) {
      status[k] <- "improper"
      improper[[k]] <- parts$improper
      next
    }
    if (nrow(und) > 0L) {
      undetermined[[k]] <- und
      if (!identical(denominator, "model_implied")) {
        status[k] <- "undetermined"
        next
      }
    }
    rho_den <- if (identical(denominator, "model_implied")) {
      parts$rho_model
    } else {
      parts$rho_sat
    }
    est <- .gy_omega_from_components(parts$thresh, parts$rho_num, rho_den)
    if (!is.finite(est)) {
      status[k] <- "non_finite"
      next
    }
    status[k] <- "proper"
    estimates[k] <- est
    r <- parts$rho_sat - parts$rho_model
    for (u in seq_len(nrow(und))) {
      r[und[u, 1L], und[u, 2L]] <- NA_real_
      r[und[u, 2L], und[u, 1L]] <- NA_real_
    }
    resid[[k]] <- r
  }
  fallbacks <- comp$polr_fallbacks_each %||% integer(m)
  ridge <- comp$ridge_each %||% integer(m)
  list(estimates = estimates, status = status, messages = messages,
       improper = improper, undetermined = undetermined, resid = resid,
       polr_fallbacks = as.integer(fallbacks), unidentified = unidentified,
       ridge = as.integer(ridge), df = df, items = items,
       denominator = denominator)
}


# .omega_c_mi_combine(a, b) / .omega_c_mi_take(p, idx)
#
# The completed data sets of two runs of .omega_c_mi_point(), as one, in
# the order made; and the sets idx of one.
.OMEGA_C_MI_PER_SET <- c("estimates", "status", "messages", "improper",
                         "undetermined", "resid", "polr_fallbacks",
                         "unidentified", "ridge")

.omega_c_mi_combine <- function(a, b) {
  out <- lapply(.OMEGA_C_MI_PER_SET, function(nm) c(a[[nm]], b[[nm]]))
  names(out) <- .OMEGA_C_MI_PER_SET
  out$df <- a$df %||% b$df
  out$items <- a$items %||% b$items
  out$denominator <- a$denominator %||% b$denominator
  out
}

.omega_c_mi_take <- function(p, idx) {
  out <- lapply(.OMEGA_C_MI_PER_SET, function(nm) p[[nm]][idx])
  names(out) <- .OMEGA_C_MI_PER_SET
  out$df <- p$df
  out$items <- p$items
  out$denominator <- p$denominator
  out
}


# .omega_c_mi_tally(status) / .omega_c_mi_why(status) /
# .omega_c_mi_pairs(undetermined)
#
# The completed data sets that returned no value, counted by reason (a
# named integer vector over the five reasons, in a fixed order), the
# same as a phrase ("2 gave an improper (Heywood) solution, 1 did not
# converge"), and the distinct undetermined pairs of a list of per-set
# pair matrices, as one two-column matrix.
.OMEGA_C_MI_REASONS <- c(
  improper = "gave an improper (Heywood) solution",
  undetermined = "left a polychoric correlation undetermined",
  error = "raised an error",
  non_converged = "did not converge",
  non_finite = "returned a non-finite value")

.omega_c_mi_tally <- function(status) {
  tab <- table(factor(status[status != "proper"],
                      levels = names(.OMEGA_C_MI_REASONS)))
  stats::setNames(as.integer(tab), names(tab))
}

.omega_c_mi_why <- function(status) {
  tab <- .omega_c_mi_tally(status)
  paste(sprintf("%d %s", tab[tab > 0L], .OMEGA_C_MI_REASONS[tab > 0L]),
        collapse = ", ")
}

.omega_c_mi_pairs <- function(undetermined) {
  pairs <- matrix(character(0), ncol = 2L)
  for (u in undetermined) pairs <- rbind(pairs, u)
  unique(pairs)
}


# .omega_c_mi_check(point, denominator)
#
# The rule for a completed data set of the observed data that returns no
# categorical omega. The analysis cannot be computed on such a set: an
# improper solution (a standardized loading above one, a Heywood case)
# above all, for which categorical omega is not defined; under either
# denominator, determined polychoric correlations that do not identify
# the single-factor loadings (.omega_c_identification()), counted with
# the undetermined sets; and under the observed denominator any
# polychoric correlation the set leaves undetermined. Rubin's (1987)
# rules presuppose an estimate on every
# completed data set, and practice treats a set without one in two ways.
# Stopping when any one fails, as Jia and Wu (2019) counted a replication
# of their simulation as converged only if the model converged on all 50
# of its completed data sets, makes the method fail by design: with a
# share p of such sets, M imputations fail with probability
# 1 - (1 - p)^M, which goes to one as von Hippel's rule asks for more.
# Setting them aside and averaging the rest is what lavaan.mi
# (Jorgensen, version 0.1-1) does by default with imputations that did
# not converge, and what Lim and Cheung (2022) did, who counted an
# analysis as failed when fewer than half of its imputations converged;
# it conditions the estimate on a property of the estimate. Both keep an
# improper solution, whose estimates exist (lavaan.mi's documentation:
# such solutions "are likely to occur due to sampling error"), but
# categorical omega has no value on one, so here it joins the sets that
# did not converge. The function sets them all aside, as it sets aside
# the bootstrap replications, counts them by reason, and reports them in
# the warning (.omega_c_mi_point_notes()) and in the printed lines; it
# stops only when fewer than half of the completed data sets return a
# value, Lim and Cheung's criterion, since the average of the rest would
# then describe a minority of the imputations, chosen by whether they
# returned one. Because the decision rests on a share, its outcome
# settles as M grows instead of tending to failure.
#
# Evidence for the help page's mixed-skew example (61 samples of its
# design, N = 1000, 20 imputations each; the scripts are not part of
# the package): no complete data fit and no pairwise deletion fit was
# improper (one pairwise fit did not converge), but 13 samples had at
# least one improper completed data set under predictive mean matching,
# with a share of 3% on average and up to 65%.
# The improper sets are those with the largest loadings, so setting them
# aside lowers the estimate: against refitting them with every loading
# bounded at one (constrained estimation), the mean of the rest was lower
# by about 0.034 times the share set aside (0.0017 at 5%, a tenth of the
# complete data estimate's sampling standard deviation of 0.0165), and
# the two treatments were equally far from the complete data estimate
# (root mean square 0.0126 and 0.0122 over the 13 samples). See the help
# page for why imputation produces improper solutions the observed data
# do not.
#
# Inputs:
#   point : a list from .omega_c_mi_point() (or combined runs).
# Returns:
#   invisible(TRUE) when at least half of the sets returned a value;
#   otherwise an error naming how many did not and why.
.omega_c_mi_check <- function(point) {
  bad <- point$status != "proper"
  m <- length(point$status)
  if (sum(!bad) >= m / 2) return(invisible(TRUE))
  tab <- .omega_c_mi_tally(point$status)
  advice <- .omega_c_mi_advice(point, tab)
  .dmar_refuse_argument(c(
    "Categorical omega could not be computed on ", sum(bad), " of the ",
    m, " completed data sets of the observed data (",
    .omega_c_mi_why(point$status), "), more than half, so an average ",
    "of the rest would describe a minority of the imputations, chosen ",
    "by whether they returned a value. ",
    paste(advice, collapse = " ")), arg = NULL)
}


# .omega_c_mi_advice(point, tab)
#
# The sentences that say, for each reason a completed data set of the
# observed data returned no value, what it was and what to do. A set is
# "undetermined" under either denominator when its determined
# correlations do not identify the loadings (`unidentified`), and under
# the observed denominator also when they do but leave a polychoric
# correlation undetermined (.omega_c_mi_point()). Only the second kind
# returns a value under the model implied denominator: its fit is the
# same, proper and identified. So the model implied denominator is
# offered only when those sets, with the ones that returned a value,
# reach half of the completed data sets, the share the call needs
# (.omega_c_mi_check()); otherwise the advice says that switching would
# not be enough, rather than send the user to the same stop.
.omega_c_mi_advice <- function(point, tab) {
  advice <- character()
  if (tab[["improper"]] > 0L) {
    advice <- c(advice, paste0(
      "The improper solutions had a standardized loading above one for ",
      paste(unique(unlist(point$improper)), collapse = " and "),
      ". Imputation produces them when the imputed values leave few ",
      "respondents in the sparse cells that determine an item's ",
      "polychoric correlations (see Missing data in ",
      "?reliability_omega_categorical); a larger sample, or reconsidering ",
      "the item, is the usual remedy."))
  }
  # Sets whose determined correlations do not identify the loadings have
  # no value under either denominator, so only the others are told that
  # the model implied denominator would keep them.
  unid <- .omega_c_mi_unidentified(point)
  if (any(unid)) {
    first <- which(unid)[1L]
    und <- point$undetermined[[first]]
    items <- point$items %||% unique(as.vector(und))
    advice <- c(advice, paste0(
      "In ", sum(unid), " of them the correlations the data determine ",
      "leave the single-factor loadings unidentified. In the first: ",
      .omega_c_unidentified_message(.omega_c_identification(items, und),
                                    und)))
  }
  only <- point$status == "undetermined" & !unid
  if (any(only)) {
    m <- length(point$status)
    would <- sum(point$status == "proper") + sum(only)
    what <- paste0(
      "The data do not determine ",
      .undetermined_message(.omega_c_mi_pairs(point$undetermined[only])),
      " in ", sum(only), " completed data ",
      if (sum(only) == 1L) "set" else "sets",
      " whose determined correlations do identify the loadings, and the ",
      "observed denominator needs every polychoric correlation. ")
    those <- if (sum(only) == 1L) "that set" else "those sets"
    advice <- c(advice, if (would >= m / 2) {
      paste0(what, "Use denominator = \"model_implied\" (Green & Yang, ",
             "2009), which needs no unrestricted polychoric correlation ",
             "and would return a value on ", those, " as well, ", would,
             " of the ", m, ".")
    } else {
      paste0(what, "Switching to denominator = \"model_implied\", which ",
             "needs no unrestricted polychoric correlation, would return a ",
             "value on ", those, " as well, but on only ", would, " of the ",
             m, " in all, still fewer than half, so it is not a remedy here.")
    })
  }
  if (tab[["error"]] > 0L) {
    first <- point$messages[!is.na(point$messages)][1L]
    advice <- c(advice, paste0("The first error was: ", first))
  }
  advice
}


# .omega_c_mi_unidentified(point)
#
# Which completed data sets were set aside because their determined
# correlations do not identify the loadings (FALSE throughout for a
# record without the element, such as a test's).
.omega_c_mi_unidentified <- function(point) {
  u <- point$unidentified
  if (is.null(u)) return(logical(length(point$status)))
  as.logical(u)
}


# .omega_c_mi_point_notes(point)
#
# The sentences of the warning when some completed data sets of the
# observed data that the point estimate would average were set aside
# (.omega_c_mi_check()); empty when none was.
.omega_c_mi_point_notes <- function(point) {
  bad <- point$status != "proper"
  if (!any(bad)) return(character())
  m <- length(point$status)
  tab <- .omega_c_mi_tally(point$status)
  notes <- sprintf(paste0(
    "%d of the %d completed data sets of the observed data %s set aside ",
    "(%s), and the point estimate averages the %d that returned a value"),
    sum(bad), m, if (sum(bad) == 1L) "was" else "were",
    .omega_c_mi_why(point$status), sum(!bad))
  if (tab[["improper"]] > 0L) {
    notes <- c(notes, paste0(
      "An improper set is one in which an error variance fell below zero ",
      "(a standardized loading past one; here for ",
      paste(unique(unlist(point$improper)), collapse = " and "),
      "), so setting these aside lowers the estimate, the more so the ",
      "more are set aside; report their number with the estimate"))
  }
  unid <- .omega_c_mi_unidentified(point)
  if (any(unid)) {
    notes <- c(notes, paste0(
      "In ", sum(unid), " of the completed data sets set aside as leaving ",
      "a polychoric correlation undetermined, the correlations left ",
      "determined do not identify the single-factor loadings (they do not ",
      "connect every item to the others and include a cycle of an odd ",
      "number of items), so categorical omega has no value there under ",
      "either denominator"))
  }
  if (any(point$status == "undetermined" & !unid)) {
    notes <- c(notes, paste0(
      "denominator = \"model_implied\" needs no unrestricted polychoric ",
      "correlation and keeps every completed data set that leaves one ",
      "undetermined, provided the correlations left determined identify ",
      "the loadings"))
  }
  if (tab[["error"]] > 0L) {
    first <- point$messages[!is.na(point$messages)][1L]
    notes <- c(notes, paste0("The first error was: ",
                             trimws(gsub("\\s+", " ", first))))
  }
  notes
}


# -----------------------------------------------------------------------------
# (E) Printed lines
# -----------------------------------------------------------------------------

# .missing_treatment_footer_lines(x)
#
# The lines a table prints under itself to name a missing data treatment
# whose assumption the reader needs beside the estimate: pairwise
# deletion (missing = "pairwise"), which is consistent only when values
# are missing completely at random, and multiple imputation (missing =
# "mi"), with its imputation model, its assumption, and the counts behind
# the point estimate and the interval. .dmar_tbl_footer()
# (R/dmar_output_helpers.R) calls it, so the console and a knitted report
# show the same lines. A table whose `missing` attribute is anything else
# (listwise deletion, full information maximum likelihood) gets no line
# here, so every result printed before 1.1.0 prints as it did. Numbers
# are formatted for display only.
#
# Inputs:
#   x : a table, read for its `missing`, `aux`, `imputation`, and
#       `bootstrap` attributes.
# Returns:
#   A character vector named by the kind of each line ("missing",
#   "imputation"); empty when there is nothing to print.
.missing_treatment_footer_lines <- function(x) {
  missing <- attr(x, "missing", exact = TRUE)
  if (!is.character(missing) || length(missing) != 1L) return(character())
  if (identical(missing, "pairwise")) {
    return(c(missing = paste0(
      "Missing data: pairwise deletion, consistent only when values are ",
      "missing completely at random, as listwise deletion is (Asparouhov ",
      "& Muth\u00e9n, 2010); it uses more of the data, not a weaker ",
      "assumption, and missing = \"mi\" is the remedy when values are ",
      "missing at random")))
  }
  if (!identical(missing, "mi")) return(character())
  mi <- attr(x, "imputation", exact = TRUE)
  if (!is.list(mi)) return(character())
  aux <- attr(x, "aux", exact = TRUE)
  model <- if (identical(mi$impute, "polr")) {
    "the proportional odds model"
  } else {
    "predictive mean matching"
  }
  bootstrap <- is.list(attr(x, "bootstrap", exact = TRUE))
  lines <- c(missing = paste0(
    "Missing data: multiple imputation by ", model, " (mice), assuming ",
    "values missing at random given the items",
    if (length(aux)) {
      paste0(" and the auxiliary ",
             if (length(aux) == 1L) "variable " else "variables ",
             paste(aux, collapse = ", "))
    } else {
      ""
    },
    if (bootstrap) {
      "; interval by bootstrap then impute (Schomaker & Heumann, 2018)"
    } else {
      ""
    }))
  fmt2 <- function(v) formatC(v, format = "f", digits = 2)
  imp <- sprintf("Imputations: %d of the observed data", mi$imputations_point)
  rule <- bootstrap && !is.na(mi$imputations_rule)
  if (rule) {
    imp <- paste0(imp, sprintf(
      " (von Hippel's rule asks for %d, the floor is %d)",
      mi$imputations_rule, mi$imputations_floor))
  }
  # Completed data sets set aside because they returned no value
  # (.omega_c_mi_check()); results made before the counts were kept
  # carry no `status`.
  if (is.character(mi$status) && any(mi$status != "proper")) {
    n_aside <- sum(mi$status != "proper")
    imp <- paste0(imp, sprintf(", %d of them set aside (%s)", n_aside,
                               .omega_c_mi_why(mi$status)))
  }
  # No interval was formed when fewer than half of the replications (or
  # fewer than two) returned a value (.omega_c_mi()): the fraction of
  # missing information is then that of the replications that did, which
  # still chose the number of imputations above and is given as such, and
  # there is no widening of an interval to give.
  bs <- attr(x, "bootstrap", exact = TRUE)
  no_interval <- is.list(bs) && is.numeric(bs$returned) &&
    length(bs$returned) == 1L &&
    (bs$returned < 2L || (isTRUE(bs$require_half) &&
                            isTRUE(bs$returned < bs$attempted / 2)) ||
       identical(bs$bootstrap_status, "insufficient_count"))
  if (rule && no_interval) {
    imp <- paste0(imp, sprintf(paste0(
      "; %d in each bootstrap replication; fraction of missing ",
      "information %s (95%% upper limit %s), estimated from the %d ",
      "replications that returned a value, a minority chosen by whether ",
      "they returned one; no interval is reported, so there is no widening ",
      "to give"),
      mi$imputations_replicate, fmt2(mi$fmi), fmt2(mi$fmi_ci_upper),
      as.integer(bs$returned)))
  } else if (rule) {
    imp <- paste0(imp, sprintf(paste0(
      "; %d in each bootstrap replication; fraction of missing ",
      "information %s (95%% upper limit %s); imputation noise widens the ",
      "interval by about %s"),
      mi$imputations_replicate, fmt2(mi$fmi), fmt2(mi$fmi_ci_upper),
      .pct(mi$widening - 1)))
  } else if (bootstrap) {
    imp <- paste0(imp, sprintf(paste0(
      "; %d in each bootstrap replication; too few replications ",
      "returned a value to estimate the fraction of missing information"),
      mi$imputations_replicate))
  } else {
    imp <- paste0(imp, paste0(
      "; the number von Hippel's (2020) rule asks for needs the fraction ",
      "of missing information, which only the bootstrap estimates"))
  }
  fallbacks <- sum(mi$polr_fallbacks)
  if (isTRUE(fallbacks > 0L)) {
    imp <- paste0(imp, sprintf(paste0(
      "; the proportional odds model fell back to a multinomial model ",
      "%d %s"), fallbacks, if (fallbacks == 1L) "time" else "times"))
  }
  c(lines, imputation = imp)
}
