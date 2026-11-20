# Every error and every warning DMAR signals on purpose is a classed
# condition made by .dmar_refuse(), .dmar_refuse_argument(), a
# .dmar_check_*() helper, or .dmar_warn(), and every class is listed on
# ?dmar_metadata. The checks read the functions of the
# installed namespace, so they hold for the package as built.

dmar_calls <- function(names) {
  ns <- asNamespace("DMAR")
  found <- list()
  walk <- function(x, owner) {
    if (is.call(x)) {
      head <- x[[1L]]
      if (is.symbol(head) && as.character(head) %in% names) {
        found[[length(found) + 1L]] <<- list(owner = owner, call = x)
      }
      # x[[i]] is read in place: an empty argument (the blank in x[i, ])
      # cannot be assigned to a variable.
      for (i in seq_along(x)) {
        if (is.call(x[[i]]) || is.pairlist(x[[i]])) walk(x[[i]], owner)
      }
    } else if (is.pairlist(x)) {
      for (i in seq_along(x)) {
        if (is.call(x[[i]])) walk(x[[i]], owner)
      }
    }
  }
  for (nm in ls(ns, all.names = TRUE)) {
    obj <- get(nm, envir = ns)
    if (is.function(obj) && !is.primitive(obj)) walk(body(obj), nm)
  }
  found
}

test_that("no function signals a plain error or warning of its own", {
  skip_on_cran()
  helpers <- c(".dmar_refuse", ".dmar_warn")
  plain <- dmar_calls(c("stop", "warning", "stopifnot"))
  plain <- plain[!vapply(plain, function(p) p$owner %in% helpers, logical(1))]
  # A condition caught and passed on unchanged keeps the class it came
  # with: stop(e), warning(w). The two conditions of a bootstrap's progress
  # callback are documented as not being refusals (?dmar_metadata).
  passed_on <- vapply(plain, function(p) {
    cl <- p$call
    length(cl) == 2L &&
      (is.symbol(cl[[2L]]) ||
         (is.call(cl[[2L]]) &&
            identical(cl[[2L]][[1L]], as.name("errorCondition"))))
  }, logical(1))
  left <- plain[!passed_on]
  expect_identical(
    vapply(left, function(p) paste0(p$owner, ": ",
                                    paste(deparse(p$call), collapse = " ")),
           character(1)),
    character(0))
})

test_that("every condition class in use is listed on ?dmar_metadata", {
  skip_on_cran()
  calls <- dmar_calls(c(".dmar_refuse", ".dmar_warn"))
  codes <- unlist(lapply(calls, function(p) {
    cl <- as.list(p$call)[-1L]
    code <- if (!is.null(names(cl)) && "code" %in% names(cl)) cl$code else
      if (length(cl) >= 2L) cl[[2L]] else NULL
    if (is.character(code)) code else NULL
  }))
  codes <- sort(unique(c(codes, "argument_invalid")))
  expect_gt(length(codes), 40L)
  # The help sources when the tests run from the source tree, the installed
  # help when they run from a built package.
  src <- file.path("..", "..")
  rd <- if (dir.exists(file.path(src, "man"))) tools::Rd_db(dir = src) else
    tools::Rd_db("DMAR")
  pages <- c("dmar_metadata.Rd", "dmar_bootstrap.Rd")
  expect_true(all(pages %in% names(rd)))
  text <- paste(unlist(lapply(rd[pages], function(p)
    paste(as.character(p), collapse = ""))), collapse = "")
  listed <- vapply(codes, function(code)
    grepl(paste0("dmar_", code), text, fixed = TRUE), logical(1))
  expect_identical(codes[!listed], character(0))
})

test_that("an argument check is a refusal that names the argument", {
  e <- tryCatch(ci_nc_t(ncp = 2, df = -1), error = function(e) e)
  expect_s3_class(e, c("dmar_argument_invalid", "dmar_refusal"), exact = FALSE)
  expect_identical(e$code, "argument_invalid")
  expect_identical(e$arg, "df")
  expect_null(conditionCall(e))
  e <- tryCatch(convert_cor_cov(matrix(1, 2, 3), sd = 1:2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "cor_mat")
})

test_that("the general warning classes carry their code", {
  # ci_nc_t() warned (class dmar_nc_t_accuracy in development builds of
  # 1.1.0) past a noncentrality of 37.62; its limits are now exact there,
  # and it does not warn (see test-ci_nc_t.R).
  expect_no_warning(ci_nc_t(ncp = 50, df = 10))
  w <- tryCatch(
    ss_aipe_mean(true_sd = 10, width = 4, true_kurtosis = 1),
    warning = function(w) w)
  expect_s3_class(w, "dmar_argument_ignored")
  expect_identical(w$code, "argument_ignored")
})

test_that("a name in ... that lavaan does not have is refused by the function, with the name as typed", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")
  model <- "f =~ x1 + x2 + x3 + x4"
  S <- matrix(0.5, 4, 4, dimnames = list(paste0("x", 1:4), paste0("x", 1:4)))
  diag(S) <- 1
  refused <- function(expr, bad) {
    err <- tryCatch(expr, error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, bad)
    for (b in bad) expect_match(conditionMessage(err), paste0("'", b, "'"), fixed = TRUE)
    expect_match(conditionMessage(err), "?lavaan::lavOptions", fixed = TRUE)
  }
  refused(ss_aipe_sem_path(model, true_cov_matrix = S, width = .3, which_path = "f=~x2", N = 50), "N")
  refused(ss_aipe_sem_path(model, Sigma = S, desired_width = .3, which_path = "f=~x2"),
          c("Sigma", "desired_width"))
  refused(ss_aipe_sem_path_sensitivity(model, supposed_cov_matrix = S, which_path = "f=~x2",
                                       width = .3, G = 2, desired_width = .3), "desired_width")
  refused(ss_aipe_rmsea_sensitivity(width = .05, model = model, true_cov_matrix = S, G = 2, N = 100), "N")
  refused(ss_aipe_composite_sem(model = model, true_cov_matrix = S, width = .3, mu = rep(0, 4)), "mu")
  refused(ss_power_composite_sem(model = model, true_cov_matrix = S, pop_model = model), "pop_model")
  # An argument without a name is refused too.
  expect_error(ss_aipe_sem_path(model, S, .3, "f=~x2", .95, NULL, FALSE, FALSE, 50),
               class = "dmar_argument_invalid")
  # An argument or an option of lavaan passes the check.
  expect_null(DMAR:::.dmar_check_lavaan_dots(std.lv = TRUE, estimator = "ML", fn = "f"))
  expect_null(DMAR:::.dmar_check_lavaan_dots(fn = "f"))
})

test_that("the shared entry checks refuse by class, with the argument and the value given", {
  skip_on_cran()
  # conf_level: one number strictly inside (0, 1).
  for (bad in list(-1, 0, 1, 1.2, NA, NA_real_, "a", TRUE, c(.9, .95), numeric(0))) {
    err <- tryCatch(ci_c(means = c(2, 4, 9, 13), s_anova = .8, c_weights = c(.5, -.5, -.5, .5),
                         n_per_group = rep(3, 4), N = 12, conf_level = bad),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "conf_level")
  }
  # The level is checked before the tail areas are derived from it.
  err <- tryCatch(ci_pvaf(F_value = 5, df_1 = 2, df_2 = 57, N = 60, conf_level = 1.2),
                  error = function(e) e)
  expect_identical(err$arg, "conf_level")
  err <- tryCatch(ci_eta_squared(F_value = 5, df_effect = 2, df_error = 57, N = 60, conf_level = 1.2),
                  error = function(e) e)
  expect_identical(err$arg, "conf_level")
  # A refusal with a message: a NULL level once gave an empty one.
  err <- tryCatch(ss_aipe_R2(true_R2 = .5, p = 5, width = .2, conf_level = NULL),
                  error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_true(nzchar(conditionMessage(err)))
  # width: the planners that searched without end for a width of 0 or below.
  for (call in list(quote(ss_aipe_sm(true_sm = 10, width = -1)),
                    quote(ss_aipe_R2(true_R2 = .5, p = 5, width = -1)),
                    quote(ss_aipe_rmsea(true_rmsea = .035, df = 50, width = -1)),
                    quote(ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = 0,
                                              n_items = 6, true_item_covariance = .6,
                                              true_residual_var = .4)))) {
    err <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "width")
  }
  # alpha_level.
  err <- tryCatch(ss_power_reg_coef(true_R2_Y_X = .4, true_R2_Y_X_without_j = .3, p = 5,
                                    desired_power = .8, alpha_level = 0),
                  error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "alpha_level")
  # A number that is missing.
  err <- tryCatch(ci_smd(smd = .5, n_1 = NA, n_2 = 10), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "n_1")
  expect_match(conditionMessage(err), "NA was given", fixed = TRUE)
  expect_match(conditionMessage(err), "ci_smd()", fixed = TRUE)
})

test_that(".dmar_check_numbers() leaves alone what it was not asked about", {
  f <- function(a, b = NULL, d = 1) {
    DMAR:::.dmar_check_numbers("f", c("a", "b", "d"), required = "a")
    "ran"
  }
  expect_identical(f(1), "ran")
  expect_identical(f(1, b = NULL), "ran")
  expect_identical(f(c(1, NA, 3)), "ran")          # a vector with a missing value is not a single NA
  expect_identical(f(TRUE, 2, Inf), "ran")
  expect_identical(f(), "ran")                     # an argument left out is not looked at
  for (bad in list(NA, NaN, "a", list(1), numeric(0))) {
    err <- tryCatch(f(bad), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "a")
  }
  expect_error(f(NULL), class = "dmar_argument_invalid")
  expect_identical(tryCatch(f(1, d = NA), error = function(e) e)$arg, "d")
})

test_that("a planning value of 0 for the effect is refused where a sample size is solved for", {
  err <- tryCatch(ss_power_c(true_contrast = 0, c_weights = c(.5, .5, -.5, -.5), true_sd = 1,
                             desired_power = .8), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "true_contrast")
  err <- tryCatch(ss_power_sc(true_std_contrast = 0, c_weights = c(.5, .5, -.5, -.5),
                              desired_power = .8), error = function(e) e)
  expect_identical(err$arg, "true_std_contrast")
  # With a size given the power at 0 is the Type I error rate.
  at0 <- ss_power_c(true_contrast = 0, c_weights = c(.5, .5, -.5, -.5), true_sd = 1, specified_n_per_group = 20)
  expect_equal(at0$value[at0$term == "actual_power"], .05, tolerance = 1e-8)
})

test_that("arg is read from a message that opens with the argument's name", {
  expect_identical(DMAR:::.dmar_leading_arg_names("'x' must be one number."), "x")
  expect_identical(DMAR:::.dmar_leading_arg_names("'a' and 'b' cannot both be given."), c("a", "b"))
  expect_identical(DMAR:::.dmar_leading_arg_names("'a', 'b', and 'c' are needed."), c("a", "b", "c"))
  expect_null(DMAR:::.dmar_leading_arg_names("The package 'lavaan' is needed."))
  err <- tryCatch(bayes_one_sample_t(mean = .5, sd = 1, N = 30, prior_location = NULL),
                  error = function(e) e)
  expect_identical(err$arg, "prior_location")
})
