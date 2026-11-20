# .dmar_check_count(): the shared check of a count of resamples,
# permutations, or replications (R/dmar_conditions.R).

check_count <- function(...) DMAR:::.dmar_check_count(...)

test_that("a whole count is returned as an integer, from a double or an integer", {
  expect_identical(check_count(1000, "B"), 1000L)
  expect_identical(check_count(1000L, "B"), 1000L)
  expect_identical(check_count(1, "n_resamples"), 1L)
  expect_identical(check_count(0, "n_resamples", min = 0L), 0L)
  expect_identical(check_count(10, "G", min = 10L), 10L)
  expect_identical(check_count(50, "G", min = 2L, max = 50), 50L)
  expect_identical(check_count(.Machine$integer.max, "B"), .Machine$integer.max)
  expect_visible(check_count(5, "B"))
})

test_that("a value that is not a single finite whole number is refused", {
  refused <- list(NULL, NA, NA_real_, NA_integer_, NaN, Inf, -Inf, "10", TRUE,
                  factor(10), list(10), c(10, 20), numeric(0), integer(0),
                  2.5, 999.5, 1 + 1e-9, 0, -1, -3L)
  for (value in refused) {
    expect_error(check_count(value, "n_resamples"),
                 "'n_resamples' must be a single whole number of at least 1.",
                 fixed = TRUE, class = "dmar_count_invalid")
  }
})

test_that("the minimum is the one stated, and the message names it", {
  expect_error(check_count(-1, "n_resamples", min = 0L),
               "'n_resamples' must be a single whole number of at least 0.",
               fixed = TRUE)
  expect_error(check_count(9, "G", min = 10L),
               "'G' must be a single whole number of at least 10.",
               fixed = TRUE)
  expect_error(check_count(1, "G", min = 2L),
               "'G' must be a single whole number of at least 2.",
               fixed = TRUE)
  expect_error(check_count(99, "B", min = 100L), "at least 100.", fixed = TRUE)
})

test_that("a stated maximum is enforced and named", {
  expect_error(check_count(51, "G", min = 2L, max = 50),
               "'G' must be a single whole number of at least 2 and at most 50.",
               fixed = TRUE, class = "dmar_count_invalid")
  expect_error(check_count(1, "G", min = 2L, max = 50),
               "'G' must be a single whole number of at least 2 and at most 50.",
               fixed = TRUE)
  expect_error(check_count(2.5, "G", min = 2L, max = 100000),
               "at least 2 and at most 100000.", fixed = TRUE)
})

test_that("a count above the largest integer is refused, with the limit named", {
  expect_error(check_count(1e10, "B"),
               "'B' must be a single whole number of at least 1 and at most 2147483647.",
               fixed = TRUE, class = "dmar_count_invalid")
  expect_error(check_count(1e10, "B", max = 1e12), "at most 2147483647.",
               fixed = TRUE)
})

test_that("the refusal is classed and carries the argument and its limits", {
  e <- tryCatch(check_count(2.5, "n_resamples", min = 0L), error = function(e) e)
  expect_identical(class(e),
                   c("dmar_count_invalid", "dmar_refusal", "error", "condition"))
  expect_identical(e$code, "count_invalid")
  expect_identical(e$arg, "n_resamples")
  expect_identical(e$min, 0L)
  expect_null(e$max)
  expect_null(conditionCall(e))

  e_max <- tryCatch(check_count(7, "G", max = 5), error = function(e) e)
  expect_equal(e_max$max, 5)
})

test_that("the check is made before anything is allocated", {
  allocated <- FALSE
  draw <- function(B) {
    B <- check_count(B, "B")
    allocated <<- TRUE
    numeric(B)
  }
  expect_error(draw(2.5), class = "dmar_count_invalid")
  expect_false(allocated)
  expect_length(draw(3), 3L)
})

test_that("degrees of freedom are not counts: a fractional value is accepted", {
  # The count check is for numbers of draws. A degrees of freedom argument
  # takes any positive value, and none of these calls goes through it.
  br <- df_barnard_rubin(between = 0.02, within = 0.10, m = 20,
                         df_complete = 48.5)
  expect_true(is.finite(br$df))
  f_chisq <- convert_F_chisq(F_value = 2.4, df_numerator = 3,
                             df_denominator = 12.5)
  expect_true(all(is.finite(f_chisq$value)))

  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  fit_t <- cfa_1(holzinger_swineford[, items], wald_df = 12.5)
  fit_z <- cfa_1(holzinger_swineford[, items])
  expect_identical(attr(fit_t, "wald_df"), 12.5)
  expect_gt(fit_t$ci_upper[1] - fit_t$ci_lower[1],
            fit_z$ci_upper[1] - fit_z$ci_lower[1])
})

test_that("desired_power is strictly inside (0, 1), by one classed check", {
  check <- DMAR:::.dmar_check_desired_power
  for (ok in c(1e-6, 0.5, 0.80, 0.999999)) expect_identical(check(ok), ok)
  for (bad in list(0, 1, -0.1, 1.5, NA_real_, Inf, c(0.8, 0.9), "0.8", NULL,
                   TRUE)) {
    err <- expect_error(check(bad), class = "dmar_desired_power_invalid")
    expect_s3_class(err, "dmar_refusal")
    expect_identical(err$arg, "desired_power")
    expect_identical(err$value, bad)
    expect_match(conditionMessage(err),
                 "'desired_power' must be a single number in (0, 1)",
                 fixed = TRUE)
  }
  expect_match(conditionMessage(tryCatch(check(1), error = identity)),
               "1 was given", fixed = TRUE)
})

test_that("every power planner with desired_power goes through the one check", {
  # A planner either makes the check itself or hands desired_power to the
  # planner that does. The bodies are read, so a new ss_power_ function
  # that checks desired_power in its own words fails here.
  ns <- asNamespace("DMAR")
  planners <- grep("^ss_power_", getNamespaceExports(ns), value = TRUE)
  with_power <- Filter(function(f) "desired_power" %in% names(formals(get(f, ns))),
                       planners)
  expect_gte(length(with_power), 25L)
  for (f in with_power) {
    body_text <- paste(deparse(body(get(f, ns))), collapse = " ")
    checks <- grepl(".dmar_check_desired_power(desired_power)", body_text,
                    fixed = TRUE)
    hands_on <- grepl("desired_power = desired_power", body_text, fixed = TRUE)
    expect_true(checks || hands_on, label = f)
  }
})
