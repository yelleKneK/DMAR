## ss_aipe_crossed_sensitivity(): the replications are ci_crossed()'s
## interval on arrays drawn as documented, the inputs are echoed, the
## seed is local, and a well specified plan delivers what it promised.
## The family-wide schema check is in test-ss_aipe_sensitivity_family.R.

val <- function(x, term) x$value[x$term == term]

test_that("each replication is ci_crossed() on an array drawn units, raters, residuals", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  res <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    specified_n_units = 15, specified_n_raters = 5, width = 0.8,
    true_grand_mean = 3, G = 6, seed = 11, filename = path)
  per_rep <- utils::read.csv(path)
  expect_named(per_rep, c("grand_mean", "ci_lower", "ci_upper", "ci_width",
                          "df", "type_I_lower", "type_I_upper"))
  expect_equal(nrow(per_rep), 6L)

  # Redraw the same arrays by hand and analyze them with the public
  # function.
  set.seed(11)
  for (g in 1:6) {
    x <- 3 + outer(rnorm(15, 0, sqrt(0.5)), rnorm(5, 0, sqrt(0.1)), "+") +
      matrix(rnorm(75, 0, sqrt(0.4)), 15, 5)
    ci <- ci_crossed(x)
    expect_equal(per_rep$grand_mean[g], val(ci, "grand_mean"), tolerance = 1e-12)
    expect_equal(per_rep$ci_lower[g], val(ci, "ci_lower"), tolerance = 1e-12)
    expect_equal(per_rep$ci_upper[g], val(ci, "ci_upper"), tolerance = 1e-12)
    expect_equal(per_rep$df[g], val(ci, "df"), tolerance = 1e-12)
  }
  expect_equal(val(res, "mean_ci_width"), mean(per_rep$ci_width),
               tolerance = 1e-12)
  expect_equal(val(res, "pct_ci_less_width"), mean(per_rep$ci_width <= 0.8))
  expect_equal(val(res, "mean_grand_mean"), mean(per_rep$grand_mean),
               tolerance = 1e-12)

  # A second call appends to the file.
  ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    specified_n_units = 15, specified_n_raters = 5, width = 0.8,
    G = 3, seed = 11, filename = path)
  expect_equal(nrow(utils::read.csv(path)), 9L)
})

test_that("the planned design is the planner's, and the inputs are echoed", {
  res <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.2, true_var_residual = 0.4,
    supposed_var_unit = 0.5, supposed_var_rater = 0.1,
    supposed_var_residual = 0.4, width = 0.4, rater_cost = 40, G = 4,
    seed = 1)
  plan <- ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, rater_cost = 40)
  expect_equal(val(res, "n_units"), val(plan, "necessary_n_units"))
  expect_equal(val(res, "n_raters"), val(plan, "necessary_n_raters"))
  expect_equal(val(res, "true_var_rater"), 0.2)
  expect_equal(val(res, "supposed_var_rater"), 0.1)
  expect_equal(val(res, "width"), 0.4)
  expect_false("assurance" %in% res$term)

  fixed <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    supposed_var_unit = 0.5, supposed_var_rater = 0.1,
    supposed_var_residual = 0.4, width = 0.4, n_raters = 30,
    assurance = 0.8, G = 4, seed = 1)
  plan <- ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, n_raters = 30,
                          assurance = 0.8)
  expect_equal(val(fixed, "n_units"), val(plan, "necessary_n_units"))
  expect_equal(val(fixed, "n_raters"), 30)
  expect_equal(val(fixed, "assurance"), 0.8)

  spec <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    specified_n_units = 40, specified_n_raters = 8, width = 0.4, G = 4,
    seed = 1)
  # With the sizes specified no planning values are used, so the three
  # supposed_ echoes are present and NA.
  supposed <- c("supposed_var_unit", "supposed_var_rater",
                "supposed_var_residual")
  expect_identical(spec$term[startsWith(spec$term, "supposed_")], supposed)
  expect_true(all(is.na(spec$value[match(supposed, spec$term)])))
  expect_equal(c(val(spec, "n_units"), val(spec, "n_raters")), c(40, 8))
})

test_that("a seed makes the call reproducible without moving the caller's stream", {
  args <- list(true_var_unit = 0.5, true_var_rater = 0.1,
               true_var_residual = 0.4, specified_n_units = 20,
               specified_n_raters = 6, width = 0.6, G = 5)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(ss_aipe_crossed_sensitivity, c(args, seed = 7))
  expect_identical(.Random.seed, before)
  b <- do.call(ss_aipe_crossed_sensitivity, c(args, seed = 7))
  expect_identical(a$value, b$value)
  # Without a seed the current stream is used, and advanced.
  set.seed(7)
  c1 <- do.call(ss_aipe_crossed_sensitivity, args)
  expect_identical(c1$value, a$value)
})

test_that("the planning values and the design are supplied one way only", {
  base <- list(true_var_unit = 0.5, true_var_rater = 0.1,
               true_var_residual = 0.4, width = 0.4, G = 2)
  expect_error(do.call(ss_aipe_crossed_sensitivity, base),
               "must specify either")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base, supposed_var_unit = 0.5)),
               "all three")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base, specified_n_units = 20)),
               "both")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base, supposed_var_unit = 0.5,
                         supposed_var_rater = 0.1,
                         supposed_var_residual = 0.4,
                         specified_n_units = 20, specified_n_raters = 5)),
               "not both")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base, specified_n_units = 20,
                         specified_n_raters = 5, n_raters = 5)),
               "'n_raters' fixes the raters")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base[-3], true_var_residual = 0,
                         specified_n_units = 20, specified_n_raters = 5)),
               "'true_var_residual'")
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(base, specified_n_units = 20,
                         specified_n_raters = 5, filename = 1)),
               "'filename' must be NULL or a single character string")
  expect_null(formals(ss_aipe_crossed_sensitivity)$filename)
  expect_null(formals(ss_aipe_crossed_sensitivity)$seed)
})

test_that("an assurance is checked with a specified design as well as a planned one", {
  spec <- list(true_var_unit = 0.5, true_var_rater = 0.1,
               true_var_residual = 0.4, specified_n_units = 20,
               specified_n_raters = 5, width = 0.4, G = 2, seed = 1)
  plan <- list(true_var_unit = 0.5, true_var_rater = 0.1,
               true_var_residual = 0.4, supposed_var_unit = 0.5,
               supposed_var_rater = 0.1, supposed_var_residual = 0.4,
               width = 0.4, G = 2, seed = 1)
  for (args in list(spec, plan)) {
    expect_error(do.call(ss_aipe_crossed_sensitivity, c(args, assurance = 7)),
                 "'assurance' must be a single number strictly between 0 and 1")
    expect_error(do.call(ss_aipe_crossed_sensitivity,
                         c(args, assurance = "0.8")),
                 "'assurance'")
  }
  # A valid assurance with a specified design is echoed, and the value
  # column stays numeric.
  res <- do.call(ss_aipe_crossed_sensitivity, c(spec, assurance = 0.8))
  expect_equal(val(res, "assurance"), 0.8)
  expect_type(res$value, "double")
})

test_that("the search limits reach the planner", {
  # Units costly and raters free: the cheapest design uses every rater the
  # search allows, and the planner's warning can be acted on here.
  args <- list(true_var_unit = 0.5, true_var_rater = 0.1,
               true_var_residual = 0.4, supposed_var_unit = 0.5,
               supposed_var_rater = 0.1, supposed_var_residual = 0.4,
               width = 0.4, unit_cost = 1, judgment_cost = 0, G = 2,
               seed = 1)
  expect_warning(res <- do.call(ss_aipe_crossed_sensitivity,
                                c(args, max_raters = 40)),
                 "search limit 'max_raters'")
  expect_equal(val(res, "n_raters"), 40)
  expect_error(do.call(ss_aipe_crossed_sensitivity,
                       c(args[setdiff(names(args), c("unit_cost",
                                                     "judgment_cost"))],
                         max_units = 20)),
               "a larger 'max_units'")
})

test_that("a well specified plan covers at the nominal rate and meets its width about half the time", {
  skip_on_cran()
  G <- 2000
  res <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    supposed_var_unit = 0.5, supposed_var_rater = 0.1,
    supposed_var_residual = 0.4, width = 0.4, G = G, seed = 2026)
  se <- sqrt(0.95 * 0.05 / G)
  expect_lt(abs(val(res, "total_type_I_error") - 0.05), 4 * se)
  # The plan targets the width at the planning values, so the realized
  # width is near it, on either side about equally often.
  expect_gt(val(res, "pct_ci_less_width"), 0.35)
  expect_lt(val(res, "pct_ci_less_width"), 0.75)
  expect_equal(val(res, "mean_ci_width"), 0.4, tolerance = 0.05)
  expect_equal(val(res, "sd_grand_mean"),
               sqrt(0.5 / val(res, "n_units") + 0.1 / val(res, "n_raters") +
                      0.4 / (val(res, "n_units") * val(res, "n_raters"))),
               tolerance = 0.1)

  # With an assurance of .80 at least about 80% of intervals meet it.
  assured <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.1, true_var_residual = 0.4,
    supposed_var_unit = 0.5, supposed_var_rater = 0.1,
    supposed_var_residual = 0.4, width = 0.4, assurance = 0.8, G = G,
    seed = 2026)
  expect_gt(val(assured, "pct_ci_less_width"), 0.8 - 4 * sqrt(0.16 / G))

  # Underestimating the rater component leaves the design short of raters.
  short <- ss_aipe_crossed_sensitivity(
    true_var_unit = 0.5, true_var_rater = 0.2, true_var_residual = 0.4,
    supposed_var_unit = 0.5, supposed_var_rater = 0.1,
    supposed_var_residual = 0.4, width = 0.4, G = 500, seed = 2026)
  expect_lt(val(short, "pct_ci_less_width"), 0.2)
})

test_that("with no rater variance the normal plan's two raters fall far short, and the t plan's do not", {
  # The help's account of quantile = "normal": with var_rater = 0 the
  # normal plan uses two raters (148 units here), the analysis keeps a
  # positive rater estimate in about a third of samples (MS_r / MS_e is F
  # on 1 and 147 df, P(F > 1) = .32), and the Satterthwaite df then fall
  # to about 2. The t plan sees the few df two raters give and uses three
  # (154 units).
  skip_on_cran()
  G <- 1000
  args <- list(true_var_unit = 0.3, true_var_rater = 0, true_var_residual = 0.6,
               supposed_var_unit = 0.3, supposed_var_rater = 0,
               supposed_var_residual = 0.6, width = 0.25, G = G,
               seed = 2026)
  nrm <- do.call(ss_aipe_crossed_sensitivity, c(args, quantile = "normal"))
  tq  <- do.call(ss_aipe_crossed_sensitivity, args)
  expect_equal(c(val(nrm, "n_units"), val(nrm, "n_raters")), c(148, 2))
  expect_equal(c(val(tq, "n_units"), val(tq, "n_raters")), c(154, 3))
  expect_lt(val(nrm, "pct_ci_less_width"), 0.40)
  expect_gt(val(nrm, "mean_ci_width"), 1.4 * 0.25)
  expect_gt(val(tq, "pct_ci_less_width"), 0.50)
  expect_lt(val(tq, "mean_ci_width"), 1.15 * 0.25)
})
