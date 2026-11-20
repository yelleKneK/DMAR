# Sizes stored as R integers (50000L, or a count from length(), nrow(), or
# sum() of a logical vector) whose product passes R's largest integer,
# 2^31 - 1: 46,341 in each of two groups, or 40,000 and 60,000. R forms
# the product of two integers as an integer, which overflows to NA with
# the warning "NAs produced by integer overflow". DMAR 1.0.x did so in
# those of the functions below that it had, which then returned NA or, in
# most of them, stopped with an error. ci_two_phase(), ss_seq_two_phase(),
# and ss_aipe_two_phase_sensitivity() are new in 1.1.0, so their overflow
# was in code that no released DMAR had. The sizes are now multiplied as
# doubles. One test covers a single specified size above 2^31 - 1, which
# DMAR 1.0.x turned into NA with as.integer().

numbers <- function(v) {
  if (is.data.frame(v)) v <- v[vapply(v, is.numeric, logical(1))]
  unlist(v, use.names = FALSE)
}

test_that("integer group sizes whose product passes 2^31 - 1 give the same results as doubles, without warning", {
  calls <- list(
    function(a, b) ci_smd(smd = .2, n_1 = a, n_2 = b),
    function(a, b) ci_smd(ncp = 3, n_1 = a, n_2 = b),
    function(a, b) ci_smd_c(smd_c = .2, n_C = a, n_E = b),
    function(a, b) ci_smd_c(ncp = 3, n_C = a, n_E = b),
    function(a, b) ci_mahalanobis(D2 = .04, n_1 = a, n_2 = b, p = 3),
    function(a, b) convert_d_r(.2, n_1 = a, n_2 = b),
    function(a, b) convert_r_d(.2, n_1 = a, n_2 = b),
    function(a, b) convert_delta_lambda(.2, a, b),
    function(a, b) convert_lambda_delta(3, a, b),
    function(a, b) cles(.2, n_1 = a, n_2 = b),
    function(a, b) nnt_from_smd(.2, n_1 = a, n_2 = b),
    function(a, b) proportion_of_superiority(.2, n_1 = a, n_2 = b),
    function(a, b) meta_smd(c(.2, .3, .1), n_1 = rep(a, 3), n_2 = rep(b, 3)),
    function(a, b) ss_power_smd(.02, n_1 = a, n_2 = b),
    function(a, b) var_smd(.2, a, b),
    function(a, b) bayes_independent_t(mean_1 = 0, sd_1 = 1, n_1 = a,
                                       mean_2 = .02, sd_2 = 1, n_2 = b)
  )
  # Equal groups of 50,000 (product 2.5e9), and unequal groups of 40,000
  # and 60,000 (product 2.4e9), not both above 46,340.
  for (sizes in list(c(50000L, 50000L), c(40000L, 60000L))) {
    for (i in seq_along(calls)) {
      expect_no_warning(
        as_integer <- suppressMessages(calls[[i]](sizes[1], sizes[2])))
      as_double <- suppressMessages(calls[[i]](as.numeric(sizes[1]),
                                               as.numeric(sizes[2])))
      expect_false(anyNA(numbers(as_integer)))
      expect_equal(numbers(as_integer), numbers(as_double), tolerance = 1e-12,
                   label = paste("call", i, "at", sizes[1], "and", sizes[2]))
    }
  }
})

test_that("ci_smd_c() converts only integer sizes: factor sizes are refused, not read as their codes", {
  # as.numeric() of a factor gives its codes, so converting every size to a
  # double turned factor(50) with levels 10, 20, ..., 50 into 5 and
  # returned the interval for 5 and 6 with no warning. The factors reach
  # the arithmetic as they are, which warns that each operator is not
  # meaningful for factors; the NA it gives is refused (class
  # dmar_argument_invalid) for the argument of ci_nc_t() it reaches first,
  # 'ncp' when smd_c is given and 'df' when ncp is given. These are the
  # warnings and refusals of the same calls at 31fa5f5, before the sizes
  # were converted.
  fC <- factor(50, levels = c(10, 20, 30, 40, 50))
  fE <- factor(60, levels = c(10, 20, 30, 40, 50, 60))
  outcome <- function(expr) {
    w <- character()
    res <- withCallingHandlers(
      tryCatch(expr, error = function(e) e),
      warning = function(m) {
        w <<- c(w, conditionMessage(m))
        invokeRestart("muffleWarning")
      })
    list(result = res, operators = gsub("[^*+-]", "", w),
         all_factor = all(grepl("not meaningful for factors", w)))
  }
  cases <- list(
    list(o = outcome(ci_smd_c(smd_c = .2, n_C = fC, n_E = fE)),
         arg = "ncp", operators = c("-", "*", "+")),
    list(o = outcome(ci_smd_c(ncp = 1, n_C = fC, n_E = fE)),
         arg = "df", operators = c("+", "*", "-")),
    list(o = outcome(ci_smd_c(smd_c = .2, n_C = fC, n_E = 60)),
         arg = "ncp", operators = c("-", "*", "+")),
    list(o = outcome(ci_smd_c(smd_c = .2, n_C = fC, n_E = 60L)),
         arg = "ncp", operators = c("-", "*", "+")),
    list(o = outcome(ci_smd_c(smd_c = .2, n_C = 50, n_E = fE)),
         arg = "ncp", operators = c("*", "+")))
  for (k in seq_along(cases)) {
    o <- cases[[k]]$o
    expect_s3_class(o$result, "dmar_argument_invalid")
    expect_identical(o$result$arg, cases[[k]]$arg, label = paste("case", k))
    expect_identical(o$operators, cases[[k]]$operators, label = paste("case", k))
    expect_true(o$all_factor, label = paste("case", k))
  }
})

test_that("randomization_test() gives the SMD interval when the product of its group counts passes 2^31 - 1", {
  # Its sizes are counts of a logical vector, always integers, so DMAR 1.0.x
  # stopped with "missing value where TRUE/FALSE needed" for any data with
  # n_1 * n_2 above 2^31 - 1, whatever the type of the scores.
  set.seed(123)
  g1 <- rnorm(46341); g2 <- rnorm(46341, .02)
  expect_no_warning(
    res <- randomization_test(group_1 = g1, group_2 = g2, n_resamples = 20L,
                              seed = 1, shift_ci = FALSE))
  v <- setNames(res$value, res$term)
  d <- smd(group_1 = g1, group_2 = g2)$value
  ci <- ci_smd(smd = d, n_1 = 46341, n_2 = 46341)
  expect_equal(unname(v[c("smd_ci_lower", "smd", "smd_ci_upper")]),
               ci$value, tolerance = 1e-10)
  expect_true(all(is.finite(v[c("cles", "cles_ci_lower", "cles_ci_upper")])))
})

test_that("plot_smd() and ss_aipe_smd_sensitivity() pass integer sizes whose product passes 2^31 - 1 to ci_smd() without warning", {
  skip_if_not_installed("ggplot2")
  # plot_smd() counts the sizes of raw data with length(), so they are
  # integers whatever the type of the scores, and ss_aipe_smd_sensitivity()
  # uses its one size for both groups. Both pass the sizes to ci_smd(), and
  # DMAR 1.0.x stopped with "missing value where TRUE/FALSE needed" after
  # the overflow warning.
  set.seed(5)
  g1 <- rnorm(46341); g2 <- rnorm(46341, .02)
  expect_no_warning(p <- plot_smd(group_1 = g1, group_2 = g2))
  # The interval is drawn as the first segment, from ci_lower to ci_upper.
  seg <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)[[1]]
  d <- smd(group_1 = g1, group_2 = g2)$value[1]
  ci <- ci_smd(smd = d, n_1 = 46341, n_2 = 46341)$value
  expect_equal(c(seg$data$x, seg$data$xend), ci[c(1, 3)], tolerance = 1e-12)

  run <- function(n)
    ss_aipe_smd_sensitivity(true_smd = .1, specified_n_per_group = n,
                            width = .03, G = 3, seed = 1)
  expect_no_warning(as_integer <- run(50000L))
  expect_true(is.finite(as_integer$value[as_integer$term == "mean_ci_width"]))
  expect_equal(numbers(as_integer), numbers(run(50000)), tolerance = 1e-12)
})

test_that("the two-phase functions give their intervals when labeled times unlabeled passes 2^31 - 1", {
  # The sizes come from length(), so they are integers whatever the user
  # gives: 2,000 labeled units with 1.1 million unlabeled (product 2.2e9)
  # gave NA standard errors, an unclassed error, and a false refusal.
  set.seed(1)
  n <- 2000
  x <- rnorm(n); y <- .7 * x + sqrt(.51) * rnorm(n); xu <- rnorm(1100000)
  expect_no_warning(res <- ci_two_phase(y, x, xu))
  v <- setNames(res$value, res$term)
  # The standard error from its definition: the residual mean square over
  # n, plus the squared slope times the proxy's variance over the whole
  # pool, over the size of the pool.
  b <- stats::cov(y, x) / stats::var(x)
  e <- stats::resid(stats::lm(y ~ x))
  M <- n + length(xu)
  se <- sqrt(sum(e^2) / (n - 2) / n + b^2 * stats::var(c(x, xu)) / M)
  expect_equal(unname(v["se"]), se, tolerance = 1e-10)
  expect_true(all(is.finite(v[c("ci_lower", "ci_upper",
                                "se_ratio_target_only")])))

  expect_no_warning(sq <- ss_seq_two_phase(y, x, xu, half_width = .05))
  ci90 <- ci_two_phase(y, x, xu, conf_level = .90)
  half <- (ci90$value[ci90$term == "ci_upper"] -
             ci90$value[ci90$term == "ci_lower"]) / 2
  expect_equal(sq$value[sq$term == "half_width_current"], half,
               tolerance = 1e-12)

  # Integer sizes in the sensitivity function give the table of the same
  # call with doubles, rather than a refusal.
  run <- function(n_lab, pool)
    ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = .7,
                                  specified_n_labeled = n_lab, n_pool = pool,
                                  width = .05, G = 5, seed = 1)
  expect_no_warning(as_integer <- run(50000L, 100000L))
  expect_equal(numbers(as_integer), numbers(run(50000, 100000)),
               tolerance = 1e-12)
})

test_that("the whole-number size check returns a double, so a specified size above 2^31 - 1 is used", {
  # as.integer() turned a size above 2^31 - 1 into NA ("NAs introduced by
  # coercion to integer range"), and made a product of two checked sizes
  # an integer product.
  expect_type(DMAR:::.check_whole_n(5L, "n", 2), "double")
  expect_identical(DMAR:::.check_whole_n(3e9, "specified_N", 5), 3e9)
  expect_no_warning(res <- ss_power_reg_coef(true_f2 = 1e-9, p = 3,
                                             specified_N = 3e9))
  v <- setNames(res$value, res$term)
  expect_equal(unname(v["specified_N"]), 3e9)
  # The noncentrality parameter is sqrt(f2 N) = sqrt(3), and at so many
  # degrees of freedom the power of the two-sided test is the normal one.
  expect_equal(unname(v["noncentral_t_parm"]), sqrt(3), tolerance = 1e-12)
  z <- stats::qnorm(.975)
  expect_equal(unname(v["actual_power"]),
               1 - stats::pnorm(z - sqrt(3)) + stats::pnorm(-z - sqrt(3)),
               tolerance = 1e-6)
})

test_that("krippendorff_alpha() forms the expected disagreement when N (N - 1) passes 2^31 - 1", {
  # 23,171 units rated twice give N = 46,342 pairable values, and
  # 46,342 * 46,341 = 2,147,534,622. N is a length, so N * (N - 1L)
  # overflowed to NA, and the call stopped with "missing value where
  # TRUE/FALSE needed" after summing the distances over all 1.07e9 pairs
  # (about 20 minutes). The pair sum is replaced here by its value from the
  # frequencies of the values, exact for the nominal metric: the number of
  # pairs of unlike values, (N^2 - sum of squared frequencies) / 2.
  local_mocked_bindings(.kripp_pair_sum = function(all_vals, level) {
    f <- tabulate(all_vals)
    (length(all_vals)^2 - sum(f^2)) / 2
  })
  set.seed(7)
  units <- 23171L
  r1 <- sample(1:3, units, TRUE)
  r2 <- ifelse(runif(units) < .7, r1, sample(1:3, units, TRUE))
  expect_no_warning(res <- krippendorff_alpha(cbind(r1, r2)))
  v <- setNames(res$value, res$term)
  N <- 2 * units
  f <- tabulate(c(r1, r2))
  D_o <- 2 * sum(r1 != r2) / N
  D_e <- (N^2 - sum(f^2)) / (N * (N - 1))
  expect_equal(unname(v["n_pairable"]), N)
  expect_equal(unname(v["D_expected"]), D_e, tolerance = 1e-12)
  expect_equal(unname(v["krippendorff_alpha"]), 1 - D_o / D_e,
               tolerance = 1e-12)
})
