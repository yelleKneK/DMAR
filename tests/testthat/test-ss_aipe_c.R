test_that("ss_aipe_c() returns a 1-row tidy data frame with necessary_n_per_group", {
  result <- ss_aipe_c(true_var_error = 4, c_weights = c(1, -.5, -.5), width = 0.5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_n_per_group")
})

test_that("ss_aipe_c() requires a larger n for larger error variance", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  n_small <- ss_aipe_c(true_var_error = 1, c_weights = c(1, -.5, -.5), width = 0.5)$value
  n_large <- ss_aipe_c(true_var_error = 9, c_weights = c(1, -.5, -.5), width = 0.5)$value
  expect_lt(n_small, n_large)
})

test_that("ss_aipe_c() floors an impossibly wide target at the admissible minimum", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # An enormous target width was previously met at a degenerate n = 1, which
  # cannot estimate the within-group error variance (df = J * (n - 1) = 0). The
  # smallest admissible design is n = 2 per group.
  r <- ss_aipe_c(true_var_error = 40, c_weights = c(1, -.5, -.5), width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})

test_that("ss_aipe_c() validates its boundary inputs", {
  expect_error(ss_aipe_c(true_var_error = 40, c_weights = c(1, -.5, -.5),
                         width = 3, conf_level = 1.2), "conf_level")
  expect_error(ss_aipe_c(true_var_error = 40, c_weights = c(1, -.5, -.5),
                         width = 0), "width")
  expect_error(ss_aipe_c(true_var_error = 40, c_weights = c(1, -.5, -.5),
                         width = Inf), "width")
  expect_error(ss_aipe_c(true_var_error = -1, c_weights = c(1, -.5, -.5),
                         width = 3), "error variance")
  expect_error(ss_aipe_c(true_var_error = 40, c_weights = c(1, -.5, -.5),
                         width = 3, assurance = 1.5), "assurance")
})

## The kurtosis argument (1.1.0). The 1.0.0 search is reproduced below as
## a reference, verbatim in its arithmetic and in its iteration cap, so that
## the default true_kurtosis = 0 can be held to identity with what 1.0.0
## returned, plans and refusals alike; with `true_kurtosis` supplied, the
## reference computes the effective degrees of freedom in the form
## nu = 2 / (2 / df + kappa / (n J)) that the help states, to check the
## algebraically equal form the function computes.
.ss_aipe_c_reference_n <- function(error_variance, c_weights, width,
                                   conf_level, assurance, kurtosis = NULL) {
  alpha <- 1 - conf_level
  J <- length(c_weights)
  sigma <- sqrt(error_variance)
  n <- (sigma^2 * 4 * (qnorm(1 - alpha / 2))^2 * sum(c_weights^2)) / width^2
  tol <- 1e-6
  dif <- tol + 1
  iter <- 0L
  while (dif > tol) {
    iter <- iter + 1L
    if (iter > 10000L)
      stop("The sample size search did not converge within 10000 iterations; the specified target may be infeasible.", call. = FALSE)
    n_p <- n
    df <- max(n, 2) * J - J
    multiplier <- if (is.null(assurance)) {
      1
    } else if (is.null(kurtosis)) {
      qchisq(assurance, df) / df
    } else {
      nu <- 2 / (2 / df + kurtosis / (max(n, 2) * J))
      qchisq(assurance, nu) / nu
    }
    n <- if (is.null(assurance)) {
      (sigma^2 * 4 * (qt(1 - alpha / 2, df))^2 * sum(c_weights^2)) / width^2
    } else {
      ((sigma^2 * 4 * (qt(1 - alpha / 2, df))^2 * sum(c_weights^2)) / width^2) * multiplier
    }
    dif <- abs(n - n_p)
  }
  max(ceiling(n), 2L)
}

test_that("ss_aipe_c() with the default kurtosis reproduces the 1.0.0 plan of its example", {
  # 114 per group is the plan 1.0.0 printed for the first example of the
  # help page.
  expect_identical(
    ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5), width = 3,
              assurance = .90)$value,
    114
  )
  expect_identical(
    ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5), width = 3,
              assurance = .90, true_kurtosis = 0)$value,
    114
  )
})

test_that("ss_aipe_c() with true_kurtosis = 0 is identical to the 1.0.0 search on a grid", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # The plan, or the refusal message where the 1.0.0 search refused (a
  # target so wide that the search straddles the floor of 2 per group can
  # alternate without converging, in 1.0.0 as now).
  outcome <- function(expr) {
    tryCatch(expr, error = function(e) paste("error:", conditionMessage(e)))
  }
  weights <- list(c(1, -1), c(1, -0.5, -0.5), c(0.5, 0.5, -0.5, -0.5),
                  c(1, 0, 0, 0, -1))
  assurances <- list(NULL, 0.3, 0.6, 0.8, 0.9, 0.99)
  for (w in weights) {
    for (a in assurances) {
      for (ev in c(0.5, 4, 40)) {
        for (wd in c(0.5, 2, 6, 100)) {
          for (cl in c(0.90, 0.95, 0.99)) {
            reference <- outcome(.ss_aipe_c_reference_n(ev, w, wd, cl, a))
            label <- paste0("true_var_error = ", ev, ", c_weights = (",
                            paste(w, collapse = ", "), "), width = ", wd,
                            ", conf_level = ", cl, ", assurance = ",
                            if (is.null(a)) "NULL" else a)
            expect_identical(
              outcome(ss_aipe_c(true_var_error = ev, c_weights = w, width = wd,
                                conf_level = cl, assurance = a)$value),
              reference, label = label
            )
            expect_identical(
              outcome(ss_aipe_c(true_var_error = ev, c_weights = w, width = wd,
                                conf_level = cl, assurance = a,
                                true_kurtosis = 0)$value),
              reference, label = label
            )
          }
        }
      }
    }
  }
})

test_that("ss_aipe_c() computes the effective degrees of freedom the help states", {
  skip_on_cran()
  for (k in c(-1.9, -1.2, -0.5, 0.5, 1, 6)) {
    for (w in list(c(1, -1), c(1, -0.5, -0.5))) {
      for (a in c(0.6, 0.8, 0.95)) {
        expect_equal(
          ss_aipe_c(true_var_error = 4, c_weights = w, width = 1,
                    assurance = a, true_kurtosis = k)$value,
          .ss_aipe_c_reference_n(4, w, 1, 0.95, a, kurtosis = k),
          label = paste0("kurtosis = ", k, ", J = ", length(w),
                         ", assurance = ", a)
        )
      }
    }
  }
})

test_that("ss_aipe_c() plans more per group as the kurtosis of the errors grows", {
  kappas <- c(-1.99, -1.2, 0, 1, 3, 10)
  n_80 <- vapply(kappas, function(k)
    ss_aipe_c(true_var_error = 1, c_weights = c(1, -1), width = 0.8,
              assurance = 0.80, true_kurtosis = k)$value, numeric(1))
  expect_true(all(diff(n_80) >= 0))
  expect_lt(n_80[1], n_80[3])
  expect_gt(n_80[6], n_80[3])
  # Below an assurance of one half the chi square quantile falls under its
  # degrees of freedom, so a more variable S^2 makes the target easier to
  # meet at that probability and the order reverses.
  n_30 <- vapply(kappas, function(k)
    ss_aipe_c(true_var_error = 1, c_weights = c(1, -1), width = 0.8,
              assurance = 0.30, true_kurtosis = k)$value, numeric(1))
  expect_true(all(diff(n_30) <= 0))
})

test_that("ss_aipe_c() leaves the expected width plan untouched by kurtosis", {
  base <- ss_aipe_c(true_var_error = 4, c_weights = c(1, -0.5, -0.5),
                    width = 0.5)$value
  for (k in c(-1.5, 1, 20)) {
    expect_identical(
      ss_aipe_c(true_var_error = 4, c_weights = c(1, -0.5, -0.5),
                width = 0.5, true_kurtosis = k)$value,
      base
    )
  }
})

test_that("ss_aipe_c() validates kurtosis", {
  for (bad in list(-2, -3, NA_real_, Inf, c(0, 1), "1", NULL)) {
    expect_error(
      ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5), width = 3,
                assurance = .90, true_kurtosis = bad),
      "kurtosis", label = deparse(bad)
    )
  }
  # Validated without an assurance too, where it has no effect, so a
  # mistyped value is not silently ignored.
  expect_error(
    ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5), width = 3,
              true_kurtosis = -5),
    "kurtosis"
  )
})

test_that("a plan made with a nonzero kurtosis records it and says so under the table", {
  skip_on_cran()
  # The object 1.0.0 returned, built the way 1.0.0 built it: with the
  # default kurtosis, or without an assurance, no attribute and no line
  # are added.
  expected <- DMAR:::.as_dmar_tbl(
    data.frame(term = "necessary_n_per_group", value = 114),
    conf_level = 0.95, subclass = "dmar_ss_aipe")
  plan <- function(...) {
    ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5), width = 3,
              ...)
  }
  expect_identical(plan(assurance = .90), expected)
  expect_identical(plan(assurance = .90, true_kurtosis = 0), expected)
  expect_identical(plan(true_kurtosis = 1), plan())
  expect_null(attr(plan(true_kurtosis = 1), "true_kurtosis"))
  expect_false(any(grepl("kurtosis",
                         utils::capture.output(print(plan(assurance = .90))))))
  expect_false("true_kurtosis" %in% names(generics::glance(plan(assurance = .90))))

  for (k in c(-1.2, 1)) {
    out <- plan(assurance = .90, true_kurtosis = k)
    # The table keeps its one row and its class.
    expect_named(out, c("term", "value"))
    expect_identical(out$term, "necessary_n_per_group")
    expect_identical(class(out), class(expected))
    expect_identical(attr(out, "true_kurtosis"), k)
    expect_output(print(out), paste0(
      "The assurance allows for errors with an excess kurtosis of ",
      format(k), " (normal errors have 0)."), fixed = TRUE)
    gl <- generics::glance(out)
    expect_identical(gl$estimate, out$value)
    expect_identical(gl$true_kurtosis, k)
  }
  # The line shows the value to three significant digits; the attribute
  # keeps it at full precision.
  out <- plan(assurance = .90, true_kurtosis = -0.9045304)
  expect_identical(attr(out, "true_kurtosis"), -0.9045304)
  expect_output(print(out), "an excess kurtosis of -0.905 (normal", fixed = TRUE)
})
