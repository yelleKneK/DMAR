test_that("dmacs() reduces to the standardized intercept difference", {
  # Equal loadings, intercepts differing by exactly d, one pooled SD: the
  # integral collapses to d^2 and dMACS must be abs(d) / s exactly.
  d <- c(0.30, -0.45, 1.25, 0)
  s <- 1.37
  res <- dmacs(lambda_reference = rep(0.82, 4), lambda_focal = rep(0.82, 4),
               nu_reference = c(2.0, 2.0, 2.0, 2.0),
               nu_focal = c(2.0, 2.0, 2.0, 2.0) - d,
               sd_pooled = s, mean_focal = -0.4, sd_focal = 1.6)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res),
                   c("item", "lambda_reference", "lambda_focal",
                     "nu_reference", "nu_focal", "sd_pooled", "dmacs"))
  expect_equal(res$dmacs, abs(d) / s, tolerance = 1e-15)
  expect_true(all(res$dmacs >= 0))
  expect_identical(res$item, paste0("item_", 1:4))
  expect_identical(attr(res, "reference"), "reference")
  expect_identical(attr(res, "focal"), "focal")
  expect_identical(attr(res, "mean_focal"), -0.4)
  expect_identical(attr(res, "sd_focal"), 1.6)
})

test_that("dmacs() is exactly zero under full invariance", {
  set.seed(113)
  lambda <- runif(6, 0.4, 1.2)
  nu <- runif(6, -2, 3)
  res <- dmacs(lambda_reference = lambda, lambda_focal = lambda,
               nu_reference = nu, nu_focal = nu,
               sd_pooled = runif(6, 0.8, 2.5),
               mean_focal = 0.35, sd_focal = 0.9,
               item_names = paste0("v", 1:6))
  expect_identical(res$dmacs, rep(0, 6))
  expect_identical(res$item, paste0("v", 1:6))
})

test_that("the closed form equals numerical integration of the definition", {
  set.seed(4831)
  for (rep in seq_len(25L)) {
    lam_r <- stats::runif(1, 0.2, 1.5)
    lam_f <- stats::runif(1, 0.2, 1.5)
    nu_r  <- stats::runif(1, -3, 3)
    nu_f  <- stats::runif(1, -3, 3)
    mu    <- stats::runif(1, -2, 2)
    sg    <- stats::runif(1, 0.3, 2.0)
    sp    <- stats::runif(1, 0.5, 3.0)

    res <- dmacs(lambda_reference = lam_r, lambda_focal = lam_f,
                 nu_reference = nu_r, nu_focal = nu_f,
                 sd_pooled = sp, mean_focal = mu, sd_focal = sg)

    integrand <- function(eta) {
      ((nu_r + lam_r * eta) - (nu_f + lam_f * eta))^2 *
        stats::dnorm(eta, mean = mu, sd = sg)
    }
    numeric_integral <- stats::integrate(
      integrand, lower = mu - 15 * sg, upper = mu + 15 * sg,
      rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
    )$value

    expect_equal(res$dmacs, sqrt(numeric_integral) / sp, tolerance = 1e-8)
  }
})

test_that("dmacs() validates its inputs and fails loudly", {
  ok <- list(lambda_reference = c(1, 1), lambda_focal = c(1, 1),
             nu_reference = c(0, 0), nu_focal = c(0, 0), sd_pooled = 1)

  expect_error(do.call(dmacs, utils::modifyList(ok, list(lambda_focal = 1))),
               "same length")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(nu_focal = c(0, 0, 0)))),
               "same length")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(sd_pooled = c(1, 1, 1)))),
               "length 1 or length 2")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(sd_pooled = 0))),
               "must be positive")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(sd_focal = 0))),
               "single positive number")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(nu_focal = c(0, NA)))),
               "must be finite")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(mean_focal = c(0, 1)))),
               "single number")
  expect_error(do.call(dmacs, utils::modifyList(ok, list(item_names = "a"))),
               "one label per item")
  expect_error(dmacs(), "Supply 'fit'")
  expect_error(dmacs(lambda_reference = c(1, 1), lambda_focal = c(1, 1)),
               "Supply 'fit'")

  # Names that disagree across groups are a misalignment, not a relabeling.
  expect_error(
    dmacs(lambda_reference = c(a = 1, b = 1), lambda_focal = c(b = 1, a = 1),
          nu_reference = c(a = 0, b = 0), nu_focal = c(a = 0, b = 1),
          sd_pooled = 1),
    "do not align"
  )

  # No silent NA: a legitimate input set returns finite, nonnegative values.
  res <- do.call(dmacs, utils::modifyList(ok, list(nu_focal = c(0.5, -0.5))))
  expect_true(all(is.finite(res$dmacs)))
})


# ---- Fitted multiple group lavaan models -----------------------------------

.dmacs_spatial_items <- function() {
  c("t1_visual_perception", "t2_cubes", "t3_paper_form_board", "t4_lozenges")
}

.dmacs_partial_fit <- function(data) {
  items <- .dmacs_spatial_items()
  model <- paste("spatial =~", paste(items, collapse = " + "))
  lavaan::cfa(model, data = data, group = "school",
              group.equal = c("loadings", "intercepts"),
              group.partial = c("t2_cubes ~ 1", "t4_lozenges ~ 1"))
}

test_that("dmacs() reproduces the fit's parameters and pooled SDs by hand", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, envir = environment())
  items <- .dmacs_spatial_items()
  fit <- .dmacs_partial_fit(holzinger_swineford)
  res <- dmacs(fit)

  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$item, items)
  expect_true(all(is.finite(res$dmacs)))
  expect_true(all(res$dmacs >= 0))

  labs <- as.character(lavaan::lavInspect(fit, "group.label"))
  expect_identical(attr(res, "reference"), labs[1])
  expect_identical(attr(res, "focal"), labs[2])

  # Hand computation from lavInspect.
  est <- lavaan::lavInspect(fit, "est")
  lam_r <- est[[1]]$lambda[items, 1]
  lam_f <- est[[2]]$lambda[items, 1]
  nu_r  <- est[[1]]$nu[items, 1]
  nu_f  <- est[[2]]$nu[items, 1]
  mu    <- est[[2]]$alpha[1, 1]
  sg    <- sqrt(est[[2]]$psi[1, 1])
  expect_equal(res$lambda_reference, unname(lam_r), tolerance = 1e-12)
  expect_equal(res$lambda_focal, unname(lam_f), tolerance = 1e-12)
  expect_equal(res$nu_reference, unname(nu_r), tolerance = 1e-12)
  expect_equal(res$nu_focal, unname(nu_f), tolerance = 1e-12)
  expect_equal(unname(attr(res, "mean_focal")), mu, tolerance = 1e-12)
  expect_equal(unname(attr(res, "sd_focal")), sg, tolerance = 1e-12)

  # Pooled SD from the raw data, computed independently.
  s <- split(holzinger_swineford[, items], holzinger_swineford$school)
  n_r <- vapply(s[[labs[1]]], function(x) sum(!is.na(x)), numeric(1))
  n_f <- vapply(s[[labs[2]]], function(x) sum(!is.na(x)), numeric(1))
  v_r <- vapply(s[[labs[1]]], stats::var, numeric(1), na.rm = TRUE)
  v_f <- vapply(s[[labs[2]]], stats::var, numeric(1), na.rm = TRUE)
  sd_pooled <- sqrt(((n_r - 1) * v_r + (n_f - 1) * v_f) / (n_r + n_f - 2))
  expect_equal(res$sd_pooled, unname(sd_pooled), tolerance = 1e-10)

  a <- unname(nu_r - nu_f)
  b <- unname(lam_r - lam_f)
  expect_equal(res$dmacs,
               sqrt(a^2 + 2 * a * b * mu + b^2 * (sg^2 + mu^2)) /
                 unname(sd_pooled),
               tolerance = 1e-12)

  # The two anchored items carry equality constraints, so their dMACS is 0.
  expect_equal(res$dmacs[c(1, 3)], c(0, 0), tolerance = 1e-10)
})

test_that("a deliberately shifted item recovers shift / pooled SD", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, envir = environment())
  items <- .dmacs_spatial_items()
  shift <- 20

  fit0 <- .dmacs_partial_fit(holzinger_swineford)
  base <- dmacs(fit0)
  labs <- as.character(lavaan::lavInspect(fit0, "group.label"))

  shifted <- holzinger_swineford
  focal_rows <- shifted$school == labs[2]
  shifted$t2_cubes[focal_rows] <- shifted$t2_cubes[focal_rows] + shift
  fit1 <- .dmacs_partial_fit(shifted)
  res <- dmacs(fit1)

  j <- match("t2_cubes", res$item)
  others <- res$dmacs[-j]

  # Adding a constant to one group's item raises only that intercept, and
  # leaves every variance alone, so the pooled SD is unchanged.
  expect_equal(res$sd_pooled, base$sd_pooled, tolerance = 1e-10)
  expect_equal(res$lambda_reference, base$lambda_reference, tolerance = 1e-6)
  expect_equal(res$lambda_focal, base$lambda_focal, tolerance = 1e-6)
  expect_equal(res$nu_focal[j], base$nu_focal[j] + shift, tolerance = 1e-5)

  # Only the intercepts were freed, so b = 0 and dMACS is exactly abs(a) / SD.
  a_base <- base$nu_reference[j] - base$nu_focal[j]
  expect_equal(res$dmacs[j], abs(a_base - shift) / res$sd_pooled[j],
               tolerance = 1e-5)

  # The recovered value is close to the constant divided by the pooled SD,
  # and its distance from that target is exactly the noninvariance the item
  # already carried before the constant was added.
  expect_equal(res$dmacs[j], shift / res$sd_pooled[j], tolerance = 0.06)
  expect_equal(abs(res$dmacs[j] - shift / res$sd_pooled[j]), base$dmacs[j],
               tolerance = 1e-5)

  # It dwarfs the untouched items, one of which (the lozenges test) carries
  # real noninvariance of its own.
  expect_gt(res$dmacs[j], 5 * max(others))
  expect_gt(max(others), 0.1)
  expect_true(all(is.finite(res$dmacs)))
})

test_that("dmacs() resolves groups and errors clearly on more than two", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, envir = environment())
  items <- .dmacs_spatial_items()
  fit <- .dmacs_partial_fit(holzinger_swineford)
  labs <- as.character(lavaan::lavInspect(fit, "group.label"))

  by_label <- dmacs(fit, reference = labs[1], focal = labs[2])
  by_index <- dmacs(fit, reference = 1, focal = 2)
  expect_equal(by_label$dmacs, by_index$dmacs, tolerance = 1e-12)
  expect_equal(by_label$dmacs, dmacs(fit)$dmacs, tolerance = 1e-12)

  # With equal loadings the index is symmetric in the two roles, and the
  # reported components swap.
  swapped <- dmacs(fit, reference = labs[2], focal = labs[1])
  expect_equal(swapped$dmacs, by_label$dmacs, tolerance = 1e-8)
  expect_equal(swapped$nu_reference, by_label$nu_focal, tolerance = 1e-12)
  expect_identical(attr(swapped, "reference"), labs[2])

  expect_error(dmacs(fit, reference = "Nowhere"), "must name one of")
  expect_error(dmacs(fit, reference = 9), "must index one of")
  expect_error(dmacs(fit, reference = 1, focal = 1), "must be different")
  expect_error(dmacs(fit, lambda_reference = 1), "not both")
  expect_error(dmacs(fit, mean_focal = 1), "read from 'fit'")
  expect_error(dmacs(fit = stats::lm(1 ~ 1)), "fitted lavaan object")

  # Four groups (school by grade) require the user to name the two compared.
  four <- holzinger_swineford
  four$school <- droplevels(interaction(four$school, four$grade, sep = "_"))
  model <- paste("spatial =~", paste(items, collapse = " + "))
  fit4 <- lavaan::cfa(model, data = four, group = "school",
                      group.equal = c("loadings", "intercepts"))
  expect_error(dmacs(fit4), "Name the two to compare")
  expect_error(dmacs(fit4, reference = 1), "Name both")
  expect_true(all(is.finite(dmacs(fit4, reference = 1, focal = 3)$dmacs)))
})

test_that("dmacs() pools from sample moments when there is no raw data", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, envir = environment())
  items <- .dmacs_spatial_items()
  model <- paste("spatial =~", paste(items, collapse = " + "))

  parts <- split(holzinger_swineford[, items], holzinger_swineford$school)
  moment_fit <- lavaan::cfa(
    model,
    sample.cov  = lapply(parts, stats::cov),
    sample.mean = lapply(parts, colMeans),
    sample.nobs = lapply(parts, nrow),
    group.equal = c("loadings", "intercepts"),
    group.partial = c("t2_cubes ~ 1", "t4_lozenges ~ 1")
  )
  from_moments <- dmacs(moment_fit)
  from_data <- dmacs(.dmacs_partial_fit(holzinger_swineford))

  # Same pooled standard deviations either way, and with the loadings
  # constrained equal the index does not depend on which group is reference.
  expect_equal(from_moments$sd_pooled, from_data$sd_pooled, tolerance = 1e-10)
  expect_equal(from_moments$dmacs, from_data$dmacs, tolerance = 1e-6)
})

test_that("dmacs() requires a multiple group model with a mean structure", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, envir = environment())
  items <- .dmacs_spatial_items()
  model <- paste("spatial =~", paste(items, collapse = " + "))

  one_group <- lavaan::cfa(model, data = holzinger_swineford)
  expect_error(dmacs(one_group), "multiple group model")

  no_means <- lavaan::cfa(model, data = holzinger_swineford, group = "school",
                          meanstructure = FALSE)
  expect_error(dmacs(no_means), "no mean structure")
})
