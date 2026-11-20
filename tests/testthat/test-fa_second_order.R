# fa_second_order(): the second-order model fit by fa_confirmatory()'s
# engine. Every number is checked against lavaan::cfa(std.lv = TRUE) on
# the same syntax, against the Schmid-Leiman algebra written out here, and
# against fa_confirmatory() and fa_bifactor() where the models coincide or
# nest.

skip_if_not_installed("lavaan")

hs_four <- list(
  spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
  verbal  = c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning"),
  speed   = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"),
  memory  = c("t14_word_recognition", "t15_number_recognition",
              "t17_object_number"))
hs_second <- list(g = names(hs_four))

# The lavaan syntax of the same model, for lavaan::cfa(std.lv = TRUE),
# which fixes the variance of every second-order factor and of every
# first-order factor's disturbance at 1.
so_syntax <- function(factors, second_order) {
  paste(c(paste0(names(factors), " =~ ",
                 vapply(factors, paste, character(1L), collapse = " + ")),
          paste0(names(second_order), " =~ ",
                 vapply(second_order, paste, character(1L),
                        collapse = " + "))),
        collapse = "\n")
}

# The indices of a lavaan second-order fit with one second-order factor,
# by hand from its estimates: the Schmid-Leiman loadings lambda_i gamma_f
# (general) and lambda_i (group), then omega, omega_H, and omega_HS for
# sums of the raw items, ECV from the standardized Schmid-Leiman loadings,
# H with every other source of an item's variance as its error, the
# maximal reliability c' Sigma^-1 c with Sigma lavaan's model implied
# covariance matrix, and PUC from the structure.
so_hand <- function(fit, factors, g = "g") {
  pe <- lavaan::parameterEstimates(fit)
  items <- unlist(factors, use.names = FALSE)
  f_of <- rep(names(factors), lengths(factors))
  ld <- pe[pe$op == "=~" & pe$rhs %in% items, ]
  lambda <- ld$est[match(items, ld$rhs)]
  gm <- pe[pe$op == "=~" & pe$lhs == g, ]
  gamma <- gm$est[match(f_of, gm$rhs)]
  psi_rows <- pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items
  psi <- pe$est[psi_rows][match(items, pe$lhs[psi_rows])]
  lg <- lambda * gamma
  ls <- lambda
  v <- lg^2 + ls^2 + psi
  L_f <- tapply(lambda, factor(f_of, names(factors)), sum)
  total <- sum(lg)^2 + sum(L_f^2) + sum(psi)
  Sigma <- lavaan::lavInspect(fit, "implied")$cov[items, items]
  H <- function(l, other) { t <- sum(l^2 / other); t / (1 + t) }
  out <- c(omega = (sum(lg)^2 + sum(L_f^2)) / total,
           omega_H = sum(lg)^2 / total,
           ECV = sum(lg^2 / v) / sum((lg^2 + ls^2) / v),
           H = H(lg, v - lg^2),
           maximal = sum(lg * solve(Sigma, lg)),
           PUC = 1 - sum(choose(lengths(factors), 2)) /
             choose(length(items), 2))
  for (f in names(factors)) {
    r <- f_of == f
    sub <- L_f[[f]]^2 * (1 + gamma[r][1L]^2) + sum(psi[r])
    c_f <- ifelse(r, lambda, 0)
    out[paste0(c("omega_", "omega_HS_", "ECV_", "H_", "maximal_"), f)] <- c(
      L_f[[f]]^2 * (1 + gamma[r][1L]^2) / sub, L_f[[f]]^2 / sub,
      sum(lg[r]^2 / v[r]) / sum((lg[r]^2 + ls[r]^2) / v[r]),
      H(ls[r], v[r] - ls[r]^2),
      sum(c_f * solve(Sigma, c_f)))
  }
  out
}

est <- function(x, term) x$estimate[match(term, x$term)]

test_that("the estimates and the fit are lavaan's second-order model, under FIML and listwise deletion", {
  skip_on_cran()
  for (m in c("fiml", "listwise")) {
    x <- fa_second_order(holzinger_swineford, hs_four, hs_second, missing = m)
    ref <- lavaan::cfa(so_syntax(hs_four, hs_second),
                       data = holzinger_swineford, std.lv = TRUE,
                       missing = if (m == "fiml") "ml" else "listwise")
    pe <- lavaan::parameterEstimates(ref)
    for (f in names(hs_four)) {
      ld <- pe[pe$op == "=~" & pe$lhs == f, ]
      j <- seq_along(hs_four[[f]])
      expect_equal(est(x, paste0("lambda_", f, "_", j)), ld$est,
                   tolerance = 1e-5)
      expect_equal(x$se[match(paste0("lambda_", f, "_", j), x$term)], ld$se,
                   tolerance = 1e-4)
      psi <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% hs_four[[f]], ]
      expect_equal(est(x, paste0("psi_", f, "_", j)), psi$est,
                   tolerance = 1e-5)
      # The disturbance's variance is fixed at 1, the identification.
      expect_identical(est(x, paste0("zeta_", f)), 1)
    }
    gm <- pe[pe$op == "=~" & pe$lhs == "g", ]
    expect_equal(est(x, paste0("gamma_g_", names(hs_four))), gm$est,
                 tolerance = 1e-5)
    expect_equal(x$se[match(paste0("gamma_g_", names(hs_four)), x$term)],
                 gm$se, tolerance = 1e-4)
    expect_identical(est(x, "phi_g"), 1)
    expect_equal(est(x, c("chi_square", "df")),
                 unname(lavaan::fitMeasures(ref, c("chisq", "df"))),
                 tolerance = 1e-5)
    expect_identical(any(startsWith(x$term, "nu_")), m == "fiml")
    expect_identical(attr(x, "missing"), m)
    # The indices are the Schmid-Leiman algebra above, from lavaan's fit.
    h <- so_hand(ref, hs_four)
    expect_equal(est(x, c("omega_g", "omega_H_g", "ECV_g", "H_g",
                          "maximal_reliability_g", "PUC_g")),
                 unname(h[c("omega", "omega_H", "ECV", "H", "maximal",
                            "PUC")]), tolerance = 1e-5)
    for (f in names(hs_four)) {
      expect_equal(est(x, paste0(c("omega_", "omega_HS_", "ECV_", "H_",
                                   "maximal_reliability_"), f)),
                   unname(h[paste0(c("omega_", "omega_HS_", "ECV_", "H_",
                                     "maximal_"), f)]), tolerance = 1e-5)
    }
  }
})

test_that("the Schmid-Leiman attribute is lambda_std gamma_std and lambda_std sqrt(1 - gamma_std^2), and a first-order factor's ECV is gamma_std^2", {
  skip_on_cran()
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       missing = "listwise")
  ref <- lavaan::cfa(so_syntax(hs_four, hs_second),
                     data = holzinger_swineford, std.lv = TRUE)
  std <- lavaan::standardizedSolution(ref)
  sl <- attr(x, "schmid_leiman")
  expect_identical(names(sl), c("item", "factor", "second_order",
                                "general_loading", "group_loading",
                                "uniqueness"))
  items <- unlist(hs_four, use.names = FALSE)
  expect_identical(sl$item, items)
  expect_identical(sl$factor, rep(names(hs_four), lengths(hs_four)))
  expect_identical(sl$second_order, rep("g", length(items)))
  lam <- std[std$op == "=~" & std$rhs %in% items, ]
  lambda_std <- lam$est.std[match(items, lam$rhs)]
  gm <- std[std$op == "=~" & std$lhs == "g", ]
  gamma_std <- gm$est.std[match(sl$factor, gm$rhs)]
  expect_equal(sl$general_loading, lambda_std * gamma_std, tolerance = 1e-6)
  expect_equal(sl$group_loading, lambda_std * sqrt(1 - gamma_std^2),
               tolerance = 1e-6)
  # The uniqueness is the standardized error variance.
  psi <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, ]
  expect_equal(sl$uniqueness, psi$est.std[match(items, psi$lhs)],
               tolerance = 1e-6)
  expect_equal(sl$general_loading^2 + sl$group_loading^2 + sl$uniqueness,
               rep(1, length(items)))
  # Every item of a first-order factor divides its common variance in the
  # same proportion, so the factor's ECV is its squared standardized
  # second-order loading (the proportionality constraint), and the ratio
  # of an item's two loadings is the same for every item of the factor.
  expect_equal(est(x, paste0("ECV_", names(hs_four))),
               gm$est.std[match(names(hs_four), gm$rhs)]^2,
               tolerance = 1e-6)
  ratio <- sl$general_loading / sl$group_loading
  expect_equal(as.vector(tapply(ratio, sl$factor, stats::sd)),
               rep(0, 4L), tolerance = 1e-10)
  # A first-order factor's omega is the omega of its sum score, as
  # fa_confirmatory() writes it, at the second-order estimates: the squared
  # sum of its loadings times the factor's variance, 1 + gamma^2.
  lambda <- est(x, paste0("lambda_verbal_", 1:3))
  gamma <- est(x, "gamma_g_verbal")
  common <- sum(lambda)^2 * (1 + gamma^2)
  expect_equal(est(x, "omega_verbal"),
               common / (common + sum(est(x, paste0("psi_verbal_", 1:3)))),
               tolerance = 1e-8)
})

test_that("the second-order factor's share of each subscale score is the value semTools::compRelSEM() gives", {
  skip_on_cran()
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       output = "measurement")
  # compRelSEM() reads a second-order fit's first-order factors as
  # composites whose true score is the second-order factor, the
  # disturbances counted as error: the share of a subscale score's variance
  # due to g, which is omega_<f> less omega_HS_<f> here. Pinned from
  # semTools::compRelSEM(fit, obs.var = FALSE) on the fit of
  # output = "fit" (semTools 0.5.9, lavaan 0.7.3, 2026-11-05), for
  # spatial, verbal, speed, and memory; they are the hand algebra of
  # so_hand() on lavaan's own fit of the model, which agreed with them to
  # 1e-15, and the live comparison is
  # tools/oracle_snippets/semtools_second_order.R.
  pinned <- c(0.4263133166, 0.2504824422, 0.2305681256, 0.1234972521)
  f <- names(hs_four)
  expect_equal(est(x, paste0("omega_", f)) - est(x, paste0("omega_HS_", f)),
               pinned, tolerance = 1e-6)
  ref <- lavaan::cfa(so_syntax(hs_four, hs_second),
                     data = holzinger_swineford, std.lv = TRUE,
                     missing = "ml")
  h <- so_hand(ref, hs_four)
  expect_equal(unname(h[paste0("omega_", f)] - h[paste0("omega_HS_", f)]),
               pinned, tolerance = 1e-6)
})

test_that("three first-order factors under one second-order factor fit as fa_confirmatory()'s correlated factors model", {
  skip_on_cran()
  nine <- hs_four[c("spatial", "verbal", "speed")]
  x <- fa_second_order(holzinger_swineford, nine, list(g = names(nine)))
  y <- fa_confirmatory(holzinger_swineford, nine)
  # The structural part is just identified: the same chi square on the
  # same degrees of freedom (85.31 on 24, as the help says).
  expect_equal(est(x, c("chi_square", "df")), est(y, c("chi_square", "df")),
               tolerance = 1e-6)
  expect_equal(round(est(x, "chi_square"), 2), 85.31)
  expect_identical(est(x, "df"), 24)
  expect_identical(attr(x, "structural_df"), 0)
  lines <- .dmar_tbl_footer(x)
  expect_match(lines[names(lines) == "second_order_identified"],
               paste0("^With exactly three first-order factors under one ",
                      "second-order factor, the second-order structure is ",
                      "just identified"))
  expect_match(lines[names(lines) == "second_order_identified"],
               "positive product and imply standardized loadings below one")
  # The standardized second-order loadings reproduce the factor
  # correlations: phi_spatial_verbal = gamma*_spatial gamma*_verbal.
  g_std <- est(x, paste0("gamma_g_", names(nine)))
  g_std <- g_std / sqrt(1 + g_std^2)
  expect_equal(g_std[1L] * g_std[2L], est(y, "phi_spatial_verbal"),
               tolerance = 1e-5)
  # Four first-order factors leave the second-order structure 2 degrees of
  # freedom, and no such line.
  z <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       output = "measurement")
  expect_identical(attr(z, "structural_df"), 2)
  expect_false("second_order_identified" %in% names(.dmar_tbl_footer(z)))
  # Two first-order factors with equal second-order loadings are just
  # identified too, and the second-order factor is turned the way its
  # first-order factors point (the optimizer reaches -0.61 here).
  two <- list(verbal = c("t5_general_information", "t6_paragraph_comprehension",
                         "t7_sentence", "t8_word_classification",
                         "t9_word_meaning"),
              memory = c("t14_word_recognition", "t15_number_recognition",
                         "t16_figure_recognition", "t17_object_number"))
  w <- fa_second_order(holzinger_swineford, two, list(g = names(two)),
                       equal_second_order_loading = TRUE)
  expect_equal(est(w, c("chi_square", "df")),
               est(fa_confirmatory(holzinger_swineford, two),
                   c("chi_square", "df")), tolerance = 1e-6)
  expect_gt(est(w, "gamma_g"), 0)
  expect_true(all(attr(w, "schmid_leiman")$general_loading > 0))
  expect_identical(attr(w, "structural_df"), 0)
  expect_identical(attr(w, "equal")$second_order_loading, "g")
  lines <- .dmar_tbl_footer(w)
  expect_match(lines[names(lines) == "second_order_identified"],
               paste0("^With two first-order factors under one second-order ",
                      "factor and their loadings equated.*whenever that ",
                      "correlation is positive"))
})

test_that("two correlated second-order factors are lavaan's model, and each one's rows describe its own items", {
  skip_on_cran()
  # Two first-order factors under each, identified through the correlation
  # of the second-order factors (.82 here).
  so <- list(a = c("speed", "memory"), b = c("spatial", "verbal"))
  x <- fa_second_order(holzinger_swineford, hs_four, so, missing = "listwise")
  ref <- lavaan::cfa(so_syntax(hs_four, so), data = holzinger_swineford,
                     std.lv = TRUE)
  pe <- lavaan::parameterEstimates(ref)
  expect_equal(est(x, c("chi_square", "df")),
               unname(lavaan::fitMeasures(ref, c("chisq", "df"))),
               tolerance = 1e-5)
  expect_equal(est(x, "phi_a_b"),
               pe$est[pe$lhs == "a" & pe$op == "~~" & pe$rhs == "b"],
               tolerance = 1e-5)
  # The second-order structure: 6 correlations, 4 loadings, and one
  # correlation of the second-order factors.
  expect_identical(attr(x, "structural_df"), 1)
  # The rows of 'a' are those of a model of its own items, and its maximal
  # reliability counts the items of 'b' through the correlation.
  h <- so_hand(ref, hs_four[so$a], g = "a")
  expect_equal(est(x, c("omega_a", "omega_H_a", "ECV_a", "H_a", "PUC_a")),
               unname(h[c("omega", "omega_H", "ECV", "H", "PUC")]),
               tolerance = 1e-5)
  items <- unlist(hs_four, use.names = FALSE)
  Sigma <- lavaan::lavInspect(ref, "implied")$cov[items, items]
  cov_a <- lavaan::lavInspect(ref, "cov.all")[items, "a"]
  expect_equal(est(x, "maximal_reliability_a"),
               sum(cov_a * solve(Sigma, cov_a)), tolerance = 1e-5)
  # The table is ordered by second-order factor, each followed by its
  # first-order factors.
  idx <- x$term[startsWith(x$term, "omega_") & !startsWith(x$term, "omega_H")]
  expect_identical(idx, paste0("omega_", c("a", "speed", "memory", "b",
                                           "spatial", "verbal")))
  expect_identical(unique(attr(x, "schmid_leiman")$second_order),
                   c("b", "a"))
  expect_true(any(grepl("each second-order factor's rows describe the sum ",
                        .dmar_tbl_footer(x))))
})

test_that("the standard errors are the delta method's, and the limits are formed on the logit scale", {
  skip_on_cran()
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       missing = "listwise")
  fit <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                         missing = "listwise", output = "fit")
  # A delta method written here over lavaan's coef() and vcov(), by label:
  # it reproduces lavaan's standard error of omega_H, and then gives the
  # maximal reliability's, which lavaan cannot define.
  theta <- lavaan::coef(fit)
  V <- lavaan::vcov(fit)
  items <- unlist(hs_four, use.names = FALSE)
  f_of <- rep(names(hs_four), lengths(hs_four))
  j_of <- unlist(lapply(lengths(hs_four), seq_len), use.names = FALSE)
  index <- function(th) {
    lambda <- th[paste0("lambda_", f_of, "_", j_of)]
    psi <- th[paste0("psi_", f_of, "_", j_of)]
    gamma <- th[paste0("gamma_g_", names(hs_four))]
    L <- matrix(0, length(items), 4L)
    L[cbind(seq_along(items), match(f_of, names(hs_four)))] <- lambda
    Sigma <- L %*% (tcrossprod(gamma) + diag(4L)) %*% t(L) + diag(psi)
    lg <- lambda * gamma[match(f_of, names(hs_four))]
    total <- sum(lg)^2 + sum(tapply(lambda, f_of, sum)^2) + sum(psi)
    c(omega_H = sum(lg)^2 / total, maximal = sum(lg * solve(Sigma, lg)))
  }
  J <- vapply(seq_along(theta), function(i) {
    h <- 1e-6 * max(1, abs(theta[i]))
    up <- dn <- theta
    up[i] <- up[i] + h
    dn[i] <- dn[i] - h
    (index(up) - index(dn)) / (2 * h)
  }, numeric(2L))
  se <- sqrt(diag(J %*% V %*% t(J)))
  expect_equal(x$se[x$term == "omega_H_g"], se[["omega_H"]], tolerance = 1e-5)
  expect_equal(x$se[x$term == "maximal_reliability_g"], se[["maximal"]],
               tolerance = 1e-5)
  for (t in c("omega_H_g", "ECV_g", "H_g", "maximal_reliability_g",
              "omega_HS_verbal", "maximal_reliability_memory")) {
    r <- x[x$term == t, ]
    eta <- stats::qlogis(r$estimate)
    d <- stats::qnorm(0.975) * r$se / (r$estimate * (1 - r$estimate))
    expect_equal(c(r$ci_lower, r$ci_upper), stats::plogis(eta + c(-d, d)),
                 tolerance = 1e-8)
  }
  # PUC depends on the structure alone: 1 - 4 * 3 / 66.
  r <- x[x$term == "PUC_g", ]
  expect_equal(r$estimate, 1 - 12 / 66)
  expect_true(all(is.na(unlist(r[c("se", "test_statistic", "p_value", "ci_lower",
                                   "ci_upper")]))))
  y <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       missing = "listwise", wald_df = 20,
                       output = "measurement")
  r <- y[y$term == "maximal_reliability_g", ]
  d <- stats::qt(0.975, 20) * r$se / (r$estimate * (1 - r$estimate))
  expect_equal(c(r$ci_lower, r$ci_upper),
               stats::plogis(stats::qlogis(r$estimate) + c(-d, d)),
               tolerance = 1e-8)
  expect_identical(attr(y, "wald_df"), 20)
})

test_that("covariance input, robust ML, and incomplete data follow fa_confirmatory()", {
  skip_on_cran()
  items <- unlist(hs_four, use.names = FALSE)
  S <- stats::cov(holzinger_swineford[items])
  s_fit <- fa_second_order(S = S, N = 301, factors = hs_four,
                           second_order = hs_second)
  raw <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                         missing = "listwise")
  keep <- seq_len(match("chi_square", raw$term) - 1L)
  expect_equal(est(s_fit, raw$term[keep]), raw$estimate[keep],
               tolerance = 1e-5)
  expect_identical(attr(s_fit, "missing"), "none")
  mlr <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                         estimator = "MLR")
  ref <- lavaan::cfa(so_syntax(hs_four, hs_second),
                     data = holzinger_swineford, std.lv = TRUE,
                     estimator = "MLR", missing = "ml")
  gm <- lavaan::parameterEstimates(ref)
  gm <- gm[gm$op == "=~" & gm$lhs == "g", ]
  expect_equal(mlr$se[match(paste0("gamma_g_", names(hs_four)), mlr$term)],
               gm$se, tolerance = 1e-4)
  expect_identical(attr(mlr, "test"), "scaled")
  hs_inc <- holzinger_swineford
  set.seed(113)
  hs_inc$t7_sentence[sample(301, 60)] <- NA
  expect_no_warning(inc <- fa_second_order(hs_inc, hs_four, hs_second))
  ref <- lavaan::cfa(so_syntax(hs_four, hs_second), data = hs_inc,
                     std.lv = TRUE, missing = "ml")
  expect_equal(est(inc, "chi_square"),
               unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
  expect_identical(attr(inc, "missing"), "fiml")
  expect_error(fa_second_order(hs_inc, hs_four, hs_second, se = "robust.sem"),
               "no cases are deleted unasked")
})

test_that("equality constraints are written with shared labels", {
  skip_on_cran()
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       missing = "listwise",
                       equal_loading = c(spatial = FALSE, verbal = TRUE,
                                         speed = FALSE, memory = FALSE),
                       equal_error = c(spatial = FALSE, verbal = FALSE,
                                       speed = FALSE, memory = TRUE),
                       equal_second_order_loading = TRUE)
  model <- paste(c(
    paste0("spatial =~ ", paste(hs_four$spatial, collapse = " + ")),
    paste0("verbal =~ a*", paste(hs_four$verbal, collapse = " + a*")),
    paste0("speed =~ ", paste(hs_four$speed, collapse = " + ")),
    paste0("memory =~ ", paste(hs_four$memory, collapse = " + ")),
    paste0(hs_four$memory, " ~~ e*", hs_four$memory),
    paste0("g =~ b*", paste(names(hs_four), collapse = " + b*"))),
    collapse = "\n")
  ref <- lavaan::cfa(model, data = holzinger_swineford, std.lv = TRUE)
  expect_equal(est(x, c("chi_square", "df")),
               unname(lavaan::fitMeasures(ref, c("chisq", "df"))),
               tolerance = 1e-5)
  expect_true(all(c("lambda_verbal", "psi_memory", "gamma_g") %in% x$term))
  expect_identical(attr(x, "equal")$loading, "verbal")
  expect_identical(attr(x, "equal")$error, "memory")
  expect_identical(attr(x, "equal")$second_order_loading, "g")
  h <- so_hand(ref, hs_four)
  expect_equal(est(x, c("omega_H_g", "H_verbal", "omega_HS_memory")),
               unname(h[c("omega_H", "H_verbal", "omega_HS_memory")]),
               tolerance = 1e-5)
  # Equal second-order loadings leave the structure 5 degrees of freedom,
  # its 6 correlations less one loading.
  expect_identical(attr(x, "structural_df"), 5)
  expect_match(.dmar_tbl_footer(x)[["header"]], "Second-order model")
  expect_true(any(grepl("equal second-order loadings of g",
                        .dmar_tbl_footer(x))))
})

test_that("structures outside the model, not identified, or with names that clash are refused by class, with the reason", {
  hs <- holzinger_swineford
  refused <- function(expr, pattern, reason) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_second_order_unidentified")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$code, "second_order_unidentified")
    expect_identical(e$reason, reason)
    expect_match(conditionMessage(e), pattern)
    invisible(e)
  }
  # Two first-order factors under one second-order factor with free
  # loadings: one correlation for two loadings.
  e <- refused(fa_second_order(hs, hs_four[1:2], list(g = c("spatial",
                                                            "verbal"))),
               "do not identify its loadings.*equal_second_order_loading",
               "unidentified")
  expect_setequal(e$factors, c("g", "spatial", "verbal"))
  one <- hs_four
  one$memory <- "t14_word_recognition"
  e <- refused(fa_second_order(hs, one, hs_second), "'memory' has one item",
               "unidentified")
  expect_identical(e$factors, "memory")
  refused(fa_second_order(hs, hs_four,
                          list(g = c("spatial", "verbal", "speed"),
                               h = "memory")),
          "'h' has one \\('memory'\\)", "unidentified")
  both <- hs_four
  both$verbal <- c(both$verbal, "t2_cubes")
  e <- refused(fa_second_order(hs, both, hs_second),
               "t2_cubes is listed under spatial and verbal", "scope")
  expect_identical(e$items, "t2_cubes")
  # An item listed twice under one first-order factor is named as a slip
  # in the list.
  twice <- hs_four
  twice$spatial <- c(twice$spatial, "t2_cubes")
  e <- refused(fa_second_order(hs, twice, hs_second),
               "t2_cubes is listed more than once under spatial", "names")
  expect_identical(e$items, "t2_cubes")
  refused(fa_second_order(hs, hs_four,
                          list(g = names(hs_four), h = c("verbal", "speed"))),
          "'verbal' is listed under g and h", "scope")
  refused(fa_second_order(hs, hs_four, list(g = c(names(hs_four), "verbal"))),
          "'verbal' is listed twice under g", "names")
  refused(fa_second_order(hs, hs_four, list(g = c("spatial", "verbal",
                                                  "speed"))),
          "'memory' is listed under none", "names")
  refused(fa_second_order(hs, hs_four, list(g = c(names(hs_four), "math"))),
          "names 'math', which is not a first-order factor", "names")
  refused(fa_second_order(hs, hs_four, list(verbal = names(hs_four))),
          "first-order and a second-order factor are both named 'verbal'",
          "names")
  refused(fa_second_order(hs, hs_four, list(t2_cubes = names(hs_four))),
          "A factor and an item are both named 't2_cubes'", "names")
  # Labels made from the names coincide: factor 'a' numbers its first
  # loading lambda_a_1, which a factor 'a_1' with equal loadings also
  # writes.
  clash <- hs_four
  names(clash)[1:2] <- c("a", "a_1")
  refused(fa_second_order(hs, clash, list(g = names(clash)),
                          equal_loading = c(a = FALSE, a_1 = TRUE,
                                            speed = FALSE, memory = FALSE)),
          "labels that coincide \\(lambda_a_1\\)", "names")
  # The checks run before anything is fit.
  expect_error(fa_second_order(hs, hs_four, c("spatial", "verbal")),
               "'second_order' must be a named list")
})

test_that("what is not available is refused with the reason", {
  hs <- holzinger_swineford
  expect_error(fa_second_order(hs, hs_four, hs_second, ordered = TRUE),
               "ordered items are not available")
  expect_error(fa_second_order(hs, hs_four, hs_second, ord = TRUE),
               "ordered items are not available")
  expect_error(fa_second_order(hs, hs_four, hs_second, estimator = "WLSMV"),
               "must be one of .ML., .MLR.\\. fa_second_order\\(\\) fits")
  expect_error(fa_second_order(hs, hs_four, hs_second,
                               missing = "two_stage"),
               "'missing' must be one of")
  expect_error(fa_second_order(hs, hs_four, hs_second,
                               wald_df = "missing_information"),
               "not available for fa_second_order\\(\\): the calibration")
  expect_error(fa_second_order(hs, hs_four, hs_second, group = "school"),
               "'group' is not available in fa_second_order\\(\\)")
  expect_error(fa_second_order(hs, hs_four, hs_second,
                               equal_second_order_loading = c(h = TRUE)),
               "named with exactly the second-order factor names")
  expect_error(fa_second_order(hs, hs_four, hs_second,
                               equal_loading = c(g = TRUE)),
               "named with exactly the first-order factor names")
})

test_that("an item the input lacks is refused with the class of the factor analysis items", {
  hs <- holzinger_swineford
  misnamed <- hs_four
  misnamed$memory[3L] <- "t17_object_numbers"
  items_refusal <- function(expr, message) {
    e <- tryCatch(expr, dmar_items_invalid = function(e) e)
    expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
    expect_identical(e$code, "items_invalid")
    expect_identical(e$reason, "unknown")
    expect_identical(e$items, "t17_object_numbers")
    expect_identical(conditionMessage(e), message)
  }
  # The engine of fa_confirmatory() refuses the name, with the message of
  # each input that DMAR 1.0.x gave without a class.
  items_refusal(fa_second_order(hs, misnamed, hs_second),
                "Items not found in the data: t17_object_numbers.")
  S <- cov(hs[unlist(hs_four)], use = "complete.obs")
  items_refusal(fa_second_order(S = S, factors = misnamed,
                                second_order = hs_second, N = nrow(hs)),
                "Items not found in the covariance matrix: t17_object_numbers.")
})

test_that("an improper solution leaves the maximal reliabilities NA, with warnings", {
  skip_on_cran()
  # A covariance matrix whose first item has a standardized loading of
  # about 1.04 on its factor (r12 r13 / r23 = .7 * .7 / .45), so that its
  # error variance is estimated below zero.
  l1 <- sqrt(0.7 * 0.7 / 0.45)
  L <- matrix(0, 9L, 3L)
  L[1:3, 1L] <- c(l1, 0.7 / l1, 0.7 / l1)
  L[4:6, 2L] <- L[7:9, 3L] <- 0.7
  Phi <- matrix(0.5, 3L, 3L)
  diag(Phi) <- 1
  R <- L %*% Phi %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("x", 1:9), paste0("x", 1:9))
  f <- list(a = paste0("x", 1:3), b = paste0("x", 4:6), c = paste0("x", 7:9))
  seen <- new.env(parent = emptyenv())
  seen$msgs <- character()
  x <- withCallingHandlers(
    fa_second_order(S = R, N = 300, factors = f,
                    second_order = list(g = names(f))),
    warning = function(w) {
      seen$msgs <- c(seen$msgs, paste(class(w)[1L], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  expect_true(any(startsWith(seen$msgs, "dmar_heywood_warning")))
  expect_true(any(grepl("Improper solution: the error variance of x1",
                        seen$msgs)))
  expect_lt(attr(x, "improper")[["x1"]], 0)
  m <- x[startsWith(x$term, "maximal_reliability_"), ]
  expect_true(all(is.na(m$estimate)))
  expect_false(is.na(est(x, "omega_H_g")))
  expect_true(any(grepl("Improper solution", .dmar_tbl_footer(x))))
  expect_length(attr(x, "improper_second_order"), 0L)
})

# Records the warnings of a call, muffled, beside its value.
with_warnings <- function(expr) {
  seen <- new.env(parent = emptyenv())
  seen$msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen$msgs <- c(seen$msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen$msgs)
}

# The model of `factors` and `second_order` with the disturbance variance
# of each first-order factor of `at` fixed at zero and its second-order
# loading at 1, for lavaan::cfa(std.lv = TRUE): the face of the boundary
# on which those standardized second-order loadings are one, whose maximum
# is the second-order model's maximum there, on one more degree of freedom
# for each loading held.
so_face_syntax <- function(factors, second_order, at) {
  paste(c(paste0(names(factors), " =~ ",
                 vapply(factors, paste, character(1L), collapse = " + ")),
          paste0(names(second_order), " =~ ",
                 vapply(second_order, function(fs) {
                   paste(ifelse(fs %in% at, paste0("1*", fs), fs),
                         collapse = " + ")
                 }, character(1L))),
          paste0(at, " ~~ 0*", at)), collapse = "\n")
}

test_that("the fit is the maximum of the likelihood where the optimizer stops short on the way to a standardized second-order loading of one", {
  skip_on_cran()
  # The 1.1.0 review's case: one second-order factor over spatial, verbal,
  # speed, and t24, t25, and t26 (the items with real missingness), by
  # FIML. In this parameterization a standardized second-order loading of
  # one is reached only as gamma runs off to infinity, and nlminb stopped
  # on the way, at a chi square of 251.77 with three second-order loadings
  # from 43 to 233; the first version of fa_second_order() reported that
  # point, and so did lavaan::cfa(). The maximum puts extra's standardized
  # loading at one and the others at .74, .66, and .57.
  four <- c(hs_four[c("spatial", "verbal", "speed")],
            list(extra = c("t24_woody_mccall", "t25_paper_form_board_r",
                           "t26_flags")))
  so <- list(g = names(four))
  r <- with_warnings(fa_second_order(holzinger_swineford, four, so))
  x <- r$value
  face <- lavaan::cfa(so_face_syntax(four, so, "extra"),
                      data = holzinger_swineford, std.lv = TRUE,
                      missing = "ml")
  expect_equal(est(x, "chi_square"),
               unname(lavaan::fitMeasures(face, "chisq")), tolerance = 1e-6)
  expect_equal(round(est(x, "chi_square"), 2), 169.77)
  expect_identical(est(x, "df"),
                   as.numeric(lavaan::fitMeasures(face, "df")) - 1)
  expect_identical(names(attr(x, "improper_second_order")), "gamma_g_extra")
  expect_identical(unname(attr(x, "improper_second_order")), 1)
  expect_true(any(grepl(paste0(
    "^Improper solution at the second order: at the maximum of the ",
    "likelihood the standardized second-order loading of extra is one"),
    r$warnings)))
  # The other second-order loadings, the loadings, and the error variances
  # are the face's, with its standard errors, the loading held at one.
  pe <- lavaan::parameterEstimates(face)
  gm <- pe[pe$op == "=~" & pe$lhs == "g" & pe$rhs != "extra", ]
  expect_equal(est(x, paste0("gamma_g_", gm$rhs)), gm$est, tolerance = 1e-4)
  expect_equal(x$se[match(paste0("gamma_g_", gm$rhs), x$term)], gm$se,
               tolerance = 1e-4)
  ld <- pe[pe$op == "=~" & pe$lhs == "speed", ]
  expect_equal(est(x, paste0("lambda_speed_", 1:3)), ld$est, tolerance = 1e-4)
  ps <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% four$extra, ]
  expect_equal(est(x, paste0("psi_extra_", 1:3)), ps$est, tolerance = 1e-4)
  expect_equal(x$se[match(paste0("psi_extra_", 1:3), x$term)], ps$se,
               tolerance = 1e-4)
  # The rows of extra's loadings and of its disturbance have no finite
  # estimate; the subscale's omega and every other row are reported.
  na_terms <- c(paste0("lambda_extra_", 1:3), "gamma_g_extra",
                paste0(c("omega_HS_", "ECV_", "H_", "maximal_reliability_"),
                       "extra"))
  expect_true(all(is.na(unlist(x[x$term %in% na_terms,
                                 c("estimate", "se", "ci_lower",
                                   "ci_upper")]))))
  fit_terms <- c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
                 "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                 "rmsea_ci_level", "srmr", "AIC", "BIC", "H0", "H1")
  keep <- x[!x$term %in% c(na_terms, fit_terms, "PUC_g") &
              !startsWith(x$term, "zeta_") & x$term != "phi_g", ]
  expect_false(anyNA(keep$estimate) || anyNA(keep$ci_lower))
  expect_false(anyNA(x$estimate[x$term %in% fit_terms]))
  # The Schmid-Leiman group loadings of extra's items are zero there.
  sl <- attr(x, "schmid_leiman")
  expect_identical(sl$group_loading[sl$factor == "extra"], rep(0, 3L))
  expect_match(.dmar_tbl_footer(x)[["second_order_improper"]],
               "standardized second-order loading at one, gamma_g_extra")
  # lavaan's own fit of the parameterization does no better.
  ref <- suppressWarnings(lavaan::cfa(so_syntax(four, so),
                                      data = holzinger_swineford,
                                      std.lv = TRUE, missing = "ml"))
  expect_lte(est(x, "chi_square"),
             unname(lavaan::fitMeasures(ref, "chisq")) + 1e-6)
  # The lavaan fit holds the loading just short of one, within 1e-4 of the
  # maximum, and says so.
  r <- with_warnings(fa_second_order(holzinger_swineford, four, so,
                                     output = "fit"))
  expect_equal(as.numeric(lavaan::logLik(r$value)),
               as.numeric(lavaan::logLik(face)), tolerance = 1e-8)
  expect_true(any(grepl("The fit holds gamma_g_extra at 10000",
                        r$warnings, fixed = TRUE)))
  # Robust ML: the scaled test of the model at its maximum, whose scaling
  # factor lavaan computes as .985 there (computed at the point without
  # the optimizer, which leaves the bound unmarked, it came out .790).
  m <- suppressWarnings(fa_second_order(holzinger_swineford, four, so,
                                        estimator = "MLR"))
  expect_lt(abs(est(x, "chi_square") / est(m, "chi_square") - 1), 0.05)
})

test_that("a standardized second-order loading at one beside a correlation of the second-order factors beyond one", {
  skip_on_cran()
  # Five first-order factors, verbal and math under one second-order
  # factor and the other three under another. The first version stopped
  # with math's second-order loading near 31 and reported its other rows;
  # at the maximum math's standardized loading is one.
  five <- c(hs_four, list(math = c("t20_deduction", "t21_numerical_puzzles",
                                   "t22_problem_reasoning",
                                   "t23_series_completion")))
  so <- list(a = c("verbal", "math"), b = c("spatial", "speed", "memory"))
  r <- with_warnings(fa_second_order(holzinger_swineford, five, so,
                                     output = "measurement"))
  x <- r$value
  expect_true(any(grepl(paste0("^Improper solution at the second order: at ",
                               "the maximum of the likelihood the ",
                               "standardized second-order loading of math ",
                               "is one"), r$warnings)))
  bad <- attr(x, "improper_second_order")
  expect_identical(names(bad), "gamma_a_math")
  expect_identical(unname(bad), 1)
  expect_true(all(is.na(est(x, paste0(c("omega_HS_", "ECV_", "H_",
                                        "maximal_reliability_"), "math")))))
  keep <- x[!x$term %in% paste0(c("omega_HS_", "ECV_", "H_",
                                  "maximal_reliability_"), "math") &
              !startsWith(x$term, "PUC_"), ]
  expect_false(anyNA(keep$estimate) || anyNA(keep$ci_lower))
  expect_match(.dmar_tbl_footer(x)[["second_order_improper"]],
               "standardized second-order loading at one, gamma_a_math")
  face <- lavaan::cfa(so_face_syntax(five, so, "math"),
                      data = holzinger_swineford, std.lv = TRUE,
                      missing = "ml")
  v <- suppressWarnings(fa_second_order(holzinger_swineford, five, so))
  expect_equal(est(v, "chi_square"),
               unname(lavaan::fitMeasures(face, "chisq")), tolerance = 1e-6)
  # The split spatial, verbal | speed, memory, math, which the first
  # version reported at a chi square of 292.29 with memory at one as well:
  # at the maximum (215.33) math alone is at one, and the correlation of
  # the second-order factors is beyond one.
  so <- list(a = c("spatial", "verbal"), b = c("speed", "memory", "math"))
  r <- with_warnings(fa_second_order(holzinger_swineford, five, so))
  x <- r$value
  expect_identical(names(attr(x, "improper_second_order")),
                   c("gamma_b_math", "phi_a_b"))
  expect_equal(round(est(x, "chi_square"), 2), 215.33)
  expect_true(any(grepl(paste0("standardized second-order loading of math ",
                               "is one.* Also, the correlation of the ",
                               "second-order factors a and b \\(1\\.1[0-9]+",
                               "\\) is beyond one"), r$warnings)))
  expect_true(all(is.na(x$estimate[startsWith(x$term,
                                              "maximal_reliability_")])))
})

# A correlation matrix of nine items, three on each of three factors with
# loadings of .7, whose correlations, .9, .9, and .6, imply a standardized
# second-order loading of sqrt(.9 * .9 / .6) = 1.16 for the first: the
# maximum over the second-order model puts it at one.
so_at_one_R <- function() {
  L <- matrix(0, 9L, 3L)
  for (j in 1:3) L[3L * j - 2:0, j] <- 0.7
  Phi <- matrix(c(1, .9, .9, .9, 1, .6, .9, .6, 1), 3L)
  R <- L %*% Phi %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("x", 1:9), paste0("x", 1:9))
  R
}
so_at_one_f <- list(a = paste0("x", 1:3), b = paste0("x", 4:6),
                    c = paste0("x", 7:9))

test_that("a fit that does not converge because the data put a standardized second-order loading at one is reported at the boundary", {
  skip_on_cran()
  # The parameterization cannot reach the maximum, and nlminb did not
  # converge (the first version stopped with "The model did not
  # converge.").
  R <- so_at_one_R()
  f <- so_at_one_f
  so <- list(g = names(f))
  r <- with_warnings(fa_second_order(S = R, N = 300, factors = f,
                                     second_order = so))
  x <- r$value
  expect_identical(names(attr(x, "improper_second_order")), "gamma_g_a")
  expect_true(any(grepl("standardized second-order loading of a is one",
                        r$warnings)))
  face <- lavaan::cfa(so_face_syntax(f, so, "a"), sample.cov = R,
                      sample.nobs = 300, std.lv = TRUE)
  expect_equal(est(x, c("chi_square", "gamma_g_b", "gamma_g_c")),
               c(unname(lavaan::fitMeasures(face, "chisq")),
                 lavaan::parameterEstimates(face)$est[
                   lavaan::parameterEstimates(face)$op == "=~" &
                     lavaan::parameterEstimates(face)$lhs == "g"][2:3]),
               tolerance = 1e-5)
  # When no parameterization converges the error says which loading the
  # first fit was heading to one with.
  local_mocked_bindings(.fa_second_order_check_fit = function(args) NULL)
  e <- tryCatch(suppressWarnings(fa_second_order(S = R, N = 300, factors = f,
                                                 second_order = so)),
                error = function(e) e)
  expect_s3_class(e, "dmar_not_converged")
  expect_match(conditionMessage(e), paste0(
    "^The model did not converge, in fa_second_order\\(\\)'s ",
    "parameterization or in a second one.*loading of a \\((0\\.99[0-9]+|",
    "0\\.9999 or more)\\) was near one"))
})

test_that("the boundary is found with the call's equality constraints held, and its standard errors are those of the model there", {
  skip_on_cran()
  # The 1.1.0 review's case with equal loadings and equal error variances
  # within every first-order factor. lavaan 0.7.3 leaves equal labels
  # unenforced in a model with a bound, so the check fit from lavaan's own
  # starting values reached a larger log likelihood at a point where
  # spatial's error variances were 31, 14, and 52, and the first version
  # kept it, could not carry it back, and stopped with "The model did not
  # converge". With equal loadings alone it reported the boundary but set
  # every standard error NA: the singularity check restricted the
  # information to the constrained directions before scaling it, and so
  # mixed the held factor's loadings (near 1e-4) with the others.
  four <- c(hs_four[c("spatial", "verbal", "speed")],
            list(extra = c("t24_woody_mccall", "t25_paper_form_board_r",
                           "t26_flags")))
  so <- list(g = names(four))
  r <- with_warnings(fa_second_order(holzinger_swineford, four, so,
                                     equal_loading = TRUE,
                                     equal_error = TRUE))
  x <- r$value
  # The face of the boundary with the same equalities: extra's items on g,
  # with one loading, and one loading and one error variance for the
  # items of each factor.
  first <- c("spatial", "verbal", "speed")
  face <- lavaan::cfa(paste(c(
    paste0(first, " =~ ", vapply(first, function(f) {
      paste0("l_", f, "*", four[[f]], collapse = " + ")
    }, character(1L))),
    paste0("g =~ spatial + verbal + speed + ",
           paste0("l_extra*", four$extra, collapse = " + ")),
    unlist(lapply(names(four), function(f) {
      paste0(four[[f]], " ~~ e_", f, "*", four[[f]])
    }))), collapse = "\n"),
    data = holzinger_swineford, std.lv = TRUE, missing = "ml")
  expect_equal(est(x, "chi_square"),
               unname(lavaan::fitMeasures(face, "chisq")), tolerance = 1e-6)
  expect_identical(est(x, "df"),
                   as.numeric(lavaan::fitMeasures(face, "df")) - 1)
  expect_identical(names(attr(x, "improper_second_order")), "gamma_g_extra")
  pe <- lavaan::parameterEstimates(face)
  gm <- pe[pe$op == "=~" & pe$lhs == "g" & pe$rhs %in% first, ]
  expect_equal(est(x, paste0("gamma_g_", gm$rhs)), gm$est, tolerance = 1e-4)
  expect_equal(x$se[match(paste0("gamma_g_", gm$rhs), x$term)], gm$se,
               tolerance = 1e-4)
  e <- pe[pe$label %in% paste0("e_", names(four)) & !duplicated(pe$label), ]
  expect_equal(est(x, sub("^e_", "psi_", e$label)), e$est, tolerance = 1e-4)
  expect_equal(x$se[match(sub("^e_", "psi_", e$label), x$term)], e$se,
               tolerance = 1e-4)
  expect_false(any(grepl("information matrix of this fit is singular",
                         r$warnings)))
})

test_that("whether the point that holds a loading at one is the maximum does not depend on the number of cases", {
  skip_on_cran()
  # The log likelihood is a sum over the cases, and so is its gap at the
  # held point: 5e-7 on this fixture at N = 300 and 5e-4 at N = 300,000.
  # The first version compared it with an absolute 1e-4 and, from about
  # N = 300,000 here, stopped with "The model did not converge".
  R <- so_at_one_R()
  f <- so_at_one_f
  so <- list(g = names(f))
  for (N in c(3e5, 3e6)) {
    x <- suppressWarnings(fa_second_order(S = R, N = N, factors = f,
                                          second_order = so))
    face <- lavaan::cfa(so_face_syntax(f, so, "a"), sample.cov = R,
                        sample.nobs = N, std.lv = TRUE)
    expect_equal(est(x, "chi_square"),
                 unname(lavaan::fitMeasures(face, "chisq")),
                 tolerance = 1e-6, info = paste("N =", N))
    expect_identical(names(attr(x, "improper_second_order")), "gamma_g_a",
                     info = paste("N =", N))
  }
})

test_that("at the boundary lavaan's warnings about the data are given, and bootstrap standard errors are said to be absent", {
  skip_on_cran()
  R <- so_at_one_R()
  f <- so_at_one_f
  so <- list(g = names(f))
  set.seed(113)
  d <- as.data.frame(MASS::mvrnorm(300, rep(0, 9L), R, empirical = TRUE))
  # Three empty cases: the first version gave lavaan's warning about them
  # on a proper fit and dropped it at the boundary.
  d_empty <- d
  d_empty[1:3, ] <- NA
  r <- with_warnings(fa_second_order(d_empty, f, so))
  expect_identical(names(attr(r$value, "improper_second_order")),
                   "gamma_g_a")
  expect_true(any(grepl("some cases are empty and will be ignored",
                        r$warnings)))
  # The fit at the held point is computed without the optimizer, and a
  # bootstrap of it would not refit the samples; the first version left
  # every standard error NA without saying why, and its warning about the
  # loading said that they held it at one.
  r <- with_warnings(fa_second_order(d, f, so, missing = "listwise",
                                     se = "bootstrap", bootstrap = 5))
  x <- r$value
  expect_identical(names(attr(x, "improper_second_order")), "gamma_g_a")
  expect_true(all(is.na(x[c("se", "test_statistic", "p_value", "ci_lower",
                            "ci_upper")])))
  expect_true(any(startsWith(r$warnings, paste0(
    "Bootstrap standard errors and intervals are not available for this ",
    "fit: with a standardized second-order loading at one"))))
  expect_false(any(grepl("the other standard errors hold", r$warnings)))
  r <- with_warnings(fa_second_order(d, f, so, missing = "listwise"))
  expect_true(any(grepl("the other standard errors hold the loading at one",
                        r$warnings)))
  # The lavaan objects of the other outputs have no standard errors to be
  # NA (se = "none"); the first version said that every one was NA there
  # too.
  for (out in c("fit", "summary", "standardized")) {
    r <- with_warnings(fa_second_order(d, f, so, missing = "listwise",
                                       se = "bootstrap", bootstrap = 5,
                                       output = out))
    expect_true(any(grepl(paste0("so the ", out, " output carries no ",
                                 "standard errors, tests, or intervals"),
                          r$warnings, fixed = TRUE)), info = out)
    expect_false(any(grepl("is NA.", r$warnings, fixed = TRUE)), info = out)
  }
})

test_that("under FIML the summary at the boundary replaces lavaan's warnings about its robust fit indices with one", {
  skip_on_cran()
  # lavaan cannot form its FIML-corrected robust CFI, TLI, RMSEA, and GFI
  # for a fit computed without the optimizer; the first version passed on
  # its three warnings (one from each of lav_fit_cfi_lavobject(),
  # lav_fit_rmsea_lavobject(), and lav_fit_gfi_lavobject()).
  R <- so_at_one_R()
  f <- so_at_one_f
  so <- list(g = names(f))
  set.seed(113)
  d <- as.data.frame(MASS::mvrnorm(300, rep(0, 9L), R, empirical = TRUE))
  d[1:5, 1] <- NA
  for (e in c("ML", "MLR")) {
    r <- with_warnings(fa_second_order(d, f, so, estimator = e,
                                       output = "summary"))
    expect_false(any(grepl("lav_fit_", r$warnings, fixed = TRUE)), info = e)
    expect_identical(sum(grepl(paste0("lavaan could not compute its robust ",
                                      "CFI, TLI, RMSEA, and GFI"),
                               r$warnings, fixed = TRUE)), 1L, info = e)
    expect_true(is.na(r$value$fit[["cfi.robust"]]), info = e)
    # The summary's chi square is the table's: the scaled statistic under
    # MLR.
    x <- suppressWarnings(fa_second_order(d, f, so, estimator = e))
    chi <- if (e == "MLR") "chisq.scaled" else "chisq"
    expect_equal(unname(r$value$fit[[chi]]), est(x, "chi_square"),
                 tolerance = 1e-6, info = e)
  }
  # Off the boundary and with listwise deletion nothing is replaced.
  r <- with_warnings(fa_second_order(d, f, so, missing = "listwise",
                                     output = "summary"))
  expect_false(any(grepl("lavaan could not compute its robust", r$warnings,
                         fixed = TRUE)))
})

test_that("a fit whose maximum cannot be checked says so", {
  skip_on_cran()
  local_mocked_bindings(.fa_second_order_check_fit = function(args) NULL)
  r <- with_warnings(fa_second_order(holzinger_swineford, hs_four, hs_second,
                                     missing = "listwise"))
  expect_true(any(grepl(paste0("^fa_second_order\\(\\) could not confirm ",
                               "that this fit is the maximum of the ",
                               "likelihood"), r$warnings)))
  ref <- lavaan::cfa(so_syntax(hs_four, hs_second),
                     data = holzinger_swineford, std.lv = TRUE)
  expect_equal(est(r$value, "chi_square"),
               unname(lavaan::fitMeasures(ref, "chisq")), tolerance = 1e-5)
})

test_that("each first-order factor is turned so that its loadings sum to a positive value, its second-order loading changing sign with it", {
  skip_on_cran()
  # The verbal items scored in reverse: the model is the same, and verbal is
  # reported with positive loadings and a negative second-order loading.
  rev_hs <- holzinger_swineford
  rev_hs[hs_four$verbal] <- -rev_hs[hs_four$verbal]
  for (m in c("fiml", "listwise")) {
    x <- fa_second_order(rev_hs, hs_four, hs_second, missing = m)
    y <- fa_second_order(holzinger_swineford, hs_four, hs_second, missing = m)
    expect_true(all(est(x, paste0("lambda_verbal_", 1:3)) > 0), info = m)
    expect_equal(est(x, "gamma_g_verbal"), -est(y, "gamma_g_verbal"),
                 tolerance = 1e-4, info = m)
    expect_equal(est(x, paste0("gamma_g_", c("spatial", "speed", "memory"))),
                 est(y, paste0("gamma_g_", c("spatial", "speed", "memory"))),
                 tolerance = 1e-4, info = m)
    expect_equal(est(x, "chi_square"), est(y, "chi_square"), tolerance = 1e-8,
                 info = m)
  }
  # When the reflection cannot be applied the fit is reported as it stands,
  # with a warning: two first-order factors with equal second-order
  # loadings, whose second-order factor the optimizer leaves at -0.61.
  two <- list(verbal = c("t5_general_information", "t6_paragraph_comprehension",
                         "t7_sentence", "t8_word_classification",
                         "t9_word_meaning"),
              memory = c("t14_word_recognition", "t15_number_recognition",
                         "t16_figure_recognition", "t17_object_number"))
  local_mocked_bindings(.fa_second_order_refit = function(...) {
    list(fit = NULL, warnings = list())
  })
  r <- with_warnings(fa_second_order(holzinger_swineford, two,
                                     list(g = names(two)),
                                     equal_second_order_loading = TRUE))
  expect_lt(est(r$value, "gamma_g"), 0)
  expect_true(any(grepl("could not be turned to the reported direction",
                        r$warnings)))
})

test_that("first-order factors never observed together are refused when the pairs that are observed do not identify the second-order loadings", {
  skip_on_cran()
  # A two-form planned missing design: no case observes a spatial item
  # with a verbal item. Under one second-order factor over spatial, verbal,
  # and speed the two observed correlations cannot identify three
  # loadings; the first version fit such data and reported the values
  # lavaan left the loadings at, with the indices built on them.
  set.seed(3)
  form <- sample(1:2, nrow(holzinger_swineford), replace = TRUE)
  d <- holzinger_swineford
  d[form == 1, hs_four$spatial] <- NA
  d[form == 2, hs_four$verbal] <- NA
  three <- hs_four[c("spatial", "verbal", "speed")]
  e <- tryCatch(fa_second_order(d, three, list(g = names(three))),
                error = function(e) e)
  expect_s3_class(e, "dmar_second_order_unidentified")
  expect_identical(e$reason, "unidentified")
  expect_setequal(e$factors, c("spatial", "verbal"))
  expect_match(conditionMessage(e), paste0(
    "^No case observes an item of 'spatial' with an item of 'verbal', and ",
    "the pairs of first-order factors the data do observe together do not ",
    "identify the second-order structure"))
  # Two under one second-order factor with equal loadings, the same.
  e <- tryCatch(fa_second_order(d, three[1:2], list(g = names(three)[1:2]),
                                equal_second_order_loading = TRUE),
                error = function(e) e)
  expect_s3_class(e, "dmar_second_order_unidentified")
  # A fourth first-order factor leaves five observed correlations for four
  # loadings, and the model is fit.
  x <- suppressWarnings(fa_second_order(d, hs_four, hs_second))
  expect_true(all(is.finite(est(x, paste0("gamma_g_", names(hs_four))))))
})

test_that("a singular or non-invertible information matrix leaves every standard error, test, and interval NA, the maximal reliability's included", {
  skip_on_cran()
  # Two second-order factors over two first-order factors each are
  # identified only through the correlation of the second-order factors;
  # in data where it is exactly zero each pair's two loadings have one
  # correlation between them. Under FIML lavaan inverts the singular
  # observed information; listwise and from a covariance matrix it cannot
  # invert the expected information.
  L <- matrix(0, 12L, 4L)
  for (j in 1:4) L[3L * j - 2:0, j] <- 0.7
  Phi <- diag(4L)
  Phi[1, 2] <- Phi[2, 1] <- Phi[3, 4] <- Phi[4, 3] <- 0.49
  R <- L %*% Phi %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("x", 1:12), paste0("x", 1:12))
  set.seed(113)
  d <- as.data.frame(MASS::mvrnorm(300, rep(0, 12L), R, empirical = TRUE))
  f <- list(A = paste0("x", 1:3), B = paste0("x", 4:6), C = paste0("x", 7:9),
            D = paste0("x", 10:12))
  so <- list(g = c("A", "B"), h = c("C", "D"))
  cols <- c("se", "test_statistic", "p_value", "ci_lower", "ci_upper")
  for (input in c("fiml", "listwise", "S")) {
    r <- with_warnings(if (input == "S") {
      fa_second_order(S = stats::cov(d), N = 300, factors = f,
                      second_order = so, output = "measurement")
    } else {
      fa_second_order(d, f, so, missing = input, output = "measurement")
    })
    x <- r$value
    expect_true(all(is.na(unlist(x[cols]))), info = input)
    expect_true(all(is.finite(x$estimate[startsWith(x$term,
                                                    "maximal_reliability_")])),
                info = input)
    expect_true(any(grepl(paste0("^(The information matrix of this fit is ",
                                 "singular|lavaan could not compute the ",
                                 "standard errors)"), r$warnings)),
                info = input)
  }
})

test_that("a second-order correlation beyond one leaves every maximal reliability NA, with a warning", {
  skip_on_cran()
  so <- list(a = c("spatial", "speed"), b = c("verbal", "memory"))
  r <- with_warnings(fa_second_order(holzinger_swineford, hs_four, so,
                                     missing = "listwise",
                                     output = "measurement"))
  x <- r$value
  expect_true(any(grepl(paste0("correlation of the second-order factors a ",
                               "and b \\(1\\.01[0-9]\\) is beyond one"),
                        r$warnings)))
  bad <- attr(x, "improper_second_order")
  expect_identical(names(bad), "phi_a_b")
  expect_gt(bad[["phi_a_b"]], 1)
  m <- x[startsWith(x$term, "maximal_reliability_"), ]
  expect_true(all(is.na(m$estimate)) && all(is.na(m$se)))
  # The other rows describe each second-order factor's own items and do not
  # read the correlation; they keep their estimates and intervals.
  keep <- x[!startsWith(x$term, "maximal_reliability_") &
              !startsWith(x$term, "PUC_"), ]
  expect_false(anyNA(keep$estimate) || anyNA(keep$ci_lower))
  expect_match(.dmar_tbl_footer(x)[["second_order_improper"]],
               "a second-order correlation beyond one, phi_a_b")
})

test_that("the likelihood ratio test against the bifactor model has as many degrees of freedom as items less first-order factors", {
  skip_on_cran()
  so_fit <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                            output = "fit")
  bf_fit <- fa_bifactor(holzinger_swineford, hs_four, output = "fit")
  lrt <- lavaan::lavTestLRT(so_fit, bf_fit)
  expect_identical(as.numeric(lrt[["Df diff"]][2L]), 12 - 4)
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second)
  y <- fa_bifactor(holzinger_swineford, hs_four)
  expect_equal(lrt[["Chisq diff"]][2L],
               est(x, "chi_square") - est(y, "chi_square"), tolerance = 1e-6)
})

test_that("the table's class, header, tidy(), glance(), the other outputs, and metadata", {
  skip_on_cran()
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second)
  expect_s3_class(x, c("dmar_fa_second_order", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(names(x), c("syntax", "term", "estimate", "se", "test_statistic",
                               "p_value", "ci_lower", "ci_upper"))
  idx <- c("omega_g", "omega_H_g", "ECV_g", "H_g", "maximal_reliability_g",
           "PUC_g", "omega_spatial", "omega_HS_spatial", "ECV_spatial",
           "H_spatial", "maximal_reliability_spatial")
  expect_identical(x$term[match("omega_g", x$term) + 0:10], idx)
  expect_identical(utils::tail(x$term, 1L), "H1")
  lines <- .dmar_tbl_footer(x)
  expect_match(lines[names(lines) == "header"][1L],
               "Second-order model: g over spatial \\(3 items\\), verbal")
  expect_true(any(grepl("Schmid-Leiman loadings", lines)))
  td <- generics::tidy(x)
  expect_false("chi_square" %in% td$term)
  expect_true(all(c("gamma_g_spatial", "maximal_reliability_g") %in% td$term))
  gl <- generics::glance(x)
  expect_identical(nrow(gl), 1L)
  expect_equal(gl$chi_square, est(x, "chi_square"))
  expect_identical(gl$missing, "fiml")
  md <- dmar_metadata(x)
  expect_identical(md$result, "dmar_fa_second_order")
  expect_identical(md$estimator, "ML")
  expect_identical(attr(x, "omega_metric"), "raw_item_sum")
  m <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       output = "measurement")
  expect_s3_class(m, "dmar_fa_second_order")
  expect_identical(m$term, x$term[x$term %in% m$term])
  expect_s3_class(attr(m, "schmid_leiman"), "data.frame")
  # The standardized solution gives the correlations of the first-order
  # factors with their second-order factor, and no index rows.
  std <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                         output = "standardized")
  expect_false(any(std$op == ":="))
  g <- std[std$op == "=~" & std$lhs == "g", ]
  gamma <- est(x, paste0("gamma_g_", names(hs_four)))
  expect_equal(g$est.std, gamma / sqrt(1 + gamma^2), tolerance = 1e-6)
  sm <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                        output = "summary")
  expect_s3_class(sm, "lavaan.summary")
})

test_that("the RMSEA interval is at rmsea_conf_level, apart from the intervals of the parameters and indices", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three second-order fits
  rows <- c("rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level")
  x <- fa_second_order(holzinger_swineford, hs_four, hs_second)
  y <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       rmsea_conf_level = 0.95)
  z <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                       conf_level = 0.99)
  fit <- fa_second_order(holzinger_swineford, hs_four, hs_second,
                         output = "fit")
  fm <- function(level) {
    unname(lavaan::fitMeasures(fit, c("rmsea.ci.lower", "rmsea.ci.upper"),
                               fm.args = list(rmsea.ci.level = level)))
  }
  expect_equal(eval(formals(fa_second_order)$rmsea_conf_level), 0.90)
  expect_equal(est(x, rows), c(fm(0.90), 0.90), tolerance = 1e-10)
  expect_equal(est(y, rows), c(fm(0.95), 0.95), tolerance = 1e-10)
  expect_identical(est(z, rows), est(x, rows))
  par <- !x$term %in% c(rows, "rmsea")
  expect_identical(y$ci_lower[par], x$ci_lower[par])
  expect_false(isTRUE(all.equal(z$ci_lower[par], x$ci_lower[par])))
  expect_true("Confidence level: 95% (90% for the RMSEA interval)" %in%
                utils::capture.output(print(x)))
})

test_that("the standardized output passes on no lavaan warning about the index rows it leaves out", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two second-order fits
  # Equal loadings or error variances give several items one label, and
  # lavaan warns that the standardized values of the defined indices that
  # use it are ambiguous; the standardized output leaves those rows out,
  # and through 1.1.0 development the warning reached the caller anyway.
  for (eq in c("equal_loading", "equal_error")) {
    args <- list(holzinger_swineford, hs_four, hs_second,
                 output = "standardized")
    args[[eq]] <- TRUE
    expect_no_warning(z <- do.call(fa_second_order, args))
    expect_false(any(z$op %in% c(":=", "==")))
    expect_true(all(is.finite(z$se[z$op == "=~"])))
  }
})
