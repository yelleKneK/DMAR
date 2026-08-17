# Reference implementations used by the tests below. They are written
# straight from the definitions in the documentation, independently of the
# vectorized internals of irt_information(), so agreement is evidence about
# the implementation rather than a restatement of it.

# Category probabilities of one graded response item, P_ik(theta).
.grm_probabilities <- function(a, b, theta) {
  p_star <- c(1, stats::pnorm(a * (theta - b)), 0)
  p_star[-length(p_star)] - p_star[-1L]
}

# Analytic derivatives P'_ik(theta) from the formula in @details.
.grm_derivatives <- function(a, b, theta) {
  density <- c(0, stats::dnorm(a * (theta - b)), 0)
  a * (density[-length(density)] - density[-1L])
}

test_that("a dichotomous item reproduces the two parameter normal ogive", {
  a <- 1.37
  b <- 0.42
  theta <- seq(-4, 4, length.out = 81)
  res <- irt_information(a = a, b = b, theta = theta)

  # Closed form for the two parameter normal ogive. The complement is taken
  # from the upper tail, since 1 - pnorm(z) underflows to zero (and the
  # ratio to Inf) for z beyond about 8.3, which is inside this grid.
  z <- a * (theta - b)
  closed_form <- a^2 * stats::dnorm(z)^2 /
    (stats::pnorm(z) * stats::pnorm(z, lower.tail = FALSE))

  expect_equal(max(abs(res$test_information - closed_form)), 0, tolerance = 1e-10)
  expect_lt(max(abs(res$test_information - closed_form)), 1e-12)

  # The peak of a dichotomous item's information is at theta = b, and the
  # value there is a^2 * phi(0)^2 / 0.25 = a^2 / (2 * pi) * 4.
  at_b <- irt_information(a = a, b = b, theta = b)
  expect_equal(at_b$test_information, a^2 * stats::dnorm(0)^2 / 0.25,
               tolerance = 1e-12)
  expect_equal(attr(at_b, "theta_max_information"), b)

  # Several discriminations at once, still against the closed form.
  for (a_i in c(0.35, 1, 2.5)) {
    z_i <- a_i * (theta + 0.8)
    expect_equal(
      irt_information(a = a_i, b = -0.8, theta = theta)$test_information,
      a_i^2 * stats::dnorm(z_i)^2 /
        (stats::pnorm(z_i) * stats::pnorm(z_i, lower.tail = FALSE)),
      tolerance = 1e-12
    )
  }
})

test_that("the analytic category derivative matches a central difference", {
  a <- 1.15
  b <- c(-1.4, -0.3, 0.6, 1.8)                 # five categories
  h <- 1e-5
  for (theta in c(-2.5, -1, -0.25, 0, 0.7, 1.9, 3)) {
    numerical <- (.grm_probabilities(a, b, theta + h) -
                    .grm_probabilities(a, b, theta - h)) / (2 * h)
    expect_equal(.grm_derivatives(a, b, theta), numerical, tolerance = 1e-6)
  }

  # The information computed from the verified pieces matches the function.
  theta <- seq(-3, 3, by = 0.25)
  by_hand <- vapply(theta, function(t) {
    p  <- .grm_probabilities(a, b, t)
    dp <- .grm_derivatives(a, b, t)
    sum(dp[p > 0]^2 / p[p > 0])
  }, numeric(1))
  res <- irt_information(a = a, b = b, item = rep("x", length(b)),
                         theta = theta)
  expect_equal(res$test_information, by_hand, tolerance = 1e-12)
})

test_that("test information is the sum of the item information functions", {
  res <- irt_information(
    a = c(x1 = 1.4, x2 = 0.9, x3 = 1.1),
    b = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
    item = c(rep("x1", 4), "x2", "x3")
  )
  item_information <- attr(res, "item_information")

  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res), c("theta", "test_information", "se"))
  expect_equal(dim(item_information), c(nrow(res), 3L))
  expect_identical(colnames(item_information), c("x1", "x2", "x3"))
  expect_identical(rownames(item_information), as.character(res$theta))
  expect_identical(attr(res, "item"), c("x1", "x2", "x3"))
  expect_identical(attr(res, "a"), c(x1 = 1.4, x2 = 0.9, x3 = 1.1))
  expect_identical(unname(attr(res, "b")), c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8))
  expect_identical(names(attr(res, "b")),
                   c(rep("x1", 4), "x2", "x3"))

  # Additivity, to machine precision.
  expect_equal(res$test_information, unname(rowSums(item_information)),
               tolerance = 1e-15)

  # Each column is what the item on its own produces.
  for (nm in colnames(item_information)) {
    b_i <- attr(res, "b")[names(attr(res, "b")) == nm]
    alone <- irt_information(a = unname(attr(res, "a")[nm]), b = unname(b_i),
                             item = rep(nm, length(b_i)), theta = res$theta)
    expect_equal(unname(item_information[, nm]), alone$test_information,
                 tolerance = 1e-15)
  }

  # Two copies of an item carry exactly twice the information.
  one <- irt_information(a = 1.2, b = c(-1, 0.5), item = c("a", "a"))
  two <- irt_information(a = c(1.2, 1.2), b = c(-1, 0.5, -1, 0.5),
                         item = c("a", "a", "z", "z"))
  expect_equal(two$test_information, 2 * one$test_information,
               tolerance = 1e-15)

  # The peak is where test information is largest on the grid.
  expect_identical(attr(res, "theta_max_information"),
                   res$theta[which.max(res$test_information)])
})

test_that("the standard error is the reciprocal square root of information", {
  res <- irt_information(a = c(1.1, 0.7), b = c(-0.4, 1.2))
  expect_identical(res$se, 1 / sqrt(res$test_information))
  # The reverse map holds to floating point accuracy (the square root and
  # its inverse are not exact reciprocals in binary arithmetic).
  expect_equal(1 / res$se^2, res$test_information, tolerance = 1e-14)
})

test_that("extreme theta stays finite, nonnegative, and free of NaN", {
  res <- irt_information(
    a = c(x1 = 2.5, x2 = 0.8),
    b = c(-2, -1, 0, 1, 0.3),
    item = c(rep("x1", 4), "x2"),
    theta = c(-40, -25, -10, 0, 10, 25, 40)
  )
  expect_false(anyNA(res$test_information))
  expect_true(all(is.finite(res$test_information)))
  expect_true(all(res$test_information >= 0))
  expect_false(anyNA(attr(res, "item_information")))
  expect_true(all(is.finite(attr(res, "item_information"))))
  expect_true(all(attr(res, "item_information") >= 0))
  expect_false(anyNA(res$se))

  # Information is still strictly positive at theta = -10 and 10, and the
  # dichotomous closed form confirms the value there for the second item.
  info_10 <- res$test_information[res$theta %in% c(-10, 10)]
  expect_true(all(info_10 > 0))
  z <- 0.8 * (c(-10, 10) - 0.3)
  expect_equal(
    unname(attr(res, "item_information")[res$theta %in% c(-10, 10), "x2"]),
    0.8^2 * stats::dnorm(z)^2 /
      (stats::pnorm(z) * stats::pnorm(z, lower.tail = FALSE)),
    tolerance = 1e-10
  )

  # Where every category probability underflows the guard returns zero
  # information, not NaN, and the standard error is correspondingly Inf.
  far <- irt_information(a = 3, b = 0, theta = c(-1e3, 1e3))
  expect_identical(far$test_information, c(0, 0))
  expect_identical(far$se, c(Inf, Inf))
})

test_that("the discrimination may be given per item or per boundary", {
  b <- c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8)
  item <- c(rep("x1", 4), "x2", "x3")
  per_row <- irt_information(a = c(rep(1.4, 4), 0.9, 1.1), b = b, item = item)
  in_order <- irt_information(a = c(1.4, 0.9, 1.1), b = b, item = item)
  by_name <- irt_information(a = c(x3 = 1.1, x1 = 1.4, x2 = 0.9), b = b,
                             item = item)

  expect_equal(per_row$test_information, in_order$test_information,
               tolerance = 1e-15)
  expect_equal(per_row$test_information, by_name$test_information,
               tolerance = 1e-15)
  expect_identical(attr(by_name, "a"), c(x1 = 1.4, x2 = 0.9, x3 = 1.1))

  # Without 'item', each boundary is its own dichotomous item.
  default_items <- irt_information(a = c(1.1, 0.7), b = c(-0.4, 1.2))
  expect_identical(attr(default_items, "item"), c("item_1", "item_2"))
  expect_equal(
    default_items$test_information,
    irt_information(a = 1.1, b = -0.4)$test_information +
      irt_information(a = 0.7, b = 1.2)$test_information,
    tolerance = 1e-15
  )
})

test_that("irt_information() accepts an irt_grm() style parameter table", {
  # The contract of irt_grm(): one row per item and category boundary with
  # columns item, factor, category, lambda, tau, a, b.
  grm <- data.frame(
    item     = c(rep("x1", 4), "x2", "x3"),
    factor   = "f1",
    category = c(1:4, 1, 1),
    lambda   = c(rep(0.81, 4), 0.67, 0.74),
    tau      = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8) *
      c(rep(1.4, 4), 0.9, 1.1),
    a        = c(rep(1.4, 4), 0.9, 1.1),
    b        = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
    stringsAsFactors = FALSE
  )
  from_grm <- irt_information(grm = grm)
  direct <- irt_information(a = c(1.4, 0.9, 1.1), b = grm$b, item = grm$item)
  expect_equal(from_grm$test_information, direct$test_information,
               tolerance = 1e-15)
  expect_identical(attr(from_grm, "item"), c("x1", "x2", "x3"))

  # The category column orders the boundaries within an item.
  shuffled <- grm[c(3, 1, 4, 2, 5, 6), ]
  expect_equal(irt_information(grm = shuffled)$test_information,
               direct$test_information, tolerance = 1e-15)

  expect_error(irt_information(a = 1, b = 0, grm = grm), "not both")
  expect_error(irt_information(grm = grm[, c("item", "a")]), "missing the column")
  expect_error(irt_information(grm = as.matrix(1)), "data.frame")
})

test_that("irt_information() consumes a fitted irt_grm() result", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 400
  a_pop <- c(1.2, 0.9, 1.5)
  b_pop <- matrix(c(-1.2, -0.2, 0.9), nrow = 3, ncol = 3, byrow = TRUE)
  trait <- stats::rnorm(n)
  responses <- vapply(seq_along(a_pop), function(i) {
    p_star <- outer(trait, b_pop[i, ],
                    function(z, b) stats::pnorm(a_pop[i] * (z - b)))
    as.integer(1 + rowSums(stats::runif(n) < p_star))
  }, integer(n))
  colnames(responses) <- paste0("item", seq_along(a_pop))
  grm <- irt_grm(as.data.frame(responses))

  info <- irt_information(grm = grm)
  expect_s3_class(info, "dmar_tbl")
  expect_identical(attr(info, "item"), colnames(responses))
  expect_true(all(is.finite(info$test_information)))
  expect_true(all(info$test_information > 0))

  # The grm path is the direct path on the same numbers.
  a_by_item <- vapply(split(grm$a, factor(grm$item, levels = attr(info, "item"))),
                      function(z) z[1L], numeric(1))
  direct <- irt_information(a = unname(a_by_item), b = grm$b,
                            item = as.character(grm$item))
  expect_equal(info$test_information, direct$test_information,
               tolerance = 1e-15)
  expect_equal(unname(attr(info, "a")), unname(a_by_item), tolerance = 1e-15)
})

test_that("irt_information() fails loudly on inadmissible input", {
  expect_error(irt_information(), "Supply the discriminations")
  expect_error(irt_information(a = 1), "'b' must be supplied")
  expect_error(irt_information(a = 0, b = 0), "greater than zero")
  expect_error(irt_information(a = -1.2, b = 0), "greater than zero")
  expect_error(irt_information(a = NA_real_, b = 0), "'a' must be finite")
  expect_error(irt_information(a = Inf, b = 0), "'a' must be finite")
  expect_error(irt_information(a = 1, b = NA_real_), "'b' must be finite")
  expect_error(irt_information(a = 1, b = "0"), "'b' must be a numeric")
  expect_error(irt_information(a = c(1, 2), b = 0), "same length")
  expect_error(irt_information(a = c(1, 2, 3), b = c(0, 1), item = c("x", "x")),
               "one element per item")
  expect_error(
    irt_information(a = c(1, 2), b = c(0, 1), item = c("x", "x")),
    "constant within an item"
  )
  expect_error(irt_information(a = 1, b = 0, item = c("x", "y")), "same length")
  expect_error(irt_information(a = 1, b = 0, item = NA_character_),
               "must not contain missing")
  expect_error(irt_information(a = 1, b = 0, theta = numeric(0)),
               "'theta' must be a numeric")
  expect_error(irt_information(a = 1, b = 0, theta = c(0, NA)),
               "'theta' must be finite")
  expect_error(irt_information(a = 1, b = 0, theta = c(0, Inf)),
               "'theta' must be finite")

  # Unsorted boundaries within an item are sorted, with a warning.
  expect_warning(
    unsorted <- irt_information(a = 1.2, b = c(0.5, -1), item = c("x", "x")),
    "ascending order"
  )
  sorted <- irt_information(a = 1.2, b = c(-1, 0.5), item = c("x", "x"))
  expect_equal(unsorted$test_information, sorted$test_information,
               tolerance = 1e-15)
  expect_identical(unname(attr(unsorted, "b")), c(-1, 0.5))
})

test_that("the information curve agrees in shape with mirt's test information", {
  # mirt works in the logistic metric, so the normal ogive discrimination
  # is scaled by 1.702 for the comparison. The logistic and normal ogive
  # information functions are proportional in shape, not equal in value,
  # so the check is on the shape of the curve and the location of its peak.
  a <- c(1.2, 0.8, 1.5)
  b <- list(c(-1.2, -0.2, 0.9), c(-0.5, 0.7), 0.25)
  theta <- seq(-3, 3, length.out = 121)

  dmar <- irt_information(
    a = rep(a, lengths(b)),
    b = unlist(b),
    item = rep(paste0("x", seq_along(a)), lengths(b)),
    theta = theta
  )

  # Pinned from mirt::testinfo (mirt 1.46.1, 2026-08-09); live comparison in
  # tools/oracle_checks.R. The pinned curve is the logistic test information
  # of these items over the same theta grid, computed from a fixed-parameter
  # mirt model with slopes 1.702 * a and intercepts -1.702 * a * b.
  mirt_info <- c(0.1598172808079183, 0.1740424378422731, 0.1894980026980674,
                 0.2062769904234244, 0.2244761794850198, 0.244195485317077,
                 0.2655371175750638, 0.2886044862961097, 0.3135008214044266,
                 0.3403274709607002, 0.3691818469229229, 0.4001549937644814,
                 0.4333287659811103, 0.4687726162568516, 0.5065400177606056,
                 0.5466645724592821, 0.589155892876556, 0.633995387253644,
                 0.6811321266086169, 0.7304790246530944, 0.7819096144267839,
                 0.83525575377809, 0.8903066287492951, 0.9468094413965895,
                 1.004472157609633, 1.062968642232934, 1.12194641589505,
                 1.181037126351262, 1.239869637911621, 1.298085413530289,
                 1.355355611045397, 1.411399061222322, 1.46600007004552,
                 1.519024823280796, 1.570435097979784, 1.620298026377655,
                 1.668790823295772, 1.716199673541524, 1.762912358693199,
                 1.809404645734167, 1.856220915638089, 1.903949927261143,
                 1.953196944247599, 2.00455366504965, 2.058567470104156,
                 2.115711436166381, 2.176356384861196, 2.240745965706417,
                 2.308975469807154, 2.380974781818163, 2.456495657010676,
                 2.535103402057638, 2.616173071043201, 2.698890465537952,
                 2.782258518916988, 2.865109980835635, 2.946127590530447,
                 3.023873004117075, 3.096825487002849, 3.16343069922726,
                 3.222158768068609, 3.271569351203345, 3.310379765153891,
                 3.337530813750451, 3.352244069612792, 3.354064357835873,
                 3.342882233788972, 3.318933280196052, 3.282773776032985,
                 3.23523523020846, 3.177362875893321, 3.110345005573357,
                 3.035440698628321, 2.953913006469162, 2.866973205664467,
                 2.775739663365656, 2.681212601987271, 2.584263986995172,
                 2.485640172911115, 2.385973971602672, 2.285802465173161,
                 2.185587083106367, 2.085733048899431, 1.986606104164294,
                 1.888545279077328, 1.791871271031227, 1.696890634819537,
                 1.603896436853817, 1.513166278529018, 1.424958672476029,
                 1.339508698665378, 1.257023719884291, 1.177679741059456,
                 1.10161879042414, 1.028947509315456, 0.9597369781467312,
                 0.8940236866201924, 0.8318114773122461, 0.7730742493694646,
                 0.7177591964801121, 0.6657903628136451, 0.6170723247555194,
                 0.571493838500145, 0.5289313287446998, 0.4892521281167657,
                 0.4523174081384262, 0.4179847691314764, 0.3861104779832954,
                 0.3565513591836595, 0.3291663564409902, 0.3038177901239302,
                 0.2803723404588406, 0.2587017885450225, 0.2386835474672279,
                 0.2202010146419751, 0.2031437744869814, 0.1874076779137673,
                 0.1728948222884133, 0.1595134525867883, 0.1471778016296408,
                 0.1358078846123116)

  # Shape agreement across the grid, and a peak in the same place. The two
  # metrics cannot support a tolerance on the values themselves: the
  # logistic information function is flatter in the tails than the normal
  # ogive one, so the ratio of the two curves drifts across the grid even
  # though they describe the same items.
  expect_gt(stats::cor(dmar$test_information, mirt_info), 0.99)
  expect_lt(abs(theta[which.max(mirt_info)] -
                  attr(dmar, "theta_max_information")), 0.2)
})

test_that("plot_irt_information() returns a ggplot for both views", {
  skip_if_not_installed("ggplot2")
  info <- irt_information(
    a = c(x1 = 1.4, x2 = 0.9, x3 = 1.1),
    b = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
    item = c(rep("x1", 4), "x2", "x3")
  )

  p_test <- plot_irt_information(info)
  p_item <- plot_irt_information(info, what = "item")
  expect_s3_class(p_test, "ggplot")
  expect_s3_class(p_item, "ggplot")

  # The plots build without warning, which exercises the scales and layers.
  expect_silent(invisible(ggplot2::ggplot_build(p_test)))
  expect_silent(invisible(ggplot2::ggplot_build(p_item)))

  # The item view carries one line per item.
  item_layer <- ggplot2::ggplot_build(p_item)$data[[1]]
  expect_identical(length(unique(item_layer$group)), 3L)

  # Options that change what is drawn.
  expect_s3_class(plot_irt_information(info, show_se = FALSE,
                                       show_peak = FALSE), "ggplot")
  expect_s3_class(plot_irt_information(info, palette = "tableau"), "ggplot")

  labs_test <- p_test$labels
  expect_true(grepl("Latent trait", labs_test$x))
  expect_true(grepl("Test information", labs_test$y))

  expect_error(plot_irt_information(mtcars), "irt_information")
  expect_error(plot_irt_information(info, what = "items"), "arg")
  expect_error(plot_irt_information(info, show_se = NA), "single logical")
})

test_that("a logistic-metric grm table is converted before the ogive formulas", {
  skip_if_not_installed("lavaan")
  # irt_grm() records its discrimination metric on the returned table.
  # irt_information() computes normal ogive information, so a logistic
  # table must be divided by 1.702 first, the conversion its own Details
  # section describes. Without that step the same fit reported in the two
  # metrics gave test information differing by the square of the constant
  # (about 2.9 times), with SE(theta) understated by about a third, and no
  # error or warning marked the difference.
  set.seed(113)
  n <- 400; J <- 4
  eta <- rnorm(n)
  d <- as.data.frame(sapply(seq_len(J), function(j)
    as.integer(cut(0.8 * eta + rnorm(n, 0, .6),
                   breaks = c(-Inf, -0.8, 0, 0.8, Inf), labels = FALSE))))
  names(d) <- paste0("y", seq_len(J))

  g_ogive <- irt_grm(d)
  g_logit <- irt_grm(d, metric = "logistic")
  expect_identical(attr(g_ogive, "metric"), "normal_ogive")
  expect_identical(attr(g_logit, "metric"), "logistic")
  expect_equal(g_logit$a / g_ogive$a, rep(1.702, nrow(g_ogive)))

  theta <- c(-1, 0, 1)
  i_ogive <- irt_information(grm = g_ogive, theta = theta)
  i_logit <- irt_information(grm = g_logit, theta = theta)
  expect_equal(i_logit$test_information, i_ogive$test_information)
  expect_equal(i_logit$se_theta, i_ogive$se_theta)

  # A hand-built table carries no metric attribute and is still read as
  # normal ogive, so nothing changes for callers who never used irt_grm().
  plain <- as.data.frame(g_ogive[, c("item", "a", "b")])
  expect_equal(irt_information(grm = plain, theta = theta)$test_information,
               i_ogive$test_information)
})

test_that("a named 'a' that misses an item is an error, not positional matching", {
  # A single stale or mistyped name used to fall through to positional
  # matching, silently reassigning every discrimination in the order the
  # items first appear.
  expect_error(
    irt_information(a = c(y1 = 1, y2 = 1, y3 = 1, TYPO = 1),
                    b = rep(c(-1, 0, 1), 4),
                    item = rep(paste0("y", 1:4), each = 3),
                    theta = 0),
    "no element is named for the item", fixed = TRUE
  )
  # A correctly named 'a' still matches by name, and an unnamed one still
  # matches positionally.
  named <- irt_information(a = c(y1 = 1.2, y2 = 0.9),
                           b = rep(c(-1, 1), 2),
                           item = rep(c("y1", "y2"), each = 2), theta = 0)
  unnamed <- irt_information(a = c(1.2, 0.9),
                             b = rep(c(-1, 1), 2),
                             item = rep(c("y1", "y2"), each = 2), theta = 0)
  expect_equal(named$test_information, unnamed$test_information)
})
