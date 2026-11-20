test_that("pairwise_within() wide-format returns documented columns", {
  set.seed(113)
  Y <- matrix(rnorm(40), 10, 4) +
       matrix(rep(c(0, 0.3, 0.6, 0.9), each = 10), 10, 4)
  colnames(Y) <- paste0("T", 1:4)
  res <- pairwise_within(Y)
  expect_setequal(colnames(res),
                  c("contrast", "mean_difference", "sd_difference",
                    "t_value", "df", "p_value", "p_adjusted",
                    "ci_lower", "ci_upper", "n_pairs"))
  expect_equal(nrow(res), 6)  # choose(4, 2)
})

test_that("pairwise_within() Holm adjusts p-values upward", {
  set.seed(113)
  Y <- matrix(rnorm(50), 10, 5)
  colnames(Y) <- paste0("T", 1:5)
  res <- pairwise_within(Y, adjust = "holm")
  expect_true(all(res$p_adjusted >= res$p_value - 1e-12))
})

test_that("pairwise_within() long-format equals wide-format", {
  set.seed(113)
  n <- 8; k <- 3
  Y <- matrix(rnorm(n * k), n, k); colnames(Y) <- paste0("T", 1:k)
  long <- data.frame(subject = factor(rep(1:n, k)),
                     time    = factor(rep(paste0("T", 1:k), each = n)),
                     y       = as.vector(Y))
  wide_res <- pairwise_within(Y)
  long_res <- pairwise_within(long, subject = "subject",
                              condition = "time", outcome = "y")
  expect_equal(sort(wide_res$mean_difference),
               sort(long_res$mean_difference), tolerance = 1e-10)
})

test_that("pairwise_within() returns a DMAR table that records its family confidence level", {
  # The 1.1.0 review: through 1.0.1 the table was a plain data.frame that
  # recorded no level, so the console printed none and results_sentence()
  # wrote "95% CI" whatever conf_level was; the first versions of 1.1.0
  # recorded, under bonferroni_ci = TRUE, the level of each interval
  # ("98.33333% CI"), where contrast_test() and the post hoc intervals
  # record the simultaneous level of the family.
  set.seed(1)
  Y <- matrix(rnorm(60), 20) + rep(c(0, .3, .6), each = 20)
  colnames(Y) <- c("A", "B", "C")
  pw <- pairwise_within(Y)
  expect_s3_class(pw, c("dmar_tbl", "data.frame"), exact = TRUE)
  expect_identical(attr(pw, "conf_level"), 0.95)
  expect_null(attr(pw, "simultaneous"))
  out <- utils::capture.output(print(pw))
  expect_true("Confidence level: 95%" %in% out)
  expect_identical(results_sentence(pw[2, ]),
                   "mean_difference = 0.55, 95% CI [0.07, 1.03], unadjusted p = 0.0272")
  p90 <- pairwise_within(Y, conf_level = 0.90)
  expect_match(results_sentence(p90[2, ]), "90% CI", fixed = TRUE)
  # Bonferroni intervals: the family level is recorded, and the line
  # beneath the table names the adjustment and the level of each interval.
  pb <- pairwise_within(Y, bonferroni_ci = TRUE)
  expect_identical(attr(pb, "conf_level"), 0.95)
  expect_equal(attr(pb, "simultaneous")$each, 1 - 0.05 / 3,
               tolerance = 1e-12)
  out_b <- utils::capture.output(print(pb))
  expect_true(paste0("Confidence level: 95%, simultaneous across the 3 ",
                     "comparisons (Bonferroni: each interval at 98.33%)")
              %in% out_b)
  expect_identical(results_sentence(pb[2, ]),
                   "mean_difference = 0.55, 95% CI [-0.05, 1.15]; no matching multiplicity-adjusted p-value is available, simultaneous across the 3 comparisons (Bonferroni: each interval at 98.33%)")
  # The limits are those of the paired t interval at 1 - .05 / 3.
  d <- Y[, "C"] - Y[, "A"]
  expect_equal(c(pb$ci_lower[2], pb$ci_upper[2]),
               mean(d) + c(-1, 1) * stats::qt(1 - 0.05 / 6, 19) *
                 stats::sd(d) / sqrt(20), tolerance = 1e-12)
  # The same line, and the same recorded level, as contrast_test() gives
  # its Bonferroni intervals.
  g <- factor(rep(c("a", "b", "c"), each = 10))
  set.seed(2)
  v <- rnorm(30) + as.numeric(g)
  ct <- contrast_test(aov(v ~ g), adjust = "bonferroni")
  expect_identical(attr(ct, "conf_level"), 0.95)
  expect_true(any(grepl(paste0("simultaneous across the 3 comparisons ",
                               "\\(Bonferroni: each interval at 98\\.33%\\)"),
                        utils::capture.output(print(ct)))))
  # Two conditions make one comparison, which no adjustment changes.
  expect_null(attr(pairwise_within(Y[, 1:2], bonferroni_ci = TRUE),
                   "simultaneous"))
})

test_that("the level of each Bonferroni interval never reads 100%", {
  # With 8 conditions (m = 28 comparisons) at a family level of .999,
  # each interval is at 1 - .001 / 28 = .99996429. Four significant
  # digits rounded that to "100", a level only an unbounded interval has,
  # beneath finite limits; the line now adds digits until the level reads
  # below 100, and the stored level keeps full precision.
  set.seed(1)
  Y8 <- matrix(rnorm(160), 20)
  pw <- pairwise_within(Y8, bonferroni_ci = TRUE, conf_level = .999)
  expect_identical(attr(pw, "simultaneous")$each, 1 - .001 / 28)
  out <- utils::capture.output(print(pw))
  expect_true(paste0("Confidence level: 99.9%, simultaneous across the 28 ",
                     "comparisons (Bonferroni: each interval at 99.996%)")
              %in% out)
  expect_false(any(grepl("at 100%", out, fixed = TRUE)))
  expect_true(all(is.finite(c(pw$ci_lower, pw$ci_upper))))
  # The caption of a knitted table is the same line.
  cap <- utils::capture.output(print(as_kable(pw, format = "pipe")))
  expect_true(any(grepl("each interval at 99.996%", cap, fixed = TRUE)))
  # 20 comparisons at .999 put each interval at .99995 exactly, which
  # four significant digits also rounded to 100.
  pw7 <- pairwise_within(Y8[, 1:7], bonferroni_ci = TRUE, conf_level = .999)
  expect_true(any(grepl("(Bonferroni: each interval at 99.995%)",
                        utils::capture.output(print(pw7)), fixed = TRUE)))
  # contrast_test(adjust = "bonferroni") writes the same line: 21 groups
  # make 210 comparisons, each at 1 - .01 / 210 at a family level of .99.
  set.seed(7)
  g <- factor(rep(sprintf("g%02d", 1:21), each = 4))
  v <- rnorm(84)
  ct <- contrast_test(aov(v ~ g), adjust = "bonferroni", conf_level = .99)
  expect_identical(attr(ct, "simultaneous")$each, 1 - .01 / 210)
  out_ct <- utils::capture.output(print(ct))
  expect_true(any(grepl(paste0("simultaneous across the 210 comparisons ",
                               "(Bonferroni: each interval at 99.995%)"),
                        out_ct, fixed = TRUE)))
  expect_false(any(grepl("at 100%", out_ct, fixed = TRUE)))
})
