# Mixed-design multivariate ANOVA (within- and between-subjects factors).
#' Mixed-Design Multivariate ANOVA With All Four Test Statistics
#'
#' Computes the multivariate analysis of variance for a mixed-design
#' (one between-subjects factor and one within-subjects factor),
#' returning Wilks's \eqn{\Lambda}, Pillai's trace, Hotelling-Lawley
#' trace, and Roy's largest root, each with the associated
#' \emph{F}-approximation, degrees of freedom, and \emph{p}-value. The
#' three effects, between-subjects (\emph{A}), within-subjects
#' (\emph{B}), and the interaction (\eqn{A \times B}), are tested
#' separately, following the multivariate approach to the split-plot
#' design of Maxwell, Delaney, and Kelley (2027, Section 14.6). Wraps
#' \code{\link[car]{Anova}} and returns the result in the tidy DMAR
#' style.
#'
#' @param data A \code{data.frame} in wide format, one row per subject,
#'   with the within-subjects measurements in separate columns and the
#'   between-subjects factor as one additional column.
#' @param within Character vector of column names holding the
#'   repeated measures values (one column per level of the within-
#'   subjects factor). Must be in the canonical level order.
#' @param between Character name of the between-subjects factor
#'   column in \code{data}.
#' @param ss_type The sum-of-squares type, passed through to
#'   \code{\link[car]{Anova}}. Accepts the integers \code{1}, \code{2},
#'   or \code{3} or the equivalent Roman-numeral strings \code{"I"},
#'   \code{"II"}, or \code{"III"}, and defaults to \code{3L} (Type III).
#'   With one between-subjects factor the type matters only for the
#'   within-subjects main effect when the groups differ in size: Type
#'   III tests the unweighted grand means of the within-subjects
#'   contrasts and Type II the weighted grand means (see
#'   \strong{Details}). The chosen type is recorded in the returned
#'   table (see \strong{Value}). Type I (\code{1} or \code{"I"}) is not
#'   available for this mixed design, because \code{\link[car]{Anova}}
#'   computes only Type II and Type III sums of squares for the
#'   multivariate repeated measures path; requesting it raises an
#'   error.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows for each
#'   of the three effects crossed with each of the four multivariate
#'   statistics, plus one trailing row recording the sum-of-squares type.
#'   Columns: \code{effect}, \code{statistic_name}, \code{statistic_value},
#'   \code{F_approx}, \code{df_1}, \code{df_2}, \code{p_value}. The final
#'   row has \code{effect == "sum_of_squares_type"} and carries the chosen
#'   type (\code{2} or \code{3}) in its numeric \code{statistic_value}; its
#'   remaining numeric columns are \code{NA}, so the \code{statistic_value}
#'   column stays numeric. The table prints through the
#'   \code{\link{dmar_tbl}} display layer: whole numbers such as the
#'   degrees of freedom without a decimal part, the statistics to 3
#'   decimal places, and \emph{p}-values to 4 decimal places (one
#'   below 0.0001 prints as \dQuote{< 0.0001}); the stored values keep
#'   full precision. \code{tidy()} gives one row per row of the table,
#'   its \code{term} joining \code{effect} and \code{statistic_name} and
#'   its \code{estimate} the \code{statistic_value}, so the final row is
#'   the term \code{"sum_of_squares_type:NA"} (its \code{statistic_name}
#'   is \code{NA}) with the type as its estimate; \code{glance()} gives
#'   one row, the number of rows (\code{n_terms}) with \code{conf_level}
#'   and \code{B_used} \code{NA}, since no interval or bootstrap is formed
#'   (see \code{\link{dmar_tidiers}}). \code{as.data.frame()} returns the
#'   table as a plain \code{data.frame}.
#'
#' @details
#' \strong{The three tests.} The multivariate approach replaces each
#' subject's \eqn{b} repeated measures with transformed variables: their
#' mean \eqn{M}, and \eqn{b - 1} contrasts among them, the \eqn{D}
#' variables. The between-subjects effect (\emph{A}) is the one-way
#' ANOVA of \eqn{M} across the \eqn{a} groups (Maxwell, Delaney, &
#' Kelley, 2027, Equation 14.40), which is also the univariate approach's
#' test of \emph{A}. The interaction (\eqn{A \times B}) asks whether the
#' groups differ in their means on the \eqn{D} variables, a one-way
#' MANOVA of the \eqn{D} variables (their Equations 14.62 and 14.63). The
#' within-subjects main effect (\emph{B}) asks whether the grand means of
#' the \eqn{D} variables are all zero, averaging over the groups (their
#' Equation 14.54). Both within-subjects tests take as their error matrix
#' the pooled within-group sums of squares and cross-products of the
#' \eqn{D} variables, and no test depends on which \eqn{b - 1} contrasts
#' are formed. A test of \emph{B} within a single group is a simple effect
#' (their Section 14.7.4), a different question from the main effect.
#'
#' \strong{Weighted or unweighted grand means.} When the groups differ in
#' size, the grand mean of a \eqn{D} variable can be the unweighted average
#' of its \eqn{a} group means or the average weighted by group size, which
#' is its mean over all \emph{N} subjects; the two agree when the groups
#' are of equal size, and the choice affects only the test of \emph{B}.
#' Maxwell, Delaney, and Kelley (2027, Section 14.6.9, p. 796, and Chapter
#' 14, note 6, p. 854) generally prefer the unweighted means, which they
#' note SAS and SPSS test by default, while allowing that the weighted means
#' suit some studies. The sum-of-squares type selects between them. Type III
#' (the default) tests the unweighted grand means: the between-subjects
#' factor is coded with sum-to-zero contrasts, whatever
#' \code{options("contrasts")} says, so that the intercept of the
#' regression of each \eqn{D} variable on the groups is the unweighted
#' mean of its group means, and the hypothesis matrix is
#' \eqn{a \tilde{n} \bar{D}_U \bar{D}_U^\top}{a n-tilde D-bar_U D-bar_U'},
#' where \eqn{\bar{D}_U}{D-bar_U} holds the unweighted grand means and
#' \eqn{\tilde{n}}{n-tilde} is the harmonic mean of the group sizes (their
#' notes 10 and 14). Type II tests the weighted grand means
#' \eqn{\bar{D}}{D-bar}, with the hypothesis matrix
#' \eqn{N \bar{D} \bar{D}^\top}{N D-bar D-bar'} (their note 11). Under R's
#' default treatment coding, the intercept would instead be the first
#' group's means, and its Type III test the simple effect of \emph{B} in
#' the first group (their Equation 14.66). With one between-subjects
#' factor, the two types give the same tests of \emph{A} and of the
#' interaction, and with equal group sizes the same test of \emph{B} as
#' well. Type I, the sequential decomposition, is not available here:
#' \code{\link[car]{Anova}} computes only Type II and Type III for the
#' multivariate repeated measures path. The type in force is reported in
#' the returned table so the analysis is self-documenting.
#'
#' \strong{The four statistics.} For an effect with \eqn{H} and
#' \eqn{E} hypothesis- and error-cross-products matrices:
#' \itemize{
#'   \item Wilks's \eqn{\Lambda = \det(E) / \det(E + H)}
#'   \item Pillai's trace \eqn{V = \mathrm{tr}(H (E + H)^{-1})}
#'   \item Hotelling-Lawley trace \eqn{T_0^2 = \mathrm{tr}(H E^{-1})}
#'   \item Roy's largest root \eqn{\theta = \lambda_1(H E^{-1})}
#' }
#'
#' \strong{When to use which.} For many data sets the four statistics are
#' likely to reach the same conclusion about statistical significance, but that
#' agreement is not guaranteed, and without more information about the
#' population than is usually available none of them can be said to be
#' best. Wilks's \eqn{\Lambda} is historically the most widely used, and
#' there is some evidence that Pillai's trace is the most robust of the
#' four (Maxwell, Delaney, & Kelley, 2027, Section 14.6.11). Roy's largest
#' root uses only the largest eigenvalue of \eqn{H E^{-1}}, the largest
#' ratio of hypothesis to error variation over linear combinations of the
#' transformed variables. The four statistics give the same
#' \emph{F} whenever \eqn{s = 1}, where \eqn{s} is the smaller of the
#' number of transformed variables and the effect's degrees of freedom for
#' each of them (Maxwell, Delaney, & Kelley, 2027, Section 14.6.11): always
#' for \emph{A} and \emph{B}, and for the interaction when either factor
#' has two levels.
#'
#' \strong{Dependency.} Requires the \pkg{car} package on CRAN.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 14, Section 14.6, on
#'   the multivariate approach to split-plot designs, and Chapter 7,
#'   Sections 7.9 and 7.10, on weighted and unweighted marginal means.)
#'
#' Rencher, A. C., & Christensen, W. F. (2012). \emph{Methods of
#'   multivariate analysis} (3rd ed.). Wiley.
#'   \doi{10.1002/9781118391686}
#'
#' @seealso \code{\link[car]{Anova}}, \code{\link[stats]{manova}},
#'   \code{\link{anova_within_two_way}}
#'
#' @examples
#' # Two groups of ten measured at three times. Both groups start at the
#' # same place and rise over time, and group B rises twice as fast, so the
#' # data generating means differ in slope as well as in level.
#' set.seed(113)
#' n_per <- 10
#' d <- data.frame(
#'   subject = factor(1:(2 * n_per)),
#'   group   = factor(rep(c("A", "B"), each = n_per)),
#'   t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
#'   t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
#'   t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
#' )
#'
#' # Rows are the between-subjects effect (labeled A), the within-subjects
#' # main effect (labeled B), and their interaction, each with all four
#' # multivariate criteria, followed by a row recording the sum-of-squares
#' # type. The multivariate tests make no sphericity assumption, which is
#' # what recommends them over the univariate repeated measures F and its
#' # epsilon corrections.
#' manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
#'
#' # With six subjects left in group A, the types differ in the B rows
#' # alone: Type III (the default) averages the two groups' changes over
#' # time equally, and Type II weights each group by its size.
#' d_unequal <- d[-(1:4), ]
#' manova_split_plot(d_unequal, within = c("t1", "t2", "t3"),
#'                   between = "group")
#' manova_split_plot(d_unequal, within = c("t1", "t2", "t3"),
#'                   between = "group", ss_type = 2)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#' @family mixed models
#'
#' @export

manova_split_plot <- function(data, within, between, ss_type = 3L) {
  if (!requireNamespace("car", quietly = TRUE))
    .dmar_refuse(c(
      "The 'car' package is required. Install it with: ",
      "install.packages('car')"), "package_required", package = "car")
  if (!is.character(within) || length(within) < 2L ||
      !all(within %in% names(data)))
    .dmar_refuse_argument(
      "'within' must be 2+ column names in 'data'.", arg = c("within", "data"))
  if (!is.character(between) || length(between) != 1L ||
      !between %in% names(data))
    .dmar_refuse_argument(
      "'between' must be a single column name in 'data'.",
      arg = c("between", "data"))
  ss_type <- .manova_split_plot_ss_type(ss_type)
  if (ss_type == 1L)
    .dmar_refuse_argument(c(
      "ss_type = 1 (Type I) is not available for a mixed-design ",
      "MANOVA: car::Anova computes only Type II or Type III ",
      "sums of squares for the multivariate repeated measures ",
      "path. Use ss_type = 2 or ss_type = 3."), arg = NULL)

  d <- data[, c(between, within), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  group <- factor(d[[between]])
  if (length(levels(group)) < 2L)
    .dmar_refuse_argument(
      "'between' factor must have at least 2 levels.", arg = "between")

  # The model is fit to local objects with fixed names, so no column name
  # in 'data' (one named "Time", "Y", or "group", or one that is not
  # syntactic) can enter the formula or the term labels car reports.
  #
  # The between-subjects factor is coded with sum-to-zero contrasts,
  # whatever options("contrasts") says. Under that coding the intercept of
  # the regression of each within-subjects contrast on the groups is the
  # unweighted mean of its group means, so the Type III test of the
  # intercept, which car reports as the within-subjects main effect, tests
  # the unweighted grand means, as Maxwell, Delaney, and Kelley (2027,
  # Section 14.6.10 and Chapter 14, note 10) do. Under R's default
  # treatment coding the intercept is the first group's mean, and its Type
  # III test is the simple effect of the within-subjects factor in the
  # first group (their Equation 14.66), which is what DMAR 1.0.x reported.
  # The Type II tests, and the tests of the between-subjects effect and the
  # interaction under either type, do not depend on the coding.
  Y   <- as.matrix(d[, within, drop = FALSE])
  fit <- stats::lm(Y ~ group, contrasts = list(group = stats::contr.sum))
  idata <- data.frame(
    Time = factor(within, levels = within)
  )
  manfit <- car::Anova(fit, idata = idata, idesign = ~ Time,
                       type = ss_type)

  # univariate = FALSE: the univariate tests and their epsilon corrections
  # are not reported here, and computing them warns whenever the
  # Huynh-Feldt estimate exceeds 1.
  smry <- summary(manfit, multivariate = TRUE, univariate = FALSE)
  mv_tests <- smry$multivariate.tests

  effect_map <- c(
    "A"           = "group",
    "B"           = "Time",
    "interaction" = "group:Time"
  )
  out_rows <- list()
  for (eff_label in names(effect_map)) {
    nm <- effect_map[[eff_label]]
    if (!nm %in% names(mv_tests)) next
    tt <- mv_tests[[nm]]
    stat_tbl <- .manova_split_plot_test_table(tt)
    if (is.null(stat_tbl)) next
    for (i in seq_len(nrow(stat_tbl))) {
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        effect          = eff_label,
        statistic_name  = stat_tbl$statistic_name[i],
        statistic_value = stat_tbl$statistic_value[i],
        F_approx        = stat_tbl$F_approx[i],
        df_1            = stat_tbl$df_1[i],
        df_2            = stat_tbl$df_2[i],
        p_value         = stat_tbl$p_value[i],
        stringsAsFactors = FALSE)
    }
  }
  if (length(out_rows) == 0L)
    .dmar_refuse_argument(
      "Could not extract multivariate tests from car::Anova output.",
      arg = NULL)

  out <- do.call(rbind, out_rows)

  # Record the sum-of-squares type as a trailing numeric row. The value
  # column stays numeric: the type (2 or 3; Type I was refused above)
  # goes in statistic_value and the remaining numeric columns are NA.
  out <- rbind(out, data.frame(
    effect          = "sum_of_squares_type",
    statistic_name  = NA_character_,
    statistic_value = as.numeric(ss_type),
    F_approx        = NA_real_,
    df_1            = NA_real_,
    df_2            = NA_real_,
    p_value         = NA_real_,
    stringsAsFactors = FALSE))

  rownames(out) <- NULL
  # The dmar_tbl display layer prints the p_value column, found by name,
  # to fixed decimals; DMAR 1.0.x returned the plain data.frame, whose
  # print method put small p-values in scientific notation. No interval is
  # formed, so no conf_level is recorded, and the type is already a row of
  # the table rather than an ss_type attribute.
  .as_dmar_tbl(out)
}

# Normalize the user-supplied sum-of-squares type to the integer
# 1, 2, or 3 that car::Anova accepts. Accepts the integers 1/2/3 or
# the Roman-numeral strings "I"/"II"/"III" (case-insensitive).
.manova_split_plot_ss_type <- function(ss_type) {
  if (length(ss_type) != 1L)
    .dmar_refuse_argument(c(
      "'ss_type' must be a single value: 1, 2, 3, ",
      "\"I\", \"II\", or \"III\"."), arg = "ss_type")
  if (is.character(ss_type)) {
    key <- toupper(trimws(ss_type))
    val <- c("I" = 1L, "II" = 2L, "III" = 3L,
             "1" = 1L, "2" = 2L, "3" = 3L)[key]
    if (is.na(val))
      .dmar_refuse_argument(c(
        "'ss_type' must be one of 1, 2, 3, ",
        "\"I\", \"II\", or \"III\"."), arg = "ss_type")
    return(unname(val))
  }
  if (is.numeric(ss_type) && ss_type %in% c(1, 2, 3))
    return(as.integer(ss_type))
  .dmar_refuse_argument(
    "'ss_type' must be one of 1, 2, 3, \"I\", \"II\", or \"III\".",
    arg = "ss_type")
}

# Compute the multivariate test-statistic table from a car
# linearHypothesis.mlm object. Mirrors the table built inside
# car's print.linearHypothesis.mlm so that the F-approximations
# match car's printed output exactly. Returns NULL if the error SSP
# matrix is singular (multivariate tests are unavailable).
.manova_split_plot_test_table <- function(x) {
  if (isTRUE(x$singular)) return(NULL)
  SSPE.qr <- qr(x$SSPE)
  eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)

  # car ships these as internal helpers (not in NAMESPACE). Pull them
  # via getFromNamespace() so the F-approximations match car's output.
  Pillai <- utils::getFromNamespace("Pillai", "car")
  Wilks  <- utils::getFromNamespace("Wilks",  "car")
  HL     <- utils::getFromNamespace("HL",     "car")
  Roy    <- utils::getFromNamespace("Roy",    "car")

  rows <- list(
    Pillai             = if ("Pillai"           %in% x$test) Pillai(eigs, x$df, x$df.residual),
    Wilks              = if ("Wilks"            %in% x$test) Wilks(eigs,  x$df, x$df.residual),
    `Hotelling-Lawley` = if ("Hotelling-Lawley" %in% x$test) HL(eigs,     x$df, x$df.residual),
    Roy                = if ("Roy"              %in% x$test) Roy(eigs,    x$df, x$df.residual)
  )
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) return(NULL)

  tab <- do.call(rbind, rows)
  colnames(tab) <- c("statistic_value", "F_approx", "df_1", "df_2")
  out <- data.frame(
    statistic_name  = rownames(tab),
    statistic_value = tab[, "statistic_value"],
    F_approx        = tab[, "F_approx"],
    df_1            = tab[, "df_1"],
    df_2            = tab[, "df_2"],
    stringsAsFactors = FALSE
  )
  out$p_value <- stats::pf(out$F_approx, out$df_1, out$df_2,
                           lower.tail = FALSE)
  rownames(out) <- NULL
  out
}
