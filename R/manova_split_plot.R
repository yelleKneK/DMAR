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
#' separately. Wraps \code{\link[car]{Anova}} and returns the result in
#' the tidy DMAR style.
#'
#' @param data A \code{data.frame} in wide format, one row per subject,
#'   with the within-subjects measurements in separate columns and the
#'   between-subjects factor as one additional column.
#' @param within Character vector of column names holding the
#'   repeated measures values (one column per level of the within-
#'   subjects factor). Must be in the canonical level order.
#' @param between Character name of the between-subjects factor
#'   column in \code{data}.
#' @param ss_type The sum-of-squares type for the between-subjects
#'   effects, passed through to \code{\link[car]{Anova}}. Accepts the
#'   integers \code{1}, \code{2}, or \code{3} or the equivalent Roman-
#'   numeral strings \code{"I"}, \code{"II"}, or \code{"III"}, and
#'   defaults to \code{3L} (Type III). The chosen type is recorded in
#'   the returned table (see \strong{Value}). Type I (\code{1} or
#'   \code{"I"}) is not available for this mixed design, because
#'   \code{\link[car]{Anova}} computes only Type II and Type III sums
#'   of squares for the multivariate repeated measures path; requesting
#'   it raises an error.
#'
#' @return A \code{data.frame} with rows for each of the three
#'   effects crossed with each of the four multivariate statistics,
#'   plus one trailing row recording the sum-of-squares type. Columns:
#'   \code{effect}, \code{statistic_name}, \code{statistic_value},
#'   \code{F_approx}, \code{df_1}, \code{df_2}, \code{p_value}. The
#'   final row has \code{effect == "sum_of_squares_type"} and carries
#'   the chosen type (\code{1}, \code{2}, or \code{3}) in its numeric
#'   \code{statistic_value}; its remaining numeric columns are
#'   \code{NA}, so the \code{statistic_value} column stays numeric.
#'
#' @details
#' \strong{The four statistics.} For an effect with \eqn{H} and
#' \eqn{E} hypothesis- and error-cross-products matrices:
#' \itemize{
#'   \item Wilks's \eqn{\Lambda = \det(E) / \det(E + H)}
#'   \item Pillai's trace \eqn{V = \mathrm{tr}(H (E + H)^{-1})}
#'   \item Hotelling-Lawley trace \eqn{T_0^2 = \mathrm{tr}(H E^{-1})}
#'   \item Roy's largest root \eqn{\theta = \lambda_1(H E^{-1})}
#' }
#'
#' \strong{When to use which.} Pillai's trace is the most robust to
#' departures from the multivariate normal / homogeneous-covariance
#' assumptions. Wilks's \eqn{\Lambda} is the most widely reported.
#' Roy's largest root is the most powerful when the alternative
#' concentrates on a single dimension. The four statistics agree
#' exactly when the effect has 1 numerator degree of freedom.
#'
#' \strong{Sum-of-squares type.} The between-subjects effects are
#' computed by \code{\link[car]{Anova}} using the sum-of-squares type
#' selected through \code{ss_type} (Type III by default). Type II
#' conditions each effect on the others that do not contain it, and
#' Type III conditions each effect on every other effect in the model;
#' for a single between-subjects factor the two coincide, and they can
#' differ once additional between-subjects terms are present. Type I,
#' the sequential decomposition, is not available here:
#' \code{\link[car]{Anova}} computes only Type II and Type III for the
#' multivariate repeated measures path. The type in force is reported
#' in the returned table so the analysis is self-documenting.
#'
#' \strong{Dependency.} Requires the \pkg{car} package on CRAN.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 7 on higher-order
#'   designs and Chapter 14.)
#'
#' Rencher, A. C., & Christensen, W. F. (2012). \emph{Methods of
#'   multivariate analysis} (4th ed.). Wiley.
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
#' # effect (labeled B), and their interaction, each with all four
#' # multivariate criteria, followed by a row recording the sum-of-squares
#' # type. The multivariate tests make no sphericity assumption, which is
#' # what recommends them over the univariate repeated measures F and its
#' # epsilon corrections.
#' manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
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
    stop("The 'car' package is required. Install it with: ",
         "install.packages('car')")
  if (!is.character(within) || length(within) < 2L ||
      !all(within %in% names(data)))
    stop("'within' must be 2+ column names in 'data'.")
  if (!is.character(between) || length(between) != 1L ||
      !between %in% names(data))
    stop("'between' must be a single column name in 'data'.")
  ss_type <- .manova_split_plot_ss_type(ss_type)
  if (ss_type == 1L)
    stop("ss_type = 1 (Type I) is not available for a mixed-design ",
         "MANOVA: car::Anova computes only Type II or Type III ",
         "sums of squares for the multivariate repeated measures ",
         "path. Use ss_type = 2 or ss_type = 3.")

  d <- data[, c(between, within), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d[[between]] <- factor(d[[between]])
  if (length(levels(d[[between]])) < 2L)
    stop("'between' factor must have at least 2 levels.")

  Y     <- as.matrix(d[, within, drop = FALSE])
  fit   <- stats::lm(stats::as.formula(
              sprintf("Y ~ %s", between)), data = d)
  idata <- data.frame(
    Time = factor(within, levels = within)
  )
  manfit <- car::Anova(fit, idata = idata, idesign = ~ Time,
                       type = ss_type)

  smry <- summary(manfit, multivariate = TRUE)
  mv_tests <- smry$multivariate.tests

  effect_map <- c(
    "A"           = between,
    "B"           = "Time",
    "interaction" = paste0(between, ":Time")
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
    stop("Could not extract multivariate tests from car::Anova output.")

  out <- do.call(rbind, out_rows)

  # Record the sum-of-squares type used for the between-subjects
  # effects as a trailing numeric row. The value column stays numeric:
  # the type (1, 2, or 3) goes in statistic_value and the remaining
  # numeric columns are NA.
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
  out
}

# Normalize the user-supplied sum-of-squares type to the integer
# 1, 2, or 3 that car::Anova accepts. Accepts the integers 1/2/3 or
# the Roman-numeral strings "I"/"II"/"III" (case-insensitive).
.manova_split_plot_ss_type <- function(ss_type) {
  if (length(ss_type) != 1L)
    stop("'ss_type' must be a single value: 1, 2, 3, ",
         "\"I\", \"II\", or \"III\".")
  if (is.character(ss_type)) {
    key <- toupper(trimws(ss_type))
    val <- c("I" = 1L, "II" = 2L, "III" = 3L,
             "1" = 1L, "2" = 2L, "3" = 3L)[key]
    if (is.na(val))
      stop("'ss_type' must be one of 1, 2, 3, ",
           "\"I\", \"II\", or \"III\".")
    return(unname(val))
  }
  if (is.numeric(ss_type) && ss_type %in% c(1, 2, 3))
    return(as.integer(ss_type))
  stop("'ss_type' must be one of 1, 2, 3, \"I\", \"II\", or \"III\".")
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
