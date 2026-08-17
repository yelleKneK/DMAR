# Mixed-effects ANOVA: classical F-ratios for fixed/random/mixed factors.
#' Mixed-Model ANOVA F-Ratios for One- and Two-Way Designs
#'
#' Computes the classical \emph{F}-ratios for one- and two-way ANOVA
#' designs when one or more factors are random rather than fixed,
#' using the expected-mean-square (EMS) rules that determine the
#' correct denominator for each test (Searle, Casella, & McCulloch,
#' 1992). The function is the closed-form alternative to fitting via
#' \code{\link[lme4]{lmer}} and is the standard treatment in
#' classical psychometrics and design-of-experiments texts (Maxwell,
#' Delaney, & Kelley, 2027, Ch. 10).
#'
#' @param data A \code{data.frame} containing the response, the
#'   factor(s), and (when both factors are crossed) the subject
#'   identifier.
#' @param outcome Character name of the response column.
#' @param factor_A Character name of factor \emph{A}.
#' @param factor_B Character name of factor \emph{B}, or \code{NULL}
#'   for one-way designs. Default \code{NULL}.
#' @param A_type One of \code{"fixed"} (default) or \code{"random"}:
#'   the EMS classification of factor \emph{A}.
#' @param B_type One of \code{"fixed"} or \code{"random"} (default):
#'   the EMS classification of factor \emph{B}. Ignored when
#'   \code{factor_B} is \code{NULL}.
#'
#' @return A \code{data.frame} with one row per testable effect.
#'   Columns: \code{effect}, \code{ss}, \code{df}, \code{ms},
#'   \code{denominator}, \code{F_value}, \code{p_value}.
#'
#' @details
#' \strong{One way design (only \code{factor_A}).}
#' \itemize{
#'   \item Both \emph{A} fixed and \emph{A} random use the same
#'         observed \emph{F}-ratio \eqn{MS_A / MS_{\mathrm{within}}};
#'         the test of "is there an effect of \emph{A}?" is identical.
#'         The interpretation differs: the random-effects test asks
#'         whether the variance component \eqn{\sigma^2_A} is zero.
#' }
#'
#' \strong{Two-way design, both fixed (Model I).}
#' \itemize{
#'   \item \eqn{F_A = MS_A / MS_{AB}} (when interaction is present in
#'         the model and treated as error) or \eqn{F_A = MS_A /
#'         MS_{\mathrm{within}}} (when interaction is pooled into
#'         error). The function uses \eqn{MS_{\mathrm{within}}} as
#'         the denominator throughout for Model I.
#' }
#'
#' \strong{Two-way design, both random (Model II).}
#' \itemize{
#'   \item \eqn{F_A = MS_A / MS_{AB}}, \eqn{F_B = MS_B / MS_{AB}},
#'         \eqn{F_{AB} = MS_{AB} / MS_{\mathrm{within}}}.
#' }
#'
#' \strong{Two-way mixed design (Model III, e.g., \emph{A} fixed,
#' \emph{B} random).}
#' \itemize{
#'   \item \eqn{F_A = MS_A / MS_{AB}} (fixed factor against the
#'         interaction with the random factor)
#'   \item \eqn{F_B = MS_B / MS_{\mathrm{within}}} (random factor
#'         against the within-cell residual)
#'   \item \eqn{F_{AB} = MS_{AB} / MS_{\mathrm{within}}}.
#' }
#'
#' \strong{Balanced data assumed.} The classical EMS rules require
#' equal cell sizes. The function errors out on unbalanced data and
#' recommends a mixed-effects fit via \code{\link[lme4]{lmer}}.
#'
#' \strong{Sums of squares.} For the balanced designs this function
#' targets, the Type I, Type II, and Type III sums of squares for each
#' effect coincide, so the decomposition is unambiguous and no
#' sums-of-squares type needs to be chosen (Maxwell, Delaney, &
#' Kelley, 2027, Ch. 7). The returned object carries a numeric
#' \code{sum_of_squares_type} attribute equal to 3, with the
#' understanding that it equals Types I and II here; it records the
#' convention without implying a choice that would matter for these
#' designs.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 7 on the sums of
#'   squares for balanced designs and Chapter 10 on random and mixed
#'   effects.)
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (1992).
#'   \emph{Variance components}. Wiley.
#'
#' @seealso \code{\link{anova_within}}, \code{\link{anova_within_two_way}},
#'   \code{\link[lme4]{lmer}}
#'
#' @examples
#' # 1. Two-way mixed design: A fixed, B random.
#' set.seed(113)
#' grid <- expand.grid(A = factor(1:3), B = factor(1:5), rep = 1:4)
#' grid$y <- with(grid, 0.8 * as.integer(A) + rep(rnorm(5, 0, 1.5), each = 1)[as.integer(B)] +
#'                       rnorm(nrow(grid), 0, 1))
#' mixed_anova(grid, outcome = "y", factor_A = "A", factor_B = "B",
#'             A_type = "fixed", B_type = "random")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#' @family mixed models
#'
#' @export

mixed_anova <- function(data, outcome, factor_A, factor_B = NULL,
                        A_type = c("fixed", "random"),
                        B_type = c("random", "fixed")) {
  A_type <- match.arg(A_type)
  B_type <- match.arg(B_type)

  if (!is.character(outcome) || !outcome %in% names(data))
    stop("'outcome' must be a column name in 'data'.")
  if (!is.character(factor_A) || !factor_A %in% names(data))
    stop("'factor_A' must be a column name in 'data'.")

  y <- data[[outcome]]
  A <- factor(data[[factor_A]])
  if (length(unique(A)) < 2L) stop("Factor A must have at least 2 levels.")

  # One way ----
  if (is.null(factor_B)) {
    fit <- stats::aov(y ~ A)
    tab <- stats::anova(fit)
    rn  <- trimws(rownames(tab))
    ms_A <- tab[rn == "A", "Mean Sq"]
    ms_w <- tab[rn == "Residuals", "Mean Sq"]
    df_A <- tab[rn == "A", "Df"]
    df_w <- tab[rn == "Residuals", "Df"]
    F_v  <- ms_A / ms_w
    p_v  <- stats::pf(F_v, df_A, df_w, lower.tail = FALSE)

    out <- data.frame(
      effect = "A",
      ss = tab[rn == "A", "Sum Sq"], df = df_A, ms = ms_A,
      denominator = "MS_within",
      F_value = F_v, p_value = p_v,
      stringsAsFactors = FALSE, row.names = NULL
    )
    # Balanced design: Type I, II, and III sums of squares coincide
    # (MDK 2027, Ch. 7); record the convention without a string row.
    attr(out, "sum_of_squares_type") <- 3
    return(out)
  }

  # Two-way ----
  if (!factor_B %in% names(data))
    stop("'factor_B' must be a column name in 'data'.")
  B <- factor(data[[factor_B]])
  if (length(unique(B)) < 2L) stop("Factor B must have at least 2 levels.")

  # Balance check on the rows aov() will actually use. aov(y ~ A * B) drops any
  # row with a missing y, A, or B, so counting the full data would let a single
  # NA in the outcome slip through an otherwise-balanced design and defeat the
  # documented "errors on unbalanced data" contract.
  ok <- stats::complete.cases(y, A, B)
  cell_n <- table(A[ok], B[ok])
  if (length(unique(as.vector(cell_n))) != 1L || cell_n[1] < 2L)
    stop("Design is unbalanced (or has <2 obs/cell). Use lme4::lmer().")

  fit <- stats::aov(y ~ A * B)
  tab <- stats::anova(fit)
  rn  <- trimws(rownames(tab))
  ms <- function(r) tab[rn == r, "Mean Sq"]
  ss <- function(r) tab[rn == r, "Sum Sq"]
  df <- function(r) tab[rn == r, "Df"]

  ms_A   <- ms("A");   ms_B  <- ms("B");   ms_AB <- ms("A:B")
  ms_w   <- ms("Residuals")
  df_A   <- df("A");   df_B  <- df("B");   df_AB <- df("A:B")
  df_w   <- df("Residuals")

  # EMS rules:
  if (A_type == "fixed" && B_type == "fixed") {        # Model I
    den_A  <- ms_w; den_B  <- ms_w; den_AB <- ms_w
    den_A_df <- df_w; den_B_df <- df_w; den_AB_df <- df_w
    den_A_lab <- "MS_within"; den_B_lab <- "MS_within"; den_AB_lab <- "MS_within"
  } else if (A_type == "random" && B_type == "random") { # Model II
    den_A <- ms_AB; den_A_df <- df_AB; den_A_lab <- "MS_AB"
    den_B <- ms_AB; den_B_df <- df_AB; den_B_lab <- "MS_AB"
    den_AB <- ms_w; den_AB_df <- df_w; den_AB_lab <- "MS_within"
  } else if (A_type == "fixed" && B_type == "random") { # Model III
    den_A <- ms_AB; den_A_df <- df_AB; den_A_lab <- "MS_AB"
    den_B <- ms_w;  den_B_df <- df_w;  den_B_lab <- "MS_within"
    den_AB <- ms_w; den_AB_df <- df_w; den_AB_lab <- "MS_within"
  } else {                                              # A random, B fixed
    den_A <- ms_w;  den_A_df <- df_w;  den_A_lab <- "MS_within"
    den_B <- ms_AB; den_B_df <- df_AB; den_B_lab <- "MS_AB"
    den_AB <- ms_w; den_AB_df <- df_w; den_AB_lab <- "MS_within"
  }

  F_A <- ms_A / den_A
  F_B <- ms_B / den_B
  F_AB <- ms_AB / den_AB
  p_A <- stats::pf(F_A, df_A, den_A_df, lower.tail = FALSE)
  p_B <- stats::pf(F_B, df_B, den_B_df, lower.tail = FALSE)
  p_AB <- stats::pf(F_AB, df_AB, den_AB_df, lower.tail = FALSE)

  out <- data.frame(
    effect = c("A", "B", "A:B"),
    ss = c(ss("A"), ss("B"), ss("A:B")),
    df = c(df_A, df_B, df_AB),
    ms = c(ms_A, ms_B, ms_AB),
    denominator = c(den_A_lab, den_B_lab, den_AB_lab),
    F_value = c(F_A, F_B, F_AB),
    p_value = c(p_A, p_B, p_AB),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # Balanced design: Type I, II, and III sums of squares coincide
  # (MDK 2027, Ch. 7); record the convention without a string row.
  attr(out, "sum_of_squares_type") <- 3
  out
}
