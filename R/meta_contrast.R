#' Contrast Among Study Effect Sizes (Rosenthal-Rubin)
#'
#' Tests a focused hypothesis about \emph{differences} among independent
#' study effect sizes by the method of Rosenthal and Rubin (1982): given
#' effects \eqn{y_i} with sampling variances \eqn{v_i} and contrast weights
#' \eqn{\lambda_i} summing to zero,
#' \deqn{z \;=\; \frac{\sum \lambda_i y_i}{\sqrt{\sum \lambda_i^2 v_i}}}
#' is referred to the standard normal. This is how a meta-analyst asks a
#' pointed moderator question (\dQuote{do the effects decline with weeks of
#' prior teacher-student contact?}) rather than the diffuse heterogeneity
#' question (\dQuote{do the effects differ at all?}). Raudenbush (1984) used
#' exactly this test for the teacher expectancy literature, with weights
#' inversely proportional to weeks of prior contact.
#'
#' @param yi Numeric vector of study effect sizes (any metric whose sampling
#'   distribution is approximately normal; standardized mean differences and
#'   Fisher's Z correlations qualify).
#' @param vi Sampling variances of \code{yi}, one per study.
#' @param weights Contrast weights, one per study. If they do not already
#'   sum to zero they are mean-centered (with a message) when
#'   \code{center = TRUE}, the convenient route for weights built from a
#'   moderator such as \code{1 / (weeks + 2)}.
#' @param center Logical: mean-center \code{weights} that do not sum to
#'   zero? Default \code{TRUE}.
#'
#' @details
#' The two-sided \emph{p}-value is reported; halve it for a directional
#' hypothesis stated in advance (Raudenbush's \eqn{z = 2.75} carried the
#' one-tailed \eqn{p = .003}). Dividing the squared contrast \eqn{z^2} by
#' the total heterogeneity statistic \eqn{Q} from \code{\link{meta_es}}
#' gives the proportion of between-study heterogeneity the contrast
#' accounts for, the meta-analytic analog of a contrast's share of the
#' between-group sum of squares.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the
#'   contrast \code{estimate} (\eqn{\sum \lambda_i y_i}), its \code{se},
#'   the \code{z} statistic, the two-sided \code{p_value}, and \code{k}.
#'
#' @references
#' Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
#'   pupil IQ as a function of the credibility of expectancy induction: A
#'   synthesis of findings from 18 experiments. \emph{Journal of
#'   Educational Psychology, 76}(1), 85--97.
#'
#' Rosenthal, R., & Rubin, D. B. (1982). Comparing effect sizes of
#'   independent studies. \emph{Psychological Bulletin, 92}(2), 500--504.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_es}} for the pooled effect and the total
#'   heterogeneity the contrast partitions; \code{\link{combine_p}} for
#'   combined significance tests; \code{\link{contrast_test}} for the
#'   single-study ANOVA analog.
#'
#' @family meta-analysis
#'
#' @keywords htest
#'
#' @examples
#' # Raudenbush (1984): do expectancy effects decline with weeks of prior
#' # teacher-student contact? Weights inversely proportional to weeks + 2,
#' # study-level data (Pellegrini & Hicks merged), d variances from the
#' # standard large-sample formula.
#' data(teacher_expectancy)
#' study <- teacher_expectancy[-c(4, 5), ]
#' d  <- append(study$d, 0.52, after = 3)
#' wk <- append(study$weeks, 0, after = 3)
#' ne <- append(study$n_experimental, 22, after = 3)
#' nc <- append(study$n_control, 22, after = 3)
#' v  <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
#' meta_contrast(d, v, weights = 1 / (wk + 2))
#' # z near 2.75: the better teachers knew their pupils, the smaller the
#' # expectancy effect (one-tailed p = .003 in the paper).
#'
#' @export
#' @importFrom stats pnorm
meta_contrast <- function(yi, vi, weights, center = TRUE) {
  k <- length(yi)
  if (!is.numeric(yi) || k < 2L || anyNA(yi)) {
    stop("'yi' must be two or more effect sizes with no missing values.",
         call. = FALSE)
  }
  if (!is.numeric(vi) || length(vi) != k || anyNA(vi) || any(vi <= 0)) {
    stop("'vi' must give a positive sampling variance for each effect size.",
         call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) != k || anyNA(weights)) {
    stop("'weights' must be a numeric contrast weight for each effect size.",
         call. = FALSE)
  }
  if (abs(sum(weights)) > 1e-8) {
    if (!center) {
      stop("'weights' must sum to zero (or set center = TRUE to mean-center ",
           "them).", call. = FALSE)
    }
    message("Contrast weights mean-centered to sum to zero.")
    weights <- weights - mean(weights)
  }
  if (all(weights == 0)) {
    stop("All contrast weights are zero after centering.", call. = FALSE)
  }

  est <- sum(weights * yi)
  se  <- sqrt(sum(weights^2 * vi))
  z   <- est / se

  .as_dmar_tbl(data.frame(
    term  = c("estimate", "se", "z", "p_value", "k"),
    value = c(est, se, z, 2 * pnorm(-abs(z)), k),
    stringsAsFactors = FALSE
  ), p_terms = "p_value")
}
