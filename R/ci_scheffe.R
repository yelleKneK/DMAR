# Scheffe-adjusted simultaneous CIs for arbitrary contrasts.
#' Scheffe-Adjusted Simultaneous Confidence Intervals for Contrasts
#'
#' Computes the Scheffe (1953, 1959) simultaneous confidence intervals
#' on user-specified contrasts among the means of a one-way design.
#' The Scheffe procedure controls the family-wise error rate for
#' \emph{any} set of contrasts, however many and however post-hoc,
#' which makes it more conservative than Tukey-Kramer or Bonferroni
#' for the specific case of all-pairwise comparisons but optimal for
#' arbitrary post-hoc contrasts.
#'
#' @param x A fitted \code{\link[stats]{lm}} or \code{\link[stats]{aov}}
#'   object with a single one-way factor predictor, or a numeric
#'   vector of observations with \code{group} supplied.
#' @param group Optional factor of group labels when \code{x} is a
#'   numeric vector.
#' @param contrasts An \eqn{a \times m} matrix or vector of contrast
#'   coefficients (rows = levels, columns = contrasts). Each column
#'   must sum to zero. If \code{NULL} (default), the function returns
#'   intervals for all \eqn{a (a - 1) / 2} pairwise contrasts.
#' @param conf_level Family-wise confidence level. Default
#'   \code{0.95}.
#'
#' @return A \code{data.frame} with one row per contrast.
#'   Columns: \code{contrast} (a printed label),
#'   \code{contrast_value}, \code{se}, \code{F_statistic},
#'   \code{lower_limit}, \code{upper_limit}, \code{p_adjusted}.
#'
#' @details
#' \strong{Critical value.} For \eqn{a} groups with \eqn{\nu} error
#' degrees of freedom, the Scheffe critical value is
#' \deqn{S \;=\; \sqrt{(a - 1) F_{1 - \alpha, a - 1, \nu}},}
#' where \eqn{F_{1 - \alpha, a - 1, \nu}} is the upper \eqn{\alpha}
#' quantile of the central \emph{F} distribution. The Scheffe
#' simultaneous CI on a contrast \eqn{\psi = \sum_i c_i \mu_i} is
#' \deqn{\hat\psi \;\pm\; S \cdot \mathit{SE}(\hat\psi),}
#' where \eqn{\mathit{SE}(\hat\psi) = \sqrt{\mathit{MS}_E \sum_i c_i^2 / n_i}}.
#'
#' \strong{Scope.} The Scheffe family-wise coverage holds for
#' \emph{any} number of contrasts, pairwise, complex, or chosen
#' after looking at the data. The trade-off is conservativeness: for
#' all-pairwise comparisons, Tukey-Kramer is uniformly more powerful.
#'
#' @references
#' Scheffe, H. (1953). A method for judging all contrasts in the
#'   analysis of variance. \emph{Biometrika, 40}(1/2), 87--104.
#'
#' Scheffe, H. (1959). \emph{The analysis of variance}. Wiley.
#'
#' @seealso \code{\link{cv_scheffe}}, \code{\link{ci_tukey_kramer}},
#'   \code{\link{ci_dunnett}}
#'
#' @examples
#' # 1. All pairwise contrasts among the six marketing panels of the
#' #    test_market data via the default:
#' fit <- lm(brand_movement ~ panel, data = test_market)
#' ci_scheffe(fit)
#'
#' # 2. A contrast chosen after inspecting the means: the two panels with
#' #    the highest brand movement (5 and 6) against the two with the
#' #    lowest (1 and 2). The Scheffe coverage holds for a contrast picked
#' #    this way, and the interval still excludes zero even though none of
#' #    the pairwise intervals above does.
#' cmat <- matrix(c(-0.5, -0.5, 0, 0, 0.5, 0.5), nrow = 6,
#'                dimnames = list(levels(test_market$panel),
#'                                "panels 5,6 - panels 1,2"))
#' ci_scheffe(fit, contrasts = cmat)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

ci_scheffe <- function(x, group = NULL, contrasts = NULL,
                       conf_level = 0.95) {
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      stop("ci_scheffe() requires a fitted lm/aov with one factor predictor.")
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
    tab <- stats::anova(x); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      stop("'x' and 'group' must be the same length.")
    g <- factor(group); y <- x
    ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]
    fit <- stats::aov(y ~ g)
    tab <- stats::anova(fit); rn <- trimws(rownames(tab))
    ms_e <- tab[rn == "Residuals", "Mean Sq"]
    df_e <- tab[rn == "Residuals", "Df"]
  } else {
    stop("Supply either a fitted lm/aov, or (x, group) vectors.")
  }

  lev <- levels(g); a <- length(lev)
  means <- tapply(y, g, mean); ns <- tapply(y, g, length)

  if (is.null(contrasts)) {
    m <- a * (a - 1L) / 2L
    cmat <- matrix(0, nrow = a, ncol = m)
    rownames(cmat) <- lev
    labs <- character(m); idx <- 0L
    for (i in seq_len(a - 1L)) for (j in (i + 1L):a) {
      idx <- idx + 1L
      cmat[i, idx] <- -1; cmat[j, idx] <- 1
      labs[idx] <- paste(lev[j], "-", lev[i])
    }
    colnames(cmat) <- labs
  } else {
    if (is.vector(contrasts)) contrasts <- matrix(contrasts, ncol = 1L,
                                                  dimnames = list(NULL, "C1"))
    if (nrow(contrasts) != a)
      stop(sprintf("Contrast matrix must have %d rows (= number of levels).", a))
    if (any(abs(colSums(contrasts)) > 1e-8))
      stop("Each contrast column must sum to zero.")
    cmat <- contrasts
    if (is.null(rownames(cmat))) rownames(cmat) <- lev
    if (is.null(colnames(cmat))) colnames(cmat) <- paste0("C", seq_len(ncol(cmat)))
  }

  S_crit <- sqrt((a - 1L) * stats::qf(conf_level, a - 1L, df_e))

  rows <- vector("list", ncol(cmat))
  for (k in seq_len(ncol(cmat))) {
    c_vec <- cmat[, k]
    psi   <- sum(c_vec * means)
    se    <- sqrt(ms_e * sum(c_vec^2 / ns))
    F_v   <- (psi / se)^2 / (a - 1L)
    p_adj <- stats::pf(F_v, a - 1L, df_e, lower.tail = FALSE)
    hw    <- S_crit * se
    rows[[k]] <- data.frame(
      contrast        = colnames(cmat)[k],
      contrast_value  = psi,
      se              = se,
      F_statistic     = F_v,
      lower_limit     = psi - hw,
      upper_limit     = psi + hw,
      p_adjusted      = p_adj,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  class(out) <- c("dmar_post_hoc_ci", class(out))
  out
}
