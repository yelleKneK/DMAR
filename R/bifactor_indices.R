#' Bifactor Model Dimensionality and Reliability Indices
#'
#' Computes the indices used to judge whether a multidimensional scale is
#' nonetheless unidimensional enough to score as a single total
#' (Rodriguez, Reise, and Haviland, 2016): the explained common variance
#' (ECV), coefficient omega and omega hierarchical (omega_H) for the
#' general factor, omega hierarchical subscale (omega_HS) for each group
#' factor, the percentage of uncontaminated correlations (PUC), and
#' coefficient H, the construct reliability (maximal reliability) of
#' Hancock and Mueller (2001). The input is a fitted bifactor model in which one general
#' factor loads on every item and each item loads on exactly one orthogonal
#' group factor.
#'
#' @param fit A fitted bifactor \pkg{lavaan} model: one general factor on
#'   all items plus orthogonal group factors, each item on one group
#'   factor. The factors must be orthogonal (the bifactor specification).
#' @param general Optional name of the general factor. When \code{NULL}
#'   (default) the general factor is detected as the one that loads on every
#'   item.
#'
#' @details
#' Let \eqn{\lambda^g_i} be item \eqn{i}'s standardized loading on the
#' general factor, \eqn{\lambda^s_i} its loading on its group factor, and
#' \eqn{\theta_i} its standardized residual variance. With orthogonal
#' factors the overall ECV is
#' \eqn{\sum_i (\lambda^g_i)^2 / \sum_i [(\lambda^g_i)^2 + (\lambda^s_i)^2]};
#' omega and omega_H share the total-score variance
#' \eqn{(\sum_i \lambda^g_i)^2 + \sum_g (\sum_{i \in g} \lambda^s_i)^2 +
#' \sum_i \theta_i} as denominator, with the general part
#' \eqn{(\sum_i \lambda^g_i)^2} in the numerator of omega_H. Each group
#' factor's omega_HS uses the analogous numerator and that subscale's own
#' total variance. PUC is one minus the share of item pairs that fall within
#' the same group factor. Coefficient H is computed from
#' \eqn{\sum \lambda^2 / (1 - \lambda^2)} on the relevant loadings (see
#' \code{\link{reliability_H}}). High ECV and PUC with a high omega_H
#' support scoring a single total; substantial subscale omega_HS argues for
#' reporting subscales as well.
#'
#' A bifactor model that is over-parameterized for the data frequently
#' yields an improper (Heywood) solution -- a standardized loading outside
#' \eqn{[-1, 1]} or a negative residual variance -- whose indices are not
#' trustworthy. When that happens the indices are still returned but a
#' warning is issued and the \code{"improper"} attribute is set to
#' \code{TRUE}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   factor (the general factor first, then each group factor) and columns
#'   \code{factor}, \code{ECV}, \code{omega}, \code{omega_H},
#'   \code{omega_HS}, \code{PUC}, and \code{H}. Quantities that do not apply
#'   to a row are \code{NA} (for example omega_H and PUC on a group factor).
#'   The \code{"improper"} attribute flags a Heywood solution.
#'
#' @references
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit,
#'   & D. Sörbom (Eds.), \emph{Structural equation modeling: Present and
#'   future} (pp. 195--216). Scientific Software International.
#'
#' Reise, S. P. (2012). The rediscovery of bifactor measurement models.
#'   \emph{Multivariate Behavioral Research, 47}(5), 667--696.
#'   \doi{10.1080/00273171.2012.715555}
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   \emph{Psychological Methods, 21}(2), 137--150.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{reliability_omega}}, \code{\link{reliability_H}}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Nine items: one general factor and three orthogonal group factors.
#' set.seed(113)
#' n <- 600
#' g <- rnorm(n); grp <- list(rnorm(n), rnorm(n), rnorm(n))
#' X <- vapply(1:9, function(i)
#'   0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
#' colnames(X) <- paste0("x", 1:9)
#' model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#'           f1 =~ x1 + x2 + x3
#'           f2 =~ x4 + x5 + x6
#'           f3 =~ x7 + x8 + x9"
#' fit <- lavaan::cfa(model, data = as.data.frame(X),
#'                    orthogonal = TRUE, std.lv = TRUE)
#' bifactor_indices(fit)
#'
#' @export
bifactor_indices <- function(fit, general = NULL) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required. Install it with ",
         "install.packages(\"lavaan\").", call. = FALSE)
  }
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a lavaan fit object.", call. = FALSE)
  }
  std <- lavaan::standardizedSolution(fit)
  lam <- std[std$op == "=~", c("lhs", "rhs", "est.std"), drop = FALSE]
  if (nrow(lam) == 0L) {
    stop("The fit contains no measurement (=~) loadings.", call. = FALSE)
  }
  items <- unique(lam$rhs)
  factors <- unique(lam$lhs)

  # Identify the general factor: the one loading on every item.
  loads_all <- vapply(factors, function(f) all(items %in% lam$rhs[lam$lhs == f]),
                      logical(1))
  if (is.null(general)) {
    if (sum(loads_all) != 1L) {
      stop("Could not identify a unique general factor (one factor loading ",
           "on every item). Name it with 'general'.", call. = FALSE)
    }
    general <- factors[loads_all]
  } else if (!general %in% factors || !all(items %in% lam$rhs[lam$lhs == general])) {
    stop("'general' must name a factor that loads on every item.", call. = FALSE)
  }
  groups <- setdiff(factors, general)
  if (!length(groups)) {
    stop("No group factors found; a bifactor model needs a general factor ",
         "plus two or more group factors.", call. = FALSE)
  }

  # Per-item general loading, group loading, group membership.
  lg <- stats::setNames(lam$est.std[lam$lhs == general][match(items, lam$rhs[lam$lhs == general])], items)
  ls <- stats::setNames(rep(NA_real_, length(items)), items)
  grp <- stats::setNames(rep(NA_character_, length(items)), items)
  for (g in groups) {
    gi <- lam$rhs[lam$lhs == g]
    ls[gi] <- lam$est.std[lam$lhs == g][match(gi, lam$rhs[lam$lhs == g])]
    grp[gi] <- g
  }
  if (anyNA(ls)) {
    stop("Every item must load on exactly one group factor.", call. = FALSE)
  }

  # Standardized residual variances (fall back to 1 - communality).
  rv <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, c("lhs", "est.std")]
  theta <- stats::setNames(1 - lg^2 - ls^2, items)
  if (nrow(rv)) theta[rv$lhs] <- rv$est.std

  improper <- any(abs(lg) > 1) || any(abs(ls) > 1) || any(theta < 0)
  if (improper) {
    warning("Improper (Heywood) solution: a standardized loading exceeds one ",
            "or a residual variance is negative. The indices are not ",
            "trustworthy until the model is respecified.", call. = FALSE)
  }

  H_of <- function(l) {                       # coefficient H from loadings
    l <- l[abs(l) < 1]
    if (!length(l)) return(NA_real_)
    t <- sum(l^2 / (1 - l^2)); t / (1 + t)
  }
  choose2 <- function(k) k * (k - 1) / 2

  # Overall (general-factor) row.
  total_common <- sum(lg)^2 + sum(vapply(groups, function(g) sum(ls[grp == g])^2, numeric(1)))
  total_var <- total_common + sum(theta)
  ecv_overall <- sum(lg^2) / (sum(lg^2) + sum(ls^2))
  k <- length(items)
  within_pairs <- sum(vapply(groups, function(g) choose2(sum(grp == g)), numeric(1)))
  puc <- if (choose2(k) > 0) 1 - within_pairs / choose2(k) else NA_real_

  rows <- list(data.frame(
    factor = general, ECV = ecv_overall,
    omega = if (total_var > 0) total_common / total_var else NA_real_,
    omega_H = if (total_var > 0) sum(lg)^2 / total_var else NA_real_,
    omega_HS = NA_real_, PUC = puc, H = H_of(lg),
    stringsAsFactors = FALSE
  ))
  for (g in groups) {
    sel <- grp == g
    lg_s <- lg[sel]; ls_s <- ls[sel]; th_s <- theta[sel]
    denom <- sum(lg_s^2) + sum(ls_s^2)
    sub_common <- sum(lg_s)^2 + sum(ls_s)^2
    sub_total <- sub_common + sum(th_s)
    rows[[length(rows) + 1L]] <- data.frame(
      factor = g,
      ECV = if (denom > 0) sum(lg_s^2) / denom else NA_real_,
      omega = if (sub_total > 0) sub_common / sub_total else NA_real_,
      omega_H = NA_real_,
      omega_HS = if (sub_total > 0) sum(ls_s)^2 / sub_total else NA_real_,
      PUC = NA_real_, H = H_of(ls_s),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  res <- .as_dmar_tbl(out)
  attr(res, "improper") <- improper
  res
}
