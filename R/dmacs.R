#' The dMACS Effect Size of Measurement Noninvariance
#'
#' Quantifies how much a violation of measurement invariance actually matters
#' for an item, rather than only whether it is statistically detectable. A
#' likelihood ratio or score test can flag a loading or intercept difference
#' that is too small to change anyone's score in a meaningful way, and with a
#' large sample size it usually will. The dMACS index of Nye and Drasgow
#' (2011) answers the size question directly: it is the expected difference
#' between the reference group's and the focal group's measurement equations
#' for that item, averaged over the focal group's latent distribution and
#' standardized by the pooled item standard deviation, so it reads on the
#' familiar standardized mean difference scale. Input is either a fitted
#' multiple group \pkg{lavaan} model (requires \pkg{lavaan}) or the loadings,
#' intercepts, and pooled standard deviations a paper reports.
#'
#' @param fit A fitted multiple group \pkg{lavaan} object carrying a mean
#'   structure, with the two groups' loadings and intercepts on a common
#'   metric (see Details). Supply either \code{fit} or the parameter vectors,
#'   never both. A fit whose optimizer did not converge
#'   (\code{lavaan::lavInspect(fit, "converged")} is \code{FALSE}) is
#'   refused with an error of class \code{dmar_not_converged}, since its
#'   estimates mark where the optimizer stopped.
#' @param reference,focal Which groups play the reference and focal roles,
#'   each given as a single group label or a single group index in
#'   \code{lavInspect(fit, "group.label")}. When both are \code{NULL}
#'   (default) and the fit has exactly two groups, the first is the reference
#'   and the second the focal; with more than two groups they must be named.
#' @param factor_loading_reference,factor_loading_focal Numeric vectors of
#'   unstandardized factor loadings (\eqn{\lambda_R} and \eqn{\lambda_F}
#'   of Details), one per item, in the reference and focal groups. Any
#'   finite values are admissible.
#' @param intercept_reference,intercept_focal Numeric vectors of
#'   unstandardized item intercepts (\eqn{\nu_R} and \eqn{\nu_F} of
#'   Details), one per item, in the reference and focal groups, in the same
#'   item order as the loadings. Any finite values are admissible.
#' @param mean_focal Focal group latent mean \eqn{\mu_F} on the common metric,
#'   a single finite number (default \code{0}, the usual identification in
#'   which the reference group's latent mean is fixed at zero). Used only on
#'   the parameter vector path; with a \code{fit} it is read from the fit.
#' @param sd_focal Focal group latent standard deviation \eqn{\sigma_F} on the
#'   common metric, a single positive number (default \code{1}). Used only on
#'   the parameter vector path; with a \code{fit} it is read from the fit.
#' @param sd_pooled Pooled observed standard deviation of each item across the
#'   two groups, either one positive number applied to every item or one per
#'   item. Required on the parameter vector path; with a \code{fit} it is
#'   computed from the group sample sizes and observed item variances.
#' @param item_names Optional character vector of item labels, one per item.
#'   Defaults to the names carried by the parameter vectors, to the indicator
#'   names in the fit, or to \code{item_1}, \code{item_2}, and so on.
#'
#' @details
#' For item \eqn{i}, let \eqn{\nu_R} and \eqn{\lambda_R} be the reference
#' group's intercept and loading (\code{intercept_reference} and
#' \code{factor_loading_reference}), let \eqn{\nu_F} and \eqn{\lambda_F} be
#' the focal group's (\code{intercept_focal} and
#' \code{factor_loading_focal}), and let the focal group's latent variable be
#' \eqn{\eta \sim N(\mu_F, \sigma_F^2)} with density \eqn{f_F}. Each group's
#' measurement equation gives an expected item score at every value of
#' \eqn{\eta}, and dMACS is the root mean squared vertical distance between
#' those two lines over the focal group's latent distribution, divided by the
#' pooled item standard deviation:
#' \deqn{d_{MACS, i} = \frac{1}{SD_i} \sqrt{\int \left[ (\nu_R + \lambda_R
#'   \eta) - (\nu_F + \lambda_F \eta) \right]^2 f_F(\eta) \, d\eta}.}
#'
#' Writing \eqn{a = \nu_R - \nu_F} for the intercept difference and
#' \eqn{b = \lambda_R - \lambda_F} for the loading difference, the integrand is
#' \eqn{(a + b\eta)^2} and the integral is the second moment of a linear
#' function of a normal variate, so it has the closed form
#' \eqn{a^2 + 2ab\mu_F + b^2(\sigma_F^2 + \mu_F^2)}. No numerical integration
#' is performed. The standardizer is the pooled observed standard deviation of
#' the item,
#' \deqn{SD_i = \sqrt{\frac{(n_R - 1)s_R^2 + (n_F - 1)s_F^2}{n_R + n_F - 2}},}
#' the same pooling used by the standardized mean difference, which puts dMACS
#' on a scale a reader of \code{\link{smd}} already understands.
#'
#' dMACS is a root mean square and so is nonnegative by construction; it
#' carries no sign and is not given one here. The direction of the violation
#' is read from the returned components: \eqn{\nu_R - \nu_F} says which group
#' is scored higher at the mean of the latent variable, and
#' \eqn{\lambda_R - \lambda_F} says in which group the item discriminates more
#' sharply. When only the intercepts differ, the integral collapses to
#' \eqn{a^2} and dMACS reduces to \eqn{|a| / SD_i}, a plain standardized
#' intercept difference.
#'
#' The index is interpretable only when the two groups' loadings and
#' intercepts are expressed on a common metric. In practice that means a
#' partial invariance model in which a set of anchor items is constrained
#' equal across groups while the suspect items are freed, and the focal
#' group's latent mean and variance are freely estimated. A configural model,
#' which sets each group's latent scale separately, does not put the groups on
#' a common metric, and dMACS computed from one is not meaningful. The usual
#' workflow is therefore \code{\link{measurement_invariance}} to locate where
#' the ladder breaks, a partial invariance refit that frees the offending
#' parameters, and then \code{dmacs()} on that refit to judge whether the
#' violation is large enough to matter.
#'
#' On the fit path each group's item variance is computed from the data in the
#' fit with the usual \eqn{n - 1} divisor, using the number of nonmissing
#' observations for that item in that group. When the model was fitted from
#' sample moments rather than raw data, the sample covariances in the fit are
#' used instead, rescaled to the \eqn{n - 1} divisor when the fit's likelihood
#' option calls for it.
#'
#' @return A wide \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   item and columns
#'   \describe{
#'     \item{\code{item}}{Item label.}
#'     \item{\code{factor_loading_reference}}{Reference group
#'       unstandardized loading.}
#'     \item{\code{factor_loading_focal}}{Focal group unstandardized
#'       loading.}
#'     \item{\code{intercept_reference}}{Reference group unstandardized
#'       intercept.}
#'     \item{\code{intercept_focal}}{Focal group unstandardized intercept.}
#'     \item{\code{sd_pooled}}{Pooled observed standard deviation of the item.}
#'     \item{\code{dmacs}}{The dMACS effect size, nonnegative.}
#'   }
#'   The returned object carries four attributes: \code{"reference"} and
#'   \code{"focal"}, the two group labels, and \code{"mean_focal"} and
#'   \code{"sd_focal"}, the focal group's latent mean and standard deviation
#'   used in the integral (named by latent variable on the fit path). A
#'   fifth, \code{"measured_terms"}, names the numeric columns for the
#'   display, so that a loading or an intercept given as a whole number
#'   prints to the decimal places of the others (see
#'   \code{\link{dmar_tbl}}).
#'
#' @references
#' Meredith, W. (1993). Measurement invariance, factor analysis and factorial
#'   invariance. \emph{Psychometrika, 58}(4), 525--543.
#'   \doi{10.1007/BF02294825}
#'
#' Millsap, R. E. (2011). \emph{Statistical approaches to measurement
#'   invariance}. Routledge.
#'
#' Nye, C. D., Bradburn, J., Olenick, J., Bialko, C., & Drasgow, F. (2019).
#'   How big are my effects? Examining the magnitude of effect sizes in studies
#'   of measurement equivalence. \emph{Organizational Research Methods, 22}(3),
#'   678--709. \doi{10.1177/1094428118761122}
#'
#' Nye, C. D., & Drasgow, F. (2011). Effect size indices for analyses of
#'   measurement equivalence: Understanding the practical importance of
#'   differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#'   966--980. \doi{10.1037/a0022955}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{measurement_invariance}} for the invariance ladder
#'   that locates a violation; \code{\link{smd}} for the standardized mean
#'   difference whose pooling and scale dMACS borrows; \code{\link{cfa_1}} for
#'   the single group measurement model.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Reported measurement equations: the two items share loadings, and the
#' # second item's intercept is 0.30 higher in the reference group. With no
#' # loading difference, dMACS is just 0.30 divided by the pooled SD.
#' dmacs(factor_loading_reference = c(0.80, 0.75),
#'       factor_loading_focal = c(0.80, 0.75),
#'       intercept_reference = c(2.00, 2.30),
#'       intercept_focal = c(2.00, 2.00),
#'       sd_pooled = c(1.20, 1.10), item_names = c("optimism", "worry"))
#'
#' # A partial invariance model for the four spatial tests at the two
#' # Holzinger and Swineford schools. The anchors are constrained equal; the
#' # cubes and lozenges tests are freed, so only those two can move.
#' data(holzinger_swineford)
#' items <- c("t1_visual_perception", "t2_cubes",
#'            "t3_paper_form_board", "t4_lozenges")
#' model <- paste("spatial =~", paste(items, collapse = " + "))
#' fit <- lavaan::cfa(model, data = holzinger_swineford, group = "school",
#'                    group.equal = c("loadings", "intercepts"),
#'                    group.partial = c("spatial =~ t2_cubes", "t2_cubes ~ 1",
#'                                      "spatial =~ t4_lozenges",
#'                                      "t4_lozenges ~ 1"))
#' dmacs(fit)
#'
#' # The broom verbs: one row per item, and the group metadata.
#' generics::tidy(dmacs(fit))
#' generics::glance(dmacs(fit))
#'
#' @export
dmacs <- function(fit = NULL, reference = NULL, focal = NULL,
                  factor_loading_reference = NULL, factor_loading_focal = NULL,
                  intercept_reference = NULL, intercept_focal = NULL,
                  mean_focal = 0, sd_focal = 1, sd_pooled = NULL,
                  item_names = NULL) {

  from_vectors <- c(!is.null(factor_loading_reference),
                    !is.null(factor_loading_focal),
                    !is.null(intercept_reference), !is.null(intercept_focal),
                    !is.null(sd_pooled))
  if (!is.null(fit) && any(from_vectors)) {
    .dmar_refuse_argument(c(
      "Supply either 'fit' or the parameter vectors ",
      "('factor_loading_reference', 'factor_loading_focal', ",
      "'intercept_reference', 'intercept_focal', 'sd_pooled'), not both."),
      arg = c("fit", "factor_loading_reference", "factor_loading_focal"))
  }
  if (!is.null(fit) && (!missing(mean_focal) || !missing(sd_focal))) {
    .dmar_refuse_argument(c(
      "'mean_focal' and 'sd_focal' are read from 'fit'. Drop them, or ",
      "supply the parameter vectors instead of 'fit'."),
      arg = c("mean_focal", "sd_focal", "fit"))
  }
  if (is.null(fit) && !all(from_vectors)) {
    .dmar_refuse_argument(c(
      "Supply 'fit', or all of 'factor_loading_reference', ",
      "'factor_loading_focal', 'intercept_reference', 'intercept_focal', ",
      "and 'sd_pooled'."),
      arg = c("fit", "factor_loading_reference", "factor_loading_focal"))
  }

  if (!is.null(fit)) {
    parts <- .dmacs_from_fit(fit, reference, focal)
    factor_loading_reference <- parts$factor_loading_reference
    factor_loading_focal     <- parts$factor_loading_focal
    intercept_reference      <- parts$intercept_reference
    intercept_focal          <- parts$intercept_focal
    sd_pooled                <- parts$sd_pooled
    mu                       <- parts$mu
    sigma                    <- parts$sigma
    items                    <- parts$items
    label_reference          <- parts$label_reference
    label_focal              <- parts$label_focal
    mean_focal               <- parts$mean_focal
    sd_focal                 <- parts$sd_focal
  } else {
    check_numeric <- function(x, nm) {
      if (!is.numeric(x) || length(x) == 0L) {
        .dmar_refuse_argument(
          sprintf("'%s' must be a nonempty numeric vector.", nm), arg = nm)
      }
      if (any(!is.finite(x))) {
        .dmar_refuse_argument(
          sprintf("'%s' must be finite (no NA, NaN, or Inf).", nm), arg = nm)
      }
      invisible(NULL)
    }
    check_numeric(factor_loading_reference, "factor_loading_reference")
    check_numeric(factor_loading_focal, "factor_loading_focal")
    check_numeric(intercept_reference, "intercept_reference")
    check_numeric(intercept_focal, "intercept_focal")
    check_numeric(sd_pooled, "sd_pooled")
    check_numeric(mean_focal, "mean_focal")
    check_numeric(sd_focal, "sd_focal")

    k <- length(factor_loading_reference)
    lens <- c(factor_loading_focal = length(factor_loading_focal),
              intercept_reference = length(intercept_reference),
              intercept_focal = length(intercept_focal))
    if (any(lens != k)) {
      .dmar_refuse_argument(c(
        "'factor_loading_reference', 'factor_loading_focal', ",
        "'intercept_reference', and 'intercept_focal' must all have the ",
        "same length. Lengths are ", k, ", ", paste(lens, collapse = ", "),
        "."),
        arg = c("factor_loading_reference", "factor_loading_focal",
                "intercept_reference"))
    }
    if (!(length(sd_pooled) %in% c(1L, k))) {
      .dmar_refuse_argument(c(
        "'sd_pooled' must have length 1 or length ", k,
        ", one per item."), arg = "sd_pooled")
    }
    sd_pooled <- rep(sd_pooled, length.out = k)
    if (any(sd_pooled <= 0)) {
      .dmar_refuse_argument("'sd_pooled' must be positive.", arg = "sd_pooled")
    }
    if (length(mean_focal) != 1L) {
      .dmar_refuse_argument(
        "'mean_focal' must be a single number.", arg = "mean_focal")
    }
    if (length(sd_focal) != 1L || sd_focal <= 0) {
      .dmar_refuse_argument(
        "'sd_focal' must be a single positive number.", arg = "sd_focal")
    }

    supplied_names <- Filter(Negate(is.null),
                             list(names(factor_loading_reference),
                                  names(factor_loading_focal),
                                  names(intercept_reference),
                                  names(intercept_focal)))
    if (length(supplied_names) > 1L &&
        !all(vapply(supplied_names[-1L], identical, logical(1L),
                    supplied_names[[1L]]))) {
      .dmar_refuse_argument(c(
        "The item names carried by the parameter vectors do not align ",
        "across groups. Put the items in the same order in every vector."),
        arg = NULL)
    }
    items <- if (length(supplied_names)) {
      supplied_names[[1L]]
    } else {
      paste0("item_", seq_len(k))
    }
    mu <- rep(mean_focal, k)
    sigma <- rep(sd_focal, k)
    label_reference <- "reference"
    label_focal <- "focal"
  }

  k <- length(factor_loading_reference)
  if (!is.null(item_names)) {
    if (!is.character(item_names) || length(item_names) != k) {
      .dmar_refuse_argument(c(
        "'item_names' must be a character vector with one label per item (",
        k, " expected)."), arg = "item_names")
    }
    items <- item_names
  }

  a <- unname(intercept_reference) - unname(intercept_focal)
  b <- unname(factor_loading_reference) - unname(factor_loading_focal)
  # E[(a + b eta)^2] for eta ~ N(mu, sigma^2); pmax() guards only against a
  # negative value arising from cancellation, never from the mathematics.
  integral <- pmax(a^2 + 2 * a * b * mu + b^2 * (sigma^2 + mu^2), 0)

  out <- data.frame(
    item                     = as.character(items),
    factor_loading_reference = unname(factor_loading_reference),
    factor_loading_focal     = unname(factor_loading_focal),
    intercept_reference      = unname(intercept_reference),
    intercept_focal          = unname(intercept_focal),
    sd_pooled                = unname(sd_pooled),
    dmacs                    = sqrt(integral) / unname(sd_pooled),
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )
  # Every numeric column is a measured value, so a loading or an intercept
  # given as a whole number (an intercept of 2 beside one of 2.3) prints to
  # the decimal places of the others.
  res <- .as_dmar_tbl(out, subclass = "dmar_dmacs",
                      measured_terms = c("factor_loading_reference",
                                         "factor_loading_focal",
                                         "intercept_reference", "intercept_focal",
                                         "sd_pooled", "dmacs"))
  attr(res, "reference") <- label_reference
  attr(res, "focal") <- label_focal
  attr(res, "mean_focal") <- mean_focal
  attr(res, "sd_focal") <- sd_focal
  res
}


# Resolve a group given as a label or an index into a group index, returning
# NA_integer_ when the user did not name the group.
.dmacs_group_index <- function(x, labels, what) {
  if (is.null(x)) return(NA_integer_)
  if (length(x) != 1L) {
    .dmar_refuse_argument(
      sprintf("'%s' must be a single group label or index.", what), arg = what)
  }
  if (is.numeric(x)) {
    i <- as.integer(x)
    if (is.na(i) || i < 1L || i > length(labels)) {
      .dmar_refuse_argument(
        sprintf("'%s' must index one of the %d groups in the fit.",
                what, length(labels)), arg = what)
    }
    return(i)
  }
  i <- match(as.character(x), labels)
  if (is.na(i)) {
    .dmar_refuse_argument(
      sprintf("'%s' must name one of the fit's groups (%s).",
              what, paste(labels, collapse = ", ")), arg = what)
  }
  i
}


# Pull the per-item loadings, intercepts, pooled standard deviations, and the
# focal group's latent mean and standard deviation out of a multiple group
# lavaan fit. Called only by dmacs().
.dmacs_from_fit <- function(fit, reference, focal) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required when 'fit' is supplied. Install it ",
      "with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }
  if (!inherits(fit, "lavaan")) {
    .dmar_refuse_argument("'fit' must be a fitted lavaan object.", arg = "fit")
  }
  # A fit whose optimizer did not converge holds the point where it
  # stopped, not a solution (.lavaan_check_converged(),
  # R/no_maximum_internals.R). Through DMAR 1.0.1 and the first builds of
  # 1.1.0 its loadings and intercepts were read as any other: the partial
  # invariance model of the example, stopped after one iteration, gave the
  # cubes and lozenges dMACS values of .029 and .151, where the converged
  # fit gives .175 and .461.
  .lavaan_check_converged(fit, "dmacs()")
  n_groups <- as.integer(lavaan::lavInspect(fit, "ngroups"))
  if (n_groups < 2L) {
    .dmar_refuse_argument(c(
      "'fit' must be a multiple group model; dMACS compares an item's ",
      "measurement equation across two groups."), arg = "fit")
  }
  labels <- as.character(lavaan::lavInspect(fit, "group.label"))

  i_ref <- .dmacs_group_index(reference, labels, "reference")
  i_foc <- .dmacs_group_index(focal, labels, "focal")
  if (is.na(i_ref) && is.na(i_foc)) {
    if (n_groups > 2L) {
      .dmar_refuse_argument(c(
        "'fit' has ", n_groups, " groups (", paste(labels, collapse = ", "),
        "). Name the two to compare with 'reference' and 'focal'."),
        arg = c("fit", "reference", "focal"))
    }
    i_ref <- 1L
    i_foc <- 2L
  } else if (is.na(i_ref) || is.na(i_foc)) {
    if (n_groups > 2L) {
      .dmar_refuse_argument(c(
        "'fit' has ", n_groups, " groups (", paste(labels, collapse = ", "),
        "). Name both 'reference' and 'focal'."),
        arg = c("fit", "reference", "focal"))
    }
    if (is.na(i_ref)) i_ref <- setdiff(1:2, i_foc) else i_foc <- setdiff(1:2, i_ref)
  }
  if (i_ref == i_foc) {
    .dmar_refuse_argument(
      "'reference' and 'focal' must be different groups.",
      arg = c("reference", "focal"))
  }

  est <- lavaan::lavInspect(fit, "est")
  par_ref <- est[[i_ref]]
  par_foc <- est[[i_foc]]
  if (is.null(par_ref$lambda) || is.null(par_foc$lambda) ||
      ncol(par_ref$lambda) == 0L) {
    .dmar_refuse_argument(
      "'fit' contains no measurement (=~) loadings.", arg = "fit")
  }
  if (is.null(par_ref$nu) || is.null(par_foc$nu)) {
    .dmar_refuse_argument(c(
      "'fit' has no mean structure, so the intercepts dMACS needs are not ",
      "estimated. Refit with meanstructure = TRUE."), arg = "fit")
  }
  lam_ref <- par_ref$lambda
  lam_foc <- par_foc$lambda
  if (!identical(dim(lam_ref), dim(lam_foc)) ||
      !identical(rownames(lam_ref), rownames(lam_foc)) ||
      !identical(colnames(lam_ref), colnames(lam_foc))) {
    .dmar_refuse_argument(c(
      "The two groups' loading matrices do not align; the items or the ",
      "latent variables differ across groups."), arg = NULL)
  }

  nonzero <- (lam_ref != 0) | (lam_foc != 0)
  n_loadings <- rowSums(nonzero)
  keep <- n_loadings > 0
  if (!any(keep)) {
    .dmar_refuse_argument(
      "No observed variable in 'fit' loads on a latent variable.", arg = "fit")
  }
  if (any(n_loadings[keep] > 1L)) {
    .dmar_refuse_argument(c(
      "Every item must load on exactly one latent variable; ",
      paste(rownames(lam_ref)[n_loadings > 1L], collapse = ", "),
      " loads on more than one. Supply the parameter vectors directly for ",
      "a cross-loading item."), arg = NULL)
  }
  items <- rownames(lam_ref)[keep]
  which_lv <- apply(nonzero[keep, , drop = FALSE], 1L,
                    function(z) which(z)[1L])
  idx <- cbind(seq_along(items), which_lv)

  factor_loading_reference <- lam_ref[keep, , drop = FALSE][idx]
  factor_loading_focal <- lam_foc[keep, , drop = FALSE][idx]
  intercept_reference <- par_ref$nu[items, 1L]
  intercept_focal <- par_foc$nu[items, 1L]

  lv_names <- colnames(lam_ref)
  alpha_focal <- if (is.null(par_foc$alpha)) {
    stats::setNames(rep(0, length(lv_names)), lv_names)
  } else {
    stats::setNames(as.numeric(par_foc$alpha[, 1L]), lv_names)
  }
  psi_focal <- stats::setNames(diag(as.matrix(par_foc$psi)), lv_names)
  if (any(psi_focal < 0)) {
    .dmar_refuse_argument(c(
      "The focal group's latent variance is negative, so its standard ",
      "deviation is undefined. Respecify the model."), arg = NULL)
  }
  sd_focal_lv <- sqrt(psi_focal)

  moments <- .dmacs_pooled_sd(fit, items, i_ref, i_foc)

  list(
    factor_loading_reference = factor_loading_reference,
    factor_loading_focal     = factor_loading_focal,
    intercept_reference      = unname(intercept_reference),
    intercept_focal          = unname(intercept_focal),
    sd_pooled                = moments,
    mu                       = unname(alpha_focal[which_lv]),
    sigma                    = unname(sd_focal_lv[which_lv]),
    items                    = items,
    label_reference          = labels[i_ref],
    label_focal              = labels[i_foc],
    mean_focal               = alpha_focal,
    sd_focal                 = sd_focal_lv
  )
}


# Pooled observed standard deviation of each item across the two groups,
# preferring the raw data in the fit and falling back to its sample moments.
.dmacs_pooled_sd <- function(fit, items, i_ref, i_foc) {
  raw <- try(lavaan::lavInspect(fit, "data"), silent = TRUE)
  have_raw <- !inherits(raw, "try-error") && is.list(raw) &&
    is.matrix(raw[[i_ref]]) && is.matrix(raw[[i_foc]]) &&
    all(items %in% colnames(raw[[i_ref]])) &&
    all(items %in% colnames(raw[[i_foc]]))

  if (have_raw) {
    moments <- lapply(c(i_ref, i_foc), function(g) {
      X <- raw[[g]][, items, drop = FALSE]
      list(n = colSums(!is.na(X)),
           s2 = apply(X, 2L, stats::var, na.rm = TRUE))
    })
    n_ref <- moments[[1L]]$n; s2_ref <- moments[[1L]]$s2
    n_foc <- moments[[2L]]$n; s2_foc <- moments[[2L]]$s2
  } else {
    samp <- lavaan::lavInspect(fit, "sampstat")
    n_all <- as.numeric(lavaan::lavInspect(fit, "nobs"))
    n_ref <- rep(n_all[i_ref], length(items))
    n_foc <- rep(n_all[i_foc], length(items))
    # lavaan divides the sample covariance by n under the normal likelihood
    # and by n - 1 under the Wishart likelihood; rescale to the n - 1 divisor.
    likelihood <- lavaan::lavInspect(fit, "options")$likelihood
    scale_ref <- if (identical(likelihood, "normal")) n_ref / (n_ref - 1) else 1
    scale_foc <- if (identical(likelihood, "normal")) n_foc / (n_foc - 1) else 1
    s2_ref <- diag(samp[[i_ref]]$cov)[items] * scale_ref
    s2_foc <- diag(samp[[i_foc]]$cov)[items] * scale_foc
  }

  if (any(n_ref + n_foc - 2 <= 0) || any(!is.finite(s2_ref)) ||
      any(!is.finite(s2_foc))) {
    .dmar_refuse_argument(c(
      "The observed item variances needed for the pooled standard ",
      "deviation could not be computed from 'fit'. Supply 'sd_pooled' ",
      "with the parameter vectors instead."), arg = c("fit", "sd_pooled"))
  }
  sd_pooled <- sqrt(((n_ref - 1) * s2_ref + (n_foc - 1) * s2_foc) /
                      (n_ref + n_foc - 2))
  if (any(sd_pooled <= 0)) {
    .dmar_refuse_argument(c(
      "An item has zero pooled variance in 'fit', so dMACS is undefined ",
      "for it."), arg = "fit")
  }
  unname(sd_pooled)
}
