#' Model Implied Covariance Matrix From a Lavaan-Specified SEM
#'
#' @description
#' Given a structural equation model written in lavaan model syntax with all
#' of its parameters fixed to their population values, compute the model
#' implied population covariance matrix \eqn{\Sigma(\theta)} of the observed
#' variables and, when the model has a mean structure, the model implied
#' population mean vector \eqn{\mu(\theta)}. This function requires
#' \pkg{lavaan} to be installed.
#'
#' This is the helper that drives the \emph{population} side of the sample
#' size planning workflow for SEM: it lets the user state a population
#' model, obtain the \eqn{\Sigma(\theta)} (and \eqn{\mu(\theta)}) those
#' fixed values imply, and then pass that population to
#' \code{\link{ss_aipe_sem_path}},
#' \code{\link{ss_aipe_sem_path_sensitivity}},
#' \code{\link{ss_aipe_rmsea_sensitivity}},
#' \code{\link{ss_power_composite_sem}}, or
#' \code{\link{ss_aipe_composite_sem}}.
#'
#' @param model A single character string giving a structural equation model
#'   in lavaan model syntax (see \code{\link[lavaan]{model.syntax}}), with
#'   every parameter fixed to its population value. A factor loading, a
#'   structural path, a variance, or a covariance is fixed by prefixing the
#'   numeric value to the variable with the \code{*} operator, for example
#'   \code{"f1 =~ 1*y1 + 0.8*y2 + 0.8*y3"} for the loadings,
#'   \code{"f2 ~ 0.5*f1"} for a structural path, and \code{"y1 ~~ 0.5*y1"}
#'   for a residual variance. A model with a mean structure (for example a
#'   latent growth curve model) also fixes every intercept and latent mean,
#'   for example \code{"t1 ~ 0*1"} and \code{"s ~ 0.3*1"}. lavaan syntax
#'   embeds the values and knows which variables are latent, so no separate
#'   parameter vector or list of latent variables is needed.
#'
#' @details
#' The function builds a non-fitted lavaan object from \code{model} with all
#' parameters held at the population values written into the syntax, and reads
#' back the model implied covariance matrix of the observed variables. Because
#' the object is created with \code{do.fit = FALSE}, no estimation is performed
#' and the placeholder sample covariance lavaan needs to construct the object
#' is never used; the returned \eqn{\Sigma(\theta)} comes entirely from the
#' fixed parameter values. The observed-variable names are taken from the model
#' syntax and fix the row and column order of the returned matrix.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{sigma_theta}}{The model implied population covariance matrix
#'     of the observed variables, with rows and columns named.}
#'   \item{\code{mu_theta}}{The model implied population mean vector of the
#'     observed variables, named, in the row order of \code{sigma_theta}. A
#'     vector of zeros when the model has no mean structure.}
#'   \item{\code{observed_vars}}{Character vector of observed variable names
#'     in the row/column order of \code{sigma_theta}.}
#' }
#'
#' @references
#' Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for targeted
#' effects in structural equation modeling: Sample size planning for narrow
#' confidence intervals. \emph{Psychological Methods, 16}(2), 127--148.
#'   \doi{10.1037/a0021764}
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' \emph{Journal of Statistical Software, 48}(2), 1--36.
#' \doi{10.18637/jss.v048.i02}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[lavaan]{sem}}, \code{\link[lavaan]{model.syntax}},
#'   \code{\link{ss_aipe_sem_path}},
#'   \code{\link{ss_aipe_sem_path_sensitivity}},
#'   \code{\link{ss_aipe_rmsea_sensitivity}}, \code{\link{covmat_from_cfa}}.
#'
#' @examples
#' # Population model with all parameters fixed to their values: two factors,
#' # three indicators each, and a structural path f2 ~ f1 of 0.5.
#' pop_model <- "
#'   f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#'   f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#'   f2 ~ 0.5*f1
#'   f1 ~~ 1*f1
#'   f2 ~~ 0.75*f2
#'   y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#'   y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#' "
#' cov_sem(pop_model)$sigma_theta
#'
#' # A population model with a mean structure: a linear latent growth curve
#' # over four waves. The intercepts and latent means are fixed too, and
#' # mu_theta carries the model implied means (5.0, 5.3, 5.6, 5.9).
#' pop_lgm <- "
#'   i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
#'   s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
#'   i ~~ 1*i
#'   s ~~ 0.2*s
#'   i ~~ -0.15*s
#'   t1 ~~ 0.5*t1; t2 ~~ 0.5*t2; t3 ~~ 0.5*t3; t4 ~~ 0.5*t4
#'   t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
#'   i ~ 5*1
#'   s ~ 0.3*1
#' "
#' cov_sem(pop_lgm)$mu_theta
#'
#' @keywords multivariate
#'
#' @export
cov_sem <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string of lavaan model syntax.",
         call. = FALSE)
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The package 'lavaan' is needed; please install it and try again.",
         call. = FALSE)
  }

  # auto.var = TRUE makes lavaan add the variances the syntax implies but
  # does not state. Without it those rows arrive fixed at zero with
  # free = 0, so the guard below sees nothing free and a model missing a
  # variance yields a silently wrong Sigma rather than an error.
  pt <- lavaan::lavaanify(model, auto.var = TRUE)
  observed_vars <- lavaan::lavNames(pt, type = "ov")
  if (length(observed_vars) == 0L) {
    stop("No observed variables found in 'model'.", call. = FALSE)
  }

  # cov_sem() returns the population covariance implied by a fully specified
  # model, so every parameter must be fixed to its population value. A free
  # parameter (including a variance or covariance left unspecified, which
  # lavaan adds as free) would otherwise be silently filled with a lavaan
  # start value, yielding a wrong Sigma. Catch that here.
  free_params <- pt[pt$free > 0L, , drop = FALSE]
  if (nrow(free_params) > 0L) {
    free_terms <- unique(paste(free_params$lhs, free_params$op, free_params$rhs))
    stop("'model' must have every parameter fixed to its population value; ",
         "these parameters are free: ", paste(free_terms, collapse = ", "),
         ". Fix each one with the 'value*term' syntax (for example '0.5*f1', ",
         "or '0*1' for an intercept), including every variance and ",
         "covariance, and every intercept and latent mean when the model ",
         "has a mean structure.", call. = FALSE)
  }

  dummy <- diag(length(observed_vars))
  dimnames(dummy) <- list(observed_vars, observed_vars)

  # A model with a mean structure needs a placeholder sample mean vector to
  # construct the (non-fitted) object, just as the identity matrix stands in
  # for the sample covariance; neither placeholder is used.
  fit_args <- list(model, sample.cov = dummy, sample.nobs = 1000L,
                   do.fit = FALSE)
  if (any(pt$op == "~1")) {
    fit_args$sample.mean <- stats::setNames(rep(0, length(observed_vars)),
                                            observed_vars)
  }
  fit <- do.call(lavaan::sem, fit_args)
  sigma_theta <- lavaan::lavInspect(fit, "cov.ov")[observed_vars,
                                                   observed_vars,
                                                   drop = FALSE]
  mu_theta <- lavaan::lavInspect(fit, "mean.ov")
  mu_theta <- if (length(mu_theta)) {
    mu_theta[observed_vars]
  } else {
    stats::setNames(rep(0, length(observed_vars)), observed_vars)
  }

  list(
    sigma_theta   = sigma_theta,
    mu_theta      = mu_theta,
    observed_vars = observed_vars
  )
}
