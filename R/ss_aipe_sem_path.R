#' Sample Size Planning for SEM Targeted Effects
#'
#' Plan sample size for structural equation models so that the confidence interval for the targeted model parameter is sufficiently narrow
#'
#' @param model A single character string giving the free analysis model in
#'   lavaan model syntax (see \code{\link[lavaan]{model.syntax}}). The target
#'   path must carry a parameter label so it can be referred to by name, for
#'   example \code{"f2 ~ b*f1"} labels the structural path \code{b}. This is
#'   the model that would be fit to the data; its parameters are free, not
#'   fixed to population values
#' @param true_cov_matrix Population covariance matrix of the observed
#'   variables planned for, with row and column names matching the observed variables in
#'   \code{model}. It is typically obtained from a fully fixed population model
#'   via \code{\link{cov_sem}}
#' @param width Desired confidence interval width for the model parameter of interest
#' @param which_path The parameter label of the targeted path, given as a
#'   character string, for example \code{"b"} for the path labeled
#'   \code{f2 ~ b*f1} in \code{model}
#' @param conf_level Confidence level (i.e., 1 - Type I error rate)
#' @param assurance The assurance that the confidence interval obtained in a particular study will be no wider than desired (must be \code{NULL} or a value between 0.50 and 1)
#' @param \dots Further arguments passed to \code{\link[lavaan]{sem}} for
#'   the fit to \code{true_cov_matrix}, from which the plan, with or without
#'   \code{assurance}, is read (for example \code{std.lv = TRUE}).
#'   Every name given here must be an argument or an option of \pkg{lavaan}
#'   (see \code{\link[lavaan]{lavOptions}}); any other name, and any
#'   argument without a name, is refused before a model is fitted, with a
#'   condition of class \code{dmar_argument_invalid} that gives the name
#'   as it was typed
#' @param detail if \code{TRUE}, additionally print the model parameter names and the observed variable names (the returned table is unchanged)
#' @param internal option to output a list for internal use (for ss_aipe_sem_path_sensitivity)
#'
#' @details
#' This function implements the sample size planning methods proposed in Lai
#' and Kelley (2011). It requires \pkg{lavaan} to be installed and uses
#' \code{\link[lavaan]{sem}} to obtain the expected information, that is the
#' asymptotic covariance matrix of the parameter estimates, by fitting the
#' free analysis model to the population covariance matrix
#' \code{true_cov_matrix} at a very large sample size. The analysis model is
#' written in lavaan model
#' syntax with the targeted path given a parameter label; see
#' \code{\link[lavaan]{model.syntax}} for the syntax and
#' \code{\link[lavaan]{sem}} for the fitting machinery. The population
#' covariance matrix \code{true_cov_matrix} is most naturally produced by
#' \code{\link{cov_sem}} from a fully fixed population model. An analysis
#' model whose likelihood has no maximum at \code{true_cov_matrix}, an
#' error variance running off below zero while the likelihood is flat or still
#' rising (see
#' \code{\link{fa_bifactor}}), stops the function with an error of class
#' \code{dmar_no_maximum} that prints the \pkg{lavaan} code of the fit,
#' since the asymptotic covariance matrix would describe where the
#' optimizer stopped rather than the model.
#'
#' When \code{assurance} is supplied, the assurance adjustment is based on a chi
#' square approximation to the sampling variability of the confidence interval
#' width and can undershoot the nominal assurance in finite samples; use
#' \code{\link{ss_aipe_sem_path_sensitivity}} to check the realized width and
#' coverage at the planned sample size.
#'
#' @return
#' A \code{data.frame} (a \code{dmar_tbl}) with \code{term} and
#' \code{value} columns whose rows are \code{necessary_N} (the planned sample
#' size), \code{path_index} (the position of the target path among the model
#' parameters), and \code{var_theta_j} (the population sampling variance of the
#' target path at the planned sample size). The returned table is the same
#' whether or not \code{detail = TRUE}. When \code{internal = TRUE} a list is
#' returned for use by \code{\link{ss_aipe_sem_path_sensitivity}}.
#'
#' @references
#' Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
#'   targeted effects in structural equation modeling: Sample size
#'   planning for narrow confidence intervals.
#'   \emph{Psychological Methods, 16}(2), 127--148. \doi{10.1037/a0021764}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#'   \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link[lavaan]{sem}}, \code{\link[lavaan]{model.syntax}},
#' \code{\link{cov_sem}}, \code{\link{ss_aipe_sem_path_sensitivity}}
#'
#' @examples
#' # Population covariance from a fully fixed model (see cov_sem()).
#' pop_model <- "
#'   f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
#'   f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
#'   f2 ~ 0.5*f1
#'   f1 ~~ 1*f1
#'   f2 ~~ 0.75*f2
#'   y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
#'   y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
#' "
#' Sigma <- cov_sem(pop_model)$sigma_theta
#'
#' # Free analysis model with the target structural path labeled "b".
#' analysis_model <- "
#'   f1 =~ y1 + y2 + y3
#'   f2 =~ y4 + y5 + y6
#'   f2 ~ b*f1
#' "
#' ss_aipe_sem_path(model = analysis_model, true_cov_matrix = Sigma,
#'                  width = 0.30, which_path = "b")
#'
#' # That sample size holds the interval to the desired width on average, so
#' # about half of the studies it plans return a wider one. Adding assurance
#' # plans for the width to be met in 90 percent of studies instead, at the
#' # cost of a larger sample size.
#' ss_aipe_sem_path(model = analysis_model, true_cov_matrix = Sigma,
#'                  width = 0.30, which_path = "b",
#'                  assurance = 0.90)
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sem_path <- function(model, true_cov_matrix, width, which_path, conf_level = .95, assurance = NULL, detail = FALSE, internal = FALSE, ...) {
  .dmar_check_lavaan_dots(..., fn = "ss_aipe_sem_path")
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  if (!requireNamespace("lavaan", quietly = TRUE)) .dmar_refuse(
    "The package 'lavaan' is needed; please install the package and try again.",
    "package_required", package = "lavaan")
  if (is.null(rownames(true_cov_matrix)) || is.null(colnames(true_cov_matrix))) {
    .dmar_refuse_argument(
      "'true_cov_matrix' must have row and column names that match the observed variables in 'model'.",
      arg = c("true_cov_matrix", "model"))
  }
  if (!is.null(assurance) &&
      (!is.numeric(assurance) || length(assurance) != 1L || assurance < 0.5 || assurance >= 1)) {
    .dmar_refuse_argument(
      "'assurance' must be NULL or a single value in [0.5, 1).",
      arg = "assurance")
  }

  N <- 1000000
  fit <- lavaan::sem(model, sample.cov = true_cov_matrix, sample.nobs = N, ...)
  # The asymptotic variance of the target is read from this fit. When its
  # likelihood has no maximum at true_cov_matrix, an error variance running off below
  # zero along a ridge with the likelihood flat or still rising
  # (.lavaan_no_maximum(), R/no_maximum_internals.R), that variance is
  # wherever lavaan's optimizer stopped, and the plan is refused with the
  # lavaan code of the fit.
  .no_max_check(fit, "ss_aipe_sem_path", syntax = model, fun = "lavaan::sem",
                args = c(list(model = model, sample.cov = true_cov_matrix,
                              sample.nobs = N), list(...)),
                where = "the analysis model fit to true_cov_matrix",
                placeholder = "the matrix given as true_cov_matrix in place of my_S")
  H <- lavaan::vcov(fit) * N
  cf <- lavaan::coef(fit)

  j <- which(names(cf) == which_path)
  if (length(j) != 1L) .dmar_refuse_argument(c(
    "The path of interest '", which_path, "' is not a labeled parameter in the model; label it in the lavaan syntax, e.g. 'f2 ~ b*f1', and pass which_path = \"b\"."),
    arg = NULL)
  h_jj <- H[j, j]
  p <- length(lavaan::lavNames(fit, "ov"))

  omega <- width
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha / 2)

  if (is.null(assurance)) {
    N <- 4 * z^2 * h_jj / omega^2
    N <- ceiling(N)
  }

  if (!is.null(assurance)) {
    # The sample size that holds the expected width, from this fit, which
    # takes the arguments given in '...'.
    N0 <- ceiling(4 * z^2 * h_jj / omega^2)
    gamma <- assurance
    signal <- TRUE
    N <- N0
    # print(cat(N, N0))
    while (signal) {
      c <- N0 * qchisq(p = gamma, df = N - 1)
      N_1 <- (1 + sqrt(1 + 4 * c)) / 2
      if (abs(N_1 - N) > 1e-5) {
        N <- N_1
      } else {
        signal <- FALSE
      }
      # print(cat(N, N_1, "\n"))
    } # end of while()
  } # end of if(!is.null(assurance))

  N <- ceiling(N)

  if (internal) {
    return(list(
      necessary_N = N, parameters = names(cf), obs_vars = lavaan::lavNames(fit, "ov"),
      path_index = j, var_theta_j = h_jj / N
    ))
  } else {
    if (detail) {
      print(paste0(c("parameters:", names(cf)), collapse = " "))
      print(paste0(c("obs_vars:", lavaan::lavNames(fit, "ov")), collapse = " "))
      return(.as_dmar_tbl(data.frame(term = c("necessary_N", "path_index", "var_theta_j"), value = c(N, j, h_jj / N)), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    } else {
      return(.as_dmar_tbl(data.frame(term = c("necessary_N", "path_index", "var_theta_j"), value = c(N, j, h_jj / N)), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
  }
}
