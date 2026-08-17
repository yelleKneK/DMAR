#' Approximate Measurement Invariance by Factor Alignment
#'
#' Estimates the group factor means and factor variances that make the
#' measurement parameters as nearly invariant as possible across groups,
#' following the alignment method of Asparouhov and Muthén (2014). The
#' invariance ladder in \code{\link{measurement_invariance}} asks whether
#' loadings and intercepts are exactly equal across groups. With many
#' groups that hypothesis is essentially never true, the ladder stalls at
#' configural invariance, and the comparison of factor means the researcher
#' wanted never happens. Alignment takes the configural solution as given
#' and searches for the group factor means and variances that concentrate
#' the noninvariance in a few parameters instead of spreading it thinly
#' over many, so that factor means stay comparable without claiming that
#' exact invariance holds. Requires \pkg{lavaan}.
#'
#' @param data A \code{data.frame} holding the items and the grouping
#'   variable.
#' @param items Character vector naming the indicator columns of
#'   \code{data}. Three or more items are required, and each must be
#'   numeric.
#' @param group Single character string naming the grouping column of
#'   \code{data}. Two or more groups are required.
#' @param model Optional \pkg{lavaan} model syntax for the measurement
#'   model, for the case where the default single factor over \code{items}
#'   is not what is wanted (a residual covariance, for example). The syntax
#'   must define exactly one latent variable and its indicators must be
#'   exactly \code{items}. The factor is standardized (mean 0, variance 1)
#'   in every group by this function, so do not set its scale in the
#'   syntax. When \code{NULL} (the default) the model
#'   \code{f =~ item1 + item2 + ...} is used.
#' @param alignment Identification rule for the metric of the factor,
#'   either \code{"fixed"} (the default) or \code{"free"}. Under
#'   \code{"fixed"} the first group's factor mean is held at 0 and its
#'   factor variance at 1. Under \code{"free"} the group factor means are
#'   constrained to average 0 and the group factor standard deviations to
#'   have a geometric mean of 1, so no group serves as the reference. The
#'   choice sets the metric the factor means and variances are reported on;
#'   see Details.
#' @param estimator Estimator passed to \pkg{lavaan} for the configural
#'   model. Defaults to \code{"ML"}; \code{"MLR"} supplies the robust
#'   corrections.
#' @param n_starts Number of starting values for the optimizer, a positive
#'   integer (default 10). The first start sets every factor mean to 0 and
#'   every factor variance to 1, the second uses a median based heuristic,
#'   and any remaining starts are random. The simplicity function has local
#'   minima, so several starts is the default rather than a refinement.
#' @param seed Optional integer seed for the random starts. Defaults to
#'   \code{NULL}, which uses the caller's random number state and leaves it
#'   alone. When an integer is supplied the seed is set locally and the
#'   caller's state is restored on exit.
#' @param epsilon Smoothing constant \eqn{\epsilon > 0} inside the
#'   component loss function, default 0.01. Smaller values push the loss
#'   closer to \eqn{\sqrt{|x|}} and make the surface harder to optimize;
#'   larger values smooth it at the cost of blurring the distinction
#'   between a few large differences and many small ones.
#'
#' @details
#' Alignment starts from the \emph{configural} solution: the multiple group
#' single factor model with every loading, intercept, and residual variance
#' free across groups and the factor standardized (mean 0 and variance 1)
#' in each group. Write \eqn{\lambda^{0}_{gi}} and \eqn{\nu^{0}_{gi}} for
#' the resulting loading and intercept of item \eqn{i} in group \eqn{g}.
#' That solution is not unique: for any group factor means \eqn{\alpha_g}
#' and factor variances \eqn{\psi_g} the reparameterized measurement
#' parameters
#' \deqn{\lambda_{gi} = \lambda^{0}_{gi} / \sqrt{\psi_g}}
#' \deqn{\nu_{gi} = \nu^{0}_{gi} - \alpha_g \lambda^{0}_{gi} / \sqrt{\psi_g}}
#' reproduce the observed means and covariances exactly and so fit the data
#' identically. Alignment picks the member of that family whose measurement
#' parameters are closest to invariant, by minimizing the total simplicity
#' function
#' \deqn{F = \sum_i \sum_{g_1 < g_2} w_{g_1 g_2}
#'            f(\lambda_{g_1 i} - \lambda_{g_2 i}) +
#'           \sum_i \sum_{g_1 < g_2} w_{g_1 g_2}
#'            f(\nu_{g_1 i} - \nu_{g_2 i})}
#' over \eqn{\alpha_g} and \eqn{\psi_g}, with weights
#' \eqn{w_{g_1 g_2} = \sqrt{N_{g_1} N_{g_2}}} and the component loss
#' function
#' \deqn{f(x) = \sqrt{\sqrt{x^2 + \epsilon}} = (x^2 + \epsilon)^{1/4}.}
#'
#' The fourth root is the substance of the method, not a technical detail.
#' A squared loss would spread a fixed amount of noninvariance evenly over
#' all the parameters, because halving two differences beats zeroing one.
#' The fourth root is concave in \eqn{|x|}, so the marginal penalty falls as
#' a difference grows: the criterion prefers a solution in which most
#' parameters agree closely and a few disagree substantially, which is
#' exactly the pattern of approximate invariance a researcher wants to find
#' and report. The constant \eqn{\epsilon} only rounds off the kink of
#' \eqn{\sqrt{|x|}} at zero so the criterion is differentiable. One
#' consequence is worth knowing: each of the \eqn{2I} terms in a group pair
#' contributes at least \eqn{\epsilon^{1/4}}, so the smallest attainable
#' value of \eqn{F} is
#' \eqn{2 I \epsilon^{1/4} \sum_{g_1 < g_2} w_{g_1 g_2}}, reached only when
#' every aligned parameter is exactly invariant. The achieved value is
#' comparable across runs on the same items, groups, and \eqn{\epsilon},
#' not across data sets.
#'
#' Two constraints identify the solution. The scale of the factor is
#' genuinely undetermined by \eqn{F}: multiplying every
#' \eqn{\sqrt{\psi_g}} and every \eqn{\alpha_g} by the same constant leaves
#' the aligned intercepts alone and shrinks the aligned loadings toward
#' each other, so without a constraint the criterion is minimized by
#' letting the factor variances run away. The location constraint plays the
#' same role in the limiting case of exact invariance, where shifting all
#' the factor means by a constant leaves \eqn{F} unchanged. Under
#' \code{alignment = "fixed"} the constraints are \eqn{\alpha_1 = 0} and
#' \eqn{\psi_1 = 1}. Under \code{alignment = "free"} they are
#' \eqn{\sum_g \alpha_g = 0} and \eqn{\prod_g \sqrt{\psi_g} = 1}, which
#' treats the groups symmetrically and is the more natural choice when no
#' group is a meaningful reference. The two rules give genuinely different
#' solutions, not a relabeling of one another, because \eqn{F} is not
#' invariant to shifting the factor means.
#'
#' The rule also sets the metric the answer is reported on, which matters
#' when the estimates are compared with anything else. The constraint
#' \eqn{\psi_1 = 1} makes the first group's factor the common metric, so a
#' reported factor mean of 0.4 says that group's factor mean is 0.4 of the
#' \emph{first group's} factor standard deviations above the first group's,
#' and a reported factor variance of 1.3 says its factor variance is 1.3
#' times the first group's. Under \code{"free"} the common metric is the
#' one in which the group factor standard deviations have a geometric mean
#' of 1, and the factor means are deviations from their own average on that
#' metric. Simulating data with known group factor means and then comparing
#' them with the estimates requires dividing the generating means by the
#' generating factor standard deviation of the reference group first.
#'
#' The minimization runs \code{\link[stats]{optim}} with the BFGS method
#' and the analytic gradient, over the factor means \eqn{\alpha_g} and the
#' log factor standard deviations \eqn{\log \sqrt{\psi_g}}, reduced to the
#' coordinates the identification rule leaves free. The optimizer is run
#' from \code{n_starts} starting values and the best solution is kept,
#' because a single start is not safe. Two things can go wrong, and they are
#' different problems.
#'
#' The first is ordinary multimodality. When several loadings and intercepts
#' are noninvariant the criterion has more than one local minimum, typically
#' within a percent or so of each other in value but at different factor
#' means. The number of distinct minima the starts reached is returned in
#' the \code{"n_optima"} attribute and the value achieved from each start in
#' \code{"simplicity_starts"}. More than one is a signal to raise
#' \code{n_starts} and to check whether the reported solution is stable.
#'
#' The second is a degenerate branch, and it is the one that bites. Sending
#' the factor variances of every group but the reference off to infinity
#' drives those groups' aligned loadings to zero, which flattens the loading
#' half of the criterion; the optimizer stops there and reports convergence
#' with factor variances of \eqn{10^{36}} or \code{Inf}. That is not an
#' estimate, and on real data a meaningful share of random starts finds it.
#' A start whose factor standard deviations or factor means leave a very
#' wide sanity range (a factor of \eqn{10^4} either way) is therefore
#' discarded and its entry in \code{"simplicity_starts"} is \code{Inf}. If
#' every start is discarded the function stops rather than return the
#' degenerate solution.
#'
#' Item level invariance is summarized with the \eqn{R^2} measure of
#' Asparouhov and Muthén (2014), which asks how much of the group to group
#' variation in an item's configural parameter is accounted for by the
#' estimated factor means and variances alone. Let \eqn{\bar{\lambda}_i}
#' and \eqn{\bar{\nu}_i} be the aligned parameters of item \eqn{i} averaged
#' over groups. The factor means and variances by themselves imply the
#' configural values \eqn{\sqrt{\psi_g} \bar{\lambda}_i} and
#' \eqn{\bar{\nu}_i + \alpha_g \bar{\lambda}_i}. Write \eqn{T_i} for the
#' sum over groups of the squared configural parameter of item \eqn{i} and
#' \eqn{E_i} for the sum over groups of its squared departure from that
#' implied value. Then
#' \deqn{R^2_i = 1 - E_i / T_i,}
#' computed separately for the loadings and for the intercepts. A value of
#' 1 for every item on the loadings is metric invariance, and on both
#' loadings and intercepts is scalar invariance. The complementary view is
#' the \code{"item_loss"} attribute, which splits the achieved simplicity
#' function into the contribution of each item, so the items carrying the
#' noninvariance can be named.
#'
#' The criterion adds up raw parameter differences, so items on very
#' different measurement scales do not contribute equally: an item whose
#' raw variance is a hundred times another's dominates the sum. When the
#' items are not already on a common metric, put them on one (standardizing
#' them is the simple choice) before aligning.
#'
#' The method is defined for two groups and is computed here for two, but
#' it has little to offer there. With two groups the standard invariance
#' ladder is tractable and interpretable, and the alignment criterion has
#' only one group pair to work with. Alignment earns its keep when the
#' number of groups makes exact invariance implausible and the ladder
#' uninformative.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   group and columns \code{group} (the group label), \code{n} (the number
#'   of cases the configural model used in that group),
#'   \code{factor_mean} (the estimated \eqn{\alpha_g}), and
#'   \code{factor_variance} (the estimated \eqn{\psi_g}). The rows follow
#'   the group ordering \pkg{lavaan} uses.
#'
#'   Attributes carry the rest of the solution:
#'   \describe{
#'     \item{\code{aligned_loadings}, \code{aligned_intercepts}}{Matrices
#'       of the aligned \eqn{\lambda_{gi}} and \eqn{\nu_{gi}}, groups in
#'       rows and items in columns.}
#'     \item{\code{configural_loadings},
#'       \code{configural_intercepts}}{Matrices of the configural
#'       \eqn{\lambda^{0}_{gi}} and \eqn{\nu^{0}_{gi}}, same layout.}
#'     \item{\code{simplicity_function}}{The achieved minimum of \eqn{F}.}
#'     \item{\code{simplicity_starts}}{The value of \eqn{F} achieved from
#'       each starting value, in start order. A start the optimizer could
#'       not complete, or one that ended in the degenerate branch described
#'       in Details, is recorded as \code{Inf} and discarded.}
#'     \item{\code{converged}}{\code{TRUE} when the optimizer reported
#'       convergence at the retained solution.}
#'     \item{\code{n_starts}}{The number of starting values used.}
#'     \item{\code{n_optima}}{The number of distinct local minima the
#'       converged starts reached.}
#'     \item{\code{alignment}}{The identification rule used,
#'       \code{"fixed"} or \code{"free"}.}
#'     \item{\code{epsilon}}{The smoothing constant used.}
#'     \item{\code{R2_loadings}, \code{R2_intercepts}}{Named numeric
#'       vectors, one entry per item, holding the per item \eqn{R^2}
#'       invariance measures defined in Details.}
#'     \item{\code{R2_total}}{The same two measures pooled over items, as a
#'       named numeric vector with elements \code{loadings} and
#'       \code{intercepts}. This is the overall effect size of approximate
#'       invariance reported by Asparouhov and Muthén (2014).}
#'     \item{\code{item_loss}}{Matrix with one row per item and columns
#'       \code{loadings}, \code{intercepts}, and \code{total}, splitting
#'       the achieved simplicity function into per item contributions.}
#'     \item{\code{fit}}{The configural \pkg{lavaan} fit object.}
#'   }
#'
#' @references
#' Asparouhov, T., & Muthén, B. (2014). Multiple-group factor analysis
#'   alignment. \emph{Structural Equation Modeling, 21}(4), 495--508.
#'   \doi{10.1080/10705511.2014.919210}
#'
#' Marsh, H. W., Guo, J., Parker, P. D., Nagengast, B., Asparouhov, T.,
#'   Muthén, B., & Dicke, T. (2018). What to do when scalar invariance
#'   fails: The extended alignment method for multi-group factor analysis
#'   comparison of latent means across many groups. \emph{Psychological
#'   Methods, 23}(3), 524--545. \doi{10.1037/met0000113}
#'
#' Muthén, B., & Asparouhov, T. (2014). IRT studies of many groups: The
#'   alignment method. \emph{Frontiers in Psychology, 5}, Article 978.
#'   \doi{10.3389/fpsyg.2014.00978}
#'
#' Muthén, B., & Asparouhov, T. (2018). Recent methods for the study of
#'   measurement invariance with many groups: Alignment and random effects.
#'   \emph{Sociological Methods & Research, 47}(4), 637--664.
#'   \doi{10.1177/0049124117701488}
#'
#' Robitzsch, A. (2025). \emph{sirt: Supplementary item response theory
#'   models}. R package version 4.2-133.
#'   \url{https://CRAN.R-project.org/package=sirt}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{measurement_invariance}} for the exact invariance
#'   ladder alignment is meant to rescue; \code{\link{cfa_1}} for the
#'   single group measurement model; \code{\link{htmt}} for discriminant
#'   validity of the same items.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords models multivariate
#'
#' @examples
#' # Five groups that differ in factor mean and factor variance, with two
#' # deliberately noninvariant measurement parameters (the loading of item
#' # 3 in group 2 and the intercept of item 5 in group 4). The first group's
#' # factor standard deviation is 1, which is the metric the default "fixed"
#' # rule reports on, so the estimates compare directly with the generating
#' # values below.
#' set.seed(113)
#' n_g <- c(400, 450, 380, 500, 420)
#' factor_mean <- c(0, 0.30, -0.50, 0.80, 0.20)
#' factor_sd   <- c(1, 1.20,  0.80, 1.10, 0.90)
#' Lambda <- matrix(0.8, nrow = 5, ncol = 6)
#' Nu     <- matrix(1.0, nrow = 5, ncol = 6)
#' Lambda[2, 3] <- 0.3
#' Nu[4, 5]     <- 1.8
#' d <- do.call(rbind, lapply(1:5, function(g) {
#'   eta <- rnorm(n_g[g], factor_mean[g], factor_sd[g])
#'   x <- sapply(1:6, function(i)
#'     Nu[g, i] + Lambda[g, i] * eta + rnorm(n_g[g], 0, 0.6))
#'   data.frame(x, cohort = paste0("cohort_", g))
#' }))
#' names(d)[1:6] <- paste0("x", 1:6)
#'
#' out <- measurement_alignment(d, items = paste0("x", 1:6),
#'                              group = "cohort", seed = 113)
#' out                                   # recovered means and variances
#' attr(out, "R2_loadings")              # item 3 stands out
#' attr(out, "R2_intercepts")            # item 5 stands out
#' attr(out, "item_loss")
#'
#' # Holzinger and Swineford's verbal tests across four groups formed by
#' # crossing school with sex. The raw tests are on very different scales,
#' # so they are standardized first (see Details).
#' data(holzinger_swineford)
#' hs <- holzinger_swineford
#' hs$school_sex <- interaction(hs$school, hs$sex, sep = ", ")
#' verbal <- c("t5_general_information", "t6_paragraph_comprehension",
#'             "t7_sentence", "t8_word_classification", "t9_word_meaning")
#' hs[verbal] <- scale(hs[verbal])
#' ma <- measurement_alignment(hs, items = verbal, group = "school_sex",
#'                             seed = 113)
#' ma
#'
#' # The broom verbs: one row per group, and the alignment summary.
#' generics::tidy(ma)
#' generics::glance(ma)
#'
#' @export
measurement_alignment <- function(data, items, group, model = NULL,
                                  alignment = c("fixed", "free"),
                                  estimator = "ML", n_starts = 10,
                                  seed = NULL, epsilon = 0.01) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for measurement_alignment(). ",
         "Install it with install.packages(\"lavaan\").", call. = FALSE)
  }
  alignment <- match.arg(alignment)
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if (!is.character(items) || length(items) < 3L ||
      anyNA(items) || !all(items %in% names(data))) {
    stop("'items' must name three or more columns of 'data'.", call. = FALSE)
  }
  if (anyDuplicated(items)) {
    stop("'items' must not repeat an item.", call. = FALSE)
  }
  if (!all(vapply(data[items], is.numeric, logical(1L)))) {
    stop("Every column named in 'items' must be numeric.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) ||
      !(group %in% names(data))) {
    stop("'group' must name one column of 'data'.", call. = FALSE)
  }
  if (length(unique(stats::na.omit(data[[group]]))) < 2L) {
    stop("'group' must have at least two groups; alignment compares ",
         "measurement parameters across groups.", call. = FALSE)
  }
  if (!is.null(model) && (!is.character(model) || length(model) != 1L)) {
    stop("'model' must be a single string of lavaan model syntax, or NULL.",
         call. = FALSE)
  }
  if (!is.numeric(n_starts) || length(n_starts) != 1L || is.na(n_starts) ||
      n_starts < 1 || n_starts != round(n_starts)) {
    stop("'n_starts' must be a single positive integer.", call. = FALSE)
  }
  n_starts <- as.integer(n_starts)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) ||
      epsilon <= 0) {
    stop("'epsilon' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         is.na(seed))) {
    stop("'seed' must be a single number, or NULL.", call. = FALSE)
  }

  # ---- configural solution ------------------------------------------------
  # All measurement parameters free by group, the factor standardized in
  # every group. std.lv = TRUE fixes each group's factor variance to 1 and,
  # with no cross-group equality constraints, lavaan fixes each group's
  # factor mean to 0; both are verified below rather than assumed.
  if (is.null(model)) {
    model <- paste("f =~", paste(items, collapse = " + "))
  }
  fit <- lavaan::cfa(model = model, data = data, group = group,
                     estimator = estimator, std.lv = TRUE,
                     meanstructure = TRUE)
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop("The configural model did not converge, so there is nothing to ",
         "align. Check the items and the group sizes.", call. = FALSE)
  }
  est <- lavaan::lavInspect(fit, "est")
  group_label <- as.character(lavaan::lavInspect(fit, "group.label"))
  n <- as.numeric(lavaan::lavInspect(fit, "nobs"))
  G <- length(group_label)
  if (G < 2L) {
    stop("The fitted configural model has fewer than two groups.",
         call. = FALSE)
  }
  if (ncol(est[[1L]]$lambda) != 1L) {
    stop("'model' must define exactly one latent variable; alignment is ",
         "defined for a single factor measurement model.", call. = FALSE)
  }
  if (!setequal(rownames(est[[1L]]$lambda), items)) {
    stop("The indicators in 'model' must be exactly the columns named in ",
         "'items'.", call. = FALSE)
  }
  # The psi/alpha check above only sees the factor variance and mean, so a
  # scale set through the loadings slips past it: a marker constraint
  # (f =~ 1*item1) or an equality label (c(a, a)*item1) leaves psi at 1 and
  # alpha at 0 in every group yet pins loadings the alignment needs free,
  # and the aligned means and variances shift materially. Both leave
  # fingerprints in the parameter table: a fixed loading row, or an
  # equality-constraint row.
  pt <- lavaan::parameterTable(fit)
  if (any(pt$op == "=~" & pt$free == 0L) || any(pt$op == "==")) {
    stop("'model' fixes or constrains a loading (a marker constraint such ",
         "as f =~ 1*item1, or an equality label). The alignment requires ",
         "every loading free in every group; the factor is standardized ",
         "(mean 0, variance 1) by this function instead. Remove the ",
         "constraint from 'model'.", call. = FALSE)
  }
  I <- length(items)
  lambda_0 <- t(vapply(est, function(z) z$lambda[items, 1L], numeric(I)))
  nu_0 <- t(vapply(est, function(z) z$nu[items, 1L], numeric(I)))
  dimnames(lambda_0) <- dimnames(nu_0) <- list(group_label, items)
  psi_g <- vapply(est, function(z) as.numeric(z$psi[1L, 1L]), numeric(1L))
  alpha_g <- vapply(est, function(z) as.numeric(z$alpha[1L, 1L]), numeric(1L))
  if (any(abs(psi_g - 1) > 1e-8) || any(abs(alpha_g) > 1e-8)) {
    stop("The configural model must standardize the factor in every group ",
         "(mean 0 and variance 1). Remove any scaling of the factor from ",
         "'model'.", call. = FALSE)
  }
  if (!all(is.finite(lambda_0)) || !all(is.finite(nu_0))) {
    stop("The configural solution contains nonfinite loadings or ",
         "intercepts, so the alignment criterion is undefined.",
         call. = FALSE)
  }

  # ---- minimize the simplicity function -----------------------------------
  if (!is.null(seed)) {
    has_old <- exists(".Random.seed", envir = globalenv())
    old <- if (has_old) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
      if (has_old) assign(".Random.seed", old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }
  machinery <- .alignment_machinery(lambda_0, nu_0, n, epsilon, alignment)
  starts <- .alignment_starts(lambda_0, nu_0, n_starts, alignment)
  runs <- lapply(starts, function(par) {
    tryCatch(
      stats::optim(par, machinery$fn, machinery$gr, method = "BFGS",
                   control = list(maxit = 1000L, reltol = 1e-13)),
      error = function(e) list(par = par, value = Inf, convergence = 99L)
    )
  })
  values <- vapply(runs, function(z) z$value, numeric(1L))
  usable <- vapply(runs, function(z) {
    all(is.finite(z$par)) && .alignment_usable(machinery$build(z$par))
  }, logical(1L))
  values[!usable] <- Inf
  if (!any(is.finite(values))) {
    stop("No starting value produced a usable alignment: the optimizer ",
         "either failed or drifted into the degenerate branch of the ",
         "criterion, where the factor variances run off and the aligned ",
         "loadings collapse to zero. Raise 'n_starts', and check that the ",
         "items are on a common scale.", call. = FALSE)
  }
  best <- runs[[which.min(values)]]
  converged <- identical(as.integer(best$convergence), 0L)
  ok <- usable & vapply(runs, function(z) {
    identical(as.integer(z$convergence), 0L)
  }, logical(1L))
  n_optima <- .alignment_n_optima(values[ok])

  # ---- aligned parameters and per item summaries --------------------------
  pars <- machinery$build(best$par)
  aligned <- machinery$aligned(pars)
  lambda <- aligned$lambda
  nu <- aligned$nu
  dimnames(lambda) <- dimnames(nu) <- list(group_label, items)

  d_lambda <- lambda[machinery$i1, , drop = FALSE] -
    lambda[machinery$i2, , drop = FALSE]
  d_nu <- nu[machinery$i1, , drop = FALSE] - nu[machinery$i2, , drop = FALSE]
  loss_lambda <- colSums(machinery$w * (d_lambda^2 + epsilon)^0.25)
  loss_nu <- colSums(machinery$w * (d_nu^2 + epsilon)^0.25)
  item_loss <- cbind(loadings = loss_lambda, intercepts = loss_nu,
                     total = loss_lambda + loss_nu)
  rownames(item_loss) <- items

  lambda_bar <- colMeans(lambda)
  nu_bar <- colMeans(nu)
  ss_res_lambda <- colSums((lambda_0 - outer(pars$sd, lambda_bar))^2)
  ss_tot_lambda <- colSums(lambda_0^2)
  ss_res_nu <- colSums((nu_0 - outer(rep(1, G), nu_bar) -
                          outer(pars$alpha, lambda_bar))^2)
  ss_tot_nu <- colSums(nu_0^2)
  R2_lambda <- ifelse(ss_tot_lambda > 0, 1 - ss_res_lambda / ss_tot_lambda,
                      NA_real_)
  R2_nu <- ifelse(ss_tot_nu > 0, 1 - ss_res_nu / ss_tot_nu, NA_real_)
  names(R2_lambda) <- names(R2_nu) <- items
  R2_total <- c(
    loadings = if (sum(ss_tot_lambda) > 0)
      1 - sum(ss_res_lambda) / sum(ss_tot_lambda) else NA_real_,
    intercepts = if (sum(ss_tot_nu) > 0)
      1 - sum(ss_res_nu) / sum(ss_tot_nu) else NA_real_
  )

  out <- data.frame(
    group           = group_label,
    n               = n,
    factor_mean     = pars$alpha,
    factor_variance = pars$sd^2,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  res <- .as_dmar_tbl(out, subclass = "dmar_measurement_alignment")
  attr(res, "aligned_loadings")      <- lambda
  attr(res, "aligned_intercepts")    <- nu
  attr(res, "configural_loadings")   <- lambda_0
  attr(res, "configural_intercepts") <- nu_0
  attr(res, "simplicity_function")   <- best$value
  attr(res, "simplicity_starts")     <- values
  attr(res, "converged")             <- converged
  attr(res, "n_starts")              <- n_starts
  attr(res, "n_optima")              <- n_optima
  attr(res, "alignment")             <- alignment
  attr(res, "epsilon")               <- epsilon
  attr(res, "R2_loadings")           <- R2_lambda
  attr(res, "R2_intercepts")         <- R2_nu
  attr(res, "R2_total")              <- R2_total
  attr(res, "item_loss")             <- item_loss
  attr(res, "fit")                   <- fit
  res
}

# Everything the alignment optimizer needs, closed over the configural
# parameter matrices so the objective and its gradient do no repeated setup.
# `lambda_0` and `nu_0` are G x I, `n` the per group sample sizes, `epsilon`
# the smoothing constant, and `alignment` the identification rule. The free
# parameters are the G - 1 factor means and the G - 1 log factor standard
# deviations that the constraint leaves free; `build()` expands them to the
# full length G vectors. Not exported; consumed only by
# measurement_alignment().
.alignment_machinery <- function(lambda_0, nu_0, n, epsilon, alignment) {
  G <- nrow(lambda_0)
  pairs <- utils::combn(G, 2L)
  i1 <- pairs[1L, ]
  i2 <- pairs[2L, ]
  n <- as.numeric(n)
  w <- sqrt(n[i1] * n[i2])
  # Signed incidence matrix, +1 for the first member of a group pair and -1
  # for the second. One matrix product with it turns per pair derivatives
  # into per group derivatives.
  D <- matrix(0, nrow = G, ncol = length(i1))
  D[cbind(i1, seq_along(i1))] <-  1
  D[cbind(i2, seq_along(i2))] <- -1

  build <- function(par) {
    k <- G - 1L
    a <- par[seq_len(k)]
    l <- par[k + seq_len(k)]
    if (identical(alignment, "fixed")) {
      list(alpha = c(0, a), sd = exp(c(0, l)))
    } else {
      list(alpha = c(-sum(a), a), sd = exp(c(-sum(l), l)))
    }
  }
  aligned <- function(pars) {
    lambda <- lambda_0 / pars$sd
    list(lambda = lambda, nu = nu_0 - pars$alpha * lambda)
  }
  fn <- function(par) {
    z <- aligned(build(par))
    dl <- z$lambda[i1, , drop = FALSE] - z$lambda[i2, , drop = FALSE]
    dn <- z$nu[i1, , drop = FALSE] - z$nu[i2, , drop = FALSE]
    sum(w * (rowSums((dl^2 + epsilon)^0.25) + rowSums((dn^2 + epsilon)^0.25)))
  }
  # Analytic gradient. With f(x) = (x^2 + epsilon)^(1/4) the per pair
  # derivatives are A = w f'(d_lambda) and B = w f'(d_nu); collapsing them
  # with D and contracting against the aligned loadings gives
  # dF/d alpha_g = -v_g and dF/d log sd_g = -u_g + alpha_g v_g, where
  # u_g = sum_i (D A)_gi lambda_gi and v_g = sum_i (D B)_gi lambda_gi.
  gr <- function(par) {
    pars <- build(par)
    z <- aligned(pars)
    dl <- z$lambda[i1, , drop = FALSE] - z$lambda[i2, , drop = FALSE]
    dn <- z$nu[i1, , drop = FALSE] - z$nu[i2, , drop = FALSE]
    A <- 0.5 * w * dl * (dl^2 + epsilon)^(-0.75)
    B <- 0.5 * w * dn * (dn^2 + epsilon)^(-0.75)
    u <- rowSums((D %*% A) * z$lambda)
    v <- rowSums((D %*% B) * z$lambda)
    g_alpha <- -v
    g_log_sd <- -u + pars$alpha * v
    if (identical(alignment, "fixed")) {
      c(g_alpha[-1L], g_log_sd[-1L])
    } else {
      c(g_alpha[-1L] - g_alpha[1L], g_log_sd[-1L] - g_log_sd[1L])
    }
  }
  list(build = build, aligned = aligned, fn = fn, gr = gr,
       w = w, i1 = i1, i2 = i2)
}

# Starting values for the alignment optimizer, in the free parameter
# coordinates build() expects. The first start is the null start (every
# factor mean 0, every factor variance 1). The second exploits
# lambda^0_gi = sd_g lambda_i and nu^0_gi = nu_i + alpha_g lambda_i under
# exact invariance: the median loading in a group is proportional to that
# group's factor standard deviation, and the median intercept moves with its
# factor mean. Later starts are random, since the surface has local minima.
# Not exported; consumed only by measurement_alignment().
.alignment_starts <- function(lambda_0, nu_0, n_starts, alignment) {
  G <- nrow(lambda_0)
  k <- G - 1L
  pack <- function(alpha, log_sd) {
    if (!identical(alignment, "fixed")) {
      alpha <- alpha - mean(alpha)
      log_sd <- log_sd - mean(log_sd)
    }
    c(alpha[-1L], log_sd[-1L])
  }
  out <- vector("list", n_starts)
  out[[1L]] <- pack(rep(0, G), rep(0, G))
  if (n_starts >= 2L) {
    sd_init <- apply(lambda_0, 1L, stats::median)
    sd_init <- sd_init / sd_init[1L]
    alpha_init <- apply(nu_0, 1L, stats::median)
    alpha_init <- alpha_init - alpha_init[1L]
    out[[2L]] <- if (all(is.finite(sd_init)) && all(sd_init > 0) &&
                     all(is.finite(alpha_init))) {
      pack(alpha_init, log(sd_init))
    } else {
      out[[1L]]
    }
  }
  for (b in seq_len(n_starts)[-seq_len(min(2L, n_starts))]) {
    out[[b]] <- pack(c(0, stats::runif(k, -1, 1)),
                     c(0, stats::runif(k, -0.5, 0.5)))
  }
  out
}

# Is a candidate solution a solution at all? The criterion has a degenerate
# branch: sending the factor variances of every group but the reference off
# to infinity drives those groups' aligned loadings to zero, which flattens
# the loading half of the criterion, so a start heading that way can stop
# with the optimizer reporting convergence and factor variances of 1e300 or
# Inf. Such a run is discarded rather than returned. The bounds are far
# outside anything a measurement model produces: a factor standard deviation
# ten thousand times the reference group's, or a factor mean ten thousand
# standard deviations away from it, is the degenerate limit, not an estimate.
# Not exported; consumed only by measurement_alignment().
.alignment_usable <- function(pars, sd_bound = 1e4, mean_bound = 1e4) {
  all(is.finite(pars$alpha)) && all(is.finite(pars$sd)) &&
    all(pars$sd > 1 / sd_bound) && all(pars$sd < sd_bound) &&
    all(abs(pars$alpha) < mean_bound)
}

# Count the distinct local minima a set of achieved criterion values
# represents, treating values within a relative tolerance of each other as
# the same minimum. Not exported; consumed only by measurement_alignment(),
# where it fills the "n_optima" attribute.
.alignment_n_optima <- function(values, tol = 1e-6) {
  v <- sort(values[is.finite(values)])
  if (!length(v)) return(0L)
  k <- 1L
  reference <- v[1L]
  for (x in v[-1L]) {
    if (x - reference > tol * max(1, abs(reference))) {
      k <- k + 1L
      reference <- x
    }
  }
  k
}
