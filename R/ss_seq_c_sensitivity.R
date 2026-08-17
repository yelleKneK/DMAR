# Monte Carlo evaluation of the sequential fixed-width contrast procedure.
#' Monte Carlo Sensitivity of the Sequential Fixed-Width Procedure
#'
#' Simulates the purely sequential fixed-width confidence interval
#' procedure of \code{\link{ss_seq_c}} under a known data generating
#' mechanism, reporting the distribution of the stopping sample size
#' and the empirical coverage of the fixed-width interval. The two
#' quantities to read are the ratio of the mean stopping size to the
#' oracle \eqn{n^*} (first-order efficiency: the ratio approaches 1 as
#' the target half-width shrinks) and the coverage (asymptotic
#' consistency: coverage approaches 1 - 2\eqn{\alpha}). The normal
#' quantile rule stops slightly early at wide targets, the
#' finite-sample undershoot anticipated by Woodroofe (1977); the
#' \emph{t} quantile rule corrects it at a small cost in sample size.
#'
#' @param c_weights The contrast weights. The weights must sum to zero
#'   with the positive weights summing to 1 and the negative weights
#'   to -1.
#' @param half_width The target half-width \eqn{h} of the
#'   100(1 - 2\eqn{\alpha})\% interval, in raw units of the response.
#' @param true_sigma The data generating error standard deviation:
#'   a single value applied to every group.
#' @param true_means Optional vector of data generating group means,
#'   aligned with \code{c_weights}. Default all zero, so the true
#'   contrast is 0.
#' @param alpha_level One-sided rate per bound; the interval is at
#'   confidence level 1 - 2\eqn{\alpha}. Default \code{0.05}.
#' @param quantile \code{"t"} (default) or \code{"normal"}; see
#'   \code{\link{ss_seq_c}}.
#' @param m0 Pilot sample size per group. Default \code{10}.
#' @param G Number of Monte Carlo replications. Default \code{1000}.
#' @param seed Optional integer seed. Default \code{NULL} (the
#'   current RNG state is used). When supplied, the caller's RNG
#'   state is restored on exit.
#'
#' @return A \code{data.frame} with rows \code{n_star} (the oracle
#'   total sample size an investigator with known \eqn{\sigma} would
#'   use), \code{mean_N}, \code{median_N}, \code{sd_N} (the stopping
#'   total across replications), \code{ratio_mean_N_n_star},
#'   \code{coverage} (the proportion of replications whose
#'   \eqn{\hat\psi_N \pm h} interval covered the true contrast),
#'   \code{se_coverage} (its simulation standard error), and the input
#'   echoes \code{half_width}, \code{true_psi} (the population contrast
#'   implied by \code{c_weights} and \code{true_means}),
#'   \code{true_sigma}, \code{alpha_level}, and \code{m0}.
#'
#' @details
#' \strong{Simulation design.} Each replication samples the groups
#' with nonzero weights in balanced fashion (one observation per
#' group per step) from normal populations with common
#' \code{true_sigma}, starting at \code{m0} per group, and stops at
#' the first step satisfying the \code{\link{ss_seq_c}} criterion with
#' the pooled variance estimate. This matches the equal-cost,
#' equal-variance case of Chattopadhyay, Bandyopadhyay, Kelley, and
#' Padalunkal (2025); unequal costs change the optimal allocation but
#' not the logic.
#'
#' \strong{The oracle.} With known \eqn{\sigma} and balanced
#' allocation over the \eqn{J_0} groups with nonzero weights, the
#' fixed-width requirement is
#' \eqn{n^* = z_{1-\alpha}^2\, \sigma^2 J_0 \sum_j c_j^2 / h^2} in
#' total. The sequential procedure spends about \eqn{n^*} without
#' knowing \eqn{\sigma}, which is its point.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of
#'   fixed-width sequential confidence intervals for the mean.
#'   \emph{The Annals of Mathematical Statistics, 36}(2), 457--462.
#'
#' Ghosh, M., Mukhopadhyay, N., & Sen, P. K. (1997). \emph{Sequential
#'   estimation}. Wiley.
#'
#' Woodroofe, M. (1977). Second order approximations for sequential
#'   point and interval estimation. \emph{The Annals of Statistics,
#'   5}(5), 984--995.
#'
#' @seealso \code{\link{ss_seq_c}}, \code{\link{ss_aipe_c}},
#'   \code{\link{ss_power_equivalence_c}}
#'
#' @examples
#' # A two-group contrast, target half-width 2.5, error SD 15.67:
#' # the t-quantile rule stops near the oracle with near-nominal
#' # coverage. (G kept small here for speed; use G = 2000 or more in
#' # earnest.)
#' ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
#'                      true_sigma = 15.67, G = 200, seed = 113)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family sequential estimation
#'
#' @export

ss_seq_c_sensitivity <- function(c_weights, half_width, true_sigma,
                                 true_means = NULL,
                                 alpha_level = 0.05,
                                 quantile = c("t", "normal"),
                                 m0 = 10, G = 1000, seed = NULL) {
  quantile <- match.arg(quantile)

  if (!is.numeric(half_width) || length(half_width) != 1L || half_width <= 0)
    stop("'half_width' must be a single positive number.")
  if (!is.numeric(true_sigma) || length(true_sigma) != 1L || true_sigma <= 0)
    stop("'true_sigma' must be a single positive number.")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be a single number in (0, 0.5).")
  if (!is.numeric(G) || length(G) != 1L || G < 2)
    stop("'G' must be a single number of at least 2.")
  if (!is.numeric(m0) || length(m0) != 1L || m0 < 2)
    stop("'m0' must be a single number of at least 2.")
  if (!identical(round(sum(c_weights), 5), 0))
    stop("The sum of the contrast weights ('c_weights') should equal zero.")
  pos <- sum(c_weights[c_weights > 0])
  neg <- sum(c_weights[c_weights < 0])
  if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
    stop("The positive weights must sum to 1 and the negative weights to ",
         "-1, so that 'half_width' is on the raw scale of the response.")

  J <- length(c_weights)
  if (is.null(true_means)) true_means <- rep(0, J)
  if (length(true_means) != J)
    stop("'true_means' must give one mean per group, aligned with ",
         "'c_weights'.")

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = globalenv()))
                rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  involved <- which(c_weights != 0)
  J0       <- length(involved)
  cw       <- c_weights[involved]
  mu       <- true_means[involved]
  true_psi <- sum(c_weights * true_means)

  z <- stats::qnorm(1 - alpha_level)
  # Oracle total under balanced allocation over the involved groups:
  # per-group requirement k* solves z^2 * sum(c_j^2) * sigma^2 / k = h^2.
  k_star <- z^2 * true_sigma^2 * sum(cw^2) / half_width^2
  n_star <- J0 * k_star

  m0   <- ceiling(m0)
  kmax <- max(m0 + 30L, ceiling(4 * k_star))

  one_rep <- function() {
    y <- matrix(stats::rnorm(kmax * J0, mean = rep(mu, each = kmax),
                             sd = true_sigma),
                nrow = kmax, ncol = J0)
    kk <- seq_len(kmax)
    # Running per-group variances via cumulative sums, then the pooled
    # estimate, without refitting at every step.
    v <- sapply(seq_len(J0), function(j) {
      (cumsum(y[, j]^2) - cumsum(y[, j])^2 / kk) / (kk - 1)
    })
    sp2 <- rowMeans(v)
    q <- if (quantile == "t") stats::qt(1 - alpha_level, df = pmax(J0 * (kk - 1), 1))
         else rep(z, kmax)
    # Criterion at per-group size k (balanced): q^2 * sum(c^2) * sp2 / k <= h^2.
    need_k <- q^2 * sum(cw^2) * sp2 / half_width^2
    ok <- which(kk >= m0 & kk >= need_k)
    K  <- if (length(ok)) ok[1] else kmax
    psi_hat <- sum(cw * colMeans(y[seq_len(K), , drop = FALSE]))
    c(N = J0 * K, cover = as.numeric(abs(psi_hat - true_psi) <= half_width))
  }

  sims <- t(replicate(G, one_rep()))
  cover <- mean(sims[, "cover"])

  out <- data.frame(
    term  = c("n_star", "mean_N", "median_N", "sd_N",
              "ratio_mean_N_n_star", "coverage", "se_coverage",
              "half_width", "true_psi", "true_sigma",
              "alpha_level", "m0"),
    value = c(n_star, mean(sims[, "N"]), stats::median(sims[, "N"]),
              stats::sd(sims[, "N"]), mean(sims[, "N"]) / n_star,
              cover, sqrt(cover * (1 - cover) / G),
              half_width, true_psi, true_sigma,
              alpha_level, m0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level)
}
