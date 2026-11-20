#' Cohen's \emph{f} Effect Size
#'
#' Computes Cohen's \emph{f}, the standard deviation of the group means
#' divided by the standard deviation within the groups, from one of three
#' complete specifications: the group means with the within-group variance
#' (and the group sizes, when the groups differ in size), eta squared, or
#' the two standard deviations.
#'
#' Cohen's \emph{f} is defined with population quantities (Cohen, 1988,
#' p. 276): \eqn{f = \sigma_m / \sigma}, where \eqn{\sigma_m} is the
#' standard deviation of the population means and \eqn{\sigma} is the
#' standard deviation within the populations. Given population values, the
#' function returns the population \emph{f}. Given sample summaries, it
#' returns the plug-in value from the same calculation. The plug-in value
#' depends on which summaries are given: from the same data, a sample eta
#' squared gives a larger value than the observed group means with the
#' pooled within-group mean square (see below). Sample inputs do not make
#' the result an unbiased estimate of the population \emph{f}: the spread of
#' sample means includes their sampling error as well as the spread of the
#' population means. Unequal group sizes determine how the group means are
#' weighted; without \code{n_per_group} the means are weighted equally.
#'
#' @param means Numeric vector of group means, one per group and at least
#'   two: population means for the population \emph{f}, or observed sample
#'   means for the plug-in value. Give with \code{var_within}.
#' @param var_within The within-group variance, a single positive number:
#'   the population variance common to the groups, or, for the plug-in
#'   value, the pooled within-group mean square (the mean square error of
#'   the one-way analysis of variance). Give with \code{means}.
#' @param n_per_group Optional, and used only with \code{means}: the group
#'   sizes, one per group, which weight each group mean by its share of the
#'   total, \eqn{n_j / N}. A single number states equal groups and gives the
#'   same result as \code{NULL}, the default, which weights the means
#'   equally (\eqn{1 / k} each, for \emph{k} groups).
#' @param eta_squared Eta squared, the proportion of the variance that the
#'   groups account for, a single number in [0, 1). Give it alone. A
#'   population value gives the population \emph{f}. A sample eta squared
#'   gives the plug-in value
#'   \eqn{\sqrt{df_1 F / df_2}}{sqrt(df_1 F / df_2)}, which is larger
#'   than the value from the observed means and the pooled mean square
#'   (see Details).
#' @param sd_means The standard deviation of the group means,
#'   \eqn{\sigma_m}, a single number of at least 0, computed with divisor
#'   \emph{k} (or with the weights \eqn{n_j / N}), not with the divisor
#'   \emph{k} - 1 that \code{stats::sd()} uses. Give with \code{sd_within}.
#' @param sd_within The within-group standard deviation, \eqn{\sigma}, a
#'   single positive number, or, for the plug-in value, the square root of
#'   the pooled within-group mean square. Give with \code{sd_means}.
#'
#' @details
#' A call gives exactly one complete specification: \code{means} with
#' \code{var_within} (and \code{n_per_group} if wanted), \code{eta_squared}
#' alone, or \code{sd_means} with \code{sd_within}. Each determines
#' \emph{f} on its own, so an argument of one specification given with an
#' argument of another is refused, as is one half of a pair given alone, or
#' \code{n_per_group} without \code{means}; the refusal has class
#' \code{dmar_argument_invalid}.
#'
#' From the means, with weights \eqn{w_j = 1 / k} by default or
#' \eqn{w_j = n_j / N} when \code{n_per_group} is given,
#' \deqn{\bar\mu = \sum_j w_j \mu_j, \qquad
#'   \sigma_m^2 = \sum_j w_j (\mu_j - \bar\mu)^2, \qquad
#'   f = \sqrt{\sigma_m^2 / \sigma^2}.}{mu_bar = sum_j w_j mu_j,
#'   sigma_m^2 = sum_j w_j (mu_j - mu_bar)^2, f = sqrt(sigma_m^2 / sigma^2).}
#' The divisor of \eqn{\sigma_m^2} is \emph{k} (or \emph{N}), because the
#' variance is that of a fixed set of means. From eta squared,
#' \eqn{f = \sqrt{\eta^2 / (1 - \eta^2)}} (Kim, 2016, p. 233), and from the
#' two standard deviations, \eqn{f = \sigma_m / \sigma}. Given population
#' values for the same populations, with the same weighting of the means,
#' the three give the same \emph{f}.
#'
#' Given sample summaries, the three agree only when the within-group
#' variance is defined the same way in each. With the observed means, the
#' group sizes, and the pooled within-group mean square,
#' \eqn{SS_{within} / (N - k)}{SS_within / (N - k)}, of a one-way design,
#' \eqn{N \sigma_m^2} is the between-groups sum of squares, so the plug-in
#' value equals \eqn{\sqrt{df_1 F / N}}{sqrt(df_1 F / N)}, with \emph{F},
#' \eqn{df_1 = k - 1}, and \eqn{df_2 = N - k} from the analysis of
#' variance. A sample eta squared,
#' \eqn{SS_{between} / SS_{total}}{SS_between / SS_total}, pairs with the
#' within-group variance \eqn{SS_{within} / N}{SS_within / N} rather than
#' with the mean square, so the value from it,
#' \eqn{\sqrt{SS_{between} / SS_{within}}}{sqrt(SS_between / SS_within)}
#' \eqn{= \sqrt{df_1 F / df_2}}{= sqrt(df_1 F / df_2)}, is larger than the
#' value from the means and the pooled mean square by the factor
#' \eqn{\sqrt{N / (N - k)}}{sqrt(N / (N - k))}. The two are equal when
#' \code{var_within} is \eqn{SS_{within} / N}{SS_within / N} (or
#' \code{sd_within} is its square root). Without \code{n_per_group}, the
#' means of unequal groups
#' are weighted equally and the value differs from both. The square root of
#' the signal-to-noise ratio of \code{\link{ci_srsnr}} is this \emph{f}
#' with the weights \eqn{n_j / N} and the pooled mean square; that function
#' gives its confidence interval from \emph{F} or from the same summaries.
#'
#' @return A one-row \code{data.frame} (printed through \code{dmar_tbl})
#'   with columns \code{term} and \code{value}; \code{term} is
#'   \code{"cohen_f"} and \code{value} is \emph{f} at full precision.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
#'   size calculation 3. Comparison of several means using one-way ANOVA.
#'   \emph{Restorative Dentistry & Endodontics, 41}(3), 231--234.
#'   \doi{10.5395/rde.2016.41.3.231}
#'
#' @seealso \code{\link{ci_srsnr}}, \code{\link{ci_snr}},
#'   \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_R2}}
#'
#' @examples
#' # Kim (2016, Table 1) plans a four-group study around means of 18.4,
#' # 22.2, 25.1, and 32.1 with a common standard deviation of 8.6, and
#' # reports f as .58. At full precision, this f plans a balanced study of
#' # 40, ten per group, for power .80 at the .05 level, as the article does.
#' f_kim <- cohen_f(means = c(18.4, 22.2, 25.1, 32.1), var_within = 8.6^2)
#' f_kim
#' ss_power_one_way_anova(a = 4, true_f = f_kim$value, desired_power = .80)
#'
#' # Cohen (1988, p. 277) places seven population means over half a
#' # within-population standard deviation, one at each end and five at the
#' # midpoint, and gives f as .13 (.267 times .5, his Equation 8.2.8).
#' cohen_f(means = c(-.25, 0, 0, 0, 0, 0, .25), var_within = 1)
#'
#' # The plug-in value from the summaries of an observed one-way design
#' # with unequal groups: the group means, the pooled within-group mean
#' # square, and the group sizes of the three arms of drinks_trial.
#' fit <- stats::lm(log_drinks ~ treatment, data = drinks_trial)
#' cohen_f(means = tapply(drinks_trial$log_drinks, drinks_trial$treatment, mean),
#'         var_within = sum(stats::residuals(fit)^2) / stats::df.residual(fit),
#'         n_per_group = as.vector(table(drinks_trial$treatment)))
#'
#' # The sample eta squared of the same fit gives a larger plug-in value,
#' # by the factor sqrt(N / (N - k)) = sqrt(88 / 85) (see Details).
#' ss <- stats::anova(fit)[["Sum Sq"]]
#' cohen_f(eta_squared = ss[1] / sum(ss))
#'
#' # From eta squared, and from the two standard deviations:
#' cohen_f(eta_squared = .10)
#' cohen_f(sd_means = 4, sd_within = 8)
#'
#' @export
cohen_f <- function(means = NULL, var_within = NULL, n_per_group = NULL,
                    eta_squared = NULL, sd_means = NULL, sd_within = NULL) {
  .dmar_check_numbers("cohen_f", c("means", "var_within", "n_per_group",
                                   "eta_squared", "sd_means", "sd_within"))
  given <- c(means = !is.null(means), var_within = !is.null(var_within),
             n_per_group = !is.null(n_per_group),
             eta_squared = !is.null(eta_squared),
             sd_means = !is.null(sd_means), sd_within = !is.null(sd_within))
  used <- c(means = any(given[c("means", "var_within", "n_per_group")]),
            eta_squared = given[["eta_squared"]],
            sd = any(given[c("sd_means", "sd_within")]))
  specifications <- c(
    "'means' with 'var_within' (and 'n_per_group' when the groups differ ",
    "in size), 'eta_squared' alone, or 'sd_means' with 'sd_within'")
  quoted <- function(x) paste0("'", x, "'", collapse = ", ")

  if (!any(used)) {
    .dmar_refuse_argument(c(
      "cohen_f() needs one complete specification of f, and none was ",
      "given. Give ", specifications, "."), arg = names(given))
  }
  if (sum(used) > 1L) {
    named <- names(given)[given]
    .dmar_refuse_argument(c(
      quoted(named), " were given together, and they belong to different ",
      "specifications of f. Each specification determines f on its own, ",
      "so values from two of them could disagree, and cohen_f() does not ",
      "choose between them. Give exactly one: ", specifications, "."),
      arg = named)
  }
  pair_missing <- function(have, need, role) {
    .dmar_refuse_argument(c(
      "'", have, "' was given without '", need, "'. ", role, " Give both, ",
      "or give one of the other specifications: ", specifications, "."),
      arg = c(have, need))
  }
  if (used[["means"]]) {
    if (!given[["means"]] && !given[["var_within"]]) {
      .dmar_refuse_argument(c(
        "'n_per_group' was given without 'means' and 'var_within'. It ",
        "weights the group means and is used only with them; give the ",
        "means and the within-group variance as well, or give ",
        "'eta_squared' alone, or 'sd_means' with 'sd_within', without ",
        "'n_per_group'."), arg = c("n_per_group", "means", "var_within"))
    }
    role <- paste0("f is the standard deviation of the group means divided ",
                   "by the within-group standard deviation, so the means ",
                   "and the within-group variance are both needed.")
    if (!given[["var_within"]]) pair_missing("means", "var_within", role)
    if (!given[["means"]]) pair_missing("var_within", "means", role)
  }
  if (used[["sd"]]) {
    role <- paste0("f is the standard deviation of the group means ",
                   "('sd_means') divided by the within-group standard ",
                   "deviation ('sd_within'), so both are needed.")
    if (!given[["sd_within"]]) pair_missing("sd_means", "sd_within", role)
    if (!given[["sd_means"]]) pair_missing("sd_within", "sd_means", role)
  }

  one_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

  if (used[["means"]]) {
    if (!is.numeric(means) || length(means) < 2L || !all(is.finite(means))) {
      .dmar_refuse_argument(c(
        "'means' must be a numeric vector of at least two finite group ",
        "means; ", .dmar_given(means), "."), arg = "means", value = means)
    }
    means <- as.vector(means)
    k <- length(means)
    if (!one_number(var_within) || var_within <= 0) {
      .dmar_refuse_argument(c(
        "'var_within' must be a single positive number, the within-group ",
        "variance; ", .dmar_given(var_within), "."),
        arg = "var_within", value = var_within)
    }
    if (is.null(n_per_group)) {
      # Equal weighting of the group means: divisor k.
      var_means <- mean((means - mean(means))^2)
    } else {
      if (!is.numeric(n_per_group) ||
          !(length(n_per_group) %in% c(1L, k))) {
        .dmar_refuse_argument(c(
          "'n_per_group' must be either a single number or a vector of the ",
          "same length as 'means' (", k, " group sizes); ",
          .dmar_given(n_per_group), "."),
          arg = c("n_per_group", "means"), value = n_per_group)
      }
      if (!all(is.finite(n_per_group)) || any(n_per_group <= 0)) {
        .dmar_refuse_argument(c(
          "All 'n_per_group' values must be positive finite numbers; they ",
          "are group sizes, which weight the group means."),
          arg = "n_per_group", value = n_per_group)
      }
      n_j <- rep_len(as.vector(n_per_group), k)
      N <- sum(n_j)
      grand_mean <- sum(n_j * means) / N
      var_means <- sum(n_j * (means - grand_mean)^2) / N
    }
    f <- sqrt(var_means / var_within)
  } else if (used[["eta_squared"]]) {
    if (!one_number(eta_squared) || eta_squared < 0 || eta_squared >= 1) {
      .dmar_refuse_argument(c(
        "'eta_squared' must be a single number in [0, 1); ",
        .dmar_given(eta_squared), "."),
        arg = "eta_squared", value = eta_squared)
    }
    f <- sqrt(eta_squared / (1 - eta_squared))
  } else {
    if (!one_number(sd_within) || sd_within <= 0) {
      .dmar_refuse_argument(c(
        "'sd_within' must be a single positive number, the within-group ",
        "standard deviation; ", .dmar_given(sd_within), "."),
        arg = "sd_within", value = sd_within)
    }
    if (!one_number(sd_means) || sd_means < 0) {
      .dmar_refuse_argument(c(
        "'sd_means' must be a single number of at least 0, the standard ",
        "deviation of the group means; ", .dmar_given(sd_means), "."),
        arg = "sd_means", value = sd_means)
    }
    f <- sd_means / sd_within
  }

  .as_dmar_tbl(data.frame(term = "cohen_f", value = f))
}
