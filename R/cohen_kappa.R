# Internal: build a k x k agreement-weight matrix from a string spec or a
# user-supplied matrix. Linear and quadratic weights treat categories as
# ordinal with equal spacing; a custom matrix accommodates arbitrary
# structures, in either of two scalings:
#
#   scaling = "agreement"    : diagonal 1, decreasing off-diagonal entries
#                              (typically in [0, 1]); used directly.
#   scaling = "disagreement" : Cohen's (1968) ratio-scaled disagreement
#                              weights v_ij (zero diagonal, larger =
#                              graver disagreement); converted to
#                              agreement weights via w = 1 - v / max(v),
#                              which leaves kappa_w unchanged (Cohen,
#                              1968, Footnote 3) and puts the matrix on
#                              the scale the Fleiss-Cohen-Everitt
#                              variance expects.
.build_kappa_weights <- function(weights, k, scaling = "agreement") {
  if (scaling == "disagreement") {
    if (!is.matrix(weights)) {
      stop("weight_scaling = \"disagreement\" requires a custom 'weights' ",
           "matrix of ratio-scaled disagreement weights; the string specs ",
           "(\"unweighted\", \"linear\", \"quadratic\") are defined on the ",
           "agreement scale.", call. = FALSE)
    }
    if (!all(dim(weights) == c(k, k))) {
      stop("Custom 'weights' matrix must be ", k, " x ", k, ".",
           call. = FALSE)
    }
    if (!is.numeric(weights) || anyNA(weights) || any(weights < 0)) {
      stop("Disagreement weights must be non-negative numbers.",
           call. = FALSE)
    }
    if (any(diag(weights) != 0)) {
      stop("Disagreement weights must be 0 on the diagonal (no ",
           "disagreement for exact agreement). Cohen (1968) notes the ",
           "zero diagonal is the convention; kappa_w is invariant to ",
           "multiplying the weights by any positive constant, so only ",
           "the ratios among the off-diagonal entries matter.",
           call. = FALSE)
    }
    v_max <- max(weights)
    if (v_max <= 0) {
      stop("Disagreement weights must include at least one positive ",
           "entry.", call. = FALSE)
    }
    return(1 - weights / v_max)
  }
  if (is.matrix(weights)) {
    if (!all(dim(weights) == c(k, k))) {
      stop("Custom 'weights' matrix must be ", k, " x ", k, ".",
           call. = FALSE)
    }
    if (any(weights > 1) || any(weights < 0)) {
      warning("Custom 'weights' matrix entries are typically in [0, 1] ",
              "(1 on the diagonal, decreasing off-diagonal). If you have ",
              "Cohen-style disagreement weights, pass weight_scaling = ",
              "\"disagreement\".", call. = FALSE)
    }
    return(weights)
  }
  spec <- match.arg(weights, c("unweighted", "linear", "quadratic"))
  if (spec == "unweighted") return(diag(k))
  i_mat <- matrix(seq_len(k), nrow = k, ncol = k)
  j_mat <- t(i_mat)
  d <- abs(i_mat - j_mat)
  if (spec == "linear")    return(1 - d / (k - 1))
  if (spec == "quadratic") return(1 - (d / (k - 1))^2)
  stop("Unknown weights specification.", call. = FALSE)
}

.weights_label <- function(weights, scaling = "agreement") {
  if (is.matrix(weights)) {
    if (scaling == "disagreement") "custom_disagreement" else "custom"
  } else {
    weights[1L]
  }
}


#' Cohen's Kappa Coefficient of Inter-Rater Agreement
#'
#' Computes Cohen's (1960) kappa coefficient of agreement between two raters
#' on a categorical variable, with optional weighting (linear, quadratic, or
#' custom) for ordinal categories, following Cohen's (1968) weighted kappa.
#' Input can be the two raters' classification vectors or a published
#' \eqn{k \times k} frequency table. A confidence interval and a Wald test
#' of \eqn{H_0\!: \kappa = 0} are returned, using the asymptotic standard
#' error of Fleiss, Cohen, and Everitt (1969).
#'
#' @param rater_1 First rater's classification vector (categorical, factor or coercible to factor).
#' @param rater_2 Second rater's classification vector (same length as \code{rater_1}; categorical, factor or coercible to factor). Together \code{rater_1} and \code{rater_2} give the two raters' classifications on the same set of subjects; pairs in which either is \code{NA} are dropped before computation.
#' @param table A \eqn{k \times k} frequency table of joint assignments
#'   (rows are rater 1's categories, columns rater 2's, in the same
#'   order), supplied instead of \code{rater_1} and \code{rater_2}. This
#'   is the form in which published agreement studies usually report
#'   their data. Row and column \code{dimnames}, when present, must
#'   match and are used as the category labels.
#' @param weights Either the string \code{"unweighted"} (the default;
#'   appropriate for nominal categories), \code{"linear"}, or
#'   \code{"quadratic"}, or a user-supplied square numeric weight matrix.
#'   A custom matrix is interpreted on the scale named by
#'   \code{weight_scaling}. Asymmetric matrices are allowed; see
#'   \emph{Details}.
#' @param weight_scaling How a custom \code{weights} matrix is scaled:
#'   \code{"agreement"} (the default; diagonal 1, decreasing off-diagonal
#'   entries, typically in \eqn{[0, 1]}) or \code{"disagreement"}
#'   (Cohen's 1968 ratio-scaled disagreement weights: zero diagonal,
#'   larger entries for graver disagreements). The two scalings give
#'   identical \eqn{\kappa_W}; see \emph{Details}.
#' @param categories Optional character vector listing the full category
#'   set in display order (used both to set the order of the confusion
#'   matrix's rows/columns and to map ordinal categories to integers
#'   \eqn{1, \ldots, k} for linear and quadratic weighting). When
#'   \code{NULL} (the default) and both raters are factors with the same
#'   level set, the factor levels are used in their own order;
#'   otherwise the sorted union of values observed in
#'   \code{rater_1} and \code{rater_2} is used. When supplied, the set
#'   must be unique, non-missing, and contain every observed rating;
#'   an omitted category raises an error, because silently dropping the
#'   corresponding ratings while still dividing by the original sample
#'   size would corrupt kappa. With \code{table} input, supply
#'   \code{categories} only when the table has no \code{dimnames}.
#' @param conf_level Confidence level for the interval (default
#'   \code{0.95}).
#' @param ci_method Interval method: \code{"wald"} (the default, the
#'   asymptotic interval from the Fleiss, Cohen, and Everitt standard
#'   error), \code{"percentile"} (bootstrap percentile), or
#'   \code{"bca"} (bootstrap bias-corrected and accelerated).
#' @param B Number of bootstrap replications when \code{ci_method} is
#'   \code{"percentile"} or \code{"bca"} (default \code{10000};
#'   ignored for \code{"wald"}). The BCa adjustment pushes the working
#'   quantiles into the tails, so reduce \code{B} for exploration, not
#'   for a reported analysis.
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator;
#'   a supplied seed is set internally and the prior state restored on
#'   exit.
#'
#' @return A one-row \code{data.frame} (class \code{dmar_tbl}) with columns
#'   \code{weights} (the form used), \code{kappa}, \code{se} (asymptotic
#'   standard error),
#'   \code{lower_limit}, \code{upper_limit}, \code{z_value}, \code{p_value}
#'   (Wald test of \eqn{H_0\!: \kappa = 0}), \code{n} (number of paired
#'   ratings), and \code{n_categories} (\eqn{k}).
#'
#'   The per-cell detail behind the coefficient travels with the result
#'   as the \code{cells} attribute, in the form of Cohen's (1968)
#'   Table 1: a \code{data.frame} with one row per cell of the
#'   confusion matrix giving \code{rater_1} and \code{rater_2} (the
#'   cell's categories), \code{observed_proportion},
#'   \code{expected_proportion} (the product of the marginal
#'   proportions, the cell's chance expectation), \code{weight} (the
#'   agreement-scale weight used in the computation), and, when
#'   \code{weight_scaling = "disagreement"}, the supplied
#'   \code{disagreement_weight}. Retrieve it with
#'   \code{attr(result, "cells")}.
#'
#' @details
#' For two raters and a confusion matrix \eqn{P} of joint proportions
#' (rater 1 \eqn{\times} rater 2), the weighted kappa is
#' \deqn{\kappa_W = \frac{p_o^{(W)} - p_e^{(W)}}{1 - p_e^{(W)}}, \quad
#'                 p_o^{(W)} = \sum_{i,j} W_{ij}\,P_{ij}, \quad
#'                 p_e^{(W)} = \sum_{i,j} W_{ij}\,p_{i.}\,p_{.j},}
#' where \eqn{p_{i.}} and \eqn{p_{.j}} are the row and column marginals.
#' For \code{weights = "unweighted"} (the diagonal of \eqn{W} is 1, off-
#' diagonal 0) this collapses to Cohen's original formulation.
#'
#' \strong{Disagreement scaling.} Cohen (1968) develops weighted kappa
#' by ratio scaling \emph{disagreement}: each cell receives a weight
#' \eqn{v_{ij} \ge 0}, zero on the agreement diagonal, with, for
#' example, a weight of 6 representing twice as much disagreement as 3.
#' The weights are part of the definition of agreement (and of any
#' hypothesis tested about it), so they must be fixed before the data
#' are collected. \eqn{\kappa_W} is invariant to multiplying the
#' \eqn{v_{ij}} by any positive constant, and a disagreement matrix is
#' related to an agreement matrix by
#' \eqn{w_{ij} = 1 - v_{ij}/v_{\max}} (Cohen, 1968, Footnote 3), which
#' is the conversion applied internally when
#' \code{weight_scaling = "disagreement"}. Either scaling therefore
#' yields the same \eqn{\kappa_W}; supply whichever is more natural.
#'
#' \strong{Asymmetric weights and validity.} Nothing in \eqn{\kappa_W}
#' requires \eqn{W_{ij} = W_{ji}}. Symmetric weights suit reliability,
#' where the two sources have equal status; asymmetric weights suit
#' validity, where one source is a criterion and the other a predictor
#' and the two directions of a confusion can carry different costs
#' (Cohen, 1968). The examples reproduce Cohen's computer-diagnosis
#' illustration.
#'
#' \strong{Standard error.} The Fleiss-Cohen-Everitt (1969) asymptotic
#' variance for weighted kappa is used:
#' \deqn{\mathrm{Var}(\hat\kappa_W) = \frac{1}{N(1 - p_e^{(W)})^2}\Bigl[\sum_{i,j} P_{ij}\bigl(W_{ij} - (\bar W_{i.} + \bar W_{.j})(1 - \hat\kappa_W)\bigr)^2 - \bigl(\hat\kappa_W - p_e^{(W)}(1 - \hat\kappa_W)\bigr)^2\Bigr],}
#' with \eqn{\bar W_{i.} = \sum_j W_{ij}\,p_{.j}} and
#' \eqn{\bar W_{.j} = \sum_i W_{ij}\,p_{i.}}. The Wald confidence interval is
#' \eqn{\hat\kappa \pm z_{1-\alpha/2}\,\widehat{\mathrm{SE}}}. Cohen's
#' (1968) own Formulas 10 and 13 for the standard error of
#' \eqn{\kappa_W} preceded this result and were superseded by it; the
#' examples reproduce his Table 1 arithmetic for the historical record
#' while the function reports the Fleiss-Cohen-Everitt interval.
#'
#' \strong{Choice of weights.} Use \code{"unweighted"} for nominal
#' categories. For ordinal categories, \code{"quadratic"} is the most
#' common choice (and mathematically equivalent to the intraclass
#' correlation under certain conditions; Fleiss & Cohen, 1973);
#' \code{"linear"} is also defensible. Cohen (1968) further shows that
#' with equal marginals and quadratic-pattern disagreement weights,
#' \eqn{\kappa_W} equals the product-moment correlation between the
#' category scores.
#'
#' \strong{Small samples and the bootstrap.} The Wald interval can have
#' poor coverage for small \eqn{N} or extreme values of
#' \eqn{\hat\kappa}; a bootstrap interval is more dependable in those
#' regimes (Blackman & Koval, 2000). With
#' \code{ci_method = "percentile"} or \code{"bca"} the subjects (the
#' rated pairs) are resampled with replacement \code{B} times, kappa is
#' recomputed on each resample with the same categories and weights, and
#' the interval is read off the bootstrap distribution: the percentile
#' interval takes the empirical quantiles, and the BCa interval adjusts
#' the quantile positions for median bias (estimated from the bootstrap
#' distribution) and for acceleration (estimated from the jackknife),
#' making it second-order accurate where the percentile interval is
#' first-order accurate (Efron & Tibshirani, 1993). \code{table} input
#' is expanded to the equivalent paired ratings and resampled the same
#' way. A resample on which kappa is undefined (chance agreement 1) is
#' dropped, and the interval is computed from the replications that
#' return a value; a single warning reports how many were dropped. The
#' \code{se}, \code{z_value}, and \code{p_value} columns keep their
#' asymptotic definitions under every \code{ci_method}; only the
#' interval changes. Bootstrap results vary from run to run; supply
#' \code{seed} for reproducibility (the RNG state is set locally and
#' the caller's state restored on exit).
#'
#' @references
#' Blackman, N. J.-M., & Koval, J. J. (2000). Interval estimation for
#'   Cohen's kappa as a measure of agreement. \emph{Statistics in
#'   Medicine, 19}(5), 723--741.
#'
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#'   \emph{Educational and Psychological Measurement, 20}(1), 37--46.
#'
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for
#'   scaled disagreement or partial credit. \emph{Psychological Bulletin,
#'   70}(4), 213--220.
#'
#' Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample
#'   standard errors of kappa and weighted kappa. \emph{Psychological
#'   Bulletin, 72}(5), 323--327.
#'
#' Fleiss, J. L., & Cohen, J. (1973). The equivalence of weighted kappa
#'   and the intraclass correlation coefficient as measures of
#'   reliability. \emph{Educational and Psychological Measurement, 33}(3),
#'   613--619.
#'
#' @examples
#' # ---------------------------------------------------------------------
#' # Cohen (1968, Table 1): two judges assign N = 200 cases to three
#' # diagnostic categories. The table ships as the diagnosis_agreement
#' # data set in the paper's own layout, Judge B in rows and Judge A in
#' # columns, carrying Cohen's per-cell disagreement weights, observed
#' # proportions, and chance-expected proportions.
#' data(diagnosis_agreement)
#' tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
#' v   <- unclass(xtabs(disagreement_weight ~ judge_b + judge_a,
#'                      data = diagnosis_agreement))
#' tab
#'
#' # Unweighted kappa, all disagreements equal (Cohen's Formula 4): .492.
#' cohen_kappa(table = tab)
#'
#' # A bootstrap interval for the same table, which expands it to its
#' # 200 paired ratings and resamples the subjects. Not run here,
#' # because 2000 refits of kappa is more than a help page should do;
#' # the call is:
#' # cohen_kappa(table = tab, ci_method = "percentile", B = 2000,
#' #             seed = 113)
#'
#' # Cohen's ratio-scaled disagreement weights: a neurosis-psychosis
#' # confusion (weight 6) is six times as grave as a personality
#' # disorder-neurosis confusion (weight 1). Weighted kappa = .348,
#' # smaller than the unweighted .492: these judges disagree less than
#' # chance expectation where it matters little and at about the chance
#' # level where it matters most.
#' res <- cohen_kappa(table = tab, weights = v,
#'                    weight_scaling = "disagreement")
#' res
#'
#' # The per-cell quantities of Cohen's Table 1 travel with the result:
#' # observed and chance-expected proportions and both weight scalings.
#' attr(res, "cells")
#'
#' # Interchanging the 6 and 1 weights reverses the story: kappa_w = .574.
#' v_swap <- v
#' v_swap[v == 6] <- 1
#' v_swap[v == 1] <- 6
#' cohen_kappa(table = tab, weights = v_swap,
#'             weight_scaling = "disagreement")
#'
#' # Cohen's own Table 1 arithmetic (his Formulas 8, 10, and 13),
#' # computed straight from the data set's per-cell columns and
#' # reproduced for the historical record. The function reports the
#' # Fleiss-Cohen-Everitt (1969) standard error, which superseded
#' # Formulas 10 and 13.
#' q_o <- with(diagnosis_agreement,
#'             sum(disagreement_weight * observed_proportion))   # = .90
#' q_c <- with(diagnosis_agreement,
#'             sum(disagreement_weight * expected_proportion))   # = 1.38
#' 1 - q_o / q_c                          # kappa_w = .348   (Formula 8)
#' v2_o <- with(diagnosis_agreement,
#'              sum(disagreement_weight^2 * observed_proportion))
#' v2_c <- with(diagnosis_agreement,
#'              sum(disagreement_weight^2 * expected_proportion))
#' sqrt((v2_o - q_o^2) / (200 * q_c^2))   # = .0901   (Formula 10)
#' sqrt((v2_c - q_c^2) / (200 * q_c^2))   # = .0916   (Formula 13)
#' (1 - q_o / q_c) + c(-1, 1) * qnorm(0.975) * 0.0901  # 95% CI [.171, .524]
#' (1 - q_o / q_c) / 0.0916               # z = 3.80, p < .001
#'
#' # Cohen's validity reinterpretation: Judge A is a diagnostic panel
#' # (the criterion), Judge B a computer diagnosis (the predictor), and
#' # the weights are asymmetric because the two directions of a
#' # confusion carry different costs. Oriented to the table's layout
#' # (rows = Judge B = computer, columns = Judge A = panel), the weights
#' # below reproduce Cohen's published quantities: sum(v * p_o) = .86,
#' # sum(v * p_c) = 1.33, kappa_w = .353, with his Formulas 10 and 13
#' # giving .0887 and .0915. The weighted kappa vignette works this
#' # example in full, including the orientation of the paper's printed
#' # weight display relative to these values.
#' v_validity <- matrix(
#'   c(0, 1, 2,
#'     1, 0, 2,
#'     4, 6, 0), nrow = 3, byrow = TRUE)
#' cohen_kappa(table = tab, weights = v_validity,
#'             weight_scaling = "disagreement")
#'
#' # ---------------------------------------------------------------------
#' # Raw rater vectors and ordinal categories (quadratic agreement
#' # weights, e.g., Likert-style severity).
#' set.seed(113)
#' x <- sample(1:5, 100, replace = TRUE)
#' y <- pmin(pmax(x + sample(-1:1, 100, replace = TRUE), 1), 5)
#' cohen_kappa(x, y, weights = "quadratic")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{diagnosis_agreement}} (Cohen's 1968 Table 1 as
#'   a data set), \code{\link{fleiss_kappa}}, \code{\link{icc}}
#'
#' @family reliability
#'
#' @keywords htest
#'
#' @export
#' @import stats
cohen_kappa <- function(rater_1 = NULL, rater_2 = NULL,
                        table      = NULL,
                        weights    = "unweighted",
                        weight_scaling = c("agreement", "disagreement"),
                        categories = NULL,
                        conf_level = 0.95,
                        ci_method  = c("wald", "percentile", "bca"),
                        B = 10000L, seed = NULL) {
  weight_scaling <- match.arg(weight_scaling)
  ci_method <- match.arg(ci_method)
  if (ci_method != "wald" &&
      (!is.numeric(B) || length(B) != 1L || B < 100)) {
    stop("'B' must be a single integer of at least 100 when ",
         "bootstrapping.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  if (!is.null(table)) {
    if (!is.null(rater_1) || !is.null(rater_2)) {
      stop("Supply either 'rater_1' and 'rater_2' or 'table', not both.",
           call. = FALSE)
    }
    tab <- as.matrix(table)
    if (!is.numeric(tab) || nrow(tab) != ncol(tab) || nrow(tab) < 2L) {
      stop("'table' must be a square numeric frequency matrix with at ",
           "least 2 categories.", call. = FALSE)
    }
    if (anyNA(tab) || any(tab < 0)) {
      stop("'table' must contain non-negative frequencies with no ",
           "missing values.", call. = FALSE)
    }
    rn <- rownames(tab)
    cn <- colnames(tab)
    if (!is.null(rn) || !is.null(cn)) {
      if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
        stop("'table' row and column names must be present together and ",
             "identical (same categories in the same order).",
             call. = FALSE)
      }
      if (!is.null(categories) &&
          !identical(as.character(categories), rn)) {
        stop("'categories' conflicts with the dimnames of 'table'; drop ",
             "one of the two.", call. = FALSE)
      }
      categories <- rn
    } else if (is.null(categories)) {
      categories <- as.character(seq_len(nrow(tab)))
    } else if (length(categories) != nrow(tab)) {
      stop("'categories' must have one label per row of 'table'.",
           call. = FALSE)
    }
    k <- nrow(tab)
    N <- sum(tab)
    if (N < 2) stop("At least 2 paired ratings are required.", call. = FALSE)
  } else {
    if (is.null(rater_1) || is.null(rater_2)) {
      stop("Supply either 'rater_1' and 'rater_2' or 'table'.",
           call. = FALSE)
    }
    if (length(rater_1) != length(rater_2)) {
      stop("'rater_1' and 'rater_2' must have the same length.",
           call. = FALSE)
    }

    ok <- !is.na(rater_1) & !is.na(rater_2)
    rater_1 <- rater_1[ok]
    rater_2 <- rater_2[ok]
    N <- length(rater_1)
    if (N < 2L) stop("At least 2 paired ratings are required.", call. = FALSE)

    observed <- sort(unique(c(as.character(rater_1), as.character(rater_2))))
    if (is.null(categories)) {
      if (is.factor(rater_1) && is.factor(rater_2) &&
          identical(levels(rater_1), levels(rater_2))) {
        # Two factors with the same level set carry their own category
        # order; respect it. Alphabetical sorting would silently
        # misalign ordinal labels (and any custom weight matrix) whose
        # natural order is not alphabetical.
        categories <- levels(rater_1)
      } else {
        categories <- observed
      }
    } else {
      # A supplied 'categories' set that omits an observed rating would silently
      # turn those ratings into NA levels, drop them from the confusion table, yet
      # still divide by the original N, corrupting kappa. Require the set to be
      # complete (and clean) so any omission fails loudly instead.
      if (anyNA(categories)) {
        stop("'categories' must not contain missing values.", call. = FALSE)
      }
      categories <- as.character(categories)
      if (anyDuplicated(categories)) {
        dup <- unique(categories[duplicated(categories)])
        stop("'categories' must be unique; duplicated: ",
             paste(dup, collapse = ", "), ".", call. = FALSE)
      }
      omitted <- setdiff(observed, categories)
      if (length(omitted) > 0L) {
        stop("'categories' omits observed rating(s): ",
             paste(omitted, collapse = ", "),
             ". The supplied 'categories' must contain every observed rating; ",
             "omitting a category would silently drop those ratings and ",
             "corrupt kappa.", call. = FALSE)
      }
    }
    rater_1 <- factor(as.character(rater_1), levels = as.character(categories))
    rater_2 <- factor(as.character(rater_2), levels = as.character(categories))
    k <- length(categories)
    if (k < 2L) stop("Need at least 2 categories.", call. = FALSE)

    tab <- base::table(rater_1, rater_2)
  }

  P     <- as.matrix(tab) / N
  p_row <- rowSums(P)
  p_col <- colSums(P)
  W     <- .build_kappa_weights(weights, k, scaling = weight_scaling)

  p_o <- sum(W * P)
  p_e <- sum(W * outer(p_row, p_col))

  if (1 - p_e <= .Machine$double.eps) {
    kappa <- NA_real_
    se    <- NA_real_
    z     <- NA_real_
    p     <- NA_real_
    lo    <- NA_real_
    hi    <- NA_real_
  } else {
    kappa <- (p_o - p_e) / (1 - p_e)
    W_row <- as.vector(W %*% p_col)        # length k: bar W_{i.}
    W_col <- as.vector(t(W) %*% p_row)     # length k: bar W_{.j}
    inner <- 0
    for (i in seq_len(k)) {
      for (j in seq_len(k)) {
        if (P[i, j] > 0) {
          term  <- W[i, j] - (W_row[i] + W_col[j]) * (1 - kappa)
          inner <- inner + P[i, j] * term^2
        }
      }
    }
    inner <- inner - (kappa - p_e * (1 - kappa))^2
    var_kappa <- inner / (N * (1 - p_e)^2)
    se <- sqrt(max(0, var_kappa))
    z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    lo <- kappa - z_crit * se
    hi <- kappa + z_crit * se
    z  <- if (se == 0) NA_real_ else kappa / se
    p  <- if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))
  }

  if (ci_method != "wald" && is.finite(kappa)) {
    # Subject-level integer codes for fast resampling; table input is
    # expanded to the equivalent paired ratings first.
    if (!is.null(table)) {
      counts <- as.vector(t(tab))
      i1 <- rep(rep(seq_len(k), each = k), counts)
      i2 <- rep(rep(seq_len(k), times = k), counts)
    } else {
      i1 <- as.integer(rater_1)
      i2 <- as.integer(rater_2)
    }
    kappa_at <- function(idx) {
      cnt <- tabulate((i1[idx] - 1L) * k + i2[idx], nbins = k * k)
      Pb  <- matrix(cnt, k, k, byrow = TRUE) / length(idx)
      po  <- sum(W * Pb)
      pe  <- sum(W * outer(rowSums(Pb), colSums(Pb)))
      if (1 - pe <= .Machine$double.eps) return(NA_real_)
      (po - pe) / (1 - pe)
    }
    if (!is.null(seed)) {
      has_old <- exists(".Random.seed", envir = globalenv())
      old_seed <- if (has_old) get(".Random.seed", envir = globalenv())
      on.exit({
        if (has_old) assign(".Random.seed", old_seed, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv()))
          rm(".Random.seed", envir = globalenv())
      }, add = TRUE)
      set.seed(seed)
    }
    B <- as.integer(B)
    boots <- vapply(seq_len(B), function(b) {
      kappa_at(sample.int(N, N, replace = TRUE))
    }, numeric(1))
    n_bad <- sum(!is.finite(boots))
    if (n_bad > 0L) {
      boots <- boots[is.finite(boots)]
      if (length(boots) < 100L) {
        stop("Only ", length(boots), " of ", B, " bootstrap replications ",
             "returned a defined kappa; the interval would not be ",
             "trustworthy.", call. = FALSE)
      }
      warning(n_bad, " of ", B, " bootstrap replications left kappa ",
              "undefined (chance agreement 1) and were dropped; the ",
              "interval is computed from the ", length(boots),
              " that did.", call. = FALSE)
    }
    alpha_2 <- (1 - conf_level) / 2
    if (ci_method == "percentile") {
      lims <- stats::quantile(boots, c(alpha_2, 1 - alpha_2),
                              names = FALSE)
    } else {
      # BCa: median bias from the bootstrap distribution, acceleration
      # from the jackknife (Efron & Tibshirani, 1993).
      z0 <- stats::qnorm(mean(boots < kappa))
      jack <- vapply(seq_len(N), function(i) {
        kappa_at(seq_len(N)[-i])
      }, numeric(1))
      jack <- jack[is.finite(jack)]
      jm <- mean(jack)
      acc <- sum((jm - jack)^3) / (6 * (sum((jm - jack)^2))^1.5)
      zq <- stats::qnorm(c(alpha_2, 1 - alpha_2))
      adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
      lims <- stats::quantile(boots, adj, names = FALSE)
    }
    lo <- lims[1L]
    hi <- lims[2L]
  }

  # Per-cell detail in the form of Cohen's (1968) Table 1: one row per
  # cell of the confusion matrix (row-major, rater 1 indexing the rows)
  # with the observed proportion, the chance-expected proportion, and
  # the weights. Attached as the "cells" attribute so the one-row
  # summary keeps its canonical shape.
  cells <- data.frame(
    rater_1 = factor(rep(categories, each = k),
                     levels = as.character(categories)),
    rater_2 = factor(rep(categories, times = k),
                     levels = as.character(categories)),
    observed_proportion = as.vector(t(P)),
    expected_proportion = as.vector(t(outer(p_row, p_col))),
    weight              = as.vector(t(W)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (weight_scaling == "disagreement") {
    cells$disagreement_weight <- as.vector(t(weights))
  }

  out <- data.frame(
    weights      = .weights_label(weights, weight_scaling),
    kappa        = kappa,
    se           = se,
    lower_limit  = lo,
    upper_limit  = hi,
    z_value      = z,
    p_value      = p,
    n            = N,
    n_categories = k,
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "cells") <- cells
  attr(out, "ci_method") <- ci_method
  if (ci_method != "wald" && is.finite(kappa)) {
    attr(out, "B_used") <- length(boots)
  }
  out
}
