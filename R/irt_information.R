#' Item and Test Information for the Graded Response Model
#'
#' Evaluates the item information functions and the test information
#' function of a graded response model on a grid of latent trait values,
#' together with the standard error of the latent trait estimate,
#' \eqn{SE(\theta) = 1 / \sqrt{I(\theta)}}. Reliability is a single number
#' that describes a scale at one place on the latent continuum; the
#' information function is the same idea expressed as a function of where
#' the respondent sits, so an item pool can be judged on where it measures
#' precisely rather than on one global summary. Because information is
#' additive across items, the curve also shows which items carry the
#' precision, and over what range, which is what makes it useful for
#' building and trimming a scale.
#'
#' @param a Discriminations, a numeric vector of positive values. Supply
#'   either one value per item (named with the item names, or in the order
#'   the items first appear in \code{item}) or one value per boundary row
#'   (constant within an item). Not used when \code{grm} is supplied.
#' @param b Boundary locations (category thresholds), a numeric vector with
#'   one element per category boundary. An item with \eqn{m} categories has
#'   \eqn{m - 1} boundaries, which the model requires to be in ascending
#'   order within the item. Not used when \code{grm} is supplied.
#' @param item Item labels, a character or factor vector the same length as
#'   \code{b} naming the item each boundary belongs to. When \code{NULL}
#'   (default), each element of \code{b} is treated as its own dichotomous
#'   item, named \code{item_1}, \code{item_2}, and so on, and \code{a} must
#'   then have the same length as \code{b}. Not used when \code{grm} is
#'   supplied.
#' @param theta Latent trait values at which to evaluate the information
#'   functions. Any finite numeric vector; the default,
#'   \code{seq(-4, 4, length.out = 81)}, covers the range in which almost
#'   all of a standard normal trait distribution falls, in steps of 0.1.
#' @param grm Optionally, the result of \code{irt_grm()}: a
#'   \code{data.frame} with one row per item and category boundary and
#'   columns \code{item}, \code{a}, and \code{b} (a \code{category} column,
#'   when present, orders the boundaries within an item). Supply this or
#'   the parameters, not both.
#'
#' @details
#' For item \emph{i} with discrimination \eqn{a_i} and ordered boundary
#' locations \eqn{b_{i1} < b_{i2} < \cdots < b_{i,m-1}} for \emph{m}
#' categories, the normal ogive graded response model of Samejima (1969)
#' defines the boundary response function
#' \deqn{P^*_{ik}(\theta) = \Phi[a_i (\theta - b_{ik})],}
#' the probability of responding \emph{above} boundary \emph{k}, that is, in
#' any category higher than the \emph{k}th, with the
#' conventions \eqn{P^*_{i0} = 1} and \eqn{P^*_{im} = 0}. The category
#' response function is the difference of adjacent boundary functions,
#' \deqn{P_{ik}(\theta) = P^*_{i,k-1}(\theta) - P^*_{ik}(\theta),}
#' and differentiating with respect to \eqn{\theta} gives
#' \deqn{P'_{ik}(\theta) = a_i \{\phi[a_i (\theta - b_{i,k-1})] -
#'   \phi[a_i (\theta - b_{ik})]\},}
#' where \eqn{\phi} is the standard normal density and the density terms
#' vanish at the two extreme categories (there is no \eqn{b_{i0}} and no
#' \eqn{b_{im}}). Item information is
#' \deqn{I_i(\theta) = \sum_{k=1}^{m} \frac{[P'_{ik}(\theta)]^2}{P_{ik}(\theta)},}
#' test information is \eqn{I(\theta) = \sum_i I_i(\theta)}, and the
#' standard error of the maximum likelihood estimate of \eqn{\theta} is
#' \eqn{SE(\theta) = 1 / \sqrt{I(\theta)}}.
#'
#' Two properties make the curve worth reading. Information is additive
#' across items, so an item's contribution can be read off directly and a
#' pool can be assembled to cover a targeted range. And the reciprocal
#' relation to the squared standard error means the peak of the curve
#' locates where the scale estimates the trait most precisely, reported
#' here as the \code{"theta_max_information"} attribute.
#'
#' For a dichotomous item the model reduces to the two parameter normal
#' ogive, whose information has the closed form
#' \deqn{I_i(\theta) = \frac{a_i^2 \phi[a_i(\theta - b_i)]^2}{
#'   \Phi[a_i(\theta - b_i)] \{1 - \Phi[a_i(\theta - b_i)]\}},}
#' which the general expression above reproduces; that identity is one of
#' the tests of this function.
#'
#' The category probabilities underflow to zero for \eqn{\theta} far from
#' every boundary, where the ratio \eqn{(P')^2 / P} would be \eqn{0/0}. A
#' category whose probability is not strictly positive contributes zero to
#' the sum, which is the limit the ratio approaches, so the returned
#' information is finite and nonnegative on any grid, however extreme, and
#' is never \code{NaN}. In the regime where \eqn{(P')^2} underflows but
#' \eqn{P} does not, the ratio is formed as
#' \eqn{\exp[2 \log |P'| - \log P]} so the contribution is kept rather than
#' flushed to zero. Where two boundaries of an item coincide, the category
#' between them has probability zero everywhere and, by the same guard,
#' contributes nothing.
#'
#' The parameters are in the normal ogive metric, which is what
#' \code{irt_grm()} returns by default. The logistic metric used by much of
#' the item response theory software scales the discrimination by
#' approximately 1.702 (Camilli, 1994); a logistic \eqn{a} is put on the
#' normal ogive scale by dividing by that constant. The two metrics give
#' information functions that are proportional in shape but not equal in
#' value, so a cross-software comparison is a comparison of curves, not of
#' numbers.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   value of \code{theta} and columns:
#'   \describe{
#'     \item{\code{theta}}{The latent trait value, as supplied.}
#'     \item{\code{test_information}}{Test information at that value, the
#'       sum of the item information functions.}
#'     \item{\code{se}}{The standard error of the latent trait estimate,
#'       \eqn{1 / \sqrt{I(\theta)}}. It is \code{Inf} where test
#'       information is zero, which is the correct statement that the
#'       items carry no information there.}
#'   }
#'   The result carries these attributes:
#'   \describe{
#'     \item{\code{"item_information"}}{A numeric matrix of item
#'       information with \code{theta} in the rows (row names are the
#'       \code{theta} values) and items in the columns (column names are
#'       the item names). Its row sums are \code{test_information}.}
#'     \item{\code{"item"}}{The item names, in the order they appear in the
#'       columns of \code{"item_information"}.}
#'     \item{\code{"a"}}{The discrimination used for each item, a numeric
#'       vector named by item.}
#'     \item{\code{"b"}}{The boundary locations used, a numeric vector in
#'       item order and, within an item, in ascending order, named by the
#'       item each boundary belongs to.}
#'     \item{\code{"theta_max_information"}}{The value of \code{theta} at
#'       which test information peaks on the supplied grid (the first such
#'       value if there are ties). It is a grid value, not the result of an
#'       optimization, so a finer \code{theta} locates the peak more
#'       sharply.}
#'   }
#'
#' @references
#' Baker, F. B., & Kim, S.-H. (2004). \emph{Item response theory:
#'   Parameter estimation techniques} (2nd ed.). Marcel Dekker.
#'
#' Camilli, G. (1994). Teacher's corner: Origin of the scaling constant
#'   \emph{d} = 1.7 in item response theory. \emph{Journal of Educational
#'   and Behavioral Statistics, 19}(3), 293--295.
#'   \doi{10.3102/10769986019003293}
#'
#' Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for
#'   psychologists}. Lawrence Erlbaum.
#'
#' Lord, F. M. (1980). \emph{Applications of item response theory to
#'   practical testing problems}. Lawrence Erlbaum.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response
#'   pattern of graded scores. \emph{Psychometrika Monograph Supplement,
#'   34}(4, Pt. 2), 1--97.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{plot_irt_information}} for the curve,
#'   \code{\link{reliability_omega}} for the single-number companion.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Three items: a five-category rating item and two dichotomous items.
#' # The discriminations are named, so they are matched to the item labels.
#' info <- irt_information(
#'   a = c(mood_1 = 1.4, mood_2 = 0.9, mood_3 = 1.1),
#'   b = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
#'   item = c(rep("mood_1", 4), "mood_2", "mood_3")
#' )
#' head(info)
#'
#' # Where does this three-item set measure most precisely?
#' attr(info, "theta_max_information")
#'
#' # Each item's contribution; the rows sum to the test information.
#' head(attr(info, "item_information"))
#'
#' # A dichotomous item matches the two parameter normal ogive closed form.
#' one <- irt_information(a = 1.5, b = 0.25, theta = c(-1, 0, 1))
#' z <- 1.5 * (c(-1, 0, 1) - 0.25)
#' 1.5^2 * dnorm(z)^2 / (pnorm(z) * (1 - pnorm(z)))
#' one$test_information
#'
#' @export
irt_information <- function(a, b = NULL, item = NULL,
                            theta = seq(-4, 4, length.out = 81),
                            grm = NULL) {

  if (!is.null(grm)) {
    if (!missing(a) || !is.null(b) || !is.null(item)) {
      stop("Supply either 'grm' or the parameters 'a', 'b', and 'item', ",
           "not both.", call. = FALSE)
    }
    if (!is.data.frame(grm)) {
      stop("'grm' must be a data.frame of graded response model ",
           "parameters, as returned by irt_grm().", call. = FALSE)
    }
    missing_cols <- setdiff(c("item", "a", "b"), names(grm))
    if (length(missing_cols)) {
      stop("'grm' is missing the column(s) ",
           paste(sQuote(missing_cols), collapse = ", "),
           ". Supply the result of irt_grm(), which has one row per item ",
           "and category boundary.", call. = FALSE)
    }
    if (!nrow(grm)) {
      stop("'grm' has no rows, so there are no items to evaluate.",
           call. = FALSE)
    }
    grm_item <- as.character(grm$item)
    if ("category" %in% names(grm)) {
      grm <- grm[order(match(grm_item, unique(grm_item)),
                       as.numeric(grm$category)), , drop = FALSE]
      grm_item <- as.character(grm$item)
    }
    item <- grm_item
    a <- as.numeric(grm$a)
    b <- as.numeric(grm$b)

    # irt_grm() records which discrimination metric its `a` column is in.
    # The information formulas below are the normal ogive ones, so a table
    # reported in the logistic metric has to be converted first, exactly
    # the division by 1.702 this function's own Details section describes.
    # Without it a logistic table was fed straight into the ogive formulas,
    # inflating test information by the square of that constant (about
    # 2.9 times) and understating SE(theta) by about a third, with no
    # error and no warning. A hand-built table carries no metric
    # attribute and is still taken at face value as normal ogive.
    grm_metric <- attr(grm, "metric")
    if (identical(grm_metric, "logistic")) {
      a <- a / 1.702
    }
  } else if (missing(a) || is.null(a)) {
    stop("Supply the discriminations 'a' and the boundary locations 'b', ",
         "or a set of graded response model parameters through 'grm'.",
         call. = FALSE)
  } else if (is.null(b)) {
    stop("'b' must be supplied: one boundary location per category ",
         "boundary, so an item with m categories contributes m - 1 ",
         "elements.", call. = FALSE)
  }

  # ---------- Validate the parameters ----------
  if (!is.numeric(a) || !length(a)) {
    stop("'a' must be a numeric vector of discriminations.", call. = FALSE)
  }
  if (!is.numeric(b) || !length(b)) {
    stop("'b' must be a numeric vector of boundary locations.", call. = FALSE)
  }
  if (any(!is.finite(a))) {
    stop("'a' must be finite; missing, NaN, and infinite discriminations ",
         "are not admissible.", call. = FALSE)
  }
  if (any(!is.finite(b))) {
    stop("'b' must be finite; missing, NaN, and infinite boundary ",
         "locations are not admissible.", call. = FALSE)
  }
  if (any(a <= 0)) {
    stop("every discrimination in 'a' must be greater than zero; ",
         "received a minimum of ", min(a), ".", call. = FALSE)
  }

  if (is.null(item)) {
    if (length(a) != length(b)) {
      stop("With 'item' not supplied, each element of 'b' is a separate ",
           "dichotomous item, so 'a' and 'b' must have the same length; ",
           "'a' has length ", length(a), " and 'b' has length ",
           length(b), ".", call. = FALSE)
    }
    item <- paste0("item_", seq_along(b))
  }
  item <- as.character(item)
  if (anyNA(item)) {
    stop("'item' must not contain missing values.", call. = FALSE)
  }
  if (length(item) != length(b)) {
    stop("'item' and 'b' must have the same length, one entry per ",
         "category boundary; 'item' has length ", length(item),
         " and 'b' has length ", length(b), ".", call. = FALSE)
  }
  items <- unique(item)
  n_items <- length(items)

  # ---------- One discrimination per item ----------
  if (!is.null(names(a)) && length(a) == n_items) {
    # Named 'a' is matched by name. A name that does not match is an error
    # rather than a silent fall-through to positional matching: one stale
    # or mistyped name would otherwise reassign every discrimination in
    # the order the items happen to appear, with no diagnostic.
    unmatched <- setdiff(items, names(a))
    if (length(unmatched)) {
      stop("'a' is named, so it is matched by name, but no element is ",
           "named for the item(s) ",
           paste(sQuote(unmatched), collapse = ", "),
           ". Correct the names, or supply 'a' unnamed to match the items ",
           "in the order they first appear.", call. = FALSE)
    }
    a_by_item <- a[items]
  } else if (length(a) == length(b)) {
    a_by_item <- vapply(items, function(it) {
      a_i <- a[item == it]
      if (any(a_i != a_i[1L])) {
        stop("The discrimination must be constant within an item, but item ",
             sQuote(it), " has more than one value of 'a'.", call. = FALSE)
      }
      a_i[1L]
    }, numeric(1L))
  } else if (length(a) == n_items) {
    a_by_item <- stats::setNames(as.numeric(a), items)
  } else {
    stop("'a' must have one element per item (", n_items,
         ") or one element per category boundary (", length(b),
         "); it has ", length(a), ".", call. = FALSE)
  }
  names(a_by_item) <- items

  # ---------- Boundaries in ascending order within an item ----------
  b_by_item <- split(b, factor(item, levels = items))
  out_of_order <- vapply(b_by_item, is.unsorted, logical(1L))
  if (any(out_of_order)) {
    warning("The boundary locations of item(s) ",
            paste(sQuote(items[out_of_order]), collapse = ", "),
            " were not in ascending order. The graded response model ",
            "requires b_i1 < b_i2 < ... < b_i,m-1, so they have been ",
            "sorted before the information was computed.", call. = FALSE)
    b_by_item <- lapply(b_by_item, sort)
  }

  # ---------- Validate the evaluation grid ----------
  if (!is.numeric(theta) || !length(theta)) {
    stop("'theta' must be a numeric vector of latent trait values.",
         call. = FALSE)
  }
  if (any(!is.finite(theta))) {
    stop("'theta' must be finite; missing, NaN, and infinite values are ",
         "not admissible.", call. = FALSE)
  }

  # ---------- Item and test information ----------
  item_information <- matrix(
    unlist(lapply(items, function(it) {
      .irt_information_item(a_by_item[[it]], b_by_item[[it]], theta)
    }), use.names = FALSE),
    nrow = length(theta), ncol = n_items,
    dimnames = list(as.character(theta), items)
  )
  test_information <- as.numeric(rowSums(item_information))
  se <- 1 / sqrt(test_information)

  out <- data.frame(
    theta = as.numeric(theta),
    test_information = test_information,
    se = se,
    stringsAsFactors = FALSE
  )
  row.names(out) <- NULL

  res <- .as_dmar_tbl(out)
  attr(res, "item_information") <- item_information
  attr(res, "item") <- items
  attr(res, "a") <- a_by_item
  attr(res, "b") <- stats::setNames(
    unlist(b_by_item, use.names = FALSE),
    rep(items, vapply(b_by_item, length, integer(1L)))
  )
  attr(res, "theta_max_information") <- theta[which.max(test_information)]
  res
}


# Information function of one graded response item, evaluated on the whole
# theta grid at once. 'b' is already sorted; 'a' is a single positive value.
# Returns a numeric vector the length of 'theta'.
.irt_information_item <- function(a, b, theta) {
  n_boundaries <- length(b)                 # m - 1 boundaries, m categories
  z <- a * outer(theta, b, "-")             # length(theta) by (m - 1)
  p_lower <- stats::pnorm(z)                # Phi(z), the boundary function
  p_upper <- stats::pnorm(z, lower.tail = FALSE)
  density <- stats::dnorm(z)

  information <- numeric(length(theta))
  for (k in seq_len(n_boundaries + 1L)) {
    if (k == 1L) {
      # P_i1 = 1 - Phi(z_1), computed in the upper tail for accuracy.
      p  <- p_upper[, 1L]
      dp <- -a * density[, 1L]
    } else if (k == n_boundaries + 1L) {
      p  <- p_lower[, n_boundaries]
      dp <- a * density[, n_boundaries]
    } else {
      lo <- k - 1L                          # z[, lo] > z[, hi]
      hi <- k
      # Take the difference in whichever tail avoids cancellation.
      p  <- ifelse(z[, lo] <= 0,
                   p_lower[, lo] - p_lower[, hi],
                   p_upper[, hi] - p_upper[, lo])
      dp <- a * (density[, lo] - density[, hi])
    }
    information <- information + .irt_info_ratio(dp, p)
  }
  information
}


# The per-category contribution (P')^2 / P, guarded so that a category
# whose probability underflows contributes zero (the limit of the ratio)
# instead of NaN, and so that a squared derivative that underflows while
# the probability does not is still recovered on the log scale.
.irt_info_ratio <- function(dp, p) {
  out <- numeric(length(p))
  usable <- is.finite(p) & p > 0 & is.finite(dp)
  underflows <- usable & abs(dp) < 1e-150
  direct <- usable & !underflows
  out[direct] <- dp[direct]^2 / p[direct]
  if (any(underflows)) {
    out[underflows] <- exp(2 * log(abs(dp[underflows])) - log(p[underflows]))
  }
  out[!is.finite(out)] <- 0
  out
}
