#' Quantify Simple Structure in a Factor Loading Matrix
#'
#' Summarizes how closely a rotated loading matrix approaches Thurstone's
#' simple structure, in which each item loads on as few factors as
#' possible so that the factors are interpretable. Three complementary
#' quantities are reported: the mean item complexity (the average number
#' of factors an item effectively loads on, one for a perfectly simple
#' item), the hyperplane proportion (the share of loadings near zero,
#' which Thurstone sought to maximize), and the counts of pure versus
#' complex items at a salience cutoff. Together they turn a visual
#' impression of a loading matrix into numbers.
#'
#' @param Lambda A numeric matrix of factor loadings, items in rows and
#'   factors in columns (for example \code{unclass(psych::fa(...)$loadings)}
#'   or a lavaan standardized loading matrix). Row names, if present, label
#'   the items.
#' @param salient Absolute loading at or above which an item is counted as
#'   loading saliently on a factor. Defaults to 0.30 (about ten percent of
#'   an item's variance), a common floor for a meaningful loading.
#' @param hyperplane Absolute loading below which a loading is treated as
#'   lying in the hyperplane (effectively zero). Defaults to 0.10.
#'
#' @details
#' Item complexity is Hofmann's complexity index, proposed in Hofmann
#' (1977) and given as Equation 1 of Hofmann (1978),
#' \eqn{c_i = (\sum_j \lambda_{ij}^2)^2 / \sum_j \lambda_{ij}^4}, which
#' equals one when an item loads on a single factor and rises toward the
#' number of factors as the loadings spread out; it is the same complexity
#' that \code{psych::fa} reports. The \code{"complexity"} attribute holds
#' these per-item values, and the \code{mean_complexity} row is their
#' arithmetic average, what Hofmann (1978) calls the total matrix
#' complexity. These are complexities, not Kaiser's (1974) simplicity
#' index; Hofmann (1978) shows that either can be derived from the other
#' at the item level, with the simplicity of item \eqn{i} in an
#' \eqn{m}-factor solution given by his Equation 3,
#' \eqn{s_i = [1/(m - 1)][(m / c_i) - 1]}. An item is \emph{pure} when
#' exactly one of its loadings is salient and \emph{complex} when more
#' than one is. The hyperplane proportion is the fraction of all loadings
#' whose absolute value is below \code{hyperplane}; a clean simple
#' structure is mostly such near-zero loadings.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   summary quantity (\code{term}, \code{value}): the number of items and
#'   factors, the mean and median item complexity, the hyperplane
#'   proportion, and the counts and proportion of pure items. The per-item
#'   complexities are attached as the \code{"complexity"} attribute (a named
#'   numeric vector), and the salience and hyperplane cutoffs as the
#'   \code{"salient"} and \code{"hyperplane"} attributes.
#'
#' @references
#' Hofmann, R. J. (1977). Indices descriptive of factor complexity.
#'   \emph{The Journal of General Psychology, 96}, 58--66.
#'
#' Hofmann, R. J. (1978). Complexity and simplicity as objective indices
#'   descriptive of factor solutions. \emph{Multivariate Behavioral
#'   Research, 13}(2), 247--250.
#'
#' Kaiser, H. F. (1974). An index of factorial simplicity.
#'   \emph{Psychometrika, 39}(1), 31--36.
#'
#' Thurstone, L. L. (1947). \emph{Multiple-factor analysis}. University of
#'   Chicago Press.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{average_variance_extracted}} and \code{\link{htmt}}
#'   for the convergent and discriminant sides of an exploratory solution.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # A nearly simple two-factor structure: six items, three per factor.
#' Lambda <- rbind(
#'   i1 = c(0.80, 0.05), i2 = c(0.75, 0.10), i3 = c(0.70, -0.05),
#'   i4 = c(0.08, 0.78), i5 = c(-0.04, 0.72), i6 = c(0.30, 0.60))
#' simple_structure(Lambda)
#' attr(simple_structure(Lambda), "complexity")
#'
#' @export
simple_structure <- function(Lambda, salient = 0.30, hyperplane = 0.10) {
  if (is.data.frame(Lambda)) Lambda <- as.matrix(Lambda)
  if (!is.matrix(Lambda) || !is.numeric(Lambda) || nrow(Lambda) < 2L ||
      ncol(Lambda) < 1L || anyNA(Lambda)) {
    stop("'Lambda' must be a numeric loading matrix (items in rows, ",
         "factors in columns) with no missing values.", call. = FALSE)
  }
  for (a in list(salient = salient, hyperplane = hyperplane)) {
    if (!is.numeric(a) || length(a) != 1L || is.na(a) || a < 0 || a > 1) {
      stop("'salient' and 'hyperplane' must each be a single number in ",
           "[0, 1].", call. = FALSE)
    }
  }
  items <- rownames(Lambda)
  if (is.null(items)) items <- paste0("item", seq_len(nrow(Lambda)))
  p <- nrow(Lambda); m <- ncol(Lambda)

  ss2 <- rowSums(Lambda^2)
  ss4 <- rowSums(Lambda^4)
  complexity <- ifelse(ss4 > 0, ss2^2 / ss4, NA_real_)
  names(complexity) <- items

  n_salient <- rowSums(abs(Lambda) >= salient)
  n_pure <- sum(n_salient == 1L)
  n_complex <- sum(n_salient > 1L)
  hyperplane_proportion <- mean(abs(Lambda) < hyperplane)

  out <- data.frame(
    term = c("items", "factors", "mean_complexity", "median_complexity",
             "hyperplane_proportion", "n_pure", "n_complex", "proportion_pure"),
    value = c(p, m, mean(complexity, na.rm = TRUE),
              stats::median(complexity, na.rm = TRUE),
              hyperplane_proportion, n_pure, n_complex, n_pure / p),
    stringsAsFactors = FALSE, row.names = NULL
  )
  res <- .as_dmar_tbl(out)
  attr(res, "complexity") <- complexity
  attr(res, "salient") <- salient
  attr(res, "hyperplane") <- hyperplane
  res
}
