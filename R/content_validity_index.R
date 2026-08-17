# Content validity from a panel of subject matter expert ratings:
# I-CVI, modified kappa, Lawshe's CVR, and the scale level S-CVI summaries.
#' Content Validity Index From Expert Ratings
#'
#' Quantifies how well a pool of candidate items covers the construct it is
#' meant to measure, using the relevance ratings of a panel of subject matter
#' experts. Validation begins before any respondent data are collected: experts
#' rate each item for relevance, and the content validity index summarizes their
#' agreement (Lynn, 1986; Polit & Beck, 2006; Bandalos, 2018). The function
#' returns the item level index (I-CVI) with an exact binomial confidence
#' interval, the chance corrected modified kappa of Polit, Beck, and Owen
#' (2007), Lawshe's (1975) content validity ratio, and the two scale level
#' summaries S-CVI/Ave and S-CVI/UA. The confidence interval is what keeps a
#' small panel from being read as more informative than it is: with five experts
#' an I-CVI of 0.80 carries an interval roughly seven tenths of the width of the scale.
#'
#' @param ratings A matrix or \code{data.frame} of expert ratings with items in
#'   rows and experts in columns, on the conventional 4 point relevance scale
#'   (1 = not relevant, 2 = somewhat relevant, 3 = quite relevant,
#'   4 = highly relevant). Every non-missing entry must be one of 1, 2, 3, or 4.
#'   Row names, when present, name the items; otherwise items are labeled
#'   \code{item_1}, \code{item_2}, and so on. At least 1 item and at least 2
#'   experts are required. \code{NA} is allowed and means that expert did not
#'   rate that item; a missing entry lowers that item's expert count and thereby
#'   the denominator of its I-CVI, kappa, and CVR, leaving other items
#'   untouched. A row with no ratings at all is an error.
#' @param relevant The rating values counted as relevant. Default
#'   \code{c(3, 4)}, the standard dichotomization of the 4 point scale into
#'   relevant (3 or 4) versus not relevant (1 or 2). Must be a subset of
#'   \code{1:4}.
#' @param essential Optional. The rating values counted as \dQuote{essential}
#'   for Lawshe's content validity ratio, when the panel answered the
#'   essential-versus-not question on a separate part of the scale. Must be a
#'   subset of \code{1:4}. When \code{NULL} (default), the content validity
#'   ratio is computed from the same dichotomization as relevance, that is from
#'   \code{relevant}.
#' @param conf_level Confidence level for the exact binomial interval on each
#'   I-CVI. Default \code{0.95}. Must be in (0, 1).
#'
#' @details
#' Let \eqn{N} be the number of experts who rated an item and \eqn{A} the
#' number of those who rated it relevant.
#'
#' \strong{Item level index.} The item level content validity index is the
#' proportion of rating experts who called the item relevant,
#' \deqn{\mathrm{I\mbox{-}CVI} = A / N.}{I-CVI = A / N.}
#' Lynn (1986) gives the conventional criteria: with five or fewer experts an
#' item is expected to reach 1.00, and with six to ten experts at least 0.78.
#'
#' \strong{Confidence interval.} The I-CVI is a binomial proportion, so the
#' interval reported here is the exact (Clopper-Pearson) interval from
#' \code{\link[stats]{binom.test}} at \code{conf_level}. Expert panels are
#' small by design and the interval states plainly how little a handful of
#' ratings pins down the population proportion. This is the accuracy in
#' parameter estimation view of content validity: if the interval is too wide
#' to act on, the remedy is more experts.
#'
#' \strong{Chance corrected agreement.} Some of the observed agreement on
#' relevance would occur if experts responded at random with probability 0.5,
#' so Polit, Beck, and Owen (2007) correct the index in the manner of a kappa.
#' The probability of chance agreement is the binomial point probability
#' \deqn{p_c = \frac{N!}{A!\,(N - A)!}\, 0.5^{N},}{p_c = [N! / (A! (N - A)!)] * 0.5^N,}
#' and the modified kappa is
#' \deqn{\kappa = \frac{\mathrm{I\mbox{-}CVI} - p_c}{1 - p_c}.}{kappa = (I-CVI - p_c) / (1 - p_c).}
#' The binomial coefficient is evaluated on the log scale with
#' \code{\link[base]{lchoose}} and then exponentiated, so a large panel does
#' not overflow the way \code{factorial()} would.
#'
#' \strong{Content validity ratio.} Lawshe (1975) asked a panel whether each
#' item measures behavior that is essential to the performance domain. With
#' \eqn{n_e} experts calling the item essential,
#' \deqn{\mathrm{CVR} = \frac{n_e - N / 2}{N / 2},}{CVR = (n_e - N / 2) / (N / 2),}
#' which equals 1 when every expert says essential, 0 when exactly half do, and
#' -1 when none do.
#'
#' \strong{Scale level summaries.} S-CVI/Ave is the mean of the I-CVIs over
#' items, the averaging approach Polit and Beck (2006) recommend reporting.
#' S-CVI/UA is the universal agreement proportion, the fraction of items whose
#' I-CVI equals 1, a stricter and considerably more conservative summary.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per item,
#'   in the order the items appear in \code{ratings}, and columns:
#'   \describe{
#'     \item{\code{item}}{The item name, from the row names of \code{ratings}
#'       when present.}
#'     \item{\code{n_experts}}{Number of experts who rated the item (missing
#'       ratings excluded).}
#'     \item{\code{n_relevant}}{Number of those experts whose rating was in
#'       \code{relevant}.}
#'     \item{\code{i_cvi}}{The item level content validity index,
#'       \code{n_relevant / n_experts}.}
#'     \item{\code{ci_lower}, \code{ci_upper}}{Limits of the exact binomial
#'       confidence interval for \code{i_cvi} at \code{conf_level}.}
#'     \item{\code{kappa}}{The modified kappa of Polit, Beck, and Owen (2007),
#'       the I-CVI adjusted for chance agreement.}
#'     \item{\code{cvr}}{Lawshe's content validity ratio, computed from
#'       \code{essential} when supplied and from \code{relevant} otherwise.}
#'   }
#'   Attributes: \code{"s_cvi_ave"}, the mean of the \code{i_cvi} column;
#'   \code{"s_cvi_ua"}, the proportion of items with \code{i_cvi} equal to 1;
#'   \code{"relevant"}, the rating values counted as relevant;
#'   \code{"essential"}, the rating values counted as essential for the content
#'   validity ratio (equal to \code{"relevant"} when \code{essential} was
#'   \code{NULL}); and \code{"conf_level"}, the confidence level used.
#'
#' @references
#' Bandalos, D. L. (2018). \emph{Measurement theory and applications for the
#'   social sciences}. Guilford Press.
#'
#' Lawshe, C. H. (1975). A quantitative approach to content validity.
#'   \emph{Personnel Psychology, 28}(4), 563--575.
#'
#' Lynn, M. R. (1986). Determination and quantification of content validity.
#'   \emph{Nursing Research, 35}(6), 382--385.
#'
#' Polit, D. F., & Beck, C. T. (2006). The content validity index: Are you sure
#'   you know what's being reported? Critique and recommendations.
#'   \emph{Research in Nursing and Health, 29}(5), 489--497.
#'   \doi{10.1002/nur.20147}
#'
#' Polit, D. F., Beck, C. T., & Owen, S. V. (2007). Is the CVI an acceptable
#'   indicator of content validity? Appraisal and recommendations.
#'   \emph{Research in Nursing and Health, 30}(4), 459--467.
#'   \doi{10.1002/nur.20199}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{gwet_ac}}, \code{\link{fleiss_kappa}} for agreement
#'   among raters on a common set of units.
#'
#' @family agreement and measurement
#'
#' @keywords multivariate
#'
#' @examples
#' # Six experts rate five candidate items on the 4 point relevance scale.
#' ratings <- rbind(
#'   item_1 = c(4, 4, 3, 4, 4, 3),
#'   item_2 = c(4, 3, 4, 4, 3, 2),
#'   item_3 = c(2, 3, 1, 2, 3, 2),
#'   item_4 = c(4, 4, 4, 4, 4, 4),
#'   item_5 = c(3, 4, 4, 3, NA, 4))
#' colnames(ratings) <- paste0("expert_", 1:6)
#' cvi <- content_validity_index(ratings)
#' cvi
#'
#' # The scale level summaries travel with the table as attributes.
#' attr(cvi, "s_cvi_ave")
#' attr(cvi, "s_cvi_ua")
#'
#' # The worked example of Polit, Beck, and Owen (2007): 6 experts, 5 of whom
#' # rate the item relevant, gives I-CVI = 0.83, p_c = 0.094, and a modified
#' # kappa of 0.816 (the paper reports 0.81, carrying its rounded I-CVI).
#' content_validity_index(matrix(c(4, 4, 3, 4, 3, 1), nrow = 1))
#'
#' # A wide interval is the point: with 5 experts an I-CVI of 0.80 is
#' # compatible with a population proportion anywhere from about 0.28 to 0.99.
#' content_validity_index(matrix(c(4, 4, 3, 4, 1), nrow = 1))
#'
#' # The broom verbs on the earlier result: one row per item, and the
#' # scale-level summary.
#' generics::tidy(cvi)
#' generics::glance(cvi)
#'
#' @export
content_validity_index <- function(ratings, relevant = c(3, 4),
                                   essential = NULL, conf_level = 0.95) {
  if (is.data.frame(ratings)) {
    item_names <- row.names(ratings)
    if (identical(item_names, as.character(seq_len(nrow(ratings)))))
      item_names <- NULL
    ratings <- as.matrix(ratings)
    if (!is.null(item_names)) row.names(ratings) <- item_names
  }
  if (!is.matrix(ratings))
    stop("'ratings' must be a matrix or data.frame with items in rows and ",
         "experts in columns.", call. = FALSE)
  if (!is.numeric(ratings))
    stop("'ratings' must be numeric: every non-missing entry is a rating on ",
         "the 4 point relevance scale (1, 2, 3, or 4).", call. = FALSE)
  if (nrow(ratings) < 1L)
    stop("'ratings' must contain at least 1 item (row).", call. = FALSE)
  if (ncol(ratings) < 2L)
    stop("'ratings' must contain at least 2 experts (columns); a content ",
         "validity index summarizes agreement across a panel.", call. = FALSE)

  scale_values <- 1:4
  observed <- ratings[!is.na(ratings)]
  if (!all(observed %in% scale_values))
    stop("Every non-missing entry of 'ratings' must be on the 4 point ",
         "relevance scale (1, 2, 3, or 4).", call. = FALSE)
  if (!length(relevant) || !is.numeric(relevant) || !all(relevant %in% scale_values))
    stop("'relevant' must be a non-empty subset of the 4 point relevance ",
         "scale (1, 2, 3, or 4).", call. = FALSE)
  if (!is.null(essential) &&
      (!length(essential) || !is.numeric(essential) ||
       !all(essential %in% scale_values)))
    stop("'essential' must be NULL or a non-empty subset of the 4 point ",
         "relevance scale (1, 2, 3, or 4).", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)

  n_experts <- rowSums(!is.na(ratings))
  if (any(n_experts == 0L))
    stop("Every item must be rated by at least one expert; rows ",
         paste(which(n_experts == 0L), collapse = ", "),
         " of 'ratings' are entirely missing.", call. = FALSE)

  essential_used <- if (is.null(essential)) relevant else essential
  n_relevant  <- rowSums(matrix(ratings %in% relevant, nrow = nrow(ratings)) &
                           !is.na(ratings))
  n_essential <- rowSums(matrix(ratings %in% essential_used,
                                nrow = nrow(ratings)) & !is.na(ratings))

  i_cvi <- n_relevant / n_experts

  # Exact (Clopper-Pearson) interval for each item's proportion relevant.
  limits <- vapply(seq_len(nrow(ratings)), function(i) {
    stats::binom.test(n_relevant[i], n_experts[i],
                      conf.level = conf_level)$conf.int[1:2]
  }, numeric(2))

  # Chance agreement is the binomial point probability of exactly n_relevant of
  # n_experts agreeing under random responding with probability 0.5. Taken on
  # the log scale so a large panel does not overflow factorial().
  p_c <- exp(lchoose(n_experts, n_relevant) + n_experts * log(0.5))
  kappa <- (i_cvi - p_c) / (1 - p_c)

  cvr <- (n_essential - n_experts / 2) / (n_experts / 2)

  item <- row.names(ratings)
  if (is.null(item)) item <- paste0("item_", seq_len(nrow(ratings)))

  out <- data.frame(
    item       = item,
    n_experts  = as.numeric(n_experts),
    n_relevant = as.numeric(n_relevant),
    i_cvi      = as.numeric(i_cvi),
    ci_lower   = limits[1, ],
    ci_upper   = limits[2, ],
    kappa      = as.numeric(kappa),
    cvr        = as.numeric(cvr),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  res <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_content_validity")
  attr(res, "s_cvi_ave") <- mean(out$i_cvi)
  attr(res, "s_cvi_ua")  <- mean(out$i_cvi == 1)
  attr(res, "relevant")  <- relevant
  attr(res, "essential") <- essential_used
  res
}
