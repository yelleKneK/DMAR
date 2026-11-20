#' Tucker's Congruence Coefficients Between the Factors of Two Solutions
#'
#' @description
#' Compares the factors of two solutions for the same items by Tucker's
#' congruence coefficient, after matching each factor of one solution with
#' a factor of the other. The two solutions can come from two samples, two
#' rotations, two estimators, or two factor counts, and each can be a
#' loading matrix, a table from \code{\link{fa_exploratory}()}, or a
#' \pkg{lavaan} fit. The factors of \code{y} are paired one to one with
#' those of \code{x} by the assignment that maximizes the total absolute
#' congruence, and a factor whose match is negative is reflected. Each
#' factor is reported with its coefficient and with the largest absolute
#' coefficient it has with any other factor of the other solution, so that
#' a close competitor is in view, and the full matrix of coefficients
#' travels with the table.
#'
#' @param x The reference solution, in one of three forms. A numeric matrix
#'   of loadings with items in rows and factors in columns (a numeric
#'   vector is one factor, and a \code{"loadings"} object, such as
#'   \code{psych::fa()$loadings} or \code{factanal()$loadings}, is used as
#'   its matrix). A table returned by \code{\link{fa_exploratory}()}, whose
#'   standardized pattern loadings and factor correlations are used. Or a
#'   single-group, single-level \pkg{lavaan} fit, such as
#'   \code{fa_confirmatory(output = "fit")}, \code{lavaan::cfa()}, or
#'   \code{lavaan::efa(output = "lavaan")}, whose standardized loadings and
#'   model implied factor correlations are used. A table of results from
#'   \code{\link{fa_confirmatory}()}, \code{\link{cfa_1}()}, or
#'   \code{\link{cfa_2}()} is not a fit and is refused with a pointer to
#'   \code{output = "fit"}. Items are named by the row names of a matrix,
#'   the names of a vector, or the item names of a fit.
#' @param y The solution compared with \code{x}, in any of the same forms.
#'   Its factors are matched to those of \code{x}.
#' @param match How the factors are paired: \code{"congruence"} (default),
#'   the one-to-one assignment that maximizes the sum of the absolute
#'   coefficients, with each matched factor of \code{y} reflected when its
#'   coefficient is negative; or \code{"order"}, column \emph{j} of
#'   \code{y} with column \emph{j} of \code{x}, signed and never reflected,
#'   for solutions whose factors are already aligned. See Details.
#' @param loadings \code{"pattern"} (default) compares columns of the
#'   pattern matrix, the loadings themselves; \code{"structure"} compares
#'   columns of the structure matrix, the pattern times the factor
#'   correlations, the correlations of the items with the factors. The
#'   structure needs the factor correlations, so it is available for
#'   \code{fa_exploratory()} tables and \pkg{lavaan} fits, not for a
#'   matrix; a matrix is used as supplied, and for a structure comparison
#'   of matrices the structure matrices themselves are passed under the
#'   default.
#'
#' @details
#' \strong{The coefficient.} For two columns of loadings \eqn{a} and
#' \eqn{b} on the same \eqn{n} items, Tucker's congruence coefficient is
#' \deqn{\phi(a, b) = \frac{\sum_{i=1}^{n} a_i b_i}
#'   {\sqrt{\sum_{i=1}^{n} a_i^2 \sum_{i=1}^{n} b_i^2}},}
#' the cosine of the angle between the two columns (Lorenzo-Seva & ten
#' Berge, 2006, eq. 1). Burt suggested it in 1948, and it became known as
#' Tucker's congruence coefficient through Tucker's report of 1951
#' (Lorenzo-Seva & ten Berge, 2006, p. 57). It is not the Pearson
#' correlation of the two columns, which would subtract their means first.
#' The coefficient is unchanged when a column is multiplied by a positive
#' constant, so it compares the pattern of the loadings rather than their
#' size. It changes when a constant is added, as it should, since
#' loadings of .4, .4, and 0 and loadings of .1, .1, and \eqn{-.3} differ
#' by a constant and are read as different factors (their coefficient is
#' .43, their correlation 1). It is unchanged when an item's sign is
#' reversed in both solutions, and it reverses its sign when one of the
#' two factors is reflected (Lorenzo-Seva & ten Berge, 2006,
#' pp. 57--58). The columns are usually those of the pattern matrix,
#' though those of the structure matrix can be compared as well (p. 57).
#' In an oblique solution the two differ, and \code{loadings} chooses
#' between them; in an orthogonal one they coincide.
#'
#' Rescaling a whole column leaves the coefficient alone, but rescaling an
#' item does not, and an item's loadings in the original metric scale with
#' its standard deviation, which differs from sample to sample. A fit is
#' therefore read in the standardized metric. A matrix is used as
#' supplied, so two matrices should be in the same metric.
#'
#' \strong{Why matching is part of the comparison.} The order in which a
#' solution lists its factors, and the sign of each, are arbitrary: a
#' rotation can list the verbal factor first in one sample and second in
#' another, and reflecting a factor (reversing the signs of its loadings
#' and of its correlations with the other factors) gives an equivalent
#' solution. Before two solutions can be compared, each factor of
#' \code{y} has to be paired with a factor of \code{x}. Under
#' \code{match = "congruence"} the pairing is the one-to-one assignment
#' that maximizes the sum of the absolute coefficients over the pairs,
#' found by the Hungarian method (Kuhn, 1955), and each matched factor of
#' \code{y} whose coefficient is negative is then reflected, so the
#' reported coefficients are not negative. The assignment is optimal over
#' the whole matrix rather than factor by factor: pairing each factor of
#' \code{x} with its best factor of \code{y} can give two factors the same
#' partner, and taking the largest coefficient first and then the largest
#' of what remains can give a lower total. Assignments whose totals agree
#' to within rounding are decided in favor of the listed order: the first
#' factor of \code{x} takes the earliest factor of \code{y} it can, then
#' the second, and so on. Under \code{match = "order"} column \emph{j} of
#' \code{y} is paired with column \emph{j} of \code{x}, the coefficients
#' keep their signs, and nothing is reflected; this suits solutions whose
#' factors are aligned by construction, such as two confirmatory fits of
#' one model, and a negative coefficient then shows a factor whose
#' direction differs.
#'
#' When the two solutions have different numbers of factors, as many pairs
#' are formed as the smaller solution has factors, and each surplus factor
#' of the larger solution has a row with no partner and no coefficient.
#' Its \code{next_largest} shows how close it comes to a factor already
#' taken, which is how a factor of one solution that splits into two in
#' the other shows itself.
#'
#' Items are aligned by name when both solutions name them, and the two
#' sets of items must then be the same; otherwise they are aligned by
#' position, and the two numbers of items must be the same.
#'
#' \strong{The matched value is the largest of several candidates.}
#' Matching chooses for each factor the partner it agrees with best, so a
#' matched coefficient is a maximum, larger on average than the
#' coefficient of a pairing fixed in advance. The \code{next_largest}
#' column gives each factor's strongest competitor, the largest absolute
#' coefficient it has with any factor of the other solution other than its
#' partner. A matched value close to its \code{next_largest} marks a
#' pairing that a small change in the data could reverse, and the
#' \code{congruence_matrix} attribute, printed under the table, shows every
#' candidate.
#'
#' \strong{Chance congruence is not zero.} Two columns of loadings
#' unrelated to each other still have a coefficient different from zero.
#' For two unrelated factors whose loadings point in random directions,
#' Korth and Tucker (1975, p. 366) give the expected squared coefficient
#' as \eqn{1/n} for \eqn{n} items (for six items a root mean square of
#' about .41). Matching raises the chance values, and rotating the two
#' solutions toward each other raises them further. Korth and Tucker gave
#' each of two unrelated patterns an orthogonal transformation chosen to
#' make the two agree as closely as possible (pp. 364--365); for two
#' factors on ten items, the first and largest coefficient after that
#' transformation was typically about .64 (their Table 4, p. 368, a mean
#' on the Fisher \emph{Z} scale converted back). Reordering and
#' reflecting factors are orthogonal transformations too, so for the same
#' two patterns no coefficient this function reports exceeds that first
#' coefficient. Columns whose loadings are all positive, as a general
#' factor's often are, agree further by construction: two unrelated
#' columns of loadings drawn independently and uniformly between 0 and 1
#' have a coefficient that approaches 3/4 as the number of items grows. A
#' coefficient is therefore read against what unrelated solutions of the
#' same shape would give, not against zero.
#'
#' \strong{Reading the coefficient.} Tucker declined in his 1951 report to
#' name a value below 1 that would count as identity of two factors until
#' experience with the coefficient had accumulated (Lorenzo-Seva & ten
#' Berge, 2006, p. 58). Lorenzo-Seva and ten Berge (2006) asked 56
#' university professors who use exploratory factor analysis in their
#' research to rate, on a five-point scale, the similarity of pairs of
#' loading columns constructed to have coefficients from .62 to .97. The
#' mean ratings tracked the coefficient closely (\emph{r} = .974 between
#' the coefficient and the mean rating), and the raters' mean judgments
#' led the authors to read values from .85 to .94 as fair similarity and
#' values above .95 as equality. The individual ratings varied widely, and
#' they depended on the number of items and on the judges' research field;
#' the authors concluded that a rule of thumb based on the judgment of a
#' single researcher is hazardous (p. 61). Their study describes how
#' analysts judged similarity, not a property of the coefficient, and the
#' table reports the coefficient and applies no such reading.
#'
#' \strong{What the function does not do.} It rotates neither solution: it
#' compares them as they were rotated. Rotating one solution toward the
#' other before computing the coefficients, as is often done
#' (Lorenzo-Seva & ten Berge, 2006, p. 57), is a separate step. It reports
#' no standard error, interval, or test for the coefficients.
#' \code{\link{procrustes_phi}()} computes the same coefficient for one
#' pair of loading vectors that the user has already paired, with a
#' permutation \emph{p}-value, and does no matching.
#'
#' @return A \code{dmar_tbl} data frame (classes
#'   \code{dmar_factor_congruence}, \code{dmar_tbl}, and
#'   \code{data.frame}) with one row per factor of \code{x}, in the order
#'   of \code{x}, followed by a row for each factor of \code{y} left
#'   without a partner. Its columns are \code{factor_x} and
#'   \code{factor_y}, the names of the paired factors (\code{NA} for a
#'   factor without a partner); \code{congruence}, the coefficient of the
#'   pair after reflection (signed under \code{match = "order"};
#'   \code{NA} without a partner); \code{reflected}, whether the factor of
#'   \code{y} was reflected to match (\code{NA} without a partner); and
#'   \code{next_largest}, the largest absolute coefficient between the
#'   row's factor (its factor of \code{x}, or its factor of \code{y} in
#'   the row of a factor of \code{y} without a partner) and any factor of
#'   the other solution other than its partner (\code{NA} when there is
#'   none). The attributes are \code{congruence_matrix}, the signed
#'   coefficients of every factor of \code{x} (rows) with every factor of
#'   \code{y} (columns), before matching or reflection;
#'   \code{loadings_aligned}, the loadings of \code{y} that were compared,
#'   with the items in the order of \code{x} and the factors in the order
#'   of their partners, reflected where matched with a reflection,
#'   followed by any factor without a partner; \code{Phi_aligned}, the
#'   factor correlations of \code{y} in that order, \eqn{D \Phi D} with
#'   \eqn{D} the diagonal matrix of the reflections, when \code{y} carries
#'   them (absent for a matrix); \code{items}, the item names in the order
#'   of \code{x} (their positions when neither solution names its items);
#'   \code{match} and \code{loadings}, the arguments; and \code{inputs}, a
#'   list naming the kind of each input (\code{"matrix"},
#'   \code{"fa_exploratory"}, or \code{"lavaan"}) and whether the items
#'   were aligned by name or by position. A \code{measured_terms} attribute
#'   names the two columns of coefficients, so the display prints a
#'   coefficient of 1 as 1.000 beside the others (see
#'   \code{\link{dmar_tbl}}). The print method states the coefficient, the
#'   counts, and the matching rule above and below the table, and prints
#'   the full matrix, to the decimal places of the table
#'   (\code{print(x, digits = 5)} shows five in both). \code{tidy()} gives
#'   one row per row of the table, with \code{term}
#'   \code{"factor_x | factor_y"}
#'   (\code{"(none)"} standing in for the missing name of a factor without
#'   a partner, as in \code{"X3 | (none)"}), \code{estimate} (the
#'   coefficient), \code{reflected}, and \code{next_largest};
#'   \code{glance()} gives \code{n_items}, \code{n_factors_x},
#'   \code{n_factors_y}, \code{n_matched}, \code{match}, \code{loadings},
#'   and \code{mean_congruence}, the mean coefficient over the matched
#'   pairs.
#'
#' @references
#' Korth, B., & Tucker, L. R. (1975). The distribution of chance
#'   congruence coefficients from simulated data. \emph{Psychometrika,
#'   40}(3), 361--372. \doi{10.1007/BF02291763}
#'
#' Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
#'   \emph{Naval Research Logistics Quarterly, 2}(1--2), 83--97.
#'   \doi{10.1002/nav.3800020109}
#'
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
#'   coefficient as a meaningful index of factor similarity.
#'   \emph{Methodology, 2}(2), 57--64. \doi{10.1027/1614-2241.2.2.57}
#'
#' @seealso \code{\link{procrustes_phi}} for the coefficient of one pair of
#'   loading vectors with a permutation \emph{p}-value,
#'   \code{\link{fa_exploratory}} and \code{\link{fa_confirmatory}} for the
#'   solutions compared, and \code{\link{dmar_tidiers}} for \code{tidy()}
#'   and \code{glance()}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @family multivariate and latent variable methods
#'
#' @examples
#' items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'            "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
#' # Congruence compares loadings alone, so the two fits leave out the
#' # standard errors, tests, and intervals (point_estimates_only = TRUE).
#' oblimin <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
#'                           seed = 113, point_estimates_only = TRUE)
#' geomin <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
#'                          rotation = "geomin", seed = 113,
#'                          point_estimates_only = TRUE)
#'
#' # The same two factors under two rotations, compared by their pattern
#' # loadings and then by their structure loadings.
#' factor_congruence(oblimin, geomin)
#' factor_congruence(oblimin, geomin, loadings = "structure")
#'
#' # A loading matrix and a copy with its items in reverse order, its
#' # factors listed the other way around, and one factor reflected:
#' # matching recovers the pairing and the sign.
#' L <- matrix(c(0.7, 0.6, 0.8, 0.1, 0.0, 0.2,
#'               0.1, 0.2, 0.0, 0.7, 0.8, 0.6), ncol = 2,
#'             dimnames = list(paste0("item_", 1:6), c("A", "B")))
#' L_copy <- L[6:1, c("B", "A")]
#' L_copy[, "A"] <- -L_copy[, "A"]
#' matched <- factor_congruence(L, L_copy)
#' matched
#' attr(matched, "loadings_aligned")
#'
#' # Paired by position instead, the coefficients keep their signs.
#' factor_congruence(L, L_copy, match = "order")
#' generics::glance(matched)
#'
#' @export
factor_congruence <- function(x, y, match = c("congruence", "order"),
                              loadings = c("pattern", "structure")) {
  match <- match.arg(match)
  loadings <- match.arg(loadings)
  sx <- .congruence_solution(x, "x", loadings)
  sy <- .congruence_solution(y, "y", loadings)
  aligned <- .congruence_align(sx, sy)
  Lx <- aligned$Lx
  Ly <- aligned$Ly
  C <- .congruence_matrix(Lx, Ly)
  kx <- nrow(C)
  ky <- ncol(C)
  fx <- rownames(C)
  fy <- colnames(C)

  # partner[i] is the factor of y paired with factor i of x, NA for none.
  if (match == "order") {
    if (kx != ky) {
      .dmar_refuse_argument(
        sprintf(paste0(
     "match = \"order\" pairs column j of 'y' with column j of 'x' and ",
     "needs the same number of factors in each; 'x' has %d and 'y' has ",
     "%d. Use match = \"congruence\" to pair them by the coefficient."),
     kx, ky), arg = c("y", "x"))
    }
    partner <- seq_len(kx)
  } else {
    partner <- .assignment_max(abs(C))
  }
  paired <- !is.na(partner)
  phi <- rep(NA_real_, kx)
  phi[paired] <- C[cbind(which(paired), partner[paired])]
  reflected <- rep(NA, kx)
  reflected[paired] <- match == "congruence" & phi[paired] < 0
  phi[which(reflected)] <- -phi[which(reflected)]

  # The strongest competitor of each row's factor: the largest absolute
  # coefficient with a factor of the other solution other than its partner.
  competitor <- function(v, drop) {
    v <- abs(v)
    if (!is.na(drop)) v <- v[-drop]
    if (length(v)) max(v) else NA_real_
  }
  next_x <- vapply(seq_len(kx), function(i) competitor(C[i, ], partner[i]),
                   numeric(1))
  alone_y <- setdiff(seq_len(ky), partner[paired])
  next_y <- vapply(alone_y, function(j) competitor(C[, j], NA), numeric(1))

  out <- data.frame(
    factor_x = c(fx, rep(NA_character_, length(alone_y))),
    factor_y = c(fy[partner], fy[alone_y]),
    congruence = c(phi, rep(NA_real_, length(alone_y))),
    reflected = c(reflected, rep(NA, length(alone_y))),
    next_largest = c(next_x, next_y),
    stringsAsFactors = FALSE)

  # y's loadings and factor correlations in the order of their partners,
  # reflected where the match was, then the factors without a partner.
  order_y <- c(partner[paired], alone_y)
  d <- rep(1, length(order_y))
  d[seq_len(sum(paired))][reflected[paired]] <- -1
  L_aligned <- sweep(Ly[, order_y, drop = FALSE], 2L, d, `*`)
  Phi_aligned <- if (!is.null(sy$Phi)) {
    d * sy$Phi[order_y, order_y, drop = FALSE] * rep(d, each = length(d))
  }

  # The coefficients are measured values, so a coefficient of 1 prints as
  # 1.000 beside 0.995, as it does in the congruence matrix beneath the
  # table (see the measured_terms of .as_dmar_tbl()).
  out <- .as_dmar_tbl(out, subclass = "dmar_factor_congruence",
                      measured_terms = c("congruence", "next_largest"))
  attr(out, "congruence_matrix") <- C
  attr(out, "loadings_aligned") <- L_aligned
  attr(out, "Phi_aligned") <- Phi_aligned
  attr(out, "items") <- aligned$items
  attr(out, "match") <- match
  attr(out, "loadings") <- loadings
  attr(out, "inputs") <- list(x = sx$source, y = sy$source,
                              items = aligned$by)
  out
}

# The lines above and below the table: what was compared, how the factors
# were paired, and the full matrix of coefficients. A subset of the
# columns, which keeps the class but loses the attributes, prints as an
# ordinary DMAR table.
#' @exportS3Method .dmar_tbl_footer dmar_factor_congruence
#' @noRd
.dmar_tbl_footer.dmar_factor_congruence <- function(
    x, digits = getOption("dmar.digits", 3L), ...) {
  C <- attr(x, "congruence_matrix", exact = TRUE)
  inputs <- attr(x, "inputs", exact = TRUE)
  if (is.null(C) || !is.list(inputs)) return(NextMethod())
  match <- attr(x, "match", exact = TRUE)
  loadings <- attr(x, "loadings", exact = TRUE)
  kx <- nrow(C)
  ky <- ncol(C)
  p <- length(attr(x, "items", exact = TRUE))
  what <- function(source, both = FALSE) {
    metric <- if (identical(source, "matrix")) "(as supplied)" else
      sprintf("(standardized %s loadings)", loadings)
    kind <- switch(source,
                   matrix = c("a loading matrix", "loading matrices"),
                   lavaan = c("a lavaan fit", "lavaan fits"),
                   c("an fa_exploratory() table", "fa_exploratory() tables"))
    paste(kind[if (both) 2L else 1L], metric)
  }
  inputs_line <- if (identical(inputs$x, inputs$y)) {
    paste("x and y:", what(inputs$x, both = TRUE))
  } else {
    sprintf("x: %s; y: %s", what(inputs$x), what(inputs$y))
  }
  lines <- .dmar_footer_add(character(0), c(
    sprintf("Tucker's congruence coefficients: %d items; x has %d %s, y has %d",
            p, kx, if (kx == 1L) "factor" else "factors", ky),
    paste0(inputs_line, if (identical(inputs$items, "position"))
      "; items aligned by position")), "header")
  lines <- .dmar_footer_add(lines, paste0(
    "(congruence: the coefficient of each pair",
    if (match == "congruence") ", after reflection" else ", signed",
    "; reflected: whether the factor of y was reflected to match; ",
    "next_largest: the largest absolute coefficient of the row's factor ",
    "with any factor of the other solution other than its partner)"),
    "table_note")
  n_pairs <- sum(!is.na(x$congruence))
  n_reflected <- sum(x$reflected %in% TRUE)
  rule <- if (match == "congruence") {
    sprintf(paste0(
      "Matching: the one-to-one assignment of y's factors to x's that ",
      "maximizes the sum of the absolute coefficients; %d %s, %d %s."),
      n_pairs, if (n_pairs == 1L) "pair" else "pairs", n_reflected,
      if (n_reflected == 1L) "factor of y reflected" else
        "factors of y reflected")
  } else {
    paste("Matching: by position, column j of y with column j of x;",
          "the coefficients keep their signs and no factor is reflected.")
  }
  alone <- c(x$factor_x[is.na(x$factor_y)], x$factor_y[is.na(x$factor_x)])
  if (length(alone)) {
    rule <- paste0(rule, sprintf(
      " Without a partner: %s of %s.", paste(alone, collapse = ", "),
      if (kx > ky) "x" else "y"))
  }
  lines <- .dmar_footer_add(lines, rule, "matching")
  shown <- matrix(.format_dmar_value(c(C), digits, whole = FALSE), kx, ky,
                  dimnames = dimnames(C))
  lines <- .dmar_footer_add(lines, paste(
    "Congruence matrix (rows: the factors of x; columns: the factors of y;",
    "before matching and reflection):"), "congruence_matrix")
  lines <- .dmar_footer_add(
    lines, utils::capture.output(print(shown, quote = FALSE, right = TRUE)),
    "congruence_matrix_table")
  c(lines, NextMethod())
}
