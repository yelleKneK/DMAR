# Formatted correlation matrix with tests and confidence intervals
#' Formatted Correlation Matrix With \emph{p}-values and Confidence Intervals
#'
#' Computes a correlation matrix along with, for every pair of variables, the
#' two-sided \emph{p}-value, a confidence interval, and the pairwise sample size,
#' and arranges them into a single annotated table. The table is a convenience
#' for inspecting the correlations, their significance, and their intervals at a
#' glance; it is not intended as a finished, publication-ready exhibit. Output
#' formats include plain text (for the console or to paste into a Word document),
#' HTML (best for Word via browser copy-paste), and LaTeX.
#'
#' @param x A \code{data.frame}, tibble, or \code{matrix}. Non-numeric columns are
#'   dropped with no warning; the remaining numeric, integer, and logical columns
#'   are used.
#' @param method The correlation method: \code{"pearson"} (default),
#'   \code{"spearman"}, or \code{"kendall"}.
#' @param conf_level Confidence level for the interval (default \code{0.95}).
#' @param listwise Logical. If \code{TRUE}, apply listwise deletion before any
#'   correlation is computed so that every pair uses the same sample. If
#'   \code{FALSE} (the default), each pair uses all of its available complete
#'   observations (\dQuote{pairwise} deletion).
#' @param stars Logical. If \code{TRUE}, append conventional significance stars
#'   (\code{*} \emph{p} < .05, \code{**} \emph{p} < .01, \code{***} \emph{p} < .001)
#'   to each correlation and add an explanatory footnote to the table. The exact
#'   \emph{p}-value is always shown as well.
#' @param decimals_r Number of decimals for correlations and confidence interval
#'   limits (default \code{2}).
#' @param decimals_p Number of decimals for \emph{p}-values (default \code{4},
#'   matching the package-wide \code{digits_p} convention used by
#'   \code{\link{dmar_tbl}}); values below \eqn{10^{-\mathrm{decimals\_p}}} are
#'   printed as \dQuote{< .0001} (with the threshold tracking \code{decimals_p}).
#' @param format One of \code{"text"} (default), \code{"html"}, or \code{"latex"}.
#'   Controls the returned/printed table. See Details.
#' @param file Optional file path. If supplied, the formatted table is written
#'   to this file. The HTML path writes a self-contained HTML document (the
#'   kable wrapped in a minimal document head that pulls in Bootstrap CSS from
#'   a CDN) so the file opens directly in a browser without any pandoc /
#'   webshot machinery. The LaTeX path writes the raw \code{tabular} fragment;
#'   embed it in a document that loads \code{\\usepackage{makecell}} and
#'   \code{\\usepackage{booktabs}}. The text path uses
#'   \code{\link[base]{writeLines}}.
#'
#' @return An object of class \code{"correlations_test"} containing the matrices
#'   \code{r}, \code{p}, \code{ci_lower}, \code{ci_upper}, and \code{n} (all
#'   \eqn{p \times p} with variable names as row/column names), plus the arguments
#'   used. When \code{format = "html"} or \code{format = "latex"} and \code{file}
#'   is \code{NULL}, a \code{\link[knitr]{kable}}
#'   object is returned instead so that the table renders inside R Markdown /
#'   Quarto documents. When \code{format = "text"}, the table is printed to the
#'   console and the raw object is returned invisibly.
#'
#' @details
#' \strong{Layout.} Each lower-triangle cell stacks four values: the correlation
#' (with optional significance stars), the two-sided \emph{p}-value, the
#' confidence interval, and the pairwise sample size. The upper triangle is left
#' blank so the same information is not repeated.
#'
#' \strong{\emph{p}-values.} Computed with \code{\link[stats]{cor.test}} using the
#' requested \code{method}. For Spearman and Kendall with ties, \code{cor.test}
#' cannot compute an exact \emph{p}-value and falls back to a normal-approximation
#' \emph{p}-value; the associated warnings are suppressed for a cleaner table.
#'
#' \strong{Confidence intervals.} All three methods use Fisher's variance-
#' stabilizing transformation, \eqn{z(r) = \mathrm{atanh}(r)}, and back-transform
#' through \eqn{\tanh(\cdot)}, but the standard error in the Fisher-\emph{z} scale
#' is selected to match the sampling distribution of the chosen correlation
#' coefficient:
#' \itemize{
#'   \item Pearson (\code{method = "pearson"}): \eqn{\mathrm{SE}(z) =
#'     1/\sqrt{n - 3}}. This is the classical Fisher (1921) interval. Under
#'     bivariate normality the coverage of this interval matches
#'     \code{\link[stats]{cor.test}}'s \code{conf.int} exactly; for departures
#'     from bivariate normality both intervals lose coverage in the same way.
#'     See Kelley (2007) and Maxwell, Delaney, & Kelley (2027, Chapter 9) for
#'     discussion and worked examples.
#'   \item Spearman (\code{method = "spearman"}): \eqn{\mathrm{SE}(z) =
#'     \sqrt{(1 + r^{2}/2) / (n - 3)}}. This is the Bonett and Wright (2000)
#'     adjustment, which uses Fisher's transformation but inflates the standard
#'     error to account for the heavier-than-Pearson tails of the Spearman
#'     sampling distribution. This is the form Bonett and Wright recommend for
#'     practical use; a plain Fisher \eqn{1/\sqrt{n-3}} standard error tends to
#'     produce intervals that are too narrow for Spearman correlations.
#'   \item Kendall (\code{method = "kendall"}): \eqn{\mathrm{SE}(z) =
#'     \sqrt{0.437 / (n - 4)}}. This is Bonett and Wright's (2000,
#'     equation 2) Fisher-\emph{z} interval for Kendall's \eqn{\tau}. The
#'     constant 0.437 is the asymptotic variance factor of Fieller,
#'     Hartley, and Pearson (1957), derived under bivariate normality and
#'     stated by Bonett and Wright as accurate for \eqn{|\tau| < .8}
#'     (the Spearman variance above is likewise stated as accurate for
#'     \eqn{|\rho_s| < .95}).
#'     Requires \eqn{n \ge 5}; for smaller pairwise samples the interval is
#'     returned as \code{NA}.
#' }
#'
#' The Fisher-\emph{z} machinery requires \eqn{|r| < 1} for the transformation
#' to be finite. When \eqn{r = \pm 1} (perfect correlation in the sample), the
#' transformed value is infinite and the interval is reported as \code{NA}; this
#' is the same convention used by \code{cor.test}.
#'
#' \strong{When to use each correlation.} Pick \code{method = "pearson"} when
#' both variables are continuous, approximately linearly related, and roughly
#' bivariate normal (or at least without heavy tails and influential outliers).
#' Pick \code{method = "spearman"} or \code{method = "kendall"} when the
#' relationship is monotone but not necessarily linear, when one or both
#' variables are ordinal, or when influential outliers would distort Pearson's
#' \emph{r}. Kendall's \eqn{\tau} is often preferred over Spearman's \eqn{\rho}
#' for small samples and for samples with many tied ranks because it has
#' better small-sample properties and a more interpretable concordance-based
#' meaning. See Maxwell, Delaney, & Kelley (2027, Chapter 9) for an extended
#' discussion of effect size choice and interval estimation.
#'
#' \strong{HTML/LaTeX output.} Built with \code{knitr::kable} (a
#' \code{Suggests} dependency). Cell content and variable names are escaped for
#' the target format so that ``\emph{p} < .001'' renders correctly and that
#' variable names containing characters such as \code{_}, \code{&}, or \code{\%}
#' do not break LaTeX compilation. LaTeX output uses \code{\\makecell}, which
#' requires \code{\\usepackage\{makecell\}} in the document preamble.
#'
#' \strong{Pasting into Word.} The cleanest path is \code{format = "html"} with
#' a \code{file} argument; the function writes a small self-contained HTML
#' document (no pandoc dependency). Open the result in a browser and copy/paste
#' the table into Word. Formatting (including the stacked-cell layout) is
#' preserved.
#'
#' @references
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
#' estimating Pearson, Kendall and Spearman correlations.
#' \emph{Psychometrika, 65}(1), 23--28. \doi{10.1007/BF02294183}
#'
#' Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). Tests for
#' rank correlation coefficients. I. \emph{Biometrika, 44}(3/4),
#' 470--481. \doi{10.1093/biomet/44.3-4.470}
#'
#' Fisher, R. A. (1915). Frequency distribution of the values of the
#' correlation coefficient in samples from an indefinitely large population.
#' \emph{Biometrika, 10}(4), 507--521. \doi{10.1093/biomet/10.4.507}
#'
#' Fisher, R. A. (1921). On the ``probable error'' of a coefficient of
#' correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical
#' Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#' experiments and analyzing data: A model comparison perspective} (4th ed.).
#' Routledge.
#'
#' @examples
#' # Worked example using four cognitive tests from the Holzinger and
#' # Swineford (1939) study (301 children in two schools). The goal of
#' # correlations_test() is to produce a formatted correlation matrix that
#' # reports, for every variable pair, the correlation, its two-sided
#' # p-value, a confidence interval on the population correlation, and the
#' # pairwise sample size. See Kelley (2007) and Maxwell, Delaney, & Kelley
#' # (2027, Chapter 9) for discussion of why effect sizes should be
#' # accompanied by confidence intervals.
#' hs_tests <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
#'                                     "t4_lozenges",
#'                                     "t6_paragraph_comprehension")]
#'
#' # Pearson correlations (the default). Each lower-triangle cell stacks r,
#' # the two-sided p-value, the 95\% confidence interval (Fisher's Z
#' # transformation; Fisher, 1915, 1921), and the pairwise N.
#' correlations_test(hs_tests)
#'
#' # Add significance stars and an explanatory footnote.
#' correlations_test(hs_tests, stars = TRUE)
#'
#' # Spearman correlations at a 99\% confidence level. The interval uses
#' # Bonett and Wright's (2000) Fisher's Z standard error
#' # sqrt((1 + r^2/2) / (n - 3)), which corrects the plain Fisher interval
#' # for the heavier tails of Spearman's sampling distribution.
#' correlations_test(hs_tests, method = "spearman", conf_level = 0.99)
#'
#' # Kendall's tau, also using Bonett and Wright's (2000) Fisher's Z standard
#' # error sqrt(0.437 / (n - 4)). Kendall is often preferred over Spearman
#' # for small samples and for samples with many tied ranks, and these
#' # integer test scores carry many ties.
#' correlations_test(hs_tests, method = "kendall")
#'
#' # Save a formatted HTML table that opens directly in a browser
#' # (then copy into Word). No pandoc required.
#' tmp_html <- tempfile(fileext = ".html")
#' correlations_test(hs_tests, stars = TRUE, format = "html", file = tmp_html)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{descriptives}}, \code{\link[stats]{cor.test}},
#'   \code{\link[stats]{cor}}, \code{\link{ci_r}}
#'
#' @keywords htest multivariate
#'
#' @family hypothesis tests
#'
#' @export
#' @import stats

correlations_test <- function(
  x,
  method      = "pearson",
  conf_level  = 0.95,
  listwise    = FALSE,
  stars       = FALSE,
  decimals_r  = 2,
  decimals_p  = 4,
  format      = "text",
  file        = NULL
) {
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  format <- match.arg(format, c("text", "html", "latex"))

  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number strictly between 0 and 1.")
  }
  if (!is.numeric(decimals_r) || decimals_r < 0 ||
      !is.numeric(decimals_p) || decimals_p < 1) {
    stop("'decimals_r' must be >= 0 and 'decimals_p' must be >= 1.")
  }

  # Normalize input to a data.frame of numeric-like columns.
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("'x' must be a data frame, tibble, or matrix.")
  keep <- vapply(x, function(col) is.numeric(col) || is.logical(col), logical(1))
  if (sum(keep) < 2) stop("'x' must contain at least two numeric-like variables.")
  x <- x[, keep, drop = FALSE]
  x[] <- lapply(x, as.numeric)

  if (isTRUE(listwise)) {
    x <- x[stats::complete.cases(x), , drop = FALSE]
    if (nrow(x) < 4) {
      stop("Fewer than 4 complete cases remain; cannot compute confidence intervals.")
    }
  }

  var_names <- names(x)
  k <- ncol(x)

  # Allocate the matrices that hold the raw results. Off-diagonal correlations
  # start at NA, not 0: a pair with too few complete cases skips its loop body
  # (below), and diag(1, k) would otherwise leave it reported as r = 0 rather
  # than missing.
  r_mat  <- matrix(NA_real_, k, k, dimnames = list(var_names, var_names)); diag(r_mat) <- 1
  p_mat  <- matrix(NA_real_,   k, k, dimnames = list(var_names, var_names))
  lo_mat <- matrix(NA_real_,   k, k, dimnames = list(var_names, var_names))
  up_mat <- matrix(NA_real_,   k, k, dimnames = list(var_names, var_names))
  n_mat  <- matrix(NA_integer_, k, k, dimnames = list(var_names, var_names))

  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  n_min  <- .cor_ci_n_min(method)

  for (i in seq_len(k - 1)) {
    for (j in (i + 1):k) {
      xi <- x[[i]]
      xj <- x[[j]]
      both <- stats::complete.cases(xi, xj)
      n_pair <- sum(both)
      n_mat[i, j] <- n_mat[j, i] <- n_pair

      if (n_pair < 4L) next  # Need 4+ complete cases for any correlation at all.

      ct <- suppressWarnings(stats::cor.test(
        xi[both], xj[both],
        method = method
      ))

      r_val <- unname(ct$estimate)
      r_mat[i, j] <- r_mat[j, i] <- r_val
      p_mat[i, j] <- p_mat[j, i] <- ct$p.value

      # Method-specific Fisher's Z confidence interval. For Pearson this is the
      # classical Fisher (1921) interval and agrees with cor.test's
      # conf.int. For Spearman and Kendall the standard error in the Fisher's Z
      # scale follows Bonett and Wright (2000); see ?correlations_test details.
      if (is.finite(r_val) && abs(r_val) < 1 && n_pair >= n_min) {
        z  <- atanh(r_val)
        se <- .cor_ci_se(r_val, n_pair, method)
        lo_mat[i, j] <- lo_mat[j, i] <- tanh(z - z_crit * se)
        up_mat[i, j] <- up_mat[j, i] <- tanh(z + z_crit * se)
      }
    }
  }

  result <- structure(
    list(
      r          = r_mat,
      p          = p_mat,
      ci_lower   = lo_mat,
      ci_upper   = up_mat,
      n          = n_mat,
      method     = method,
      conf_level = conf_level,
      listwise   = listwise,
      stars      = stars,
      decimals_r = decimals_r,
      decimals_p = decimals_p,
      format     = format
    ),
    class = "correlations_test"
  )

  # Build the formatted table in the requested format and either save, print,
  # or return it.
  out <- switch(format,
    text  = .cor_test_text(result),
    html  = .cor_test_html(result),
    latex = .cor_test_latex(result)
  )

  if (!is.null(file)) {
    if (format == "text") {
      writeLines(out, con = file)
    } else if (format == "html") {
      .write_html_doc(out, file, title = .caption(result))
    } else {
      writeLines(as.character(out), con = file)
    }
    return(invisible(result))
  }

  # No file: return the object. Auto-printing at the REPL handles the
  # text formatting via print.correlations_test (which calls cat()), so we
  # do not eagerly cat() here. Assignment in scripts and in test code
  # therefore stays silent, which is the standard R convention.
  if (format == "text") {
    return(result)
  }

  # HTML or LaTeX without a file: return the kable so it renders in R Markdown.
  out
}


# Write a self-contained HTML document wrapping the supplied kable HTML so the
# file opens directly in a browser with reasonable styling. Avoids
# any HTML-writing helper's pandoc or webshot dependency.
.write_html_doc <- function(kable_obj, file, title = "Correlations") {
  body <- as.character(kable_obj)
  bootstrap_css <- "https://cdn.jsdelivr.net/npm/bootstrap@4.6.2/dist/css/bootstrap.min.css"
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"UTF-8\">",
    sprintf("<title>%s</title>", .escape_html(title)),
    sprintf("<link rel=\"stylesheet\" href=\"%s\">", bootstrap_css),
    "<style>body{padding:24px;font-family:serif;}</style>",
    "</head>",
    "<body>",
    body,
    "</body>",
    "</html>"
  )
  writeLines(html, con = file)
}


# ----- Internal statistical helpers ----------------------------------------

# Minimum pairwise sample size that yields a well-defined Fisher's Z CI for the
# requested correlation method. Pearson and Spearman need n >= 4 so that
# n - 3 > 0; Kendall uses (n - 4) in the denominator and so requires n >= 5
# (Bonett & Wright, 2000).
.cor_ci_n_min <- function(method) {
  if (identical(method, "kendall")) 5L else 4L
}

# Standard error of the Fisher's Z transformed correlation, on the Z scale.
# Pearson:  SE(z) = 1 / sqrt(n - 3)                        (Fisher, 1915, 1921)
# Spearman: SE(z) = sqrt((1 + r^2/2) / (n - 3))            (Bonett & Wright, 2000)
# Kendall:  SE(z) = sqrt(0.437 / (n - 4))                  (Bonett & Wright, 2000)
.cor_ci_se <- function(r, n, method) {
  switch(method,
    pearson  = 1 / sqrt(n - 3),
    spearman = sqrt((1 + r^2 / 2) / (n - 3)),
    kendall  = sqrt(0.437 / (n - 4))
  )
}


# ----- Internal formatting helpers ------------------------------------------

# HTML escape for cell content. Order matters: ampersand must be first so we
# don't re-encode the entity replacements that follow.
.escape_html <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  s <- gsub(">", "&gt;",  s, fixed = TRUE)
  s
}

# Minimal LaTeX text-mode escape covering the special characters that can
# appear in either p-value strings ("<", from "p < .001") or in user-supplied
# variable names (e.g., "_", "&", "%"). The backslash is handled via a
# sentinel: its replacement "\textbackslash{}" contains braces, so escaping it
# literally in the first step would let the later "{" / "}" substitutions
# re-escape those braces into "\{\}", typesetting a stray "{}" after the
# backslash. We stash the backslash as a placeholder that cannot appear in the
# input, escape the braces of the genuine "{"/"}" characters, then expand the
# placeholder to "\textbackslash{}" last.
.escape_latex <- function(s) {
  bslash <- "\001BACKSLASH\001"
  s <- gsub("\\", bslash,               s, fixed = TRUE)
  s <- gsub("{",  "\\{",                s, fixed = TRUE)
  s <- gsub("}",  "\\}",                s, fixed = TRUE)
  s <- gsub("&",  "\\&",                s, fixed = TRUE)
  s <- gsub("%",  "\\%",                s, fixed = TRUE)
  s <- gsub("$",  "\\$",                s, fixed = TRUE)
  s <- gsub("#",  "\\#",                s, fixed = TRUE)
  s <- gsub("_",  "\\_",                s, fixed = TRUE)
  s <- gsub("~",  "\\textasciitilde{}", s, fixed = TRUE)
  s <- gsub("^",  "\\textasciicircum{}",s, fixed = TRUE)
  s <- gsub("<",  "\\textless{}",       s, fixed = TRUE)
  s <- gsub(">",  "\\textgreater{}",    s, fixed = TRUE)
  s <- gsub(bslash, "\\textbackslash{}", s, fixed = TRUE)
  s
}

# Format a correlation / CI-limit value: fixed number of decimals, APA-style
# leading-zero stripped for |x| < 1.
.format_r <- function(x, decimals) {
  if (is.na(x) || !is.finite(x)) return("")
  s <- sprintf(paste0("%.", decimals, "f"), x)
  sub("^(-?)0\\.", "\\1.", s)
}

# Format a p-value for display, with "< .001"-style substitution for values
# below the precision threshold.
.format_p_display <- function(p, decimals) {
  if (is.na(p)) return("")
  threshold <- 10^(-decimals)
  fmt <- paste0("%.", decimals, "f")
  if (p < threshold) {
    return(paste0("p < ", sub("^0\\.", ".", sprintf(fmt, threshold))))
  }
  paste0("p = ", sub("^0\\.", ".", sprintf(fmt, p)))
}

.stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Return the four display lines (r, p, CI, N) for cell (i, j). Upper-triangle
# and diagonal cells are handled explicitly.
.cell_lines <- function(x, i, j) {
  if (j > i)  return(c("", "", "", ""))
  if (j == i) return(c("-", "", "", ""))  # diagonal marker

  r_val <- x$r[i, j]
  p_val <- x$p[i, j]

  r_str <- .format_r(r_val, x$decimals_r)
  if (isTRUE(x$stars)) r_str <- paste0(r_str, .stars(p_val))

  p_str <- .format_p_display(p_val, x$decimals_p)

  lo <- x$ci_lower[i, j]
  up <- x$ci_upper[i, j]
  ci_str <- if (is.na(lo) || is.na(up)) "" else {
    sprintf("[%s, %s]",
            .format_r(lo, x$decimals_r),
            .format_r(up, x$decimals_r))
  }

  n_val <- x$n[i, j]
  n_str <- if (is.na(n_val)) "" else sprintf("N = %d", n_val)

  c(r_str, p_str, ci_str, n_str)
}

.caption <- function(x) {
  sprintf("Correlations (%s, %d%% CI)",
          tools::toTitleCase(x$method),
          round(x$conf_level * 100))
}

# Text format: fixed-width grid, four physical rows per variable-row.
.cor_test_text <- function(x) {
  var_names <- rownames(x$r)
  k <- length(var_names)

  label_width <- max(nchar(var_names), 12L) + 2L
  col_width   <- max(14L, max(nchar(var_names)) + 4L)

  pad <- function(s, w) formatC(s, width = w, flag = "-")

  header <- pad("", label_width)
  for (j in seq_len(k)) header <- paste0(header, pad(var_names[j], col_width))

  lines <- c(
    .caption(x),
    "",
    header,
    strrep("-", nchar(header))
  )

  for (i in seq_len(k)) {
    r_row  <- pad(var_names[i], label_width)
    p_row  <- pad("",            label_width)
    ci_row <- pad("",            label_width)
    n_row  <- pad("",            label_width)

    for (j in seq_len(k)) {
      cl <- .cell_lines(x, i, j)
      r_row  <- paste0(r_row,  pad(cl[1], col_width))
      p_row  <- paste0(p_row,  pad(cl[2], col_width))
      ci_row <- paste0(ci_row, pad(cl[3], col_width))
      n_row  <- paste0(n_row,  pad(cl[4], col_width))
    }

    lines <- c(lines, r_row, p_row, ci_row, n_row, "")
  }

  if (isTRUE(x$stars)) {
    lines <- c(lines, "Note. * p < .05, ** p < .01, *** p < .001.")
  }

  lines
}

# HTML format via knitr::kable alone. Each cell is a single HTML string with <br>
# line breaks. Cell content and variable names are HTML-escaped before the
# <br> separators are added so that text like "p < .001" survives the raw
# (escape = FALSE) insertion into the table.
.cor_test_html <- function(x) {
  var_names <- rownames(x$r)
  k <- length(var_names)

  cell_mat <- matrix("", nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      cl <- .cell_lines(x, i, j)
      cl <- cl[nzchar(cl)]
      cl <- vapply(cl, .escape_html, character(1))
      cell_mat[i, j] <- paste(cl, collapse = "<br>")
    }
  }

  escaped_names <- vapply(var_names, .escape_html, character(1))
  df <- data.frame(Variable = escaped_names, cell_mat,
                   stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- c("", escaped_names)

  tbl <- knitr::kable(df, format = "html", escape = FALSE, row.names = FALSE,
                      caption = .caption(x), align = c("l", rep("l", k)),
                      table.attr = paste0(
                        'class="table table-striped table-bordered ',
                        'table-condensed" ',
                        'style="width: auto; margin-left: 0;"'))

  if (isTRUE(x$stars)) {
    # The note travels inside the table markup as a caption-styled footer
    # row spanning every column, so the output remains one self-contained
    # knitr::kable object with no styling package required.
    note <- paste0(
      "<tr><td colspan=\"", k + 1L, "\" style=\"border: none; ",
      "font-size: 90%; padding-top: 4px;\"><em>Note.</em> ",
      "* p &lt; .05; ** p &lt; .01; *** p &lt; .001.</td></tr>")
    tbl <- sub("</tbody>", paste0(note, "</tbody>"), tbl, fixed = TRUE)
    tbl <- structure(tbl, format = "html", class = "knitr_kable")
  }

  tbl
}

# LaTeX format via knitr::kable alone. Each cell is wrapped in \makecell so line
# breaks work inside a tabular cell. Cell content and variable names are
# LaTeX-text-mode escaped before being placed in the table; this protects
# against silent corruption from "<" in p-value strings (which would
# otherwise typeset as the Spanish inverted exclamation), from "_" or "&"
# in user-supplied variable names, and similar.
.cor_test_latex <- function(x) {
  var_names <- rownames(x$r)
  k <- length(var_names)

  cell_mat <- matrix("", nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      cl <- .cell_lines(x, i, j)
      cl <- cl[nzchar(cl)]
      if (length(cl) == 0L) {
        cell_mat[i, j] <- ""
      } else {
        cl <- vapply(cl, .escape_latex, character(1))
        cell_mat[i, j] <- sprintf("\\makecell[l]{%s}",
                                  paste(cl, collapse = " \\\\ "))
      }
    }
  }

  escaped_names <- vapply(var_names, .escape_latex, character(1))
  df <- data.frame(Variable = escaped_names, cell_mat,
                   stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- c("", escaped_names)

  tbl <- knitr::kable(df, format = "latex", escape = FALSE, row.names = FALSE,
                      caption = .caption(x), booktabs = TRUE,
                      align = c("l", rep("l", k)))

  if (isTRUE(x$stars)) {
    # The note is appended after the tabular as ordinary LaTeX, keeping
    # the output a plain knitr::kable object with no styling package.
    note <- paste0(
      "\n\\vspace{2pt}\n{\\footnotesize \\emph{Note.} ",
      "* p \\textless{} .05; ** p \\textless{} .01; ",
      "*** p \\textless{} .001.}\n")
    tbl <- structure(paste0(tbl, note), format = "latex",
                     class = "knitr_kable")
  }

  tbl
}


#' @export
print.correlations_test <- function(x, ...) {
  fmt <- if (!is.null(x$format)) x$format else "text"
  if (fmt == "text") {
    cat(.cor_test_text(x), sep = "\n")
  } else if (fmt == "html") {
    print(.cor_test_html(x), ...)
  } else if (fmt == "latex") {
    print(.cor_test_latex(x), ...)
  }
  invisible(x)
}
