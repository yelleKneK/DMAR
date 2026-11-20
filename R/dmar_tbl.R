#' Printing for DMAR Result Tables
#'
#' Most \pkg{DMAR} estimation and testing functions return a tidy
#' \code{data.frame} with a \code{term} column and one or more numeric
#' columns. Because a single numeric column often holds quantities on
#' very different scales (for example, whole-number degrees of freedom
#' alongside an \emph{F} statistic, an effect size, and a small
#' \emph{p}-value), the base \code{\link[base]{print.data.frame}}
#' method formats the whole column with one common format and is easily
#' pushed into scientific notation with many trailing digits. The
#' \code{dmar_tbl} class supplies \code{print} and \code{format}
#' methods that format each value on its own terms: whole numbers (such
#' as degrees of freedom and sample sizes) print without a decimal
#' part, other values print to a fixed number of decimal places, and no
#' value is shown in scientific notation.
#'
#' The stored numeric values are never rounded; only their display
#' changes, so downstream arithmetic on the returned object (confidence
#' interval widths, further calculations) uses full precision.
#'
#' The same formatting applies to every \code{dmar_tbl}, whether the
#' table is long (a \code{term} column beside a single \code{value} or
#' \code{estimate} column) or wide (a leading label column such as
#' \code{term}, \code{effect}, or \code{sample_type} beside several
#' typed numeric columns). Each numeric column is formatted on its own
#' terms, so the shape of the table does not matter.
#'
#' Display precision is a number of decimal places, three by default,
#' set by the \code{digits} argument or, globally, by
#' \code{options(dmar.digits = )}. (The option is dot-named because that
#' is the R convention for package options, for example
#' \code{dplyr.width} and \code{knitr.table.format}; it is not a function
#' or argument name and so is outside the package's snake_case rule.) A
#' value that is not a whole number prints to that many decimal places
#' with its trailing zeros, so every such value in a column shows the same
#' number of decimal places: 13.018 prints as 13.018, 0.56 as 0.560,
#' 233.4107 as 233.411, and 1.5 as 1.500. A value too small for those
#' decimals prints as 0.000 rather than in scientific notation (0.000123
#' and 3.35e-52 both print as 0.000), and a value that rounds to zero
#' prints as 0.000 whatever its sign, never as -0.000.
#'
#' A whole number prints without a decimal part, as a count (degrees of
#' freedom, a sample size, a number of replications) does. A measured
#' quantity can land on a whole number as well, at the bound of its range
#' or where it is fixed, and it then prints to the decimal places of the
#' values of its kind beside it. The display recognizes a measured value
#' in four ways.
#' \itemize{
#'   \item \emph{A bound and the estimate it belongs to.} A row or column
#'     named \code{ci_lower} or \code{ci_upper}, or with a name ending
#'     in \code{_lower} or \code{_upper}, is a bound: a confidence limit
#'     (\code{tau2_ci_lower}, \code{rmsea_ci_upper},
#'     \code{cles_ci_lower}, \code{smd_ci_lower}), and in the same way an
#'     equivalence bound (\code{smd_lower}, \code{contrast_upper}) or a
#'     one-sided test statistic (\code{t_upper}). The estimate an interval belongs to is
#'     the row or column its limits are named for (\code{tau2},
#'     \code{rmsea}), the row of a long table that its \code{ci_lower}
#'     and \code{ci_upper} rows belong to (the estimate
#'     \code{\link{results_sentence}} reports), the column
#'     \code{\link[generics]{tidy}} takes as the estimate of a wide table
#'     with limit columns, and a row whose own limits sit beside it in
#'     limit columns (the \code{ci_lower} and \code{ci_upper} of a
#'     coefficient table). So a lower limit of exactly 0 reads 0.000
#'     beside an upper limit of 0.083, and a loading fixed at 1 reads
#'     1.000 beside its interval.
#'   \item \emph{A measured kind, wherever it appears.} A fit index (CFI,
#'     TLI, RMSEA, SRMR, and the others, with their differences between
#'     models); a power; a proportion, probability, or fraction (the
#'     \code{pct_} rows of a sensitivity analysis, a coverage, an
#'     assurance, a fraction of missing information); a mean, median,
#'     standard deviation, variance, or higher moment (the \code{mean_},
#'     \code{median_}, and \code{sd_} rows of a sensitivity analysis); a
#'     test statistic, a standard error, or a noncentrality parameter; the
#'     delta and the weights of a model comparison; and the heterogeneity
#'     indices of a meta-analysis. So a CFI of 1 reads 1.000 beside a TLI
#'     of 1.012, a power of 1 reads 1.000 beside a power of 0.840, and a
#'     standard deviation of the interval widths of 0 reads 0.000 beside
#'     their mean of 1.176. A name that also names a count, such as
#'     \code{delta_df} or \code{median_N}, is a count.
#'   \item \emph{A column of measured values.} A column that holds values
#'     that are not whole is measured throughout, so an \emph{F} of 4
#'     reads 4.000 beside 4.270 and a standard error of 0 reads 0.000
#'     beside 0.050. The exception is a column named for a count, with
#'     \code{n}, \code{N}, \code{df}, or \code{B} as its name or as a part
#'     of its name (\code{df_2}, \code{N_effective}), which counts its
#'     whole values even when others are not whole (degrees of freedom
#'     adjusted for sphericity read 2 beside 1.830), and a row that records
#'     the table's sum-of-squares type. The \code{value} or
#'     \code{estimate} column beside a \code{term} column is not read this
#'     way, since it can hold a different quantity on each row.
#'   \item \emph{A row or column the producing function names.} A function
#'     names the rows or columns that hold measured values of its own kind
#'     through a \code{measured_terms} attribute (the congruence
#'     coefficients of \code{\link{factor_congruence}}, the loadings and
#'     intercepts of \code{\link{dmacs}}).
#' }
#' Elsewhere a whole number prints without a decimal part. That covers the
#' counts, the codes and indicators (the sum-of-squares type, whether an
#' equivalence test concluded equivalence), and a value the call was given
#' as a whole number, since nothing in a table says which of its other
#' whole numbers are counts: in a long table the \code{value} column holds
#' a sample size of 301 on one row and, on another, a planning value given
#' as 15, which prints as 15.
#'
#' \emph{p}-values are shown to a number of decimal places of their own
#' (four by default, set by \code{digits_p}), which is the conventional
#' way to report them. A
#' \emph{p}-value smaller than the smallest magnitude those decimals
#' can represent prints as \dQuote{< 0.0001} instead of rounding to
#' \code{0.0000}. A column is treated as holding \emph{p}-values when
#' it is named \code{p_value}, \code{p.value}, \code{p_adjusted} (the
#' multiplicity-adjusted case, as in \code{\link{ci_dunnett}}), or
#' \code{p_chi_square} (the exact-fit test of a fitted model, as in
#' \code{\link{measurement_invariance}}); in a long-format
#' table whose quantities share a single \code{value} column, the rows
#' to format this way are named by the producing function through a
#' \code{p_terms} attribute. The limits of an interval for a
#' \emph{p}-value, rows named for it with \code{_ci_lower} and
#' \code{_ci_upper} (the Monte Carlo interval of the bootstrap
#' \emph{p}-value of \code{\link{sem_exact_fit_bootstrap}}), print to the
#' same decimal places, without the floor, since a limit of 0 is a limit
#' and not a \emph{p}-value below the smallest that can be shown.
#'
#' Information criteria such as AIC and BIC, and log-likelihoods, keep a
#' number of decimal places of their own, \code{digits_fixed} (three by
#' default), which \code{digits} does not change, so a model comparison
#' difference of a few points is not rounded away when the rest of a
#' table is shown to fewer decimals; such a value prints its decimals
#' even when it is a whole number. The producing function names these
#' rows through a \code{fixed_terms} attribute.
#'
#' A table that carries a \code{bootstrap} attribute (the reliability
#' family's bootstrap intervals; see \code{\link{reliability_alpha}})
#' prints one more line stating how many replications were attempted,
#' how many returned a value, and why the rest were set aside, so the
#' replications the interval conditions on are always on the page.
#'
#' A table that carries a \code{denominator} attribute (the omega
#' coefficients of \code{\link{reliability_omega}},
#' \code{\link{reliability_omega_categorical}}, and
#' \code{\link{bifactor_indices}}) prints a line naming how the total
#' variance in the denominator was estimated, and a table that carries a
#' \code{residual_cor} attribute prints the largest residual correlation
#' of the single-factor fit with the item pair it belongs to (or, for a
#' three-item fit, says that the saturated model leaves no residual
#' correlation to read), so the coefficient on the page is identified
#' without reference to a default and the local dependence beside it is
#' in view. The \pkg{knitr} methods (\code{\link{knit_print.dmar_tbl}}
#' and \code{\link{as_kable}}) write the same lines with a knitted
#' table, as they do the lines a result adds of its own (for example
#' which version of a fit measure \code{\link{sem_fit_measures}}
#' reported and why), so a report shows what the console shows.
#'
#' A sample size plan made for incomplete data (\code{\link{ss_aipe_reliability}}
#' with \code{fmi} above 0) prints a line naming the fraction of missing
#' information it assumed and the sample size complete data would need,
#' and a sensitivity analysis that planned for or simulated missing data
#' (\code{\link{ss_aipe_reliability_sensitivity}} with \code{fmi} or
#' \code{missing_prop} above 0) prints a line naming both, so the
#' assumptions behind the number stay beside it. In the same way, an
#' assurance plan of \code{\link{ss_aipe_c}} made with a nonzero
#' \code{true_kurtosis} prints a line naming the excess kurtosis of the errors
#' its assurance allows for.
#'
#' To see more precision than the display shows, raise \code{digits}
#' (\code{print(x, digits = 8)} shows eight decimal places) or read the
#' columns directly, since the stored values are never rounded:
#' \code{x$value} or \code{x[["p_value"]]} returns the numbers at full
#' precision.
#'
#' @section Using the result in your own code: You do not need to know
#'   anything about S3 classes to use a \code{dmar_tbl}. It is an
#'   ordinary \code{data.frame} with a print method, so everything you
#'   already do with a data frame works: \code{x$value} pulls the
#'   numeric column, \code{x[x$term == "smd", ]} selects a row, and the
#'   full-precision numbers are right there for any further calculation.
#'   Three common needs:
#'   \itemize{
#'     \item \emph{Read one number.} Index it like any data frame, for
#'       example \code{x$value[x$term == "ci_upper"]}. The display
#'       rounds; the stored value does not, so this returns the number
#'       at full precision.
#'     \item \emph{See more (or fewer) decimal places.} Use
#'       \code{print(x, digits = 6)} for a single table, or
#'       \code{options(dmar.digits = 6)} for the rest of the session.
#'     \item \emph{Hand the result to other tools.} \code{tidy(x)}
#'       returns a one-row-per-term table in DMAR's own column
#'       vocabulary (\code{term}, \code{estimate}, \code{se},
#'       \code{statistic}, \code{p_value}, \code{ci_lower},
#'       \code{ci_upper}, and so on), which is the convenient
#'       \dQuote{wide} view for plotting or joining; \code{glance(x)}
#'       returns a one-row model-level summary. Both come from the
#'       \pkg{generics} package and need no extra setup. For a
#'       single-estimand result such as \code{\link{ci_smd}} the table is
#'       already one row, so \code{glance()} coincides with \code{tidy()}
#'       (there are no extra model-level statistics to report); for a
#'       multi-row result such as \code{\link{mlmr}} they differ.
#'   }
#'
#' @section Stacking tables with rbind(): A \code{dmar_tbl} carries
#'   facts about its numbers as attributes: the confidence level, which
#'   rows are \emph{p}-values, and whatever else its function recorded
#'   (\code{\link{dmar_metadata}} lists them). The lines printed beneath a
#'   table are written from those attributes, so a stacked table can keep
#'   them only when they are true of every row. \code{rbind()} of DMAR
#'   tables therefore returns
#'   \itemize{
#'     \item a \code{dmar_tbl} of the same class and with the same
#'       attributes, when every table stacked is a \code{dmar_tbl} and
#'       all of them have identical classes and attributes (two calls
#'       of one function at one confidence level, for example);
#'     \item an ordinary \code{data.frame} with the same rows and
#'       columns and none of the DMAR attributes, in every other case:
#'       tables at different confidence levels, tables from different
#'       functions, or a DMAR table stacked with a plain data frame. No
#'       warning is given. The numbers are all there; what is dropped is
#'       the class and the attributes, so that no printed line states a
#'       confidence level (or anything else) that does not hold for every
#'       row.
#'   }
#'   To compare results across confidence levels, stack the
#'   \code{tidy()} rows instead, which carry the level in a
#'   \code{conf_level} column:
#'   \code{rbind(tidy(x_90), tidy(x_95))}. \code{NULL} arguments are
#'   ignored, as \code{rbind()} ignores them. \code{dplyr::bind_rows()}
#'   does not call this method, so what is said here does not apply to it;
#'   use \code{rbind()}, or stack \code{tidy()} rows.
#'
#' @param x A \code{dmar_tbl} object (a \code{data.frame} returned
#'   by a DMAR function).
#' @param digits Number of decimal places for values that are not whole
#'   numbers, and for measured values that are (see Details): a single
#'   whole number of 0 or more. Defaults to
#'   \code{getOption("dmar.digits", 3L)}. A class whose display lays out
#'   cells of its own, or writes numbers in the lines beneath the table,
#'   uses the same \code{digits} there.
#' @param digits_p Number of decimal places for \emph{p}-values.
#'   Defaults to 4. A \emph{p}-value below \code{10^(-digits_p)} prints
#'   as \dQuote{< 0.0001} (with the threshold tracking \code{digits_p}).
#' @param digits_fixed Number of decimal places for \code{fixed_terms}
#'   rows (information criteria such as AIC and BIC, and log-likelihoods).
#'   Defaults to 3.
#' @param ... For \code{print()}, additional arguments passed to
#'   \code{\link[base]{print.data.frame}}. For \code{rbind()}, the tables
#'   to stack, and any of the named arguments of
#'   \code{\link[base]{rbind.data.frame}} (\code{make.row.names},
#'   \code{stringsAsFactors}, \code{factor.exclude}).
#' @param deparse.level For \code{rbind()}: passed to
#'   \code{\link[base]{rbind.data.frame}}, which does not use it for data
#'   frames. It keeps its base R spelling because \code{rbind()} supplies
#'   it by that name.
#'
#' @return \code{print.dmar_tbl} returns \code{x} invisibly.
#'   \code{format.dmar_tbl} returns a \code{data.frame} whose numeric
#'   columns have been formatted to character for display.
#'   \code{rbind()} of DMAR tables returns a \code{dmar_tbl} when every
#'   table stacked has the same classes and attributes, and an ordinary
#'   \code{data.frame} otherwise (see the section on stacking).
#'
#' @seealso \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} for the wide one-row-per-term and the
#'   one-row summary views. For a gentle, non-technical tour of how to
#'   read and use DMAR result tables, see the \dQuote{Reading DMAR result
#'   tables} vignette: \code{vignette("dmar_output", package = "DMAR")}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Every DMAR estimation function returns a table that prints this way.
#' x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
#' x                       # rounded for reading, to three decimal places
#'
#' # The stored numbers keep full precision; only the display rounds.
#' x$value[x$term == "smd"]
#' print(x, digits = 8)    # ask the display for more decimal places
#'
#' # Pull a single number out, exactly as you would from a data frame.
#' x$value[x$term == "ci_upper"]
#'
#' # The broom verbs give the programmer-friendly wide and summary views.
#' generics::tidy(x)
#' generics::glance(x)
#'
#' # The same display rules apply to a wide table, one with several typed
#' # columns, for example an effect size with its confidence interval per
#' # effect.
#' ci_eta_squared(aov(iq_8 ~ treatment, data = pygmalion))
#'
#' # Stacking two tables computed at different confidence levels gives an
#' # ordinary data frame, because no one confidence level describes every
#' # row. The tidy() rows keep each level in a column.
#' x_90 <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.90)
#' class(rbind(x, x_90))
#' rbind(generics::tidy(x), generics::tidy(x_90))
#'
#' @keywords internal
#' @name dmar_tbl
NULL

# Element-wise display formatter shared by format.dmar_tbl / print.dmar_tbl
# and by the displays that lay out cells of their own (fa_exploratory(),
# factor_count(), factor_congruence(), factor_scores()). A whole number
# prints without a decimal part where `whole` is TRUE, as a count (degrees
# of freedom, a sample size) reads; every other finite value prints to
# `digits` decimal places with its trailing zeros (13.018 reads "13.018",
# 0.56 reads "0.560"), never in scientific notation (0.000123 and 3.35e-52
# read "0.000"), and a value that rounds to zero prints without a minus
# sign ("0.000", not "-0.000"). `whole` is recycled along `v`, so a caller
# can type each element: FALSE for a measured quantity that can land on a
# whole number (a lower limit of 0, a CFI of 1), which then prints as
# "0.000" or "1.000" beside the other values of its kind. What counts as
# whole is .dmar_is_whole(). NA and NaN return NA, as a missing value; Inf
# and -Inf print as such. Through 1.0.x `digits` counted significant
# figures and format() chose fixed or scientific notation per value. Not
# exported.
.format_dmar_value <- function(v, digits = 3L, whole = TRUE) {
  digits <- .dmar_decimals(digits)
  v <- as.double(v)
  whole <- rep_len(as.logical(whole), length(v))
  out <- rep(NA_character_, length(v))
  infinite <- is.infinite(v)
  out[infinite] <- ifelse(v[infinite] > 0, "Inf", "-Inf")
  finite <- is.finite(v)
  counted <- finite & whole & .dmar_is_whole(v)
  out[counted] <- vapply(round(v[counted]), format, character(1L),
                         scientific = FALSE, trim = TRUE)
  measured <- finite & !counted
  out[measured] <- .dmar_fixed_text(v[measured], digits)
  out
}

# Whether each value is a whole number, up to the rounding error of the
# arithmetic that produced it: a count computed in floating point
# (22.000000000000004, 1 - 1e-16) is whole, and a value 1e-9 away from a
# whole number is not. The tolerance is relative to the value, so it
# applies only near a nonzero whole number: zero is whole only when it is
# exactly zero, and a tiny nonzero value (a relative weight of 3.35e-52, a
# difference of -6.66e-16) is a measured value that rounds to 0.000, which
# through 1.0.x printed as itself (3.35e-52) rather than as 0. NA for NA
# and NaN. Not exported.
.dmar_is_whole <- function(v) {
  r <- round(v)
  v == r | (r != 0 & abs(v - r) <= 8 * .Machine$double.eps * abs(v))
}

# `digits` checked as a number of decimal places: a single whole number of
# 0 or more. The display layer's `digits`, and options(dmar.digits = ),
# count decimal places from 1.1.0 on. Not exported.
.dmar_decimals <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != round(digits)) {
    .dmar_refuse_argument(c(
      "'digits' must be a single whole number of 0 or more, the number ",
      "of decimal places shown (the default is ",
      "getOption(\"dmar.digits\", 3))."), arg = "digits")
  }
  as.integer(digits)
}

# The names DMAR gives a bound: the limits of an interval, ci_lower and
# ci_upper, and a name ending in _lower or _upper (tau2_ci_lower,
# rmsea_ci_upper, cles_ci_lower, smd_ci_lower). The same names mark the
# other bounds DMAR reports, which are measured values too: the
# equivalence bounds (smd_lower and smd_upper of equivalence_smd(),
# contrast_lower and contrast_upper of equivalence_c(), cor_lower and
# cor_upper of equivalence_r()) and the one-sided test statistics t_lower
# and t_upper. Not exported.
.dmar_limit_name <- function(nm) {
  nm %in% c("ci_lower", "ci_upper") |
    grepl("._(lower|upper)$", nm)
}

# The estimate a limit named for it belongs to: "tau2" for tau2_ci_lower,
# "rmsea" for rmsea_ci_upper, "cles" for cles_ci_lower; "" for the
# generic names (ci_lower, ci_upper), whose
# estimate the table identifies otherwise. Not exported.
.dmar_limit_stem <- function(nm) {
  stem <- sub("_(ci_)?(lower|upper)$", "", nm)
  stem[stem == nm | nm %in% c("ci_lower", "ci_upper")] <- ""
  stem
}

# The names of counts: a sample size (n or N alone or as a part of the
# name: n_1, n_per_group, necessary_N, N_effective, n1), degrees of freedom
# (df alone or as a part: df_2, df_error, delta_df, df_unbounded), and
# the number of bootstrap replications (B, B_used). A column so named
# keeps the whole-number rule even when some of its values are not whole
# (degrees of freedom adjusted for sphericity, an effective sample size),
# and a measured kind (.dmar_measured_name()) never includes one, so the
# delta_df of measurement_invariance() is a count. Not exported.
.dmar_count_name <- function(nm) {
  grepl("(^|_)(n|N|df)(_|$)|^(n|N)[0-9]+$|^df|^B(_|$)", as.character(nm))
}

# The names of quantities that are measured wherever a DMAR table shows
# them, so that a whole value among them prints to the decimals of the
# rest:
#   - the fit indices of a structural equation or factor model and their
#     differences between models (a CFI of 1, an RMSEA of 0, delta_cfi);
#   - a power (power, actual_power, composite_power, power_a);
#   - a proportion, probability, or fraction (the pct_ rows and
#     total_type_I_error of a sensitivity analysis, the
#     pct_ci_less_width_ rows of a composite plan, a coverage, an
#     assurance, the fraction of missing information, the probabilities
#     p_tied, p_y1_greater, and p_y1_less of cliff_delta());
#   - a mean, median, standard deviation, variance, or higher moment (the
#     mean_, median_, and sd_ rows of a sensitivity analysis, an sd_ci_width
#     of 0 among them, and the moments of a noncentral distribution);
#   - a test statistic or a standard error (F_value, t_value, z_value,
#     t, z, F, chi_square, statistic, se, and names built on se such as
#     regression_se), and a noncentrality parameter (ncp, noncentrality,
#     noncentral_);
#   - the delta, the relative weights, and the log evidence of a model
#     comparison (compare_sem(), factor_count());
#   - and the heterogeneity indices of a meta-analysis (tau, tau2, I2, H2,
#     Q).
# A name that is also a count (.dmar_count_name()) is not one of them, so
# median_N and delta_df are counts. Not exported.
.dmar_measured_name <- function(nm) {
  nm <- as.character(nm)
  kind <- grepl(paste0(
    "(^|_)(cfi|tli|nnfi|rni|nfi|ifi|rfi|gfi|agfi|pgfi|pnfi|mfi|rmsea|srmr|",
    "crmr|gamma_hat)(_|$)",
    "|(^|_)power(_|$)",
    "|^(pct|prop|prob|proportion)(_|$)|(^|_)(coverage|assurance)(_|$)",
    "|^(total_type_I_error|fmi)$",
    "|^p_(tied|y1_greater|y1_less)$",
    "|^(mean|median|sd|var|variance|skewness|kurtosis|excess_kurtosis)(_|$)",
    "|^true_kurtosis$",
    "|^(statistic|test_statistic|statistic_value|chi_square|delta_chi_square|t|z|F)$",
    "|^(t|z|F|q)_(value|statistic)$|(^|_)se(_|$)",
    "|^(ncp|noncentrality)$|^noncentral_",
    "|^(delta|relative_weight|criterion_delta|criterion_weight)$",
    "|(^|_)log_evidence",
    "|^(tau|tau2|I2|H2|Q)$"), nm)
  kind & !.dmar_count_name(nm)
}

# Which cells of each numeric column hold a measured value whose whole
# number prints to the decimals of the rest, rather than as a count
# (degrees of freedom, a sample size), which is whole. Nothing in a table
# says in general which of its whole numbers are counts, so a cell is
# measured when the table shows it is, in four ways, read the same way in
# both shapes of table.
#   - A bound and the estimate it belongs to. Columns: the limits
#     (.dmar_limit_name()), the column a limit is named for
#     (eta_squared_partial for eta_squared_partial_ci_lower), and, beside
#     generic limit columns, the column tidy() takes as the estimate
#     (`estimate` or `value`, else the first numeric column that is not a
#     limit, standard error, test statistic, degrees of freedom, p-value,
#     or count). Rows of a long table's primary column (`primary`, the
#     value or estimate column beside `term`): the limit rows that are not
#     p-values (the p_lower and p_upper of an equivalence test are
#     p_terms), the row a limit is named for (tau2 for tau2_ci_lower), the
#     estimate of ci_lower and ci_upper rows (the row
#     results_sentence() and tidy() take, from the ci_estimate attribute
#     or else the first row that is not a limit), and a row whose own
#     limits sit in limit columns beside it (a coefficient table's ci_lower
#     and ci_upper).
#   - A measured kind (.dmar_measured_name()), a row or a column.
#   - A column, other than a long table's primary column, that holds a
#     value that is not whole (.dmar_is_whole()), unless it is named for a
#     count (.dmar_count_name()): an F of 4 beside 4.270, a standard error
#     of 0 beside 0.050, a relative weight of 1 beside 0.000. A row that
#     records the table's sum-of-squares type (the foot of
#     manova_split_plot()'s statistic_value column) holds a code, and its
#     cells are never measured.
#   - The rows or columns the producing function names in its
#     measured_terms attribute (factor_congruence()'s coefficients).
# Returns a named list of logical vectors, one per numeric column of `x`
# (a plain data frame that keeps the table's attributes); `num` flags its
# numeric columns. Not exported.
.dmar_measured_cells <- function(x, num, primary) {
  nms <- names(x)[num]
  n <- nrow(x)
  long <- !is.na(primary) && "term" %in% names(x)
  own <- attr(x, "measured_terms", exact = TRUE)
  limit_cols <- .results_sentence_limit_cols(nms)
  lim <- .dmar_limit_name(nms)
  est <- nms %in% setdiff(.dmar_limit_stem(nms[lim]), "")
  if (!is.null(limit_cols)) {
    est_col <- intersect(c("estimate", "value"), nms)
    if (!length(est_col)) {
      est_col <- nms[!grepl(.results_sentence_not_estimate, nms)]
    }
    est <- est | nms %in% est_col[1L]
  }
  typed <- vapply(nms, function(nm) {
    v <- x[[nm]]
    any(is.finite(v) & !.dmar_is_whole(v))
  }, logical(1L)) & !.dmar_count_name(nms)
  measured <- lim | est | typed | .dmar_measured_name(nms) | nms %in% own
  labels <- names(x)[!num]
  record <- rep(FALSE, n)
  for (nm in labels) {
    record <- record | as.character(x[[nm]]) %in% "sum_of_squares_type"
  }
  out <- stats::setNames(lapply(seq_along(nms), function(j) {
    rep(measured[j] && !(long && identical(nms[j], primary)), n) & !record
  }), nms)
  if (long) {
    term <- as.character(x$term)
    p_rows <- term %in% attr(x, "p_terms", exact = TRUE)
    rows_lim <- .dmar_limit_name(term) & !p_rows
    rows_est <- term %in% setdiff(.dmar_limit_stem(term[rows_lim]), "")
    if (all(c("ci_lower", "ci_upper") %in% term)) {
      rows_est <- rows_est | term %in% .results_sentence_long_estimate(x)
    }
    if (!is.null(limit_cols)) {
      rows_est <- rows_est | !is.na(x[[limit_cols[1L]]]) |
        !is.na(x[[limit_cols[2L]]])
    }
    out[[primary]] <- rows_lim | rows_est | .dmar_measured_name(term) |
      term %in% own
  }
  out
}

# Fixed-decimal p-value formatter shared by format.dmar_tbl. p-values read
# best to a fixed number of decimal places of their own (the convention is
# four), so they never slip into scientific notation. A value below the
# smallest representable magnitude (10^(-digits_p)) prints as "< 0.0001"
# rather than rounding to "0.0000". Not exported.
.format_dmar_pvalue <- function(p, digits_p = 4L) {
  thresh <- 10^(-digits_p)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = digits_p))
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < thresh) return(floor_label)
    formatC(x, format = "f", digits = digits_p)
  }, character(1L))
}

# Fixed-decimal formatter for the rows named by a table's fixed_terms
# attribute, which keep a number of decimal places of their own
# (digits_fixed), whole numbers included, whatever `digits` is. Information
# criteria (AIC, BIC) and log-likelihoods are the motivating case: a model
# comparison difference of a few points must survive a display of the rest
# of the table to fewer decimals. Through 1.0.x, when the rest of the table
# printed to 3 significant figures, an AIC of 2284.830 would have read
# "2280". Unlike the p-value formatter there is no floor and no scientific
# fallback; a value that rounds to zero prints without a minus sign. Not
# exported; consumed via the fixed_terms attribute by format.dmar_tbl.
.format_dmar_fixed <- function(v, digits_fixed = 3L) {
  out <- rep(NA_character_, length(v))
  infinite <- is.infinite(v)
  out[infinite] <- ifelse(v[infinite] > 0, "Inf", "-Inf")
  finite <- is.finite(v)
  out[finite] <- .dmar_fixed_text(v[finite], digits_fixed)
  out
}

# Numbers `v` written to `digits` decimal places with their trailing zeros
# (formatC(), which follows getOption("OutDec")), never in scientific
# notation, with the minus sign dropped from a value that rounds to zero:
# formatC() writes -0.0004 as "-0.000", which reads as a quantity below
# zero. The test is on the digits written, so it holds for any decimal
# mark. An infinite value is written "Inf" or "-Inf", and a missing one
# "NA" (formatC() pads them to the width of the others). Not exported.
.dmar_fixed_text <- function(v, digits) {
  shown <- formatC(v, format = "f", digits = digits)
  finite <- is.finite(v)
  shown[!finite] <- ifelse(is.na(v[!finite]), "NA",
                           ifelse(v[!finite] > 0, "Inf", "-Inf"))
  zero <- finite & startsWith(shown, "-") & !grepl("[1-9]", shown)
  shown[zero] <- substring(shown[zero], 2L)
  shown
}

# Tag a result table as a dmar_tbl so it prints via print.dmar_tbl
# (counts without a decimal part, other values to `digits` decimal places,
# p-values to digits_p decimals). The rule is uniform:
# every tidy-returning DMAR function routes its result through this helper,
# whether the table is long (a `term` column plus a single `value` /
# `estimate` column) or wide (a leading label column such as `term`,
# `effect`, or `sample_type` plus one or more numeric columns). The format
# method works on every numeric column, so the shape of the table does not
# matter.
#
# `conf_level`, when supplied, is recorded for the print footer, and
# `tail_areas`, the unequal tail areas of an interval formed from
# alpha_lower and alpha_upper (.tail_areas()), beside it. `p_terms`,
# when supplied, names the rows of the long-format primary column that hold
# p-values (a wide table needs no p_terms for a column named p_value,
# p.value, p_adjusted, or p_chi_square, as the format method finds those by
# name). `fixed_terms`, when supplied, names
# rows of the primary column to print to digits_fixed decimal places
# whatever `digits` is (information criteria such as AIC / BIC).
# `measured_terms`, when supplied, names the rows of the primary column or
# the columns that hold measured values of the function's own kind, whose
# whole values print to `digits` decimal places rather than as counts (see
# .dmar_measured_cells()).
# The dmar_tbl class is inserted just before `data.frame`, so any leading
# subclass (e.g., dmar_ci_long, dmar_ss_power) keeps dispatching its own
# tidy() / glance() methods while print falls through to print.dmar_tbl.
# Idempotent. Not exported; this is the one-line opt-in every tidy-returning
# DMAR function calls on its result just before returning it.
.as_dmar_tbl <- function(out, conf_level = NULL, p_terms = NULL,
                         fixed_terms = NULL, ss_type = NULL,
                         subclass = NULL, tail_areas = NULL,
                         measured_terms = NULL) {
  if (!is.null(conf_level))  attr(out, "conf_level")  <- conf_level
  if (!is.null(tail_areas))  attr(out, "tail_areas")  <- tail_areas
  if (!is.null(p_terms))     attr(out, "p_terms")     <- p_terms
  if (!is.null(fixed_terms)) attr(out, "fixed_terms") <- fixed_terms
  if (!is.null(ss_type))     attr(out, "ss_type")     <- ss_type
  if (!is.null(measured_terms)) attr(out, "measured_terms") <- measured_terms
  cls <- class(out)
  if (!("dmar_tbl" %in% cls))
    class(out) <- c(setdiff(cls, "data.frame"), "dmar_tbl", "data.frame")
  # A leading dispatch subclass (e.g. dmar_ss_aipe, dmar_ss_power) rides in
  # front of dmar_tbl so family tidiers dispatch while printing falls
  # through to print.dmar_tbl; idempotent like the dmar_tbl insertion.
  if (!is.null(subclass) && !(subclass %in% class(out)))
    class(out) <- c(subclass, class(out))
  out
}

#' @rdname dmar_tbl
#' @export
format.dmar_tbl <- function(x, digits = getOption("dmar.digits", 3L),
                            digits_p = 4L, digits_fixed = 3L, ...) {
  raw <- x
  class(raw) <- "data.frame"
  out <- raw
  num <- vapply(out, is.numeric, logical(1L))

  # The "primary" numeric column carries the per-term quantities in a long
  # table: `value`, or `estimate` for coefficient-style tables. The p_terms
  # and fixed_terms attributes name rows of this column to reformat.
  primary <- intersect(c("value", "estimate"), names(out)[num])
  primary <- if (length(primary)) primary[1L] else NA_character_

  # Every value to `digits` decimal places, a whole number as a count,
  # except that a measured value keeps its decimals when whole (see
  # .dmar_measured_cells()).
  measured <- .dmar_measured_cells(raw, num, primary)
  for (nm in names(out)[num]) {
    out[[nm]] <- .format_dmar_value(raw[[nm]], digits = digits,
                                    whole = !measured[[nm]])
  }

  # p-values get their own fixed decimals (digits_p), overriding the format
  # above. They are detected two ways: a dedicated column (p_value /
  # p.value / p_adjusted / p_chi_square) in a wide table, or named rows of
  # the primary column in a long table (the producing function lists them in
  # a p_terms attribute). The names are package-wide conventions, so a
  # function that reports an exact-fit p-value as p_chi_square, as
  # measurement_invariance() does, formats it the same way as the p_value
  # column beside it.
  for (nm in intersect(c("p_value", "p.value", "p_adjusted", "p_chi_square"),
                       names(out)[num]))
    out[[nm]] <- .format_dmar_pvalue(raw[[nm]], digits_p = digits_p)
  p_terms <- attr(x, "p_terms")
  if (!is.null(p_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% p_terms
    out[[primary]][idx] <- .format_dmar_pvalue(raw[[primary]][idx],
                                               digits_p = digits_p)
    # The limits of an interval for a p-value, rows named for it with
    # _ci_lower and _ci_upper (the Wilson interval of the bootstrap p-value
    # of sem_exact_fit_bootstrap()), on the p-value's own decimals, so a
    # p-value of 0.0200 is not given an interval of 0.004 to 0.105. There
    # is no floor: a limit of 0 is a limit, not a p-value too small to
    # show. Before, they printed to `digits` like any other row.
    idx <- raw$term %in% paste0(rep(p_terms, each = 2L),
                                c("_ci_lower", "_ci_upper"))
    out[[primary]][idx] <- .format_dmar_fixed(raw[[primary]][idx],
                                              digits_fixed = digits_p)
  }

  # fixed_terms get their own number of decimal places (digits_fixed),
  # whole numbers included, whatever `digits` is: information criteria such
  # as AIC / BIC, where a model comparison difference of a few points must
  # not round away when the rest of the table is shown to fewer decimals.
  fixed_terms <- attr(x, "fixed_terms")
  if (!is.null(fixed_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% fixed_terms
    out[[primary]][idx] <- .format_dmar_fixed(raw[[primary]][idx],
                                              digits_fixed = digits_fixed)
  }
  out
}

#' @rdname dmar_tbl
#' @export
print.dmar_tbl <- function(x, digits = getOption("dmar.digits", 3L),
                           digits_p = 4L, digits_fixed = 3L, ...) {
  # `digits` reaches the display and the lines beneath the table as well,
  # so a class that writes numbers there (the congruence matrix beneath a
  # factor_congruence() table) shows them to the decimal places asked for.
  .dmar_tbl_write(x, digits = digits, digits_p = digits_p,
                  digits_fixed = digits_fixed, print_args = list(...),
                  view_args = list(digits = digits))
}

# rbind() of DMAR tables. rbind.data.frame() gives the stack the class and
# the attributes of its first data frame, so without this method two
# ci_smd() tables at 90 and 95 percent stack into a table of class
# dmar_ci_smd whose footer says 90 percent for all six rows. The method
# stacks with rbind.data.frame() and then keeps the class and attributes
# only when every table stacked is a dmar_tbl and all agree in every
# attribute but names and row.names (class included); otherwise it returns
# a plain data.frame with those attributes removed. NULL arguments are
# set aside, as rbind() sets them aside, and so are the named arguments
# of rbind.data.frame() that arrive through the dots. No condition is
# signaled.
#
# R dispatches rbind() on all of its arguments: this method is called when
# it is the only rbind method found among them (a data frame argument with
# no method of its own class counts as rbind.data.frame only when no other
# method is found), which is what makes rbind(<dmar_tbl>, <data.frame>) and
# rbind(<data.frame>, <dmar_tbl>) arrive here. Inside R/, a function that
# builds its table by stacking does so on plain data frames and classes
# the result last with .as_dmar_tbl(); stacking a plain data frame onto a
# table that is already classed would return a plain data frame.
#' @rdname dmar_tbl
#' @export
rbind.dmar_tbl <- function(..., deparse.level = 1) {
  out <- rbind.data.frame(..., deparse.level = deparse.level)
  args <- list(...)
  nms <- names(args)
  if (!is.null(nms))
    args <- args[!(nms %in% c("make.row.names", "stringsAsFactors",
                              "factor.exclude"))]
  args <- args[!vapply(args, is.null, logical(1))]
  facts <- function(a) {
    at <- attributes(a)
    at <- at[setdiff(names(at), c("names", "row.names"))]
    at[order(names(at))]
  }
  first <- facts(args[[1L]])
  same <- all(vapply(args, function(a)
    inherits(a, "dmar_tbl") && identical(facts(a), first), logical(1)))
  if (same) return(out)
  for (k in setdiff(names(attributes(out)), c("names", "row.names")))
    attr(out, k) <- NULL
  class(out) <- "data.frame"
  out
}

# The console engine of print.dmar_tbl() and of the print methods that add
# arguments of their own (print.dmar_factor_scores()'s `n`). It prints the
# lines that introduce the table, the table every display shows
# (.dmar_tbl_display()), and the lines beneath it (.dmar_tbl_footer()), a
# blank line before each block; both generics live in
# R/dmar_output_helpers.R and are the ones knit_print.dmar_tbl() and
# as_kable.dmar_tbl() read, so the console and a knitted report show the
# same table and the same lines. `view_args` go to both generics (a
# class's own display arguments, and the `digits` of print(), so that a
# display that lays out cells of its own and a footer that writes numbers
# follow it); `print_args` go to print.data.frame(),
# where row.names and right default to FALSE for the look but can be
# overridden, e.g. print(x, row.names = TRUE). Returns x invisibly.
.dmar_tbl_write <- function(x, digits = getOption("dmar.digits", 3L),
                            digits_p = 4L, digits_fixed = 3L,
                            print_args = list(), view_args = list()) {
  shown <- do.call(.dmar_tbl_display, c(list(x), view_args))
  lines <- do.call(.dmar_tbl_footer, c(list(x), view_args))
  disp <- format(shown, digits = digits, digits_p = digits_p,
                 digits_fixed = digits_fixed)
  if (is.null(print_args[["row.names"]])) print_args[["row.names"]] <- FALSE
  if (is.null(print_args[["right"]]))     print_args[["right"]]     <- FALSE
  parts <- .dmar_footer_parts(lines)
  if (length(parts$header)) {
    cat(paste0(unname(parts$header), "\n"), "\n", sep = "")
  }
  do.call(print.data.frame, c(list(disp), print_args))
  .dmar_footer_cat(parts$body)
  invisible(x)
}
