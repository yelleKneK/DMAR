# Omega squared, total and partial (effect size for ANOVA).
#' Omega Squared (Effect Size for ANOVA)
#'
#' Computes the sample omega squared (\eqn{\omega^2}), Hays' (1994)
#' bias-corrected estimator of the proportion of the population variance
#' of the dependent variable that a fixed effect accounts for, as two
#' estimands: the \emph{total} \eqn{\omega^2}, a proportion of the total
#' variance, and the \emph{partial} \eqn{\omega^2}, a proportion of the
#' variance due to the effect and error alone. With
#' \code{estimand = "all"}, the default, the function reports every
#' estimand its inputs determine; \code{"total"} or \code{"partial"}
#' reports one. In a one-way design the two are the same number.
#'
#' The function takes a fitted model (\code{\link[stats]{aov}},
#' \code{\link[stats]{lm}}, or an \code{aovlist} from \code{aov()} with an
#' \code{Error()} term), the sums of squares of a published ANOVA table,
#' or an \emph{F} statistic with its degrees of freedom and \emph{N}. The
#' confidence interval is computed by \code{\link{ci_omega_squared}}.
#'
#' @param object Optional. A fitted model of class \code{\link[stats]{aov}},
#'   \code{\link[stats]{lm}}, or \code{aovlist} (a multi-stratum fit such
#'   as \code{aov(y ~ A + Error(subject/A), data = d)}). The function
#'   returns one row per effect, with the stratum named for an
#'   \code{aovlist}. Give either \code{object} or summary statistics, not
#'   both.
#' @param estimand Which estimand to report: \code{"all"} (the default,
#'   every estimand the inputs determine, \code{NA} for the others),
#'   \code{"total"}, or \code{"partial"}. An estimand asked for by name
#'   that the inputs do not determine is refused (see Details).
#' @param SS_effect The sum of squares of the effect, or a vector of the
#'   sums of squares of several effects of one analysis (named for the
#'   effects, which label the rows). Given with \code{df_effect} and the
#'   error mean square.
#' @param SS_error The error (residual) sum of squares, given with
#'   \code{df_error} as an alternative to \code{MS_error}.
#' @param SS_total The corrected total sum of squares,
#'   \eqn{\sum (y - \bar{y})^2} (see the Details of
#'   \code{\link{eta_squared}}). With the sums of squares, it determines
#'   the total \eqn{\omega^2}.
#' @param MS_error The error (residual) mean square. Give either
#'   \code{MS_error} or \code{SS_error} with \code{df_error}.
#' @param F_value The observed \emph{F} statistic of the effect, or a
#'   vector of the \emph{F} statistics of several effects tested against
#'   the same error term. With \code{df_effect}, \code{df_error}, and
#'   \code{N}, it determines the partial \eqn{\omega^2} only.
#' @param df_effect The numerator degrees of freedom of the effect (one
#'   value for each effect, or one for all). Required with \code{F_value}
#'   and with the sums of squares.
#' @param df_error The error (residual) degrees of freedom. Required with
#'   \code{F_value} and with \code{SS_error}.
#' @param N The total sample size, the number of observations. With the
#'   sums of squares or \code{F_value}, it determines the partial
#'   \eqn{\omega^2}. For an \emph{F} from a design with an \code{Error()}
#'   stratum, give the size of the effect's stratum instead (see Details):
#'   the number of subjects for an effect tested between subjects, and for
#'   an effect tested within subjects the stratum's total degrees of
#'   freedom, which is \code{df_effect + df_error} when the effect is the
#'   only one tested in its stratum. \code{N} must be at least
#'   \code{df_effect + df_error}.
#'
#' @return A \code{data.frame} (class \code{dmar_anova_es}, printed through
#'   \code{\link{dmar_tbl}}) with one row per effect and the columns
#'   \code{effect}, \code{omega_squared_total} and
#'   \code{omega_squared_partial} (only the one asked for when
#'   \code{estimand} is \code{"total"} or \code{"partial"}),
#'   \code{stratum} (for an \code{aovlist} fit), \code{F_value},
#'   \code{df_effect}, \code{df_error}, and \code{N}. An effect given
#'   through the summary interfaces without names is \code{"overall"} (or
#'   \code{effect_1}, \code{effect_2}, and so on, for several), and a
#'   quantity the inputs do not give is \code{NA}. \code{N} is the number
#'   of observations (for the summary interfaces, \code{N} as given); the
#'   partial value of an effect of an \code{aovlist} fit is computed with
#'   the size of its own stratum (see Details), so that its row's \emph{F},
#'   degrees of freedom, and that size, given to the \emph{F} interface,
#'   reproduce it. The \code{estimand} attribute records the estimand asked
#'   for, and the \code{estimand_notes} attribute names each estimand
#'   column that is \code{NA} with the reason, which the line beneath the
#'   printed table states. \code{\link[generics]{tidy}()} gives one row per
#'   effect and estimand (see \code{\link{dmar_tidiers}}).
#'
#' @details
#' \strong{Total and partial.} For a fixed effect with sum of squares
#' \eqn{\mathit{SS}_{\text{effect}}} on \eqn{df_{\text{effect}}} degrees of
#' freedom, error mean square \eqn{\mathit{MS}_{\text{error}}}, corrected
#' total sum of squares \eqn{\mathit{SS}_{\text{total}}}, and \emph{N}
#' observations,
#' \deqn{\hat{\omega}^2_{\text{total}} = \frac{\mathit{SS}_{\text{effect}} - df_{\text{effect}} \mathit{MS}_{\text{error}}}{\mathit{SS}_{\text{total}} + \mathit{MS}_{\text{error}}},}
#' \deqn{\hat{\omega}^2_{\text{partial}} = \frac{\mathit{SS}_{\text{effect}} - df_{\text{effect}} \mathit{MS}_{\text{error}}}{\mathit{SS}_{\text{effect}} + (N - df_{\text{effect}}) \mathit{MS}_{\text{error}}} = \frac{df_{\text{effect}} (F - 1)}{df_{\text{effect}} (F - 1) + N}}
#' (Olejnik & Algina, 2000, p. 266 and Table 10, p. 268; Maxwell,
#' Delaney, & Kelley, 2027, Eqs. 3.97, 3.98, and 7.44 to 7.46). The total
#' value estimates the effect's share of the total population variance,
#' all other effects and the error included; the partial value estimates
#' its share of the variance due to the effect and the error alone, the
#' variance that remains when the other effects are held constant. Which
#' answers the question depends on whether the design's other factors
#' would vary where the effect is studied (see the Details of
#' \code{\link{eta_squared}}). In a one-way design
#' \eqn{N - df_{\text{effect}} = df_{\text{error}} + 1}, so the two
#' denominators agree and the estimands are the same number. Some sources
#' print the total's denominator as
#' \eqn{\mathit{SS}_{\text{total}} - \mathit{MS}_{\text{error}}} (Kim,
#' 2016, Table 4, which gives .326 for an effect whose total
#' \eqn{\hat{\omega}^2} is .323); the denominator of Hays' estimator adds
#' the error mean square, as here.
#'
#' \strong{Bias correction and negative estimates.} The sample
#' \eqn{\eta^2} is the proportion of \emph{sample} variance accounted for
#' and overestimates its population value. \eqn{\hat{\omega}^2} subtracts
#' \eqn{df_{\text{effect}} \mathit{MS}_{\text{error}}} from the effect's
#' sum of squares and adds to the denominator, which gives an estimator of
#' the population proportion with much smaller bias (Hays, 1994; Olejnik &
#' Algina, 2000; Kelley, 2007). Its numerator is negative whenever
#' \emph{F} < 1, and a negative estimate of either estimand is reported
#' as 0, since \eqn{\omega^2 \ge 0} by definition (Olejnik & Algina, 2000,
#' p. 262; Maxwell, Delaney, & Kelley, 2027, Section 3.9.2).
#'
#' \strong{What the inputs determine.} A fitted model with one error
#' stratum determines both estimands of every effect. Sums of squares
#' determine the total when \code{SS_total} is given and the partial when
#' \code{N} is given, with \code{df_effect} and the error mean square in
#' either case. An \emph{F} statistic with its degrees of freedom and
#' \emph{N} determines the partial value only: the total also depends on
#' the sums of squares of the design's other effects, and the function
#' cannot tell from \emph{F} whether there are any, so the \emph{F}
#' interface reports the total as \code{NA} even for a one-way design,
#' where it equals the partial value. With \code{estimand = "all"} an
#' estimand the inputs do not determine is \code{NA}, with the reason in
#' the \code{estimand_notes} attribute and beneath the printed table;
#' asked for by name, it is refused with an error of class
#' \code{dmar_estimand_undetermined} (see \code{\link{dmar_metadata}}).
#' Sums of squares and \emph{F} are two ways of giving the same effect,
#' and a call that gives both, or a fitted model with summary statistics,
#' or both \code{MS_error} and \code{SS_error}, is refused.
#'
#' \strong{Designs.} The total \eqn{\omega^2} is defined for fixed effects
#' in designs whose factors are all between subjects, with one error term
#' (Olejnik & Algina, 2000, p. 266). For an \code{aovlist} fit, which has
#' an error term for each stratum, the total is \code{NA}; designs with
#' random factors and designs with factors within subjects have
#' \eqn{\omega^2} formulas of their own (Olejnik & Algina, 2000, Tables
#' 11 and 12; Olejnik & Algina, 2003), which this function does not
#' compute. An interaction or a covariate of an \code{lm} fit is an
#' effect like any other, with its own sum of squares and degrees of
#' freedom. The sums of squares of a fitted model are the sequential (Type
#' I) ones of \code{anova()}; see the Details of \code{\link{eta_squared}}
#' on the corrected total, on Types II and III, and on the coding of
#' factors.
#'
#' \strong{Within-subjects designs.} An effect of an \code{aovlist} fit is
#' tested in its own error stratum, and its partial estimate is computed
#' with the size of that stratum in place of \emph{N}, the size with which
#' \code{\link{ci_eta_squared}} and \code{\link{ci_omega_squared}} convert
#' the limits of its interval: the stratum's total degrees of freedom for
#' an effect tested within subjects, and the number of subjects for an
#' effect in the between-subjects stratum of a mixed design (see the
#' Details of \code{\link{ci_eta_squared}}). Since
#' \eqn{df_{\text{effect}} (F - 1)} estimates the noncentrality parameter
#' \eqn{\lambda} of the effect's test, \eqn{\hat{\omega}^2_{\text{partial}}}
#' then estimates \eqn{\lambda/(\lambda + M)}, with \emph{M} the stratum's
#' size, which is the value the sample partial \eqn{\eta^2} converges to as
#' subjects are added and the quantity whose interval
#' \code{ci_eta_squared(estimand = "partial")} reports for the row.
#' Computed with the number of observations \emph{N}, the estimate of a
#' within-subjects effect estimates \eqn{\lambda/(\lambda + N)}, a smaller
#' quantity, which the sample partial \eqn{\eta^2} does not estimate. That
#' estimate is the partial \eqn{\omega^2} that Olejnik and Algina (2003, p. 443) compute
#' for a repeated measures factor, whose \emph{N} is the total number of
#' observations, and the \emph{F} interface given the number of rows as
#' \code{N} reproduces it. For the three IQ occasions of
#' \code{\link{pygmalion}}, \emph{F}(2, 618) with 930 rows gives .150,
#' below the partial \eqn{\eta^2} of .212 and its interval [.157, .265];
#' with the stratum's 2 + 618 = 620, the estimate is .209.
#'
#' @references
#' Bargmann, R. E. (1970). Interpretation and use of a generalized
#'   discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
#'   Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), \emph{Essays in
#'   probability and statistics} (pp. 35--60). University of North
#'   Carolina Press.
#'
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'   \doi{10.1177/001316448004000309}
#'
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Keppel, G. (1991). \emph{Design and analysis: A researcher's handbook}
#'   (3rd ed.). Prentice Hall.
#'
#' Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
#'   size calculation 3. Comparison of several means using one-way ANOVA.
#'   \emph{Restorative Dentistry & Endodontics, 41}(3), 231--234.
#'   \doi{10.5395/rde.2016.41.3.231}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\omega^2} and Chapter 7
#'   on total and partial \eqn{\omega^2} in factorial designs.)
#'
#' Olejnik, S., & Algina, J. (2000). Measures of effect size for comparative
#'   studies: Applications, interpretations, and limitations.
#'   \emph{Contemporary Educational Psychology, 25}(3), 241--286.
#'   \doi{10.1006/ceps.2000.1040}
#'
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
#'   statistics: Measures of effect size for some common research designs.
#'   \emph{Psychological Methods, 8}(4), 434--447.
#'   \doi{10.1037/1082-989X.8.4.434}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and
#'   contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Venables, W. (1975). Calculation of confidence intervals for
#'   noncentrality parameters. \emph{Journal of the Royal Statistical
#'   Society, Series B, 37}(3), 406--412.
#'   \doi{10.1111/j.2517-6161.1975.tb01554.x}
#'
#' @seealso \code{\link{ci_omega_squared}} for the confidence interval,
#'   \code{\link{eta_squared}}, \code{\link{var_omega_squared}},
#'   \code{\link{ci_pvaf}}
#'
#' @examples
#' # 1. The sums of squares of a published two-way table: Kim (2016,
#' #    Table 3), light by material, 2 by 4 with 10 per cell (N = 80).
#' omega_squared(SS_effect = c(LIGHT = 34.716, MATERIAL = 1996.998,
#'                             "LIGHT:MATERIAL" = 1570.389),
#'               df_effect = c(1, 3, 3), SS_error = 2261.612,
#'               df_error = 72, SS_total = 5863.715, N = 80)
#'
#' # 2. An F statistic determines the partial value only. Bargmann's
#' #    (1970) example, used in turn by Venables (1975), Fleishman (1980),
#' #    and Steiger (2004): a one-way ANOVA with five groups of 11 and an
#' #    observed F of 11.221 on 4 and 50 df. In a one-way design the
#' #    partial value is the total.
#' omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
#'               estimand = "partial")
#'
#' # 3. One-way ANOVA from a fitted model (depression_bdi: three
#' #    treatment arms, 10 per arm, N = 30); the two estimands agree.
#' fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
#' omega_squared(fit_one)
#'
#' # 4. Two-factor ANOVA for the manipulated expectancy treatment and the
#' #    measured grade classification (pygmalion data, unequal cell
#' #    sizes, N = 310). The treatment by grade interaction is weak here
#' #    (F = 1.19), so the additive model is used.
#' fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
#' omega_squared(fit_additive)
#'
#' # 5. The same partial values with their confidence intervals.
#' ci_omega_squared(fit_additive, estimand = "partial")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

omega_squared <- function(
  object    = NULL,
  estimand  = c("all", "total", "partial"),
  SS_effect = NULL,
  SS_error  = NULL,
  SS_total  = NULL,
  MS_error  = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL,
  N         = NULL
) {
  estimand <- .anova_es_estimand(estimand)
  spec <- .anova_es_spec("omega_squared", object = object,
                         SS_effect = SS_effect, SS_error = SS_error,
                         SS_total = SS_total, MS_error = MS_error,
                         F_value = F_value, df_effect = df_effect,
                         df_error = df_error, N = N)
  .anova_es_table("omega_squared", "omega_squared", spec,
                  .anova_es_omega(spec), estimand)
}
