# Eta squared, total and partial (effect size for ANOVA).
#' Eta Squared (Effect Size for ANOVA)
#'
#' Computes the sample eta squared (\eqn{\eta^2}) of each effect of an
#' analysis of variance, the proportion of the variance of the dependent
#' variable that the effect accounts for, as two estimands: the
#' \emph{total} \eqn{\eta^2}, the effect's sum of squares as a proportion of
#' the corrected total sum of squares, and the \emph{partial}
#' \eqn{\eta^2}, its sum of squares as a proportion of its own sum of
#' squares and the error sum of squares. With \code{estimand = "all"}, the
#' default, the function reports every estimand its inputs determine;
#' \code{"total"} or \code{"partial"} reports one. In a one-way design the
#' two are the same number.
#'
#' The function takes a fitted model (\code{\link[stats]{aov}},
#' \code{\link[stats]{lm}}, or an \code{aovlist} from \code{aov()} with an
#' \code{Error()} term), the sums of squares of a published ANOVA table,
#' or an \emph{F} statistic with its degrees of freedom. The confidence
#' interval is computed by \code{\link{ci_eta_squared}}.
#'
#' @param object Optional. A fitted model of class \code{\link[stats]{aov}},
#'   \code{\link[stats]{lm}}, or \code{aovlist} (a multi-stratum fit such
#'   as \code{aov(y ~ A + Error(subject/A), data = d)}). The function
#'   returns one row per effect (every row of \code{anova(object)}, or of
#'   each stratum of \code{summary(object)}, other than \code{Residuals}),
#'   with the stratum named for an \code{aovlist}. Give either
#'   \code{object} or summary statistics, not both.
#' @param estimand Which estimand to report: \code{"all"} (the default,
#'   every estimand the inputs determine, \code{NA} for the others),
#'   \code{"total"}, or \code{"partial"}. An estimand asked for by name
#'   that the inputs do not determine is refused (see Details).
#' @param SS_effect The sum of squares of the effect, or a vector of the
#'   sums of squares of several effects of one analysis (named for the
#'   effects, which label the rows).
#' @param SS_error The error (residual) sum of squares. With
#'   \code{SS_effect}, it determines the partial \eqn{\eta^2}.
#' @param SS_total The corrected total sum of squares,
#'   \eqn{\sum (y - \bar{y})^2} (see Details). With \code{SS_effect}, it
#'   determines the total \eqn{\eta^2}.
#' @param F_value The observed \emph{F} statistic of the effect, or a
#'   vector of the \emph{F} statistics of several effects tested against
#'   the same error term. With \code{df_effect} and \code{df_error}, it
#'   determines the partial \eqn{\eta^2} only.
#' @param df_effect The numerator degrees of freedom of the effect (one
#'   value for each effect, or one for all). Required with
#'   \code{F_value}; with sums of squares it is optional and, with
#'   \code{df_error} and \code{SS_error}, fills the \code{F_value} column.
#' @param df_error The error (residual) degrees of freedom. Required with
#'   \code{F_value}; optional with sums of squares.
#'
#' @return A \code{data.frame} (class \code{dmar_anova_es}, printed through
#'   \code{\link{dmar_tbl}}) with one row per effect and the columns
#'   \code{effect}, \code{eta_squared_total} and \code{eta_squared_partial}
#'   (only the one asked for when \code{estimand} is \code{"total"} or
#'   \code{"partial"}), \code{stratum} (for an \code{aovlist} fit),
#'   \code{F_value}, \code{df_effect}, and \code{df_error}. An effect given
#'   through the summary interfaces without names is \code{"overall"} (or
#'   \code{effect_1}, \code{effect_2}, and so on, for several), and a
#'   quantity the inputs do not give (the \emph{F} statistic of sums of
#'   squares given without degrees of freedom) is \code{NA}. The
#'   \code{estimand} attribute records the estimand asked for, and the
#'   \code{estimand_notes} attribute names each estimand column that is
#'   \code{NA} with the reason, which the line beneath the printed table
#'   states. \code{\link[generics]{tidy}()} gives one row per effect and
#'   estimand (see \code{\link{dmar_tidiers}}).
#'
#' @details
#' \strong{Total and partial.} For an effect with sum of squares
#' \eqn{\mathit{SS}_{\text{effect}}}, error sum of squares
#' \eqn{\mathit{SS}_{\text{error}}}, and corrected total sum of squares
#' \eqn{\mathit{SS}_{\text{total}}},
#' \deqn{\hat{\eta}^2_{\text{total}} = \frac{\mathit{SS}_{\text{effect}}}{\mathit{SS}_{\text{total}}}, \qquad
#'   \hat{\eta}^2_{\text{partial}} = \frac{\mathit{SS}_{\text{effect}}}{\mathit{SS}_{\text{effect}} + \mathit{SS}_{\text{error}}} = \frac{df_{\text{effect}} F}{df_{\text{effect}} F + df_{\text{error}}}}
#' (Olejnik & Algina, 2000, Table 10, p. 268; Maxwell, Delaney, & Kelley,
#' 2027, Eqs. 7.40 to 7.42). The total value is a proportion of all the
#' variation in the data, and so depends on the other factors of the
#' design: adding a factor with a large effect lowers it. The partial
#' value leaves the other effects out of its denominator, so the partial
#' values of the effects of one design do not share a denominator and can
#' sum to more than 1 (Olejnik & Algina, 2000, pp. 268--269). Maxwell,
#' Delaney, and Kelley (2027, Chapter 7) frame the choice by the other
#' factors: when another factor is extrinsic to the question (a study of
#' the effect could as well hold it constant, as a study of drugs could
#' hold biofeedback constant) the partial value answers it, and when the
#' factor is intrinsic (it varies wherever the effect is studied, as gender
#' does) the total value does. In a one-way design \eqn{\mathit{SS}_{\text{total}} =
#' \mathit{SS}_{\text{effect}} + \mathit{SS}_{\text{error}}}, and the two
#' estimands are the same number.
#'
#' \strong{What the inputs determine.} A fitted model determines both
#' estimands of every effect, and so do sums of squares that include both
#' \code{SS_error} and \code{SS_total}. \code{SS_effect} with
#' \code{SS_total} alone determines the total, and with \code{SS_error}
#' alone the partial. An \emph{F} statistic and its degrees of freedom
#' determine the partial value only: the total of an effect also depends
#' on the sums of squares of the design's other effects (from the
#' \emph{F} statistics of every effect it can be written as
#' \eqn{df_{\text{effect}} F / (\sum df F + df_{\text{error}})}, Maxwell,
#' Delaney, & Kelley, 2027, Eq. 7.43, which a fitted model provides). The
#' function cannot tell from \emph{F} whether the design has other
#' effects, so the \emph{F} interface reports the total as \code{NA} even
#' for a one-way design, where it equals the partial value. With
#' \code{estimand = "all"} an estimand the inputs do not determine is
#' \code{NA}, with the reason in the \code{estimand_notes} attribute and
#' beneath the printed table; asked for by name, it is refused with an
#' error of class \code{dmar_estimand_undetermined} (see
#' \code{\link{dmar_metadata}}). Sums of squares and \emph{F} are two ways
#' of giving the same effect, and a call that gives both, or a fitted
#' model with summary statistics, is refused, as is a set of sums of
#' squares in which \code{SS_effect + SS_error} exceeds \code{SS_total}.
#'
#' \strong{The corrected total.} \code{SS_total} is the corrected total
#' sum of squares, \eqn{\sum (y - \bar{y})^2}, the sum of every sum of
#' squares of an ANOVA table with an intercept, the error included. Some
#' programs also print an uncorrected total, \eqn{\sum y^2}, which
#' includes the grand mean: Kim (2016, Table 3) lists both, a
#' \dQuote{Total} of 53758.357 and a \dQuote{Corrected total} of 5863.715,
#' and the second is \code{SS_total}. A fitted model without an intercept
#' does not partition the corrected total, and its total \eqn{\eta^2} is
#' \code{NA}.
#'
#' \strong{Sums of squares of a fitted model.} The sums of squares and
#' \emph{F} statistics of a fitted model are those of \code{anova()} (one
#' stratum) or \code{summary()} (an \code{aovlist}), which are Type I
#' (sequential): each effect is adjusted for the effects before it in the
#' formula. Sequential sums of squares add up to the corrected total, so
#' the total \eqn{\eta^2} of every effect and the error's share sum to 1.
#' In a balanced design Types I, II, and III agree; in an unbalanced one
#' they differ, and the order of the terms changes the Type I values. For
#' Type II or III sums of squares, compute them with, for example,
#' \code{car::Anova(object, type = 3)}, and give the effect's sum of
#' squares with the error and corrected total sums of squares to the
#' sums-of-squares interface (Type III sums of squares need not add up to
#' the total in an unbalanced design, so the total values of the effects
#' then need not sum with the error's share to 1).
#' Type III tests depend on how the factors are coded. The Type III
#' main-effect tests compare unweighted marginal means only when the
#' factors are coded with contrasts that sum to zero over their levels,
#' such as \code{\link[stats]{contr.sum}} (\code{contr.helmert} and
#' \code{contr.poly} also qualify). Under R's default treatment contrasts,
#' the Type III test of a factor that enters an interaction compares its
#' levels only at the reference levels of the factors it interacts with,
#' which is seldom the hypothesis intended. Set the coding before the
#' model is fitted, either in the call
#' (\code{lm(y ~ A * B, data = d, contrasts = list(A = contr.sum,
#' B = contr.sum))}) or for the session
#' (\code{options(contrasts = c("contr.sum", "contr.poly"))}). The help
#' page of \code{car::Anova()} gives the same warning, and notes that when
#' a model has numeric covariates and interactions, the Type III test of a
#' factor compares the groups where the covariates are 0; Type II tests do
#' not depend on the coding. \code{\link{factorial_anova}} applies
#' sum-to-zero contrasts itself and reports partial \eqn{\eta^2} and
#' partial \eqn{\omega^2} with their confidence intervals for
#' between-subjects factorial designs.
#'
#' \strong{Interactions and covariates.} An interaction is an effect like
#' any other, with its own sum of squares and degrees of freedom, and so
#' is a covariate of an \code{lm} fit: each is a row, its partial value
#' computed against the error, its total against the corrected total.
#'
#' \strong{Within-subjects designs.} For an \code{aovlist} fit, each
#' effect's partial \eqn{\eta^2} is computed against the error of its own
#' stratum, so a within-subjects factor's partial value leaves the
#' differences among subjects out of its denominator. Its total value
#' divides by the corrected total over all strata (the sum of every sum of
#' squares that \code{summary(object)} reports), which includes the
#' differences among subjects: in a design with one within-subjects factor
#' the total value is the factor's sum of squares over the sum of those
#' for subjects, the factor, and their interaction (Bakeman, 2005, p. 380),
#' and is smaller than the partial value. Neither is comparable with the
#' same effect in a between-subjects design; see
#' \code{\link{eta_squared_generalized}} for an index that is (Olejnik &
#' Algina, 2003; Bakeman, 2005). Other model classes (\code{lmerMod},
#' \code{lme}) are not taken; give their sums of squares or \emph{F}
#' statistics to the summary interfaces.
#'
#' @references
#' Bakeman, R. (2005). Recommended effect size statistics for repeated
#'   measures designs. \emph{Behavior Research Methods, 37}(3), 379--384.
#'   \doi{10.3758/BF03192707}
#'
#' Bargmann, R. E. (1970). Interpretation and use of a generalized
#'   discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
#'   Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), \emph{Essays in
#'   probability and statistics} (pp. 35--60). University of North
#'   Carolina Press.
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'   \doi{10.1177/001316448004000309}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
#'   size calculation 3. Comparison of several means using one-way ANOVA.
#'   \emph{Restorative Dentistry & Endodontics, 41}(3), 231--234.
#'   \doi{10.5395/rde.2016.41.3.231}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
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
#'   intervals and tests of close fit in the analysis of variance and contrast
#'   analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Venables, W. (1975). Calculation of confidence intervals for
#'   noncentrality parameters. \emph{Journal of the Royal Statistical
#'   Society, Series B, 37}(3), 406--412.
#'   \doi{10.1111/j.2517-6161.1975.tb01554.x}
#'
#' @seealso \code{\link{ci_eta_squared}} for the confidence interval,
#'   \code{\link{omega_squared}}, \code{\link{eta_squared_generalized}},
#'   \code{\link{ci_pvaf}}
#'
#' @examples
#' # 1. The sums of squares of a published two-way table: Kim (2016,
#' #    Table 3), light by material, 2 by 4 with 10 per cell (N = 80).
#' #    SS_total is the corrected total, 5863.715.
#' eta_squared(SS_effect = c(LIGHT = 34.716, MATERIAL = 1996.998,
#'                           "LIGHT:MATERIAL" = 1570.389),
#'             SS_error = 2261.612, SS_total = 5863.715,
#'             df_effect = c(1, 3, 3), df_error = 72)
#'
#' # 2. An F statistic and its degrees of freedom determine the partial
#' #    value only. Bargmann's (1970) example, used in turn by Venables
#' #    (1975), Fleishman (1980), and Steiger (2004): a one-way ANOVA with
#' #    five groups of 11 and an observed F of 11.221 on 4 and 50 df. In a
#' #    one-way design the partial value is the total.
#' eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
#'             estimand = "partial")
#'
#' # 3. One-way ANOVA from a fitted model (depression_bdi: three
#' #    treatment arms, 10 per arm, N = 30); the two estimands agree.
#' fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
#' eta_squared(fit_one)
#'
#' # 4. Factorial ANOVA (pygmalion: expectancy treatment by grade, 2 by 6
#' #    with unequal cell sizes, N = 310). The treatment is manipulated;
#' #    grade is a measured classification of the pupils.
#' fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
#' eta_squared(fit_factorial)
#'
#' # 5. Within-subjects ANOVA: 20 simulated subjects, each measured at Pre,
#' #    Mid, and Post. The total value includes the differences among
#' #    subjects in its denominator; the partial value does not.
#' set.seed(113)
#' n <- 20
#' rm_data <- data.frame(
#'   subject = factor(rep(seq_len(n), each = 3)),
#'   time    = factor(rep(c("Pre", "Mid", "Post"), n),
#'                    levels = c("Pre", "Mid", "Post")),
#'   y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
#'             0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
#' )
#' fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
#' eta_squared(fit_rm)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family effect size estimates
#'
#' @export
#' @import stats

eta_squared <- function(
  object    = NULL,
  estimand  = c("all", "total", "partial"),
  SS_effect = NULL,
  SS_error  = NULL,
  SS_total  = NULL,
  F_value   = NULL,
  df_effect = NULL,
  df_error  = NULL
) {
  estimand <- .anova_es_estimand(estimand)
  spec <- .anova_es_spec("eta_squared", object = object,
                         SS_effect = SS_effect, SS_error = SS_error,
                         SS_total = SS_total, F_value = F_value,
                         df_effect = df_effect, df_error = df_error)
  .anova_es_table("eta_squared", "eta_squared", spec, .anova_es_eta(spec),
                  estimand)
}
