#' Pygmalion in the Classroom Teacher-Expectancy Data
#'
#' The teacher-expectancy data from Rosenthal and Jacobson's (1968)
#' \emph{Pygmalion in the Classroom}, the study that introduced the
#' "Pygmalion effect": the hypothesis that a teacher's expectations can
#' become a self-fulfilling prophecy for a pupil's intellectual growth.
#' Intelligence-test scores were obtained for \emph{N} = 310 elementary
#' school children in grades 1 through 6, of whom \emph{n} = 64 were
#' randomly designated to their teachers as likely "intellectual
#' bloomers" while the remaining \emph{n} = 246 served as controls.
#' The data set is a classic benchmark for the analysis of covariance
#' (ANCOVA) and, in particular, for ANCOVA with \emph{heterogeneity of
#' regression}: it is the running example for that topic in Maxwell,
#' Delaney, and Kelley, \emph{Designing Experiments and Analyzing Data:
#' A Model Comparison Perspective} (Routledge), where it appears as a
#' Chapter 9 example (and as a Chapter 3 exercise).
#'
#' @format A data frame with 310 observations on 6 variables.
#' \describe{
#'   \item{\code{grade}}{Grade in school at the start of the study,
#'     an integer from 1 to 6.}
#'   \item{\code{treatment}}{Factor with levels \code{Control}
#'     (reference, \emph{n} = 246) and \code{Bloomer} (\emph{n} = 64).
#'     The \code{Bloomer} children were a randomly selected ~20\% of
#'     each classroom whose teachers were told, on the basis of a
#'     fictitious test purportedly predicting intellectual blooming,
#'     that they were likely to show unusual gains during the year; the
#'     \code{Control} children were not singled out. In the AMCP source
#'     this variable is coded \code{1} = Bloomer, \code{0} = Control.}
#'   \item{\code{iq_pre}}{Pretest total IQ, measured before the
#'     expectancy manipulation. The covariate in the analysis of
#'     covariance.}
#'   \item{\code{iq_4}}{Total IQ at an intermediate follow-up
#'     assessment.}
#'   \item{\code{iq_8}}{Total IQ at the end-of-study follow-up
#'     assessment. This is the dependent variable in the book's
#'     Chapter 9 analysis of covariance.}
#'   \item{\code{iq_gain}}{Total IQ change from pretest to the
#'     end-of-study assessment, equal to \code{iq_8 - iq_pre}.}
#' }
#'
#' @details
#' \strong{The study.} Robert Rosenthal (Harvard University) and Lenore
#' Jacobson (principal of an elementary school in South San Francisco
#' referred to as "Oak School") set out to test experimentally whether
#' teacher expectations influence pupil achievement. At the start of the
#' school year all children were given a standardized test of general
#' ability, described to teachers as the "Harvard Test of Inflected
#' Acquisition," a test said to identify children poised for an
#' intellectual growth spurt. In reality the instrument was Flanagan's
#' Tests of General Ability (TOGA) and the children identified as likely
#' "bloomers" were chosen \emph{at random}, about one in five per
#' classroom. The only experimental manipulation was the expectation
#' planted in the teachers' minds. Children were re-tested over the
#' following year(s), and the question was whether the randomly labeled
#' bloomers would out-gain their controls in measured IQ. Rosenthal and
#' Jacobson reported that they did, most strongly in the earliest
#' grades, and interpreted the difference as evidence that teacher
#' expectations operate as a self-fulfilling prophecy. The study became
#' one of the most famous and most debated experiments in the social
#' sciences; subsequent critiques (e.g., Thorndike, 1968) questioned the
#' reliability of the TOGA at the extremes of the score range for the
#' youngest children, which is itself part of why the data are
#' instructive for teaching careful analysis.
#'
#' \strong{Why it is a benchmark for heterogeneity of regression.} A
#' standard ANCOVA adjusts the group comparison for the pretest
#' covariate under the assumption that the regression of the outcome on
#' the covariate has the \emph{same} slope in every group (homogeneity
#' of regression). In these data that assumption is questionable: the
#' within-group regression of \code{iq_8} on \code{iq_pre} is steeper
#' for the bloomers than for the controls, so the estimated treatment
#' effect depends on the covariate value at which it is evaluated. This
#' makes the data an ideal teaching example for (a) testing the
#' homogeneity-of-regression assumption, (b) interpreting a
#' treatment-by-covariate interaction, and (c) estimating the treatment
#' effect, and its sampling variance, \emph{at chosen covariate values}
#' rather than only at the grand mean.
#'
#' \strong{Reproducible quantities.} Fitting the separate-slopes model
#' \code{lm(iq_8 ~ iq_pre * treatment)} gives a within-group slope of
#' \eqn{0.77799} for the controls and \eqn{0.96894} for the bloomers
#' (reported as \eqn{0.96895} in \code{MBESS::var.ete}, a fifth-decimal
#' rounding difference). The pooled within-group residual variance is
#' \eqn{\hat\sigma^2 = 175.3251} on 306 degrees of freedom, and the
#' sample variance of the covariate is \eqn{348.91}. These are exactly
#' the inputs used in the worked example for the variance of the
#' estimated treatment effect at selected covariate values under
#' heterogeneity of regression (Li, McLouth, and Delaney; see
#' \code{MBESS::var.ete}).
#'
#' \strong{Relationship to the AMCP package.} The same numeric data
#' ship with the book's data companion, the \pkg{AMCP} package, as
#' \code{chapter_9_exercise_15} and \code{chapter_9_extension_exercise_3}
#' (with \code{IQGain}) and \code{chapter_3_exercise_22} (without it).
#' The version here renames the columns to DMAR's descriptive
#' snake_case style and labels the experimental condition as a factor;
#' no measured value has been altered. See \code{data-raw/pygmalion.R}
#' for the construction script and its verification checks.
#'
#' @author Ken Kelley
#'
#' @source
#' Rosenthal, R., & Jacobson, L. (1968). \emph{Pygmalion in the
#' classroom: Teacher expectation and pupils' intellectual development}.
#' Holt, Rinehart and Winston.
#'
#' Distributed with the \pkg{AMCP} data companion to Maxwell, Delaney,
#' and Kelley (see References) as \code{chapter_9_exercise_15}.
#'
#' @references
#' Rosenthal, R., & Jacobson, L. (1968). \emph{Pygmalion in the
#'   classroom: Teacher expectation and pupils' intellectual development}.
#'   Holt, Rinehart and Winston.
#'
#' Rosenthal, R., & Jacobson, L. (1968). Pygmalion in the classroom.
#' \emph{The Urban Review, 3}(1), 16--20.
#' \doi{10.1007/BF02322211}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#' \emph{Designing experiments and analyzing data: A model comparison
#' perspective} (4th ed.). Routledge. (Heterogeneity-of-regression
#' ANCOVA example, Chapter 9.)
#'
#' Thorndike, R. L. (1968). Review of \emph{Pygmalion in the
#' Classroom}. \emph{American Educational Research Journal, 5}(4),
#' 708--711.
#'
#' @seealso \code{\link{ancova}} for an ANCOVA that returns adjusted
#'   means, effect size confidence intervals, and a
#'   homogeneity-of-regression test.
#'
#' @examples
#' data(pygmalion)
#' str(pygmalion)
#'
#' # Design: pupils per condition within each grade.
#' table(pygmalion$treatment, pygmalion$grade)
#'
#' # ---- Heterogeneity-of-regression ANCOVA (book Chapter 9) ----
#' # Separate IQ8-on-IQpre slopes for the two conditions.
#' fit_het <- lm(iq_8 ~ iq_pre * treatment, data = pygmalion)
#' coef(fit_het)
#' # Control slope = 0.778; the interaction (0.191) gives the
#' # steeper Bloomer slope of 0.969.
#'
#' # The treatment-by-covariate interaction is the
#' # heterogeneity-of-regression test (1 df): compare the additive
#' # ANCOVA model to the separate-slopes model.
#' fit_add <- lm(iq_8 ~ iq_pre + treatment, data = pygmalion)
#' anova(fit_add, fit_het)
#'
#' # Pooled within-group residual variance (175.3251) and the
#' # covariate variance (348.91), as used by MBESS::var.ete.
#' sum(residuals(fit_het)^2) / fit_het$df.residual
#' var(pygmalion$iq_pre)
#'
#' # ---- DMAR's ANCOVA, with the homogeneity-of-regression check ----
#' ancova(pygmalion, outcome = "iq_8", treatment = "treatment",
#'        covariates = "iq_pre")
#' @keywords datasets
"pygmalion"
