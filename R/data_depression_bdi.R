# The Chapter 9 depression pre-post data of Maxwell, Delaney, and Kelley.
#' Depression Treatment Study With a Pretest Covariate
#'
#' The hypothetical three-group depression study that runs through the
#' analysis of covariance development of Maxwell, Delaney, and Kelley
#' (2027, Chapter 9, Table 9.7). Thirty depressive individuals are
#' randomly assigned, ten per group, to a selective serotonin reuptake
#' inhibitor (SSRI), a placebo, or a wait list control. The Beck
#' Depression Inventory (BDI) is administered before the study begins
#' and again at its end, giving a pretest that serves as the covariate
#' and a posttest that serves as the outcome.
#'
#' This is the worked example behind the ANCOVA contrast pages: the
#' pretest group means are 17, 17.7, and 17.4; the within-groups sum of
#' squares of the pretest is 752.5; the ANCOVA error variance is about
#' 29; and the covariate-adjusted posttest means are approximately 7.5,
#' 12, and 14 for the SSRI, placebo, and wait list groups. The examples
#' of \code{\link{ci_c_ancova}} and \code{\link{ci_sc_ancova}} quote
#' those summary values, and a test recomputes each of them from these
#' data so the printed numbers cannot drift.
#'
#' @format A \code{data.frame} with 30 rows and 3 columns:
#' \describe{
#'   \item{condition}{Factor with levels \code{ssri}, \code{placebo},
#'     and \code{wait_list}: the randomly assigned treatment.}
#'   \item{bdi_pre}{Beck Depression Inventory score before the study.}
#'   \item{bdi_post}{Beck Depression Inventory score at the end of the
#'     study.}
#' }
#'
#' @details
#' The data are hypothetical, constructed for the book; higher BDI
#' scores mean more severe depressive symptoms, so a treatment that
#' works pulls the posttest down. The same numeric values ship in the
#' book's data companion, the \pkg{AMCP} package, as
#' \code{chapter_9_table_7}; DMAR carries them directly so its ANCOVA
#' examples and tests need no package beyond this one.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9 on ANCOVA.)
#'
#' @examples
#' data(depression_bdi)
#'
#' # The fingerprints the ANCOVA contrast pages quote.
#' tapply(depression_bdi$bdi_pre, depression_bdi$condition, mean)
#' fit <- lm(bdi_post ~ bdi_pre + condition, data = depression_bdi)
#' anova(fit)
#'
#' # The covariate-adjusted group means at the pretest grand mean.
#' ancova(depression_bdi, outcome = "bdi_post", treatment = "condition",
#'        covariates = "bdi_pre")
#'
#' @keywords datasets
"depression_bdi"
