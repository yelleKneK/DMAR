#' Holzinger and Swineford (1939) Factor Analysis Study
#'
#' The complete data set from Holzinger and Swineford's (1939)
#' \emph{A study in factor analysis: The stability of a bi-factor
#' solution}. Scores on 26 ability tests for 301 seventh and eighth
#' grade pupils at two Chicago elementary schools, Pasteur
#' (\emph{n} = 156) and Grant-White (\emph{n} = 145). The data have
#' been used over the subsequent decades as one of the most-cited
#' benchmarks in factor analysis, confirmatory factor analysis,
#' structural equation modeling, and reliability research.
#'
#' @format A data frame with 301 observations and 34 variables.
#' \describe{
#'   \item{\code{id}}{Case identifier as in the original monograph.
#'     The numbering is not strictly consecutive; a small number of
#'     cases were dropped during data preparation and the original
#'     numbering was preserved (hence the visible skips).}
#'   \item{\code{sex}}{Factor with levels \code{Female} and \code{Male}.}
#'   \item{\code{grade}}{Grade in school, 7 or 8.}
#'   \item{\code{age}}{Age in completed years, ignoring months past
#'     the most recent birthday.}
#'   \item{\code{month_since_birthday}}{Completed months since the
#'     most recent birthday.}
#'   \item{\code{age_months}}{Age in completed months, computed as
#'     \code{12 * age + month_since_birthday}.}
#'   \item{\code{age_years}}{Age in fractional years, computed as
#'     \code{age + month_since_birthday / 12}.}
#'   \item{\code{school}}{Factor with levels \code{Grant-White} and
#'     \code{Pasteur}, naming the two Chicago elementary schools from
#'     which pupils were drawn.}
#'   \item{\code{t1_visual_perception}}{Visual perception test
#'     (spatial).}
#'   \item{\code{t2_cubes}}{Cubes test (spatial).}
#'   \item{\code{t3_paper_form_board}}{Paper form board test (spatial).}
#'   \item{\code{t4_lozenges}}{Lozenges test (spatial).}
#'   \item{\code{t5_general_information}}{General information test
#'     (verbal).}
#'   \item{\code{t6_paragraph_comprehension}}{Paragraph comprehension
#'     test (verbal).}
#'   \item{\code{t7_sentence}}{Sentence completion test (verbal).}
#'   \item{\code{t8_word_classification}}{Word classification test
#'     (verbal).}
#'   \item{\code{t9_word_meaning}}{Word meaning test (verbal).}
#'   \item{\code{t10_addition}}{Addition test (mental speed).}
#'   \item{\code{t11_code}}{Code test (mental speed).}
#'   \item{\code{t12_counting_groups_of_dots}}{Counting groups of dots
#'     test (mental speed).}
#'   \item{\code{t13_straight_and_curved_capitals}}{Straight and curved
#'     capitals test (mental speed).}
#'   \item{\code{t14_word_recognition}}{Word recognition test (memory).}
#'   \item{\code{t15_number_recognition}}{Number recognition test
#'     (memory).}
#'   \item{\code{t16_figure_recognition}}{Figure recognition test
#'     (memory).}
#'   \item{\code{t17_object_number}}{Object-number test (memory).}
#'   \item{\code{t18_number_figure}}{Number-figure test (memory).}
#'   \item{\code{t19_figure_word}}{Figure-word test (memory).}
#'   \item{\code{t20_deduction}}{Deduction test (reasoning).}
#'   \item{\code{t21_numerical_puzzles}}{Numerical puzzles test
#'     (reasoning).}
#'   \item{\code{t22_problem_reasoning}}{Problem reasoning test
#'     (reasoning).}
#'   \item{\code{t23_series_completion}}{Series completion test
#'     (reasoning).}
#'   \item{\code{t24_woody_mccall}}{Woody-McCall mixed fundamentals,
#'     form I (arithmetic).}
#'   \item{\code{t25_paper_form_board_r}}{Revised paper form board,
#'     administered only to the Grant-White pupils as an experimental
#'     substitute for \code{t3_paper_form_board}. \code{NA} for the
#'     156 Pasteur pupils.}
#'   \item{\code{t26_flags}}{Flags test, administered only to the
#'     Grant-White pupils as an experimental substitute for
#'     \code{t4_lozenges}. \code{NA} for the 156 Pasteur pupils.}
#' }
#'
#' @details
#' Karl John Holzinger (1893 to 1954) was a quantitative psychologist
#' at the University of Chicago and one of the central figures in the
#' first generation of factor analysis. He spent the 1922 to 1923
#' academic year working with Charles Spearman at University College
#' London, absorbing Spearman's two-factor theory of intelligence, and
#' later developed the \emph{bi-factor model} as an extension of that
#' theory. The bi-factor model posits a single general intelligence
#' factor that runs through all tests, plus several group factors
#' that capture residual correlation among substantively related
#' subgroups of tests. It is widely regarded as a precursor of modern
#' hierarchical and orthogonal-bifactor models in psychometrics.
#' Frances Swineford was Holzinger's research collaborator at the
#' University of Chicago and a coauthor on much of his applied work.
#'
#' The 1939 monograph reports a study of pupils in seventh and eighth
#' grade classrooms at two Chicago elementary schools, Pasteur and
#' Grant-White. Two schools were used deliberately, so that the
#' stability of a bi-factor solution could be assessed by fitting the
#' same model in each school and comparing the results. The 26 tests
#' were designed to span five hypothesized ability domains:
#' \itemize{
#'   \item \strong{Spatial}: tests 1 to 4 (visual perception, cubes,
#'     paper form board, lozenges), with tests 25 (a revised paper
#'     form board) and 26 (flags) administered only to the
#'     Grant-White sample as experimental substitutes for tests 3
#'     and 4.
#'   \item \strong{Verbal}: tests 5 to 9 (general information,
#'     paragraph comprehension, sentence completion, word
#'     classification, word meaning).
#'   \item \strong{Mental speed}: tests 10 to 13 (addition, code,
#'     counting groups of dots, straight and curved capitals).
#'   \item \strong{Memory}: tests 14 to 19 (word recognition, number
#'     recognition, figure recognition, object-number, number-figure,
#'     figure-word).
#'   \item \strong{Reasoning and arithmetic}: tests 20 to 24
#'     (deduction, numerical puzzles, problem reasoning, series
#'     completion, Woody-McCall mixed fundamentals).
#' }
#'
#' Holzinger and Swineford concluded that the bi-factor solution was
#' reasonably stable across the two schools, supporting the
#' substantive interpretation of a general factor together with group
#' factors.
#'
#' The data have far outlived their original purpose. Jöreskog (1969)
#' used a 9-test subset drawn from the Grant-White sample
#' (\emph{n} = 145) to introduce confirmatory maximum likelihood
#' factor analysis; that 9-test subset is the version most modern
#' confirmatory factor analysis tutorials use and is shipped in
#' per-item-rescaled form as \code{HolzingerSwineford1939} in the
#' \pkg{lavaan} package. The complete 26-test data shipped here
#' support a wider range of analyses, including comparisons of the
#' spatial, verbal, speed, memory, and reasoning ability blocks and
#' multiple group analyses across the Pasteur and Grant-White
#' schools.
#'
#' The values in \code{holzinger_swineford} are the corrected
#' version of the data, identical on all 26 test cells to
#' \code{MBESS::HS} from MBESS version 4.9.3 onward and to
#' \code{psychTools::holzinger.raw}. An older version of the data,
#' with approximately 53 cell values that were later corrected,
#' continues to circulate as \code{HS.data} in the \pkg{sem} package
#' and as \code{HS.ability.data} in the \pkg{OpenMx} package; both
#' are byte-identical snapshots taken from \pkg{MBESS} version 4.6.0
#' prior to the correction. The corrections are concentrated on the
#' memory and reasoning tests, with the largest cluster on
#' \code{t20_deduction} (15 cells, including a number of sign flips
#' that reflect a corrected guessing-penalty adjustment).
#'
#' @author Ken Kelley
#'
#' @source
#' Holzinger, K. J., and Swineford, F. (1939). \emph{A study in factor
#' analysis: The stability of a bi-factor solution} (Supplementary
#' Educational Monographs, No. 48). University of Chicago Press.
#'
#' @references
#' Holzinger, K. J., and Swineford, F. (1939). \emph{A study in factor
#' analysis: The stability of a bi-factor solution} (Supplementary
#' Educational Monographs, No. 48). University of Chicago Press.
#'
#' Jöreskog, K. G. (1969). A general approach to confirmatory maximum
#' likelihood factor analysis. \emph{Psychometrika, 34}, 183--202.
#'
#' Holzinger, K. J. (1944). A simple method of factor analysis.
#' \emph{Psychometrika, 9}, 257--262.
#'
#' @examples
#' data(holzinger_swineford)
#' str(holzinger_swineford)
#'
#' # School and grade breakdown.
#' table(holzinger_swineford$school, holzinger_swineford$grade)
#'
#' # The 9 test subset drawn by Jöreskog (1969) from the Grant-White
#' # sample, which became the modern confirmatory factor analysis
#' # benchmark.
#' joreskog_subset <- subset(
#'   holzinger_swineford,
#'   school == "Grant-White",
#'   select = c(t1_visual_perception, t2_cubes, t4_lozenges,
#'              t6_paragraph_comprehension, t7_sentence,
#'              t9_word_meaning, t10_addition,
#'              t12_counting_groups_of_dots,
#'              t13_straight_and_curved_capitals)
#' )
#' dim(joreskog_subset)
#'
#' @keywords datasets
"holzinger_swineford"
