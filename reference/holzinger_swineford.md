# Holzinger and Swineford (1939) Factor Analysis Study

The complete data set from Holzinger and Swineford's (1939) *A study in
factor analysis: The stability of a bi-factor solution*. Scores on 26
ability tests for 301 seventh and eighth grade pupils at two Chicago
elementary schools, Pasteur (*n* = 156) and Grant-White (*n* = 145). The
data have been used over the subsequent decades as one of the most-cited
benchmarks in factor analysis, confirmatory factor analysis, structural
equation modeling, and reliability research.

## Usage

``` r
holzinger_swineford
```

## Format

A data frame with 301 observations and 34 variables.

- `id`:

  Case identifier as in the original monograph. The numbering is not
  strictly consecutive; a small number of cases were dropped during data
  preparation and the original numbering was preserved (hence the
  visible skips).

- `sex`:

  Factor with levels `Female` and `Male`.

- `grade`:

  Grade in school, 7 or 8.

- `age`:

  Age in completed years, ignoring months past the most recent birthday.

- `month_since_birthday`:

  Completed months since the most recent birthday.

- `age_months`:

  Age in completed months, computed as
  `12 * age + month_since_birthday`.

- `age_years`:

  Age in fractional years, computed as
  `age + month_since_birthday / 12`.

- `school`:

  Factor with levels `Grant-White` and `Pasteur`, naming the two Chicago
  elementary schools from which pupils were drawn.

- `t1_visual_perception`:

  Visual perception test (spatial).

- `t2_cubes`:

  Cubes test (spatial).

- `t3_paper_form_board`:

  Paper form board test (spatial).

- `t4_lozenges`:

  Lozenges test (spatial).

- `t5_general_information`:

  General information test (verbal).

- `t6_paragraph_comprehension`:

  Paragraph comprehension test (verbal).

- `t7_sentence`:

  Sentence completion test (verbal).

- `t8_word_classification`:

  Word classification test (verbal).

- `t9_word_meaning`:

  Word meaning test (verbal).

- `t10_addition`:

  Addition test (mental speed).

- `t11_code`:

  Code test (mental speed).

- `t12_counting_groups_of_dots`:

  Counting groups of dots test (mental speed).

- `t13_straight_and_curved_capitals`:

  Straight and curved capitals test (mental speed).

- `t14_word_recognition`:

  Word recognition test (memory).

- `t15_number_recognition`:

  Number recognition test (memory).

- `t16_figure_recognition`:

  Figure recognition test (memory).

- `t17_object_number`:

  Object-number test (memory).

- `t18_number_figure`:

  Number-figure test (memory).

- `t19_figure_word`:

  Figure-word test (memory).

- `t20_deduction`:

  Deduction test (reasoning).

- `t21_numerical_puzzles`:

  Numerical puzzles test (reasoning).

- `t22_problem_reasoning`:

  Problem reasoning test (reasoning).

- `t23_series_completion`:

  Series completion test (reasoning).

- `t24_woody_mccall`:

  Woody-McCall mixed fundamentals, form I (arithmetic).

- `t25_paper_form_board_r`:

  Revised paper form board, administered only to the Grant-White pupils
  as an experimental substitute for `t3_paper_form_board`. `NA` for the
  156 Pasteur pupils.

- `t26_flags`:

  Flags test, administered only to the Grant-White pupils as an
  experimental substitute for `t4_lozenges`. `NA` for the 156 Pasteur
  pupils.

## Source

Holzinger, K. J., and Swineford, F. (1939). *A study in factor analysis:
The stability of a bi-factor solution* (Supplementary Educational
Monographs, No. 48). University of Chicago Press.

## Details

Karl John Holzinger (1893 to 1954) was a quantitative psychologist at
the University of Chicago and one of the central figures in the first
generation of factor analysis. He spent the 1922 to 1923 academic year
working with Charles Spearman at University College London, absorbing
Spearman's two-factor theory of intelligence, and later developed the
*bi-factor model* as an extension of that theory. The bi-factor model
posits a single general intelligence factor that runs through all tests,
plus several group factors that capture residual correlation among
substantively related subgroups of tests. It is widely regarded as a
precursor of modern hierarchical and orthogonal-bifactor models in
psychometrics. Frances Swineford was Holzinger's research collaborator
at the University of Chicago and a coauthor on much of his applied work.

The 1939 monograph reports a study of pupils in seventh and eighth grade
classrooms at two Chicago elementary schools, Pasteur and Grant-White.
Two schools were used deliberately, so that the stability of a bi-factor
solution could be assessed by fitting the same model in each school and
comparing the results. The 26 tests were designed to span five
hypothesized ability domains:

- **Spatial**: tests 1 to 4 (visual perception, cubes, paper form board,
  lozenges), with tests 25 (a revised paper form board) and 26 (flags)
  administered only to the Grant-White sample as experimental
  substitutes for tests 3 and 4.

- **Verbal**: tests 5 to 9 (general information, paragraph
  comprehension, sentence completion, word classification, word
  meaning).

- **Mental speed**: tests 10 to 13 (addition, code, counting groups of
  dots, straight and curved capitals).

- **Memory**: tests 14 to 19 (word recognition, number recognition,
  figure recognition, object-number, number-figure, figure-word).

- **Reasoning and arithmetic**: tests 20 to 24 (deduction, numerical
  puzzles, problem reasoning, series completion, Woody-McCall mixed
  fundamentals).

Holzinger and Swineford concluded that the bi-factor solution was
reasonably stable across the two schools, supporting the substantive
interpretation of a general factor together with group factors.

The data have far outlived their original purpose. Jöreskog (1969) used
a 9-test subset drawn from the Grant-White sample (*n* = 145) to
introduce confirmatory maximum likelihood factor analysis; that 9-test
subset is the version most modern confirmatory factor analysis tutorials
use and is shipped in per-item-rescaled form as `HolzingerSwineford1939`
in the lavaan package. The complete 26-test data shipped here support a
wider range of analyses, including comparisons of the spatial, verbal,
speed, memory, and reasoning ability blocks and multiple group analyses
across the Pasteur and Grant-White schools.

The values in `holzinger_swineford` are the corrected version of the
data, identical on all 26 test cells to `MBESS::HS` from MBESS version
4.9.3 onward and to `psychTools::holzinger.raw`. An older version of the
data, with approximately 53 cell values that were later corrected,
continues to circulate as `HS.data` in the sem package and as
`HS.ability.data` in the OpenMx package; both are byte-identical
snapshots taken from MBESS version 4.6.0 prior to the correction. The
corrections are concentrated on the memory and reasoning tests, with the
largest cluster on `t20_deduction` (15 cells, including a number of sign
flips that reflect a corrected guessing-penalty adjustment).

## References

Holzinger, K. J., and Swineford, F. (1939). *A study in factor analysis:
The stability of a bi-factor solution* (Supplementary Educational
Monographs, No. 48). University of Chicago Press.

Jöreskog, K. G. (1969). A general approach to confirmatory maximum
likelihood factor analysis. *Psychometrika, 34*, 183–202.

Holzinger, K. J. (1944). A simple method of factor analysis.
*Psychometrika, 9*, 257–262.

## Author

Ken Kelley

## Examples

``` r
data(holzinger_swineford)
str(holzinger_swineford)
#> 'data.frame':    301 obs. of  34 variables:
#>  $ id                              : int  1 2 3 4 5 6 7 8 9 11 ...
#>  $ sex                             : Factor w/ 2 levels "Female","Male": 2 1 1 2 1 1 2 1 1 1 ...
#>  $ grade                           : int  7 7 7 7 7 7 7 7 7 7 ...
#>  $ age                             : int  13 13 13 13 12 14 12 12 13 12 ...
#>  $ month_since_birthday            : int  1 7 1 2 2 1 1 2 0 5 ...
#>  $ age_months                      : int  157 163 157 158 146 169 145 146 156 149 ...
#>  $ age_years                       : num  13.1 13.6 13.1 13.2 12.2 ...
#>  $ school                          : Factor w/ 2 levels "Grant-White",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ t1_visual_perception            : int  20 32 27 32 29 32 17 34 27 21 ...
#>  $ t2_cubes                        : int  31 21 21 31 19 20 24 25 23 21 ...
#>  $ t3_paper_form_board             : int  12 12 12 16 12 11 12 13 11 10 ...
#>  $ t4_lozenges                     : int  3 17 15 24 7 18 8 15 12 6 ...
#>  $ t5_general_information          : int  40 34 20 42 37 31 40 29 29 33 ...
#>  $ t6_paragraph_comprehension      : int  7 5 3 8 8 3 10 11 8 8 ...
#>  $ t7_sentence                     : int  23 12 7 18 16 12 24 17 23 20 ...
#>  $ t8_word_classification          : int  22 22 12 21 25 25 32 25 19 25 ...
#>  $ t9_word_meaning                 : int  9 9 3 17 18 6 20 9 19 18 ...
#>  $ t10_addition                    : int  78 87 75 69 85 100 108 78 104 95 ...
#>  $ t11_code                        : int  74 84 49 65 63 92 65 80 52 74 ...
#>  $ t12_counting_groups_of_dots     : int  115 125 78 106 126 133 124 103 93 91 ...
#>  $ t13_straight_and_curved_capitals: int  229 285 159 175 213 270 175 132 265 157 ...
#>  $ t14_word_recognition            : int  170 184 170 181 187 164 121 184 184 175 ...
#>  $ t15_number_recognition          : int  86 85 85 80 99 84 71 95 91 92 ...
#>  $ t16_figure_recognition          : int  96 100 95 91 104 104 78 106 105 100 ...
#>  $ t17_object_number               : int  6 12 1 5 15 6 4 11 18 5 ...
#>  $ t18_number_figure               : int  9 12 5 3 14 6 3 13 6 8 ...
#>  $ t19_figure_word                 : int  16 10 6 10 14 14 5 9 11 11 ...
#>  $ t20_deduction                   : int  3 -3 -3 -2 29 9 18 15 12 33 ...
#>  $ t21_numerical_puzzles           : int  14 13 9 10 15 2 10 9 15 8 ...
#>  $ t22_problem_reasoning           : int  34 21 18 22 19 16 19 22 18 25 ...
#>  $ t23_series_completion           : int  5 1 7 6 4 10 3 18 17 8 ...
#>  $ t24_woody_mccall                : int  24 12 20 19 20 22 15 24 18 16 ...
#>  $ t25_paper_form_board_r          : int  NA NA NA NA NA NA NA NA NA NA ...
#>  $ t26_flags                       : int  NA NA NA NA NA NA NA NA NA NA ...

# School and grade breakdown.
table(holzinger_swineford$school, holzinger_swineford$grade)
#>              
#>                7  8
#>   Grant-White 79 66
#>   Pasteur     78 78

# The 9 test subset drawn by Jöreskog (1969) from the Grant-White
# sample, which became the modern confirmatory factor analysis
# benchmark.
joreskog_subset <- subset(
  holzinger_swineford,
  school == "Grant-White",
  select = c(t1_visual_perception, t2_cubes, t4_lozenges,
             t6_paragraph_comprehension, t7_sentence,
             t9_word_meaning, t10_addition,
             t12_counting_groups_of_dots,
             t13_straight_and_curved_capitals)
)
dim(joreskog_subset)
#> [1] 145   9
```
