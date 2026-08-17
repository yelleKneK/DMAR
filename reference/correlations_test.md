# Formatted Correlation Matrix With *p*-values and Confidence Intervals

Computes a correlation matrix along with, for every pair of variables,
the two-sided *p*-value, a confidence interval, and the pairwise sample
size, and arranges them into a single annotated table. The table is a
convenience for inspecting the correlations, their significance, and
their intervals at a glance; it is not intended as a finished,
publication-ready exhibit. Output formats include plain text (for the
console or to paste into a Word document), HTML (best for Word via
browser copy-paste), and LaTeX.

## Usage

``` r
correlations_test(
  x,
  method = "pearson",
  conf_level = 0.95,
  listwise = FALSE,
  stars = FALSE,
  decimals_r = 2,
  decimals_p = 4,
  format = "text",
  file = NULL
)
```

## Arguments

- x:

  A `data.frame`, tibble, or `matrix`. Non-numeric columns are dropped
  with no warning; the remaining numeric, integer, and logical columns
  are used.

- method:

  The correlation method: `"pearson"` (default), `"spearman"`, or
  `"kendall"`.

- conf_level:

  Confidence level for the interval (default `0.95`).

- listwise:

  Logical. If `TRUE`, apply listwise deletion before any correlation is
  computed so that every pair uses the same sample. If `FALSE` (the
  default), each pair uses all of its available complete observations
  (“pairwise” deletion).

- stars:

  Logical. If `TRUE`, append conventional significance stars (`*` *p* \<
  .05, `**` *p* \< .01, `***` *p* \< .001) to each correlation and add
  an explanatory footnote to the table. The exact *p*-value is always
  shown as well.

- decimals_r:

  Number of decimals for correlations and confidence interval limits
  (default `2`).

- decimals_p:

  Number of decimals for *p*-values (default `4`, matching the
  package-wide `digits_p` convention used by
  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md));
  values below \\10^{-\mathrm{decimals\\p}}\\ are printed as “\< .0001”
  (with the threshold tracking `decimals_p`).

- format:

  One of `"text"` (default), `"html"`, or `"latex"`. Controls the
  returned/printed table. See Details.

- file:

  Optional file path. If supplied, the formatted table is written to
  this file. The HTML path writes a self-contained HTML document (the
  kable wrapped in a minimal document head that pulls in Bootstrap CSS
  from a CDN) so the file opens directly in a browser without any pandoc
  / webshot machinery. The LaTeX path writes the raw `tabular` fragment;
  embed it in a document that loads `\usepackage{makecell}` and
  `\usepackage{booktabs}`. The text path uses
  [`writeLines`](https://rdrr.io/r/base/writeLines.html).

## Value

An object of class `"correlations_test"` containing the matrices `r`,
`p`, `ci_lower`, `ci_upper`, and `n` (all \\p \times p\\ with variable
names as row/column names), plus the arguments used. When
`format = "html"` or `format = "latex"` and `file` is `NULL`, a
[`kable`](https://rdrr.io/pkg/knitr/man/kable.html) object is returned
instead so that the table renders inside R Markdown / Quarto documents.
When `format = "text"`, the table is printed to the console and the raw
object is returned invisibly.

## Details

**Layout.** Each lower-triangle cell stacks four values: the correlation
(with optional significance stars), the two-sided *p*-value, the
confidence interval, and the pairwise sample size. The upper triangle is
left blank so the same information is not repeated.

***p*-values.** Computed with
[`cor.test`](https://rdrr.io/r/stats/cor.test.html) using the requested
`method`. For Spearman and Kendall with ties, `cor.test` cannot compute
an exact *p*-value and falls back to a normal-approximation *p*-value;
the associated warnings are suppressed for a cleaner table.

**Confidence intervals.** All three methods use Fisher's variance-
stabilizing transformation, \\z(r) = \mathrm{atanh}(r)\\, and
back-transform through \\\tanh(\cdot)\\, but the standard error in the
Fisher-*z* scale is selected to match the sampling distribution of the
chosen correlation coefficient:

- Pearson (`method = "pearson"`): \\\mathrm{SE}(z) = 1/\sqrt{n - 3}\\.
  This is the classical Fisher (1921) interval. Under bivariate
  normality the coverage of this interval matches
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html)'s `conf.int`
  exactly; for departures from bivariate normality both intervals lose
  coverage in the same way. See Kelley (2007) and Maxwell, Delaney, &
  Kelley (2027, Chapter 9) for discussion and worked examples.

- Spearman (`method = "spearman"`): \\\mathrm{SE}(z) = \sqrt{(1 +
  r^{2}/2) / (n - 3)}\\. This is the Bonett and Wright (2000)
  adjustment, which uses Fisher's transformation but inflates the
  standard error to account for the heavier-than-Pearson tails of the
  Spearman sampling distribution. This is the form Bonett and Wright
  recommend for practical use; a plain Fisher \\1/\sqrt{n-3}\\ standard
  error tends to produce intervals that are too narrow for Spearman
  correlations.

- Kendall (`method = "kendall"`): \\\mathrm{SE}(z) = \sqrt{0.437 / (n -
  4)}\\. This is Bonett and Wright's (2000, equation 2) Fisher-*z*
  interval for Kendall's \\\tau\\. The constant 0.437 is the asymptotic
  variance factor of Fieller, Hartley, and Pearson (1957), derived under
  bivariate normality and stated by Bonett and Wright as accurate for
  \\\|\tau\| \< .8\\ (the Spearman variance above is likewise stated as
  accurate for \\\|\rho_s\| \< .95\\). Requires \\n \ge 5\\; for smaller
  pairwise samples the interval is returned as `NA`.

The Fisher-*z* machinery requires \\\|r\| \< 1\\ for the transformation
to be finite. When \\r = \pm 1\\ (perfect correlation in the sample),
the transformed value is infinite and the interval is reported as `NA`;
this is the same convention used by `cor.test`.

**When to use each correlation.** Pick `method = "pearson"` when both
variables are continuous, approximately linearly related, and roughly
bivariate normal (or at least without heavy tails and influential
outliers). Pick `method = "spearman"` or `method = "kendall"` when the
relationship is monotone but not necessarily linear, when one or both
variables are ordinal, or when influential outliers would distort
Pearson's *r*. Kendall's \\\tau\\ is often preferred over Spearman's
\\\rho\\ for small samples and for samples with many tied ranks because
it has better small-sample properties and a more interpretable
concordance-based meaning. See Maxwell, Delaney, & Kelley (2027, Chapter
9) for an extended discussion of effect size choice and interval
estimation.

**HTML/LaTeX output.** Built with
[`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html) (a `Suggests`
dependency). Cell content and variable names are escaped for the target
format so that “*p* \< .001” renders correctly and that variable names
containing characters such as `_`, `&`, or `%` do not break LaTeX
compilation. LaTeX output uses `\makecell`, which requires
`\usepackage{makecell}` in the document preamble.

**Pasting into Word.** The cleanest path is `format = "html"` with a
`file` argument; the function writes a small self-contained HTML
document (no pandoc dependency). Open the result in a browser and
copy/paste the table into Word. Formatting (including the stacked-cell
layout) is preserved.

## References

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). Tests for rank
correlation coefficients. I. *Biometrika, 44*(3/4), 470–481.
[doi:10.1093/biomet/44.3-4.470](https://doi.org/10.1093/biomet/44.3-4.470)

Fisher, R. A. (1915). Frequency distribution of the values of the
correlation coefficient in samples from an indefinitely large
population. *Biometrika, 10*(4), 507–521.
[doi:10.1093/biomet/10.4.507](https://doi.org/10.1093/biomet/10.4.507)

Fisher, R. A. (1921). On the “probable error” of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md),
[`cor.test`](https://rdrr.io/r/stats/cor.test.html),
[`cor`](https://rdrr.io/r/stats/cor.html),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Worked example using four cognitive tests from the Holzinger and
# Swineford (1939) study (301 children in two schools). The goal of
# correlations_test() is to produce a formatted correlation matrix that
# reports, for every variable pair, the correlation, its two-sided
# p-value, a confidence interval on the population correlation, and the
# pairwise sample size. See Kelley (2007) and Maxwell, Delaney, & Kelley
# (2027, Chapter 9) for discussion of why effect sizes should be
# accompanied by confidence intervals.
hs_tests <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                                    "t4_lozenges",
                                    "t6_paragraph_comprehension")]

# Pearson correlations (the default). Each lower-triangle cell stacks r,
# the two-sided p-value, the 95\% confidence interval (Fisher's Z
# transformation; Fisher, 1915, 1921), and the pairwise N.
correlations_test(hs_tests)
#> Correlations (Pearson, 95% CI)
#> 
#>                             t1_visual_perception          t2_cubes                      t4_lozenges                   t6_paragraph_comprehension    
#> ----------------------------------------------------------------------------------------------------------------------------------------------------
#> t1_visual_perception        -                                                                                                                       
#>                                                                                                                                                     
#>                                                                                                                                                     
#>                                                                                                                                                     
#> 
#> t2_cubes                    .30                           -                                                                                         
#>                             p < .0001                                                                                                               
#>                             [.19, .40]                                                                                                              
#>                             N = 301                                                                                                                 
#> 
#> t4_lozenges                 .44                           .34                           -                                                           
#>                             p < .0001                     p < .0001                                                                                 
#>                             [.34, .53]                    [.24, .44]                                                                                
#>                             N = 301                       N = 301                                                                                   
#> 
#> t6_paragraph_comprehension  .37                           .15                           .16                           -                             
#>                             p < .0001                     p = .0079                     p = .0058                                                   
#>                             [.27, .47]                    [.04, .26]                    [.05, .27]                                                  
#>                             N = 301                       N = 301                       N = 301                                                     
#> 

# Add significance stars and an explanatory footnote.
correlations_test(hs_tests, stars = TRUE)
#> Correlations (Pearson, 95% CI)
#> 
#>                             t1_visual_perception          t2_cubes                      t4_lozenges                   t6_paragraph_comprehension    
#> ----------------------------------------------------------------------------------------------------------------------------------------------------
#> t1_visual_perception        -                                                                                                                       
#>                                                                                                                                                     
#>                                                                                                                                                     
#>                                                                                                                                                     
#> 
#> t2_cubes                    .30***                        -                                                                                         
#>                             p < .0001                                                                                                               
#>                             [.19, .40]                                                                                                              
#>                             N = 301                                                                                                                 
#> 
#> t4_lozenges                 .44***                        .34***                        -                                                           
#>                             p < .0001                     p < .0001                                                                                 
#>                             [.34, .53]                    [.24, .44]                                                                                
#>                             N = 301                       N = 301                                                                                   
#> 
#> t6_paragraph_comprehension  .37***                        .15**                         .16**                         -                             
#>                             p < .0001                     p = .0079                     p = .0058                                                   
#>                             [.27, .47]                    [.04, .26]                    [.05, .27]                                                  
#>                             N = 301                       N = 301                       N = 301                                                     
#> 
#> Note. * p < .05, ** p < .01, *** p < .001.

# Spearman correlations at a 99\% confidence level. The interval uses
# Bonett and Wright's (2000) Fisher's Z standard error
# sqrt((1 + r^2/2) / (n - 3)), which corrects the plain Fisher interval
# for the heavier tails of Spearman's sampling distribution.
correlations_test(hs_tests, method = "spearman", conf_level = 0.99)
#> Correlations (Spearman, 99% CI)
#> 
#>                             t1_visual_perception          t2_cubes                      t4_lozenges                   t6_paragraph_comprehension    
#> ----------------------------------------------------------------------------------------------------------------------------------------------------
#> t1_visual_perception        -                                                                                                                       
#>                                                                                                                                                     
#>                                                                                                                                                     
#>                                                                                                                                                     
#> 
#> t2_cubes                    .28                           -                                                                                         
#>                             p < .0001                                                                                                               
#>                             [.14, .42]                                                                                                              
#>                             N = 301                                                                                                                 
#> 
#> t4_lozenges                 .44                           .32                           -                                                           
#>                             p < .0001                     p < .0001                                                                                 
#>                             [.31, .56]                    [.18, .45]                                                                                
#>                             N = 301                       N = 301                                                                                   
#> 
#> t6_paragraph_comprehension  .38                           .19                           .15                           -                             
#>                             p < .0001                     p = .0007                     p = .0105                                                   
#>                             [.24, .50]                    [.05, .33]                    [-.00, .29]                                                 
#>                             N = 301                       N = 301                       N = 301                                                     
#> 

# Kendall's tau, also using Bonett and Wright's (2000) Fisher's Z standard
# error sqrt(0.437 / (n - 4)). Kendall is often preferred over Spearman
# for small samples and for samples with many tied ranks, and these
# integer test scores carry many ties.
correlations_test(hs_tests, method = "kendall")
#> Correlations (Kendall, 95% CI)
#> 
#>                             t1_visual_perception          t2_cubes                      t4_lozenges                   t6_paragraph_comprehension    
#> ----------------------------------------------------------------------------------------------------------------------------------------------------
#> t1_visual_perception        -                                                                                                                       
#>                                                                                                                                                     
#>                                                                                                                                                     
#>                                                                                                                                                     
#> 
#> t2_cubes                    .21                           -                                                                                         
#>                             p < .0001                                                                                                               
#>                             [.13, .28]                                                                                                              
#>                             N = 301                                                                                                                 
#> 
#> t4_lozenges                 .32                           .23                           -                                                           
#>                             p < .0001                     p < .0001                                                                                 
#>                             [.25, .38]                    [.15, .30]                                                                                
#>                             N = 301                       N = 301                                                                                   
#> 
#> t6_paragraph_comprehension  .28                           .14                           .10                           -                             
#>                             p < .0001                     p = .0007                     p = .0104                                                   
#>                             [.21, .34]                    [.07, .21]                    [.03, .18]                                                  
#>                             N = 301                       N = 301                       N = 301                                                     
#> 

# Save a formatted HTML table that opens directly in a browser
# (then copy into Word). No pandoc required.
tmp_html <- tempfile(fileext = ".html")
correlations_test(hs_tests, stars = TRUE, format = "html", file = tmp_html)
```
