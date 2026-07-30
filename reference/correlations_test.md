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
[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)/[`kableExtra`](https://rdrr.io/pkg/kableExtra/man/kableExtra-package.html)
object is returned instead so that the table renders inside R Markdown /
Quarto documents. When `format = "text"`, the table is printed to the
console and the raw object is returned invisibly.

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
  This is the classical Fisher (1915, 1921) interval. Under bivariate
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
  4)}\\. This is Bonett and Wright's (2000) Fisher-*z* interval for
  Kendall's \\\tau\\. The constant 0.437 is derived under bivariate
  normality (Bonett & Wright, 2000, equation 4). Requires \\n \ge 5\\;
  for smaller pairwise samples the interval is returned as `NA`.

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

**HTML/LaTeX output.** Uses the kableExtra package (a `Suggests`
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
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_r.md)

Other hypothesis tests:
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`dunnett_ci()`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Worked example using the built-in 'attitude' data set (employee survey of
# 30 departments on 7 ordinal rating items). The goal of correlations_test()
# is to produce a formatted correlation matrix that reports, for every
# variable pair, the correlation, its two-sided p-value, a confidence
# interval on the population correlation, and the pairwise sample size. See
# Kelley (2007) and Maxwell, Delaney, & Kelley (2027, Chapter 9) for
# discussion of why effect sizes should be accompanied by confidence
# intervals.

# Pearson correlations (the default). Each lower-triangle cell stacks r,
# the two-sided p-value, the 95\% confidence interval (Fisher's z
# transformation; Fisher, 1915, 1921), and the pairwise N.
correlations_test(attitude)
#> Correlations (Pearson, 95% CI)
#> 
#>               rating        complaints    privileges    learning      raises        critical      advance       
#> ----------------------------------------------------------------------------------------------------------------
#> rating        -                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#> 
#> complaints    .83           -                                                                                   
#>               p < .0001                                                                                         
#>               [.66, .91]                                                                                        
#>               N = 30                                                                                            
#> 
#> privileges    .43           .56           -                                                                     
#>               p = .0189     p = .0013                                                                           
#>               [.08, .68]    [.25, .76]                                                                          
#>               N = 30        N = 30                                                                              
#> 
#> learning      .62           .60           .49           -                                                       
#>               p = .0002     p = .0005     p = .0056                                                             
#>               [.34, .80]    [.30, .79]    [.16, .72]                                                            
#>               N = 30        N = 30        N = 30                                                                
#> 
#> raises        .59           .67           .45           .64           -                                         
#>               p = .0006     p < .0001     p = .0136     p = .0001                                               
#>               [.29, .78]    [.41, .83]    [.10, .69]    [.36, .81]                                              
#>               N = 30        N = 30        N = 30        N = 30                                                  
#> 
#> critical      .16           .19           .15           .12           .38           -                           
#>               p = .4091     p = .3205     p = .4375     p = .5417     p = .0401                                 
#>               [-.22, .49]   [-.19, .51]   [-.22, .48]   [-.25, .46]   [.02, .65]                                
#>               N = 30        N = 30        N = 30        N = 30        N = 30                                    
#> 
#> advance       .16           .22           .34           .53           .57           .28           -             
#>               p = .4132     p = .2328     p = .0633     p = .0025     p = .0009     p = .1292                   
#>               [-.22, .49]   [-.15, .54]   [-.02, .63]   [.21, .75]    [.27, .77]    [-.09, .58]                 
#>               N = 30        N = 30        N = 30        N = 30        N = 30        N = 30                      
#> 

# Add significance stars and an explanatory footnote.
correlations_test(attitude, stars = TRUE)
#> Correlations (Pearson, 95% CI)
#> 
#>               rating        complaints    privileges    learning      raises        critical      advance       
#> ----------------------------------------------------------------------------------------------------------------
#> rating        -                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#> 
#> complaints    .83***        -                                                                                   
#>               p < .0001                                                                                         
#>               [.66, .91]                                                                                        
#>               N = 30                                                                                            
#> 
#> privileges    .43*          .56**         -                                                                     
#>               p = .0189     p = .0013                                                                           
#>               [.08, .68]    [.25, .76]                                                                          
#>               N = 30        N = 30                                                                              
#> 
#> learning      .62***        .60***        .49**         -                                                       
#>               p = .0002     p = .0005     p = .0056                                                             
#>               [.34, .80]    [.30, .79]    [.16, .72]                                                            
#>               N = 30        N = 30        N = 30                                                                
#> 
#> raises        .59***        .67***        .45*          .64***        -                                         
#>               p = .0006     p < .0001     p = .0136     p = .0001                                               
#>               [.29, .78]    [.41, .83]    [.10, .69]    [.36, .81]                                              
#>               N = 30        N = 30        N = 30        N = 30                                                  
#> 
#> critical      .16           .19           .15           .12           .38*          -                           
#>               p = .4091     p = .3205     p = .4375     p = .5417     p = .0401                                 
#>               [-.22, .49]   [-.19, .51]   [-.22, .48]   [-.25, .46]   [.02, .65]                                
#>               N = 30        N = 30        N = 30        N = 30        N = 30                                    
#> 
#> advance       .16           .22           .34           .53**         .57***        .28           -             
#>               p = .4132     p = .2328     p = .0633     p = .0025     p = .0009     p = .1292                   
#>               [-.22, .49]   [-.15, .54]   [-.02, .63]   [.21, .75]    [.27, .77]    [-.09, .58]                 
#>               N = 30        N = 30        N = 30        N = 30        N = 30        N = 30                      
#> 
#> Note. * p < .05, ** p < .01, *** p < .001.

# Spearman correlations at a 99\% confidence level. The interval uses
# Bonett and Wright's (2000) Fisher-z standard error
# sqrt((1 + r^2/2) / (n - 3)), which corrects the plain Fisher interval
# for the heavier tails of Spearman's sampling distribution.
correlations_test(attitude, method = "spearman", conf_level = 0.99)
#> Correlations (Spearman, 99% CI)
#> 
#>               rating        complaints    privileges    learning      raises        critical      advance       
#> ----------------------------------------------------------------------------------------------------------------
#> rating        -                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#> 
#> complaints    .83           -                                                                                   
#>               p < .0001                                                                                         
#>               [.55, .94]                                                                                        
#>               N = 30                                                                                            
#> 
#> privileges    .48           .53           -                                                                     
#>               p = .0067     p = .0029                                                                           
#>               [.00, .78]    [.05, .80]                                                                          
#>               N = 30        N = 30                                                                              
#> 
#> learning      .62           .58           .51           -                                                       
#>               p = .0003     p = .0008     p = .0041                                                             
#>               [.18, .85]    [.13, .83]    [.03, .80]                                                            
#>               N = 30        N = 30        N = 30                                                                
#> 
#> raises        .60           .65           .46           .62           -                                         
#>               p = .0005     p = .0001     p = .0115     p = .0002                                               
#>               [.15, .84]    [.22, .87]    [-.03, .77]   [.18, .85]                                              
#>               N = 30        N = 30        N = 30        N = 30                                                  
#> 
#> critical      .05           .11           .11           .13           .29           -                           
#>               p = .8003     p = .5541     p = .5473     p = .4970     p = .1245                                 
#>               [-.42, .50]   [-.37, .54]   [-.36, .55]   [-.35, .56]   [-.21, .66]                               
#>               N = 30        N = 30        N = 30        N = 30        N = 30                                    
#> 
#> advance       .20           .22           .34           .54           .49           .26           -             
#>               p = .2821     p = .2338     p = .0680     p = .0021     p = .0061     p = .1736                   
#>               [-.29, .61]   [-.27, .62]   [-.16, .70]   [.07, .81]    [.01, .79]    [-.24, .64]                 
#>               N = 30        N = 30        N = 30        N = 30        N = 30        N = 30                      
#> 

# Kendall's tau, also using Bonett and Wright's (2000) Fisher-z standard
# error sqrt(0.437 / (n - 4)). Kendall is often preferred over Spearman
# for small samples and for samples with many tied ranks.
correlations_test(attitude, method = "kendall")
#> Correlations (Kendall, 95% CI)
#> 
#>               rating        complaints    privileges    learning      raises        critical      advance       
#> ----------------------------------------------------------------------------------------------------------------
#> rating        -                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#> 
#> complaints    .65           -                                                                                   
#>               p < .0001                                                                                         
#>               [.49, .78]                                                                                        
#>               N = 30                                                                                            
#> 
#> privileges    .36           .39           -                                                                     
#>               p = .0066     p = .0032                                                                           
#>               [.12, .56]    [.15, .58]                                                                          
#>               N = 30        N = 30                                                                              
#> 
#> learning      .45           .45           .36           -                                                       
#>               p = .0006     p = .0006     p = .0066                                                             
#>               [.23, .63]    [.23, .63]    [.12, .55]                                                            
#>               N = 30        N = 30        N = 30                                                                
#> 
#> raises        .44           .49           .34           .51           -                                         
#>               p = .0007     p = .0002     p = .0105     p < .0001                                               
#>               [.22, .62]    [.27, .66]    [.09, .54]    [.30, .67]                                              
#>               N = 30        N = 30        N = 30        N = 30                                                  
#> 
#> critical      .03           .08           .09           .08           .20           -                           
#>               p = .8158     p = .5545     p = .5075     p = .5190     p = .1277                                 
#>               [-.22, .28]   [-.17, .32]   [-.17, .33]   [-.17, .33]   [-.05, .43]                               
#>               N = 30        N = 30        N = 30        N = 30        N = 30                                    
#> 
#> advance       .13           .17           .27           .39           .39           .18           -             
#>               p = .3158     p = .1911     p = .0377     p = .0029     p = .0028     p = .1727                   
#>               [-.12, .37]   [-.08, .40]   [.03, .49]    [.16, .58]    [.16, .59]    [-.07, .41]                 
#>               N = 30        N = 30        N = 30        N = 30        N = 30        N = 30                      
#> 

# Save a formatted HTML table that opens directly in a browser
# (then copy into Word). No pandoc required.
tmp_html <- tempfile(fileext = ".html")
correlations_test(attitude, stars = TRUE, format = "html", file = tmp_html)
```
