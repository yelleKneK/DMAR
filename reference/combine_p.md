# Combine Independent P-Values Across Studies

Combines the one-tailed *p*-values from \\k\\ independent tests of a
common directional hypothesis into a single test, by any of the four
classical methods that Raudenbush (1984) applied to the teacher
expectancy experiments: Fisher's (1938) “adding logs” chi square,
Edgington's (1972) “adding *p*s”, the Mosteller and Bush (1954) “adding
*Z*s” (Stouffer) method, and a weighted adding-Zs variant (weights are
typically the studies' degrees of freedom). Combined significance tests
answer the narrow question “is there an effect in at least some
studies?”; they do not estimate its size. Pair them with
[`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md) or
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) for
estimation, which is almost always the more informative summary.

## Usage

``` r
combine_p(
  p,
  method = c("fisher", "edgington", "stouffer", "stouffer_weighted"),
  weights = NULL
)
```

## Arguments

- p:

  Numeric vector of one-tailed *p*-values, each in (0, 1), oriented so
  that small values support the common directional hypothesis.

- method:

  Character vector naming the methods to compute: any of `"fisher"`,
  `"edgington"`, `"stouffer"`, `"stouffer_weighted"`; the default
  computes all four (the `"stouffer_weighted"` row appears only when
  `weights` is supplied).

- weights:

  Optional non-negative weights for `"stouffer_weighted"`, one per
  study; degrees of freedom are the conventional choice (Mosteller &
  Bush, 1954).

## Value

A `data.frame` (class `dmar_tbl`) with, per requested method, its
statistic row(s) and a one-tailed `<method>_p` row, plus a final `k`
row. The *p* rows print to fixed decimals via the `p_terms` attribute.

## Details

Fisher's statistic is \\-2 \sum \log p_i\\, distributed chi square with
\\2k\\ degrees of freedom under the joint null. Edgington's statistic is
the plain sum \\\sum p_i\\, referred to a normal approximation with mean
\\k/2\\ and variance \\k/12\\ (accurate for \\k \ge 10\\; for smaller
\\k\\ it is conservative in the tails). The Stouffer statistic is \\\sum
z_i / \sqrt{k}\\ with \\z_i = \Phi^{-1}(1 - p_i)\\, and the weighted
variant is \\\sum w_i z_i / \sqrt{\sum w_i^2}\\. All four are reported
with one-tailed combined *p*-values, matching the directional inputs.

Methods can disagree, and the disagreement is informative: Rosenthal
(1978) notes there is no uniformly best test, and Raudenbush (1984)
found three of the four rejecting the null while the df-weighted variant
did not, an early warning that large studies were finding smaller
effects.

## References

Edgington, E. S. (1972). An additive method for combining probability
values from independent experiments. *The Journal of Psychology, 80*(2),
351–363.

Fisher, R. A. (1938). *Statistical methods for research workers* (7th
ed.). Oliver & Boyd.

Mosteller, F., & Bush, R. R. (1954). Selected quantitative techniques.
In G. Lindzey (Ed.), *Handbook of social psychology* (Vol. 1).
Addison-Wesley.

Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
pupil IQ as a function of the credibility of expectancy induction: A
synthesis of findings from 18 experiments. *Journal of Educational
Psychology, 76*(1), 85–97.

Rosenthal, R. (1978). Combining results of independent studies.
*Psychological Bulletin, 85*(1), 185–193.

## See also

[`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md) and
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) for
estimating the pooled effect rather than only testing it;
[`meta_contrast`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md)
for differences among study effects;
[`teacher_expectancy`](https://yelleknek.github.io/DMAR/reference/teacher_expectancy.md)
for the data behind the examples.

Other meta-analysis:
[`meta_contrast()`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md),
[`meta_es()`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
[`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
[`meta_smd()`](https://yelleknek.github.io/DMAR/reference/meta_smd.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Raudenbush (1984), Table 2: the four combined tests over the 18
# teacher expectancy studies (Pellegrini & Hicks at its study-level
# values), weighting the Z method by degrees of freedom.
data(teacher_expectancy)
study <- teacher_expectancy[-c(4, 5), ]
p18  <- append(study$p_one_tailed, .010, after = 3)
df18 <- append(study$n_experimental + study$n_control - 2, 42, after = 3)
combine_p(p18, weights = df18)
#>  term                value 
#>  fisher_chi_square   62.2  
#>  fisher_df           36    
#>  fisher_p            0.0043
#>  edgington_sum_p     7     
#>  edgington_p         0.0509
#>  stouffer_z          2.2   
#>  stouffer_p          0.0139
#>  stouffer_weighted_z 0.87  
#>  stouffer_weighted_p 0.1922
#>  k                   18    
# Fisher chi square 62.17 on 36 df; Edgington sum near 6.9; Stouffer
# z near 2.2; and the df-weighted z under 1: the large studies disagree.
```
