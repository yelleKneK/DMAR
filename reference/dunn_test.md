# Provides Dunn's Rank-Sum Test of All Pairwise Differences Following a Kruskal–Wallis Test

Provides Dunn's Rank-Sum Test of All Pairwise Differences Following a
Kruskal–Wallis Test

## Usage

``` r
dunn_test(x, group = NULL, method = "holm")
```

## Arguments

- x:

  Either (a) a numeric vector of the outcome, in which case `group` must
  also be supplied, or (b) a fitted
  [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a single factor
  predictor, from which the outcome and grouping factor are taken (the
  model itself is not used, since the test is distribution free).

- group:

  When `x` is a vector, a factor (or coercible to factor) giving the
  group membership of each observation.

- method:

  The multiplicity adjustment applied to the pairwise *p*-values, passed
  to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"holm"`. Use `"none"` to obtain the unadjusted values.

## Value

A `data.frame` with one row per pairwise comparison and columns
`contrast`, `mean_rank_difference`, `se`, `z_statistic`, `p_value`, and
`p_adjusted`. The table prints through the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
display layer.

## Details

**Which Dunn.** Olive Jean Dunn published two different multiple
comparison procedures, and both are called “Dunn's” in the literature.
This function implements the *nonparametric* one of Dunn (1964), which
compares mean ranks. It is not the Bonferroni procedure of Dunn (1961),
which Maxwell, Delaney, and Kelley (2027, Chapter 5) call Dunn's
procedure and which reaches DMAR through the `method` arguments of
[`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)
and [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). The two are
unrelated apart from their author, and the `method` argument here can
apply the 1961 procedure to the 1964 procedure's *p*-values.

**What it does.** The Kruskal–Wallis test asks whether *any* of the
groups differ; it does not say which. Dunn's test is its pairwise
follow-up. All \\N\\ observations are ranked together, with tied values
receiving their average rank. Writing \\\bar R_g\\ for the mean rank of
group \\g\\, each pair is compared with \$\$z = \frac{\bar R_g - \bar
R_h}{\sqrt{\left\[\frac{N(N+1)}{12} - \frac{\sum_i (t_i^3 -
t_i)}{12(N-1)}\right\]\left(\frac{1}{n_g} + \frac{1}{n_h}\right)}},\$\$
where \\t_i\\ is the number of observations tied at the \\i\\th distinct
value; the second term in the brackets is the tie correction and
vanishes when there are no ties. The statistic is referred to the
standard normal distribution, and the resulting *p*-values are adjusted
for multiplicity by `method`.

The pooled ranking is what makes this the right follow-up. Dunn's test
ranks across all groups at once and uses the variance of the ranks
implied by the Kruskal–Wallis null, so it is consistent with the omnibus
test that preceded it. Running a separate Mann–Whitney test on each pair
instead re-ranks the data within every pair, which answers a different
question for every comparison and is not coherent with the omnibus
result.

**When to use it.** Use Dunn's test when the outcome is ordinal, or when
it is continuous but the normality or homogeneity assumptions behind
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
and
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
are untenable and a rank-based analysis is preferred to a
transformation; the design is between subjects; and a significant
Kruskal–Wallis test leaves the question of which groups differ. It is
the rank analogue of Tukey's HSD, in the sense of covering all
\\a(a-1)/2\\ pairs.

**What it does not tell you.** The test compares mean *ranks*, not
medians. A significant pair means one group's observations tend to be
larger, that is, stochastic dominance; it does not by itself license a
statement about medians unless the group distributions have the same
shape. It also returns no confidence interval on any quantity in the
original units, which is a real cost: where a parametric procedure is
defensible,
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
or
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
reports intervals on the mean difference, and Maxwell, Delaney, and
Kelley (2027) emphasize interval estimation over test decisions
throughout. For a distribution free effect size with a confidence
interval, see
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md).

## References

Dunn, O. J. (1964). Multiple comparisons using rank sums.
*Technometrics, 6*(3), 241–252.
[doi:10.1080/00401706.1964.10490181](https://doi.org/10.1080/00401706.1964.10490181)

Dunn, O. J. (1961). Multiple comparisons among means. *Journal of the
American Statistical Association, 56*(293), 52–64.
[doi:10.2307/2282330](https://doi.org/10.2307/2282330) (The Bonferroni
procedure; not what this function computes.)

Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion
variance analysis. *Journal of the American Statistical Association,
47*(260), 583–621.
[doi:10.2307/2280779](https://doi.org/10.2307/2280779)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html) for the
omnibus test this follows,
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
and
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
for the parametric all-pairs procedures,
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
for a distribution free effect size with a confidence interval, and
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for the `method`
choices.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The omnibus question first: do the three treatment arms of the
# drinks_trial data differ at all? The raw drinks_per_week outcome is
# markedly right-skewed, which is no obstacle here: ranks are invariant
# to a monotone transformation, so the raw and log scales give
# identical results.
kruskal.test(drinks_per_week ~ treatment, data = drinks_trial)
#> 
#>  Kruskal-Wallis rank sum test
#> 
#> data:  drinks_per_week by treatment
#> Kruskal-Wallis chi-squared = 7.2491, df = 2, p-value = 0.02666
#> 

# Then which pairs, on the same pooled ranking the omnibus test used.
dunn_test(drinks_trial$drinks_per_week, drinks_trial$treatment)
#>  contrast                    mean_rank_difference se   z_statistic p_value
#>  CRA - Standard              -13                  5.99 -2.18       0.0294 
#>  CRA + Disulfiram - Standard -16.2                7    -2.32       0.0205 
#>  CRA + Disulfiram - CRA      -3.18                7.18 -0.442      0.6583 
#>  p_adjusted
#>  0.0615    
#>  0.0615    
#>  0.6583    

# Without a multiplicity adjustment (rarely what you want).
dunn_test(drinks_trial$drinks_per_week, drinks_trial$treatment,
          method = "none")
#>  contrast                    mean_rank_difference se   z_statistic p_value
#>  CRA - Standard              -13                  5.99 -2.18       0.0294 
#>  CRA + Disulfiram - Standard -16.2                7    -2.32       0.0205 
#>  CRA + Disulfiram - CRA      -3.18                7.18 -0.442      0.6583 
#>  p_adjusted
#>  0.0294    
#>  0.0205    
#>  0.6583    
```
