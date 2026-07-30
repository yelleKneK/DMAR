# Provides the Bonferroni-Adjusted Critical Value for an *F* Test of One of Several Contrasts

Provides the Bonferroni-Adjusted Critical Value for an *F* Test of One
of Several Contrasts

## Usage

``` r
cv_bonferroni_f(
  alpha_level = 0.05,
  df_denominator,
  n_comparisons,
  df_numerator = 1,
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  The family-wise Type I error rate (i.e., the rate for the set of
  `n_comparisons` tests taken together). Default `0.05`, the level
  Maxwell, Delaney, and Kelley (2027) tabulate.

- df_denominator:

  The denominator (error) degrees of freedom (a positive number). In a
  one-way design with \\N\\ observations and \\a\\ groups this is \\N -
  a\\.

- n_comparisons:

  The number of comparisons in the family, \\C\\ (a positive integer).

- df_numerator:

  The numerator degrees of freedom. Default `1`, because a contrast
  carries a single degree of freedom, which is the case the Appendix
  table covers.

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value in a output style (a `data.frame` following
the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)). When
`verbose = TRUE` the `area_greater` column reports the per-comparison
error rate \\\alpha/C\\ that the adjustment spends on each test.

## Details

The Bonferroni adjustment tests each of \\C\\ contrasts at \\\alpha/C\\
rather than \\\alpha\\, which holds the family-wise error rate at or
below \\\alpha\\ whatever the contrasts are and however they are
correlated. The critical value is therefore an ordinary upper-tail *F*
quantile read at the smaller per-comparison rate,
\$\$F\_{\alpha/C;\\\mathrm{df\_{num}},\\\mathrm{df\_{den}}},\$\$ which
is what [`cv_f`](https://yelleknek.github.io/DMAR/reference/cv_f.md)
would return if handed `alpha / C`. This function exists because the
adjustment is worth naming: the whole of it is the division, and seeing
\\\alpha/C\\ reported back in `area_greater` is the point. Maxwell,
Delaney, and Kelley (2027) tabulate these values for one numerator
degree of freedom and a family-wise alpha of .05 in their Appendix Table
A.3.

The default `df_numerator = 1` covers the case the table addresses and
the one that arises in practice, since a contrast among means is a
single-degree-of-freedom question. Supply a larger value to Bonferroni
adjust a family of multiple-degree-of-freedom model comparisons.

The procedure is often called Dunn's, after Dunn (1961), who first
applied the Bonferroni inequality to multiple contrasts. It is not the
rank-sum procedure of
[`dunn_test`](https://yelleknek.github.io/DMAR/reference/dunn_test.md),
which the same author published three years later.

**When something else is better.** Bonferroni makes no use of the
structure of the family, so a procedure built for a particular structure
beats it there:
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
is more powerful for all pairwise comparisons, and
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
is more powerful for comparing several treatments to one control.
Bonferroni's advantage is generality; it applies to any set of contrasts
chosen in advance, and it can beat Tukey's method when only a few of the
pairwise comparisons were planned. Because it is conservative, a
step-down variant such as Holm's is uniformly more powerful while
controlling the same rate, and is available through the `method`
argument of
[`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)
and [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## References

Dunn, O. J. (1961). Multiple comparisons among means. *Journal of the
American Statistical Association, 56*(293), 52–64.
[doi:10.2307/2282330](https://doi.org/10.2307/2282330)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem;
Appendix Table A.3 reports these critical values.)

## See also

[`cv_f`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)

Other critical values:
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five planned contrasts with 20 error degrees of freedom, holding the
# family-wise error rate at .05. The area_greater column shows the .01
# per-comparison rate the adjustment spends on each test.
cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5)
#>  term     value area_less area_greater
#>  upper_cv 8.1   0.99      0.01        

# It is the F critical value read at alpha / C.
cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 5,
                verbose = FALSE)$value
#> [1] 8.095958
cv_f(alpha_level = .05 / 5, df_numerator = 1, df_denominator = 20,
     verbose = FALSE)$value[2]
#> [1] 8.095958

# With one comparison there is nothing to adjust.
cv_bonferroni_f(alpha_level = .05, df_denominator = 20, n_comparisons = 1,
                verbose = FALSE)$value
#> [1] 4.351244
```
