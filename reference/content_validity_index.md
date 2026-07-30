# Content Validity Index From Expert Ratings

Quantifies how well a pool of candidate items covers the construct it is
meant to measure, using the relevance ratings of a panel of subject
matter experts. Validation begins before any respondent data are
collected: experts rate each item for relevance, and the content
validity index summarizes their agreement (Lynn, 1986; Polit & Beck,
2006; Bandalos, 2018). The function returns the item level index (I-CVI)
with an exact binomial confidence interval, the chance corrected
modified kappa of Polit, Beck, and Owen (2007), Lawshe's (1975) content
validity ratio, and the two scale level summaries S-CVI/Ave and
S-CVI/UA. The confidence interval is what keeps a small panel from being
read as more informative than it is: with five experts an I-CVI of 0.80
carries an interval roughly seven tenths of the width of the scale.

## Usage

``` r
content_validity_index(
  ratings,
  relevant = c(3, 4),
  essential = NULL,
  conf_level = 0.95
)
```

## Arguments

- ratings:

  A matrix or `data.frame` of expert ratings with items in rows and
  experts in columns, on the conventional 4 point relevance scale (1 =
  not relevant, 2 = somewhat relevant, 3 = quite relevant, 4 = highly
  relevant). Every non-missing entry must be one of 1, 2, 3, or 4. Row
  names, when present, name the items; otherwise items are labeled
  `item_1`, `item_2`, and so on. At least 1 item and at least 2 experts
  are required. `NA` is allowed and means that expert did not rate that
  item; a missing entry lowers that item's expert count and thereby the
  denominator of its I-CVI, kappa, and CVR, leaving other items
  untouched. A row with no ratings at all is an error.

- relevant:

  The rating values counted as relevant. Default `c(3, 4)`, the standard
  dichotomization of the 4 point scale into relevant (3 or 4) versus not
  relevant (1 or 2). Must be a subset of `1:4`.

- essential:

  Optional. The rating values counted as “essential” for Lawshe's
  content validity ratio, when the panel answered the
  essential-versus-not question on a separate part of the scale. Must be
  a subset of `1:4`. When `NULL` (default), the content validity ratio
  is computed from the same dichotomization as relevance, that is from
  `relevant`.

- conf_level:

  Confidence level for the exact binomial interval on each I-CVI.
  Default `0.95`. Must be in (0, 1).

## Value

A `data.frame` (class `dmar_tbl`) with one row per item, in the order
the items appear in `ratings`, and columns:

- `item`:

  The item name, from the row names of `ratings` when present.

- `n_experts`:

  Number of experts who rated the item (missing ratings excluded).

- `n_relevant`:

  Number of those experts whose rating was in `relevant`.

- `i_cvi`:

  The item level content validity index, `n_relevant / n_experts`.

- `ci_lower`, `ci_upper`:

  Limits of the exact binomial confidence interval for `i_cvi` at
  `conf_level`.

- `kappa`:

  The modified kappa of Polit, Beck, and Owen (2007), the I-CVI adjusted
  for chance agreement.

- `cvr`:

  Lawshe's content validity ratio, computed from `essential` when
  supplied and from `relevant` otherwise.

Attributes: `"s_cvi_ave"`, the mean of the `i_cvi` column; `"s_cvi_ua"`,
the proportion of items with `i_cvi` equal to 1; `"relevant"`, the
rating values counted as relevant; `"essential"`, the rating values
counted as essential for the content validity ratio (equal to
`"relevant"` when `essential` was `NULL`); and `"conf_level"`, the
confidence level used.

## Details

Let \\N\\ be the number of experts who rated an item and \\A\\ the
number of those who rated it relevant.

**Item level index.** The item level content validity index is the
proportion of rating experts who called the item relevant,
\$\$\mathrm{I\mbox{-}CVI} = A / N.\$\$ Lynn (1986) gives the
conventional criteria: with five or fewer experts an item is expected to
reach 1.00, and with six to ten experts at least 0.78.

**Confidence interval.** The I-CVI is a binomial proportion, so the
interval reported here is the exact (Clopper-Pearson) interval from
[`binom.test`](https://rdrr.io/r/stats/binom.test.html) at `conf_level`.
Expert panels are small by design and the interval states plainly how
little a handful of ratings pins down the population proportion. This is
the accuracy in parameter estimation view of content validity: if the
interval is too wide to act on, the remedy is more experts.

**Chance corrected agreement.** Some of the observed agreement on
relevance would occur if experts responded at random with probability
0.5, so Polit, Beck, and Owen (2007) correct the index in the manner of
a kappa. The probability of chance agreement is the binomial point
probability \$\$p_c = \frac{N!}{A!\\(N - A)!}\\ 0.5^{N},\$\$ and the
modified kappa is \$\$\kappa = \frac{\mathrm{I\mbox{-}CVI} - p_c}{1 -
p_c}.\$\$ The binomial coefficient is evaluated on the log scale with
[`lchoose`](https://rdrr.io/r/base/Special.html) and then exponentiated,
so a large panel does not overflow the way
[`factorial()`](https://rdrr.io/r/base/Special.html) would.

**Content validity ratio.** Lawshe (1975) asked a panel whether each
item measures behavior that is essential to the performance domain. With
\\n_e\\ experts calling the item essential, \$\$\mathrm{CVR} =
\frac{n_e - N / 2}{N / 2},\$\$ which equals 1 when every expert says
essential, 0 when exactly half do, and -1 when none do.

**Scale level summaries.** S-CVI/Ave is the mean of the I-CVIs over
items, the averaging approach Polit and Beck (2006) recommend reporting.
S-CVI/UA is the universal agreement proportion, the fraction of items
whose I-CVI equals 1, a stricter and considerably more conservative
summary.

## References

Bandalos, D. L. (2018). *Measurement theory and applications for the
social sciences*. Guilford Press.

Lawshe, C. H. (1975). A quantitative approach to content validity.
*Personnel Psychology, 28*(4), 563–575.

Lynn, M. R. (1986). Determination and quantification of content
validity. *Nursing Research, 35*(6), 382–385.

Polit, D. F., & Beck, C. T. (2006). The content validity index: Are you
sure you know what's being reported? Critique and recommendations.
*Research in Nursing and Health, 29*(5), 489–497.
[doi:10.1002/nur.20147](https://doi.org/10.1002/nur.20147)

Polit, D. F., Beck, C. T., & Owen, S. V. (2007). Is the CVI an
acceptable indicator of content validity? Appraisal and recommendations.
*Research in Nursing and Health, 30*(4), 459–467.
[doi:10.1002/nur.20199](https://doi.org/10.1002/nur.20199)

## See also

[`gwet_ac`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
for agreement among raters on a common set of units.

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`loa()`](https://yelleknek.github.io/DMAR/reference/loa.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Six experts rate five candidate items on the 4 point relevance scale.
ratings <- rbind(
  item_1 = c(4, 4, 3, 4, 4, 3),
  item_2 = c(4, 3, 4, 4, 3, 2),
  item_3 = c(2, 3, 1, 2, 3, 2),
  item_4 = c(4, 4, 4, 4, 4, 4),
  item_5 = c(3, 4, 4, 3, NA, 4))
colnames(ratings) <- paste0("expert_", 1:6)
cvi <- content_validity_index(ratings)
cvi
#>  item   n_experts n_relevant i_cvi ci_lower ci_upper kappa cvr   
#>  item_1 6         6          1     0.541    1        1     1     
#>  item_2 6         5          0.833 0.359    0.996    0.816 0.667 
#>  item_3 6         2          0.333 0.0433   0.777    0.129 -0.333
#>  item_4 6         6          1     0.541    1        1     1     
#>  item_5 5         5          1     0.478    1        1     1     
#> 
#> Confidence level: 95%

# The scale level summaries travel with the table as attributes.
attr(cvi, "s_cvi_ave")
#> [1] 0.8333333
attr(cvi, "s_cvi_ua")
#> [1] 0.6

# The worked example of Polit, Beck, and Owen (2007): 6 experts, 5 of whom
# rate the item relevant, gives I-CVI = 0.83, p_c = 0.094, and a modified
# kappa of 0.816 (the paper reports 0.81, carrying its rounded I-CVI).
content_validity_index(matrix(c(4, 4, 3, 4, 3, 1), nrow = 1))
#>  item   n_experts n_relevant i_cvi ci_lower ci_upper kappa cvr  
#>  item_1 6         5          0.833 0.359    0.996    0.816 0.667
#> 
#> Confidence level: 95%

# A wide interval is the point: with 5 experts an I-CVI of 0.80 is
# compatible with a population proportion anywhere from about 0.28 to 0.99.
content_validity_index(matrix(c(4, 4, 3, 4, 1), nrow = 1))
#>  item   n_experts n_relevant i_cvi ci_lower ci_upper kappa cvr
#>  item_1 5         4          0.8   0.284    0.995    0.763 0.6
#> 
#> Confidence level: 95%

# The broom verbs on the earlier result: one row per item, and the
# scale-level summary.
generics::tidy(cvi)
#>     term  estimate   ci_lower  ci_upper     kappa        cvr n_experts
#> 1 item_1 1.0000000 0.54074187 1.0000000 1.0000000  1.0000000         6
#> 2 item_2 0.8333333 0.35876542 0.9957893 0.8160920  0.6666667         6
#> 3 item_3 0.3333333 0.04327187 0.7772219 0.1292517 -0.3333333         6
#> 4 item_4 1.0000000 0.54074187 1.0000000 1.0000000  1.0000000         6
#> 5 item_5 1.0000000 0.47817625 1.0000000 1.0000000  1.0000000         5
#>   n_relevant
#> 1          6
#> 2          5
#> 3          2
#> 4          6
#> 5          5
generics::glance(cvi)
#>   s_cvi_ave s_cvi_ua n_items n_experts conf_level
#> 1 0.8333333      0.6       5         6       0.95
```
