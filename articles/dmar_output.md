# Reading DMAR Result Tables

## One Kind of Result, Everywhere

Almost every DMAR function hands back the same kind of object: a tidy
table with a `term` column that names each quantity and a `value` column
that holds it. You do not need to learn a new object for each analysis,
and you do not need to know anything about R’s class systems to use what
comes back. If you can work with a `data.frame`, you can work with a
DMAR result.

Here is a confidence interval for a standardized mean difference
(Cohen’s *d*), computed from summary statistics you might read out of a
paper:

``` r

x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
x
```

| term        | value |
|:------------|:------|
| lower_limit | 0.101 |
| smd         | 0.5   |
| upper_limit | 0.897 |

Confidence level: 95%

The point estimate sits on the `smd` row; the lower and upper limits of
the 95% interval are on their own rows. The confidence level is printed
beneath the table so you never have to guess what interval you are
looking at.

## What You See Is Rounded; What Is Stored Is Exact

The printed table is a *display*. DMAR rounds numbers so the table is
easy to read, but it never rounds the numbers it stores. Pull a value
out and you get full precision:

``` r

x$value[x$term == "smd"]
#> [1] 0.5
x$value[x$term == "upper_limit"]
#> [1] 0.8969414
```

This matters whenever you do further arithmetic. The width of the
interval, for example, uses the stored numbers, not the three digits you
saw on screen:

``` r

lo <- x$value[x$term == "lower_limit"]
hi <- x$value[x$term == "upper_limit"]
hi - lo
#> [1] 0.7963557
```

If you want to *see* more (or fewer) digits, ask the display for them.
This changes only what is shown, not what is stored:

``` r

print(x, digits = 8)
#>  term        value     
#>  lower_limit 0.10058571
#>  smd         0.5       
#>  upper_limit 0.89694143
#> 
#> Confidence level: 95%
```

To change the default for the rest of your session, set the option once:

``` r

options(dmar.digits = 4)
x
```

| term        | value  |
|:------------|:-------|
| lower_limit | 0.1006 |
| smd         | 0.5    |
| upper_limit | 0.8969 |

Confidence level: 95%

``` r

options(dmar.digits = 3)  # back to the default
```

## Getting a Single Number Out

Because the result is an ordinary data frame, you select from it the
usual ways. Any of these returns the upper confidence limit:

``` r

x$value[x$term == "upper_limit"]
#> [1] 0.8969414
x[x$term == "upper_limit", "value"]
#> [1] 0.8969414
subset(x, term == "upper_limit")$value
#> [1] 0.8969414
```

That is the whole trick: rows are named by `term`, numbers live in
`value`, and you index them like any data frame.

## The “Wide” View and the One-Row Summary

Sometimes you want the result spread across columns instead of down
rows, for example to add a row to a results table or to feed a plot. The
broom verbs `tidy()` and `glance()` do this, and they ship with DMAR
through the lightweight **generics** package, so nothing extra needs to
be installed.

`tidy()` returns one row per term with tidy, predictable column names:

``` r

generics::tidy(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
```

`glance()` returns a one-row, model-level summary. For a result that is
already a single estimate and its interval, such as
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md), the
summary is the same one row, because there are no extra model-level
statistics to add:

``` r

generics::glance(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
```

The two verbs come apart when there is more to say at the model level. A
regression fit by
[`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md), for
instance, gives one `tidy()` row per coefficient and a `glance()` row of
fit statistics (`R2`, AIC, BIC, and so on). The rule is the same
everywhere, including the wide tables whose rows are items, construct
pairs, ladder rungs, or groups
([`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
and the rest): `tidy()` is the per-term table, `glance()` is the
one-line summary. The columns use the same names as every other DMAR
surface (`estimate`, `se`, `p_value`, `ci_lower`, `ci_upper`), so
nothing needs translating between the function’s own table and its
tidied view.

This is a deliberate design choice. A DMAR function always returns the
same shape no matter how you call it, so a script that reads `x$value`
keeps working. When you want a different shape, you ask for it with a
verb (`tidy()` or `glance()`); the function’s own output never changes
shape underneath you.

## Wide Tables Read the Same Way

Some functions are naturally wide already: each row is one term, and
several typed columns describe it. The display rules are identical, and
each column is formatted on its own terms. Here every effect in a
two-factor design gets its own row, with a partial effect size, its
confidence limits, the *F* statistic, the degrees of freedom, and the
sample size:

``` r

ci_eta_squared(aov(len ~ supp * factor(dose), data = ToothGrowth))
```

| effect | eta_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|
| supp | 0.224 | 0.0531 | 0.377 | 15.6 | 1 | 54 | 60 |
| factor(dose) | 0.773 | 0.638 | 0.823 | 92 | 2 | 54 | 60 |
| supp:factor(dose) | 0.132 | 0.00132 | 0.274 | 4.11 | 2 | 54 | 60 |

Notice that the degrees of freedom and the sample size print as whole
numbers with no decimal point, while the effect sizes and the *F*
statistics print to three significant figures. You did not have to
configure any of that.

## p-Values and Information Criteria

Two kinds of numbers get special, conventional treatment so they read
the way researchers expect:

- **p-values** print to four decimal places, and a p-value too small to
  show at that precision prints as `< 0.0001` rather than rounding to
  `0.0000`.
- **Information criteria** (AIC, BIC) and log-likelihoods print to a
  fixed number of decimal places rather than to significant figures, so
  a model comparison difference of a few points is never rounded away.

As always, the stored values keep full precision; only the display is
shaped to the convention. The reference page
[`?dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
documents every display rule and the arguments (`digits`, `digits_p`,
`digits_fixed`) that control them.

## In One Paragraph

Every DMAR result is a tidy data frame: rows named by `term`, numbers in
`value`, printed with sensible rounding but stored at full precision.
Read a number by indexing the way you always have; see more digits with
`print(x, digits = )`; get a wide row with `tidy()` or a one-line
summary with `glance()`. That uniformity is the point. Learn it once and
it holds for the whole package.
