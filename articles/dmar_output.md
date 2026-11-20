# Reading DMAR Result Tables

## One Kind of Result

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

| term     | value |
|:---------|:------|
| ci_lower | 0.101 |
| smd      | 0.500 |
| ci_upper | 0.897 |

Confidence level: 95%

The point estimate sits on the `smd` row; the lower and upper limits of
the 95% interval are on their own rows. The confidence level is printed
beneath the table so you never have to guess what interval you are
looking at.

## What You See Is Rounded; What Is Stored Is Exact

The printed table is a *display*. DMAR rounds the numbers it prints: a
value that is not a whole number shows three decimal places, trailing
zeros included, so every such value in a column shows the same number of
them, and a whole number such as a sample size shows none. It never
rounds the numbers it stores. Pull a value out and you get full
precision:

``` r

x$value[x$term == "smd"]
#> [1] 0.5
x$value[x$term == "ci_upper"]
#> [1] 0.8969414
```

This matters whenever you do further arithmetic. The width of the
interval, for example, uses the stored numbers, not the three decimal
places you saw on screen:

``` r

lo <- x$value[x$term == "ci_lower"]
hi <- x$value[x$term == "ci_upper"]
hi - lo
#> [1] 0.7963557
```

If you want to *see* more (or fewer) decimal places, ask the display for
them. This changes only what is shown, not what is stored:

``` r

print(x, digits = 8)
#>  term     value     
#>  ci_lower 0.10058571
#>  smd      0.50000000
#>  ci_upper 0.89694143
#> 
#> Confidence level: 95%
```

To change the default for the rest of your session, set the option once:

``` r

options(dmar.digits = 4)
x
```

| term     | value  |
|:---------|:-------|
| ci_lower | 0.1006 |
| smd      | 0.5000 |
| ci_upper | 0.8969 |

Confidence level: 95%

``` r

options(dmar.digits = 3)  # back to the default
```

## Getting a Single Number Out

Because the result is an ordinary data frame, you select from it the
usual ways. Any of these returns the upper confidence limit:

``` r

x$value[x$term == "ci_upper"]
#> [1] 0.8969414
x[x$term == "ci_upper", "value"]
#> [1] 0.8969414
subset(x, term == "ci_upper")$value
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

To compare intervals at two confidence levels, stack the `tidy()` rows,
which carry the level in a column;
[`rbind()`](https://rdrr.io/r/base/cbind.html) of two DMAR tables at
different levels returns an ordinary data frame, because no one level
describes every row.
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
does not call the method.

## Wide Tables Read the Same Way

Some functions are naturally wide already: each row is one term, and
several typed columns describe it. The display rules are identical, and
each column is formatted on its own terms. Here every effect in a
two-factor design gets its own row, with a partial effect size, its
confidence limits, the *F* statistic, the degrees of freedom, and the
sample size (the partial value, asked for with `estimand = "partial"`;
by default the table also reports the total value, see
[`?eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)).
The design is from `pygmalion`, a teacher expectancy experiment that
ships with DMAR
([`?pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md)
describes the study and its source): the 310 pupils are classified by
expectancy condition (bloomer or control) and by grade (1 through 6),
and the outcome is each pupil’s change in total IQ from the pretest to
the end of the study.

``` r

ci_eta_squared(aov(iq_gain ~ treatment * factor(grade), data = pygmalion),
               estimand = "partial")
#> Warning: The observed F_value is below the alpha_lower critical value of the
#> central F-distribution, so the lower confidence limit on partial eta squared is
#> 0.
```

| effect | eta_squared_partial | eta_squared_partial_ci_lower | eta_squared_partial_ci_upper | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|
| treatment | 0.014 | 0.000 | 0.050 | 4.273 | 1 | 298 | 310 |
| factor(grade) | 0.125 | 0.051 | 0.181 | 8.486 | 5 | 298 | 310 |
| treatment:factor(grade) | 0.030 | 0.000 | 0.061 | 1.861 | 5 | 298 | 310 |

Confidence level: 95%

Notice that the degrees of freedom and the sample size print as whole
numbers with no decimal point, while the effect sizes and the *F*
statistics print to three decimal places. A confidence limit is never a
count, so a limit that is exactly 0 prints as 0.000, beside the other
limits of its column. The same holds for other measured values that land
on a whole number, such as a CFI or a power of exactly 1, which prints
as 1.000;
[`?dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
says how the display tells a measured value from a count. You did not
have to configure any of that.

The warning printed with the table accounts for the lower limits of
exactly 0: when an effect’s *F* statistic is below the critical value
that sets the lower limit, that limit is 0. The rows come from R’s
sequential (Type I) sums of squares, each term adjusted only for the
terms before it in the formula. Because the cells hold unequal numbers
of pupils, Type II or III sums of squares would give different rows for
the earlier terms;
[`?ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
describes how to supply those *F* statistics instead.

## p-Values and Information Criteria

Two kinds of numbers get special, conventional treatment so they read
the way researchers expect:

- **p-values** print to four decimal places, and a p-value too small to
  show at that precision prints as `< 0.0001` rather than rounding to
  `0.0000`.
- **Information criteria** (AIC, BIC) and log-likelihoods print to a
  number of decimal places of their own (three by default), which
  `digits` does not change, so a model comparison difference of a few
  points is never rounded away.

As always, the stored values keep full precision; only the display is
shaped to the convention. The reference page
[`?dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
documents every display rule and the arguments (`digits`, `digits_p`,
`digits_fixed`) that control them.

## In One Paragraph

Most DMAR results are tidy data frames: rows named by `term`, numbers in
`value`, printed rounded but stored at full precision. Read a number by
indexing the way you always have; see more decimal places with
`print(x, digits = )`; get a wide row with `tidy()` or a one-line
summary with `glance()`. The regression functions for incomplete data
([`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md),
[`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md))
return a fitted model object instead, and the ANOVA effect size
functions
([`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
and their kin) return one row per effect; `tidy()` and `glance()` apply
to these as well.
