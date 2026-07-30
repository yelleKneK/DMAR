# Two-Factor Within-Subjects ANOVA With Sphericity Adjustments

Computes the full two-factor within-subjects ANOVA table, main effects
*A*, *B*, and the \\A \times B\\ interaction with each effect tested
against its own residual stratum and with Greenhouse-Geisser,
Huynh-Feldt, and lower-bound sphericity adjustments applied per effect.
Returns a tidy long-form `data.frame` that composes with the rest of
DMAR.

## Usage

``` r
anova_within_two_way(data, outcome, factor_A, factor_B, subject)
```

## Arguments

- data:

  A `data.frame` in long format, with one row per subject-by-cell
  observation.

- outcome:

  Character name of the response column in `data`.

- factor_A:

  Character name of the first within-subjects factor.

- factor_B:

  Character name of the second within-subjects factor.

- subject:

  Character name of the subject-id column.

## Value

A `data.frame` with rows for each of the three effects (`A`, `B`, `A:B`)
crossed with each sphericity adjustment (`none`, `Greenhouse-Geisser`,
`Huynh-Feldt`, `lower_bound`). Columns: `effect`, `adjustment`,
`F_value`, `df_1`, `df_2`, `p_value`, `epsilon`, `partial_eta_squared`.

## Details

**Design.** The two within-subjects factors *A* (with \\a\\ levels) and
*B* (with \\b\\ levels) are fully crossed; every subject contributes \\a
\cdot b\\ observations. Each of the three fixed effects is tested
against its own subject-by-effect residual stratum:

- \\F_A = MS_A / MS\_{A:S}\\

- \\F_B = MS_B / MS\_{B:S}\\

- \\F\_{A:B} = MS\_{A:B} / MS\_{A:B:S}\\

**Sphericity.** Each effect's univariate *F*-ratio assumes sphericity of
its corresponding subject-by-effect residual covariance matrix. Three
adjustments are reported per effect: Greenhouse-Geisser (Greenhouse &
Geisser, 1959), Huynh-Feldt (Huynh & Feldt, 1976), and the lower-bound
\\\epsilon = 1 / (df - 1)\\.

**Per-effect partial \\\eta^2\\.** Computed as \\SS\_\mathrm{effect} /
(SS\_\mathrm{effect} + SS\_\mathrm{effect,\\ error})\\ using the
appropriate subject-by-effect residual sum of squares.

**Balanced data assumed.** The implementation assumes a fully balanced
design (every subject observed once in every cell). When the design is
unbalanced, the function errors and recommends a mixed-effects fit via
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

**Sums of squares are unambiguous here.** Because the design is
balanced, the within-subjects factors are orthogonal and the Type I,
Type II, and Type III sums of squares for each effect coincide. A Type
toggle is therefore not meaningful, and the reported decomposition is
unambiguous: the sum of squares attributed to each effect does not
depend on the order in which terms enter the model (Maxwell, Delaney, &
Kelley, 2027, Chapter 12).

## References

Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of
profile data. *Psychometrika, 24*(2), 95–112.

Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for
degrees of freedom from sample data in randomized block and split-plot
designs. *Journal of Educational Statistics, 1*(1), 69–82.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 12.)

## See also

[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`epsilon_corrections`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md)

Other within-subjects analysis:
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced 3x2 within-subjects design simulated for illustration.
set.seed(113)
n_sub <- 12
grid  <- expand.grid(subject = factor(1:n_sub),
                     A = factor(c("a1", "a2", "a3")),
                     B = factor(c("b1", "b2")))
grid$y <- with(grid,
  3 * (A == "a2") + 1 * (A == "a3") +
  2 * (B == "b2") + 1 * ((A == "a3") & (B == "b2")) +
  rnorm(nrow(grid), 0, 1) +
  rep(rnorm(n_sub, 0, 1.5), times = 6))
anova_within_two_way(grid, outcome = "y", factor_A = "A",
                     factor_B = "B", subject = "subject")
#>    effect         adjustment    F_value     df_1     df_2      p_value
#> 1       A               none 60.2728298 2.000000 22.00000 1.183465e-09
#> 2       A Greenhouse-Geisser 60.2728298 1.830122 20.13135 5.323032e-09
#> 3       A        Huynh-Feldt 60.2728298 2.000000 22.00000 1.183465e-09
#> 4       A        lower_bound 60.2728298 1.000000 11.00000 8.678406e-06
#> 5       B               none 85.6959259 1.000000 11.00000 1.589730e-06
#> 6       B Greenhouse-Geisser 85.6959259 1.000000 11.00000 1.589730e-06
#> 7       B        Huynh-Feldt 85.6959259 1.000000 11.00000 1.589730e-06
#> 8       B        lower_bound 85.6959259 1.000000 11.00000 1.589730e-06
#> 9     A:B               none  0.9731229 2.000000 22.00000 3.935843e-01
#> 10    A:B Greenhouse-Geisser  0.9731229 1.848186 20.33005 3.885934e-01
#> 11    A:B        Huynh-Feldt  0.9731229 2.000000 22.00000 3.935843e-01
#> 12    A:B        lower_bound  0.9731229 1.000000 11.00000 3.451042e-01
#>      epsilon partial_eta_squared
#> 1         NA          0.84566349
#> 2  0.9150612          0.84566349
#> 3  1.0000000          0.84566349
#> 4  0.5000000          0.84566349
#> 5         NA          0.88624133
#> 6  1.0000000          0.88624133
#> 7  1.0000000          0.88624133
#> 8  1.0000000          0.88624133
#> 9         NA          0.08127561
#> 10 0.9240930          0.08127561
#> 11 1.0000000          0.08127561
#> 12 0.5000000          0.08127561
```
