# Convert Between an *F* Value and a Chi Square Value

Converts an observed *F* value into a chi square value, and a chi square
value into an *F* value. Both functions return the converted statistic
itself, a single number, not a *p*-value.

Two conversions are available, selected by `df_denominator`.

- **Scaling (the default).** With `df_denominator = Inf`,
  `convert_F_chisq()` returns `df_numerator * F_value` and
  `convert_chisq_F()` returns `chi_square / df`. This is the standard
  conversion between the two test statistics and needs no denominator
  degrees of freedom.

- **Probability matching.** With a finite `df_denominator`,
  `convert_F_chisq()` returns the chi square value that has the same
  upper-tail probability (the same *p*-value) as the *F* value, and
  `convert_chisq_F()` returns the *F* value with the same upper-tail
  probability as the chi square value.

## Usage

``` r
convert_F_chisq(F_value, df_numerator, df_denominator = Inf)

convert_chisq_F(chi_square, df, df_denominator = Inf)
```

## Arguments

- F_value:

  Observed *F* value. Must be nonnegative.

- df_numerator:

  Numerator degrees of freedom of the *F*, which is also the degrees of
  freedom of the chi square.

- df_denominator:

  Denominator degrees of freedom of the *F*, the degrees of freedom on
  which the error variance is estimated. The default, `Inf`, treats the
  error variance as known and gives the scaling conversion (\\\chi^2 =
  \nu_1 F\\); a finite value gives the probability-matching conversion.
  See *Details*.

- chi_square:

  Observed chi square value. Must be nonnegative.

- df:

  Degrees of freedom of the chi square, which becomes the numerator
  degrees of freedom of the *F*.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` is
`"chi_square_from_F"` for `convert_F_chisq` and `"F_from_chi_square"`
for `convert_chisq_F`, and `value` is the converted statistic (a chi
square value or an *F* value, respectively; never a *p*-value).

## Details

**Why there are two conversions.** An *F* statistic is the ratio of two
independent chi squares, each divided by its degrees of freedom,
\$\$F(\nu_1, \nu_2) =
\frac{\chi^2\_{\nu_1}/\nu_1}{\chi^2\_{\nu_2}/\nu_2},\$\$ where the
denominator is the estimated error variance scaled to have a mean of 1.
A chi square is what the numerator becomes when that error variance is
known rather than estimated. This is the entire difference between the
two, and it is why the conversion depends on how the error variance is
treated.

**Scaling: `df_denominator = Inf`.** As the denominator degrees of
freedom grow, the estimated error variance converges to the true one and
\\\nu_1 F \to \chi^2(\nu_1)\\. The default therefore treats the error
variance as known and returns exactly \$\$\chi^2 = \nu_1 \\ F, \qquad F
= \chi^2 / \nu_1,\$\$ with \\\nu_1\\ the numerator degrees of freedom
(`df_numerator`, which is also the degrees of freedom of the chi
square). This is the value an *F* table prints in its
infinite-denominator row, and it is the usual conversion between a Wald
*F* and a Wald chi square. It involves no probabilities and needs no
`df_denominator`.

**Probability matching: finite `df_denominator`.** When the error
variance is estimated on \\\nu_2\\ degrees of freedom, the scaling above
runs high, because \\\nu_1 F\\ is more dispersed than \\\chi^2(\nu_1)\\.
Supplying `df_denominator` returns instead the chi square value at the
same upper-tail probability. Writing `pf` and `qchisq` for R's
distribution and quantile functions, the computation is exactly

    p <- pf(F_value, df_numerator, df_denominator, lower.tail = FALSE)
    chi_square <- qchisq(p, df_numerator, lower.tail = FALSE)

and `convert_chisq_F()` composes the same two functions in the other
order. The upper tail is used so the *p*-value is represented accurately
for large statistics (the lower-tail probability rounds to 1 in double
precision by about \\F = 500\\ at small \\\nu_2\\, which would send
[`qchisq()`](https://rdrr.io/r/stats/Chisquare.html) to infinity; the
upper-tail value stays accurate past \\F = 10^{20}\\). No logarithms are
involved and the returned value is the chi square itself, not the
*p*-value used to find it. The map is strictly increasing and therefore
one to one, and `convert_chisq_F()` is its exact inverse.

**How much the two differ.** At \\\nu_1 = 3\\ and \\F = 2.75\\,
probability matching gives \\\chi^2 = 6.290\\ at \\\nu_2 = 10\\,
\\7.711\\ at \\\nu_2 = 50\\, and \\8.220\\ at \\\nu_2 = 1000\\,
approaching the scaling value \\\nu_1 F = 8.25\\ only as \\\nu_2\\
grows. Scaling and probability matching agree in the limit and diverge
as \\\nu_2\\ shrinks; for a small \\\nu_2\\, the scaled value
understates the *p*-value (with \\\nu_1 = 3, \nu_2 = 5\\, a result whose
true *p* is .05 reads as .001 if \\\nu_1 F\\ is referred to \\\chi^2\\).

**Noncentrality.** Both conversions are defined by the central
distributions and preserve the *p*-value; neither transports a
noncentrality parameter (a noncentral *F* is not carried to a noncentral
chi square with the same \\\lambda\\, except in the \\\nu_2 \to \infty\\
limit). For noncentral work use
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
and
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md).

**Special case.** With \\\nu_1 = 1\\ this is the squared form of the
relation between *t* and *z*, since \\F(1, \nu) = t(\nu)^2\\ and
\\\chi^2(1) = z^2\\.

## References

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions* (2nd ed., Vol. 2). Wiley.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`cv_f`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_chisq`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)

Other parameterization conversions:
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Scaling (the default): chi square = df_numerator * F.
convert_F_chisq(2.75, df_numerator = 3)
#>  term              value
#>  chi_square_from_F 8.25 

# The inverse: F = chi square / df.
convert_chisq_F(8.25, df = 3)
#>  term              value
#>  F_from_chi_square 2.75 

# Probability matching at 3 and 50 degrees of freedom: the chi square
# with the same p-value as the F.
convert_F_chisq(2.75, df_numerator = 3, df_denominator = 50)
#>  term              value
#>  chi_square_from_F 7.71 

# The p-value is preserved, which is what probability matching means.
f <- 2.75
x <- convert_F_chisq(f, df_numerator = 3, df_denominator = 50)$value
c(p_from_F = pf(f, 3, 50, lower.tail = FALSE),
  p_from_chi_square = pchisq(x, 3, lower.tail = FALSE))
#>          p_from_F p_from_chi_square 
#>        0.05238213        0.05238213 
```
