# Moments of the Noncentral *F* Distribution

Returns the mean, variance, standard deviation, skewness, and excess
kurtosis of a noncentral *F* distribution with `df_1` numerator and
`df_2` denominator degrees of freedom and noncentrality parameter `ncp`,
in one table. A moment that does not exist for the given `df_2` is
returned as `NA`.

## Usage

``` r
moments_nc_F(df_1, df_2, ncp = 0)
```

## Arguments

- df_1:

  Numerator degrees of freedom, a single positive number.

- df_2:

  Denominator degrees of freedom, a single positive number.

- ncp:

  Noncentrality parameter \\\lambda\\, a single non-negative number.
  Defaults to 0, the central *F*.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
`mean`, `variance`, `sd`, `skewness`, and `excess_kurtosis` (any of
which may be `NA` when `df_2` is too small), followed by the `df_1`,
`df_2`, and `ncp` that produced them.

## Details

The noncentral *F* is the reference distribution of the *F* statistic
when an effect is present, so its moments describe the sampling behavior
of \\R^2\\, eta squared, and the omnibus *F* test under the alternative.
A central *F* (`ncp = 0`) is the special case.

**The distribution.** The noncentral *F* is \\F =
(X_1/\nu_1)/(X_2/\nu_2)\\ with \\X_1 \sim \chi^2\_{\nu_1}(\lambda)\\ a
noncentral chi square and \\X_2 \sim \chi^2\_{\nu_2}\\ independent of
it. Its raw moments are \$\$\mathrm{E}\[F^k\] =
\Bigl(\tfrac{\nu_2}{\nu_1}\Bigr)^k \mathrm{E}\[X_1^k\]\\
\prod\_{i=1}^{k}\frac{1}{\nu_2 - 2i}, \qquad \nu_2 \> 2k,\$\$ where the
noncentral chi square moments \\\mathrm{E}\[X_1^k\]\\ follow from its
cumulants \\\kappa_j = 2^{j-1}(j-1)!\\(\nu_1 + j\lambda)\\ (Johnson,
Kotz, & Balakrishnan, 1995).

**What is returned.** The mean is \\\nu_2(\nu_1 +
\lambda)/\[\nu_1(\nu_2 - 2)\]\\ and the variance is
\\2(\nu_2/\nu_1)^2\[(\nu_1 + \lambda)^2 + (\nu_1 + 2\lambda)(\nu_2 -
2)\] / \[(\nu_2 - 2)^2(\nu_2 - 4)\]\\. The skewness is the third central
moment over the variance to the power 3/2, and the excess kurtosis is
the fourth central moment over the squared variance, minus 3 (so that a
normal distribution has 0). The mean exists for \\\nu_2 \> 2\\, the
variance for \\\nu_2 \> 4\\, the skewness for \\\nu_2 \> 6\\, and the
excess kurtosis for \\\nu_2 \> 8\\; a moment whose condition is not met
is returned as `NA`, and so is the standard deviation when the variance
is.

**How the moments are computed.** A central moment formed as a
difference of raw moments loses digits when the degrees of freedom are
large, because the raw moments agree in their leading digits: for
\\\nu_1 = \nu_2 = 10^6\\ and \\\lambda = 0\\ that form gives an excess
kurtosis of \\-5.29 \times 10^{-5}\\, where the value is \\6.60006
\times 10^{-5}\\. The function therefore works with the two independent
factors of \\F\\. Write \\F = (\nu_2/\nu_1) X Y\\ with \\X = X_1\\ and
\\Y = 1/X_2\\, and let \\U = X/\mathrm{E}\[X\] - 1\\ and \\V =
Y/\mathrm{E}\[Y\] - 1\\, so that \\F/\mathrm{E}\[F\] - 1 = U + V + UV\\
with \\U\\ and \\V\\ independent and of mean 0. The central moments of
\\F/\mathrm{E}\[F\]\\ are then sums of products of the central moments
of \\U\\ and \\V\\, in which every term is nonnegative and nothing
cancels. With \\n = \nu_2\\, \$\$u_2 = \kappa_2/\kappa_1^2, \quad u_3 =
\kappa_3/\kappa_1^3, \quad c_u = \kappa_4/\kappa_1^4, \quad u_4 = c_u +
3u_2^2,\$\$ \$\$v_2 = \frac{2}{n - 4}, \quad v_3 = \frac{16}{(n - 4)(n -
6)}, \quad c_v = \frac{48(5n - 22)}{(n - 4)^2(n - 6)(n - 8)}, \quad v_4
= c_v + 3v_2^2,\$\$ the second and third central moments and the fourth
cumulant of \\U\\ and of \\V\\ (those of \\V\\ from \\\mathrm{E}\[Y^k\]
= \prod\_{i=1}^{k} 1/(n - 2i)\\). Then \$\$m_2 = u_2 + v_2 + u_2 v_2,
\qquad \mathrm{Var}(F) = \mathrm{E}\[F\]^2\\ m_2,\$\$ \$\$m_3 = u_3 +
v_3 + u_3 v_3 + 3u_3 v_2 + 3u_2 v_3 + 6u_2 v_2, \qquad \mathrm{skewness}
= m_3 / m_2^{3/2},\$\$ and the excess kurtosis is \\e_4 / m_2^2\\, where
\$\$e_4 = c_u + c_v + c_u c_v + 3c_u v_2^2 + 3u_2^2 c_v + 6u_2^2 v_2^2 +
4u_4 v_3 + 4u_3 v_4 + 6c_u v_2 + 12u_2^2 v_2 + 6u_2 c_v + 12u_2 v_2^2 +
12u_3 v_2 + 12u_2 v_3 + 12u_3 v_3\$\$ is the fourth central moment of
\\U + V + UV\\ minus \\3m_2^2\\, with the terms that cancel removed
algebraically. These expressions are the same quantities as the
differences of raw moments above, derived for DMAR; they are not taken
from the reference. Checked against exact rational arithmetic in 14
cases with degrees of freedom up to \\10^{12}\\ and noncentrality up to
\\10^6\\, the four moments are right to a relative error of \\5 \times
10^{-16}\\. At ordinary degrees of freedom (200 or fewer) the values
agree with those of the difference of raw moments to a relative
difference of \\10^{-11}\\.

## References

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions* (Vol. 2, 2nd ed., Chapter 30). Wiley.

## See also

[`moments_nc_t`](https://yelleknek.github.io/DMAR/reference/moments_nc_t.md)
for the noncentral *t*;
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) for
the noncentral *F* confidence limits used in effect size intervals;
[`df`](https://rdrr.io/r/stats/Fdist.html) for the density.

Other noncentral distribution moments:
[`moments_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/moments_nc_chisq.md),
[`moments_nc_t()`](https://yelleknek.github.io/DMAR/reference/moments_nc_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A noncentral F with 3 and 40 df and noncentrality 8.
moments_nc_F(df_1 = 3, df_2 = 40, ncp = 8)
#>  term            value
#>  mean            3.860
#>  variance        5.766
#>  sd              2.401
#>  skewness        1.335
#>  excess_kurtosis 3.121
#>  df_1            3    
#>  df_2            40   
#>  ncp             8.000

# ncp = 0 is the central F: mean df_2 / (df_2 - 2).
moments_nc_F(df_1 = 3, df_2 = 40)
#>  term            value
#>  mean            1.053
#>  variance        0.841
#>  sd              0.917
#>  skewness        1.980
#>  excess_kurtosis 6.625
#>  df_1            3    
#>  df_2            40   
#>  ncp             0.000

# The variance is undefined for four or fewer denominator df.
moments_nc_F(df_1 = 2, df_2 = 4, ncp = 5)
#>  term            value
#>  mean            7.000
#>  variance        <NA> 
#>  sd              <NA> 
#>  skewness        <NA> 
#>  excess_kurtosis <NA> 
#>  df_1            2    
#>  df_2            4    
#>  ncp             5.000
```
