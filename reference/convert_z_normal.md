# Convert a Standard Normal *z* Value to the Corresponding Value on a Normal Distribution

This function maps a value on the standard normal distribution (the
*z*-distribution, with mean 0 and variance 1) to the equivalent point on
a normal distribution with arbitrary mean and standard deviation,
\\N(mean, sd^2)\\.

## Usage

``` r
convert_z_normal(z, mean = 0, sd = 1)
```

## Arguments

- z:

  A value on the standard normal distribution (with mean 0 and variance
  1).

- mean:

  The mean of the target normal distribution.

- sd:

  The standard deviation of the target normal distribution.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` is
`"value_from_z"` and `value` is the point on \\N(mean, sd^2)\\ that lies
at the same percentile as `z` does on the standard normal distribution.

## Details

The conversion is `value = mean + z * sd`, which places the returned
value at the same percentile of \\N(mean, sd^2)\\ that `z` occupies on
the standard normal distribution. Equivalently,
`value = qnorm(pnorm(z), mean, sd)`. With the defaults (`mean = 0`,
`sd = 1`) the value is returned unchanged, since the target distribution
is then the standard normal distribution itself.

## See also

[`cv_z`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A z value of 1.96 on the standard normal distribution maps to the
# corresponding point on a normal distribution with mean 100 and sd 15.
convert_z_normal(z = 1.96, mean = 100, sd = 15)
#>  term         value
#>  value_from_z 129  

# With the default standard normal target, the value is returned unchanged.
convert_z_normal(z = 1.96)
#>  term         value
#>  value_from_z 1.96 
```
