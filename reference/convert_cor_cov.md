# Correlation Matrix to Covariance Matrix Conversion

Rescales a correlation matrix into the covariance matrix implied by a
set of standard deviations, the inverse of the standardization that
produces a correlation matrix from a covariance matrix. Useful when a
published article reports correlations and standard deviations but an
analysis needs the covariances.

## Usage

``` r
convert_cor_cov(cor_mat, sd, discrepancy = 1e-05)
```

## Arguments

- cor_mat:

  The correlation matrix to be converted

- sd:

  A vector that contains the standard deviations of the variables in the
  correlation matrix

- discrepancy:

  A small nonnegative tolerance (near 0; default `1e-5`). A value on the
  main diagonal of the correlation matrix is treated as equal to 1 when
  it is within `discrepancy` of 1, that is, when \\\|d - 1\| \le\\
  `discrepancy`

## Value

A square numeric matrix giving the covariance matrix implied by the
supplied correlation matrix and standard deviations, with the same row /
column names as `cor_mat`.

## Details

The correlation matrix to convert can be either symmetric or triangular.
The covariance matrix returned is always a symmetric matrix.

## Note

The correlation matrix input should be a square matrix, and the length
of `sd` should be equal to the number of variables in the correlation
matrix (i.e., the number of rows/columns). Sometimes the correlation
matrix input may not have exactly 1's on the main diagonal, due to,
e.g., rounding; `discrepancy` specifies the allowable discrepancy so
that the function still considers the input as a correlation matrix and
can proceed (but the function does not change the numbers on the main
diagonal).

## See also

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
Cor.Mat <- rbind(c(1.0000, 0.8254, 0.4261, 0.6237, 0.5901, 0.1564, 0.1551),
              c(0.8254, 1.0000, 0.5583, 0.5967, 0.6692, 0.1877, 0.2246),
              c(0.4261, 0.5583, 1.0000, 0.4933, 0.4455, 0.1472, 0.3433),
              c(0.6237, 0.5967, 0.4933, 1.0000, 0.6403, 0.1160, 0.5316),
              c(0.5901, 0.6692, 0.4455, 0.6403, 1.0000, 0.3769, 0.5742),
              c(0.1564, 0.1877, 0.1472, 0.1160, 0.3769, 1.0000, 0.2833),
              c(0.1551, 0.2246, 0.3433, 0.5316, 0.5742, 0.2833, 1.0000))
colnames(Cor.Mat) <- rownames(Cor.Mat) <- c("rating", "complaints", "privileges",
"learning", "raises", "critical", "advance")

SDs <- c(12.172562, 13.314757, 12.235430, 11.737013, 10.397226, 9.894908, 10.288706)
convert_cor_cov(cor_mat=Cor.Mat, sd=SDs)
#>               rating complaints privileges  learning    raises critical
#> rating     148.17127  133.77646   63.46186  89.10772  74.68357 18.83781
#> complaints 133.77646  177.28275   90.95365  93.24958  92.64173 24.72916
#> privileges  63.46186   90.95365  149.70575  70.84153  56.67407 17.82128
#> learning    89.10772   93.24958   70.84153 137.75747  78.13733 13.47185
#> raises      74.68357   92.64173   56.67407  78.13733 108.10231 38.77532
#> critical    18.83781   24.72916   17.82128  13.47185  38.77532 97.90920
#> advance     19.42471   30.76832   43.21692  64.19531  61.42447 28.84158
#>              advance
#> rating      19.42471
#> complaints  30.76832
#> privileges  43.21692
#> learning    64.19531
#> raises      61.42447
#> critical    28.84158
#> advance    105.85747
```
