# Controlled Test-Market Experiment (Bryant & Bruvold, 1980)

The controlled test-market experiment of Bryant and Bruvold (1980), used
to illustrate multiple-comparison procedures in the analysis of
covariance (ANCOVA) when the covariate is *random*. A company compared
\\k = 6\\ marketing strategies (“panels”) for a brand, randomly
assigning them to retail outlets within \\s = 4\\ blocks of outlets that
were homogeneous in size, locality, and ownership (a randomized complete
block design, one outlet per panel-by-block cell). During the experiment
a concomitant variable – the remaining category movement in each outlet
– becomes available; it cannot be controlled by the experimenter and is
best modeled as a random covariate. Adjusting brand movement for this
covariate sharply reduces unexplained error, permitting far finer
comparison of the panels than the raw outcome allows.

## Usage

``` r
test_market
```

## Format

A data frame with 24 observations (6 panels \\\times\\ 4 blocks) on 4
variables.

- `panel`:

  Factor with levels `1`–`6`: the marketing strategy (treatment)
  randomly assigned to the outlet. Different panels entail different
  methods of packaging, displaying, or pricing.

- `block`:

  Factor with levels `1`–`4`: the block of retail outlets, grouped to be
  homogeneous in size, locality, ownership, and other considerations
  that influence brand movement.

- `brand_movement`:

  Test-brand movement during the test period, in hundreds of statistical
  cases. The dependent variable (\\y\\ in the source).

- `category_movement`:

  Remaining category movement, the random concomitant variable
  (covariate; \\x\\ in the source). It is not identically distributed
  across blocks, which is precisely the setting Bryant and Bruvold's
  grouped-covariate extension was designed for.

## Source

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880 (Table 1).
[doi:10.2307/2287175](https://doi.org/10.2307/2287175)

## Details

**Why it is a benchmark for ANCOVA multiple comparisons.** The model
fitted by Bryant and Bruvold (their Eq. 3.1) is a randomized-block
ANCOVA, \$\$y\_{ij} = \theta_i + \beta_j + (x\_{ij} - \delta_j)\\ u +
e\_{ij},\$\$ with \\\theta_i\\ the \\i\\th panel (adjusted) mean,
\\\beta_j\\ the \\j\\th block effect, and \\u\\ the within-cell
covariate slope. The point of the example is that the studentized range
of the adjusted panel means does *not* follow the ordinary Tukey
distribution, because the covariate is random and its adjustment must be
estimated; the correct reference distribution is the Bryant–Paulson
generalized studentized range
([`bryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)).

**Reproducible quantities.** Fitting
`lm(brand_movement ~ panel + block + category_movement)` gives a
covariate slope of \\0.4079\\ and an error mean square of \\0.01326\\ on
\\\nu = 14\\ degrees of freedom, with adjusted panel means \\3.595,
3.619, 4.102, 4.515, 4.618, 4.876\\ – exactly the values reported in the
paper. With \\q\_{.05;\\1,6,14} = 4.83\\
([`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)),
every pairwise simultaneous 95% interval has half-width \\0.139\\ and a
critical difference of \\0.278\\. Had the covariate not been measured,
the error mean square would have been \\0.2368\\ – roughly eighteen
times larger – and the intervals about four times wider. See
`data-raw/test_market.R` for the construction script and its
verification checks.

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
for the simultaneous intervals this data set illustrates,
[`bryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
for the critical values, and
[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for an
ANCOVA fit.

## Author

Ken Kelley

## Examples

``` r
data(test_market)
str(test_market)
#> 'data.frame':    24 obs. of  4 variables:
#>  $ panel            : Factor w/ 6 levels "1","2","3","4",..: 1 1 1 1 2 2 2 2 3 3 ...
#>  $ block            : Factor w/ 4 levels "1","2","3","4": 1 2 3 4 1 2 3 4 1 2 ...
#>  $ brand_movement   : num  2.98 3.29 4.33 3.48 2.99 3.36 3.99 4.18 4.49 4.25 ...
#>  $ category_movement: num  9.19 13.78 11.2 13.88 9.61 ...

# Reproduce the published ANCOVA (slope 0.4079, error MS 0.01326, df 14).
fit <- lm(brand_movement ~ panel + block + category_movement,
          data = test_market)
coef(fit)["category_movement"]
#> category_movement 
#>         0.4078841 
sum(residuals(fit)^2) / fit$df.residual
#> [1] 0.01325859

# Adjusted panel means at the covariate grand mean.
xbar <- mean(test_market$category_movement)
adj <- vapply(levels(test_market$panel), function(p) {
  nd <- data.frame(panel = factor(p, levels = levels(test_market$panel)),
                   block = factor(1:4, levels = levels(test_market$block)),
                   category_movement = xbar)
  mean(predict(fit, nd))
}, numeric(1))
adj  # 3.595 3.619 4.102 4.515 4.618 4.876
#>        1        2        3        4        5        6 
#> 3.595119 3.619463 4.101669 4.515084 4.617822 4.875843 

# Bryant-Paulson simultaneous 95% intervals (s = 4 blocks => n = 4, df 14).
ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326),
               n = 4, num_covariates = 1, df = 14)
#>  contrast          estimate lower_limit upper_limit
#>  group_1 - group_2 -0.0243  -0.302      0.254      
#>  group_1 - group_3 -0.507   -0.785      -0.228     
#>  group_1 - group_4 -0.92    -1.2        -0.642     
#>  group_1 - group_5 -1.02    -1.3        -0.745     
#>  group_1 - group_6 -1.28    -1.56       -1         
#>  group_2 - group_3 -0.482   -0.76       -0.204     
#>  group_2 - group_4 -0.896   -1.17       -0.618     
#>  group_2 - group_5 -0.998   -1.28       -0.72      
#>  group_2 - group_6 -1.26    -1.53       -0.978     
#>  group_3 - group_4 -0.413   -0.691      -0.135     
#>  group_3 - group_5 -0.516   -0.794      -0.238     
#>  group_3 - group_6 -0.774   -1.05       -0.496     
#>  group_4 - group_5 -0.103   -0.381      0.175      
#>  group_4 - group_6 -0.361   -0.639      -0.0827    
#>  group_5 - group_6 -0.258   -0.536      0.0201     
#> 
#> Confidence level: 95%
```
