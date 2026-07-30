# Bayesian Independent-Samples *t* Analysis

The Bayesian counterpart of the two-sample (pooled-variance) *t* test.
It reports the posterior of the standardized mean difference \\\delta =
(\mu_1 - \mu_2)/\sigma\\ under the default Jeffreys-Zellner-Siow (JZS)
prior, a Cauchy prior on \\\delta\\ with Jeffreys priors on the nuisance
parameters (the variance), summarized by its median, mean, a credible
interval, and the probability statement \\P(\delta \> 0 \mid
\mathrm{data})\\. The JZS default Bayes factor (Rouder, Speckman, Sun,
Morey, & Iverson, 2009) is also reported. See
[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md)
for the model, the package's interpretive stance, and the computational
details (exact quadrature, no Monte Carlo error).

## Usage

``` r
bayes_independent_t(x, y, prior_scale = sqrt(2)/2, conf_level = 0.95)
```

## Arguments

- x, y:

  Numeric vectors: the observations of the two independent groups
  (\\\delta\\ is positive when `x` runs higher).

- prior_scale:

  The way to adjust the prior. It is the scale (width) of the Cauchy
  prior on the standardized effect \\\delta\\ (the JZS prior), so a user
  who wants a more or less informative prior sets `prior_scale`: larger
  values say larger effects are plausible a priori, smaller values
  concentrate the prior near zero. The default \\\sqrt{2}/2 \approx
  0.707\\ is the JZS “medium” prior. Fully custom or subjective priors
  beyond the Cauchy family are not supported by the BayesFactor engine.

- conf_level:

  Probability mass of the credible interval.

## Value

A `data.frame` (class `dmar_tbl`) with the same rows as
[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
plus `n_1` and `n_2`.

## Details

The likelihood of the pooled-variance *t* statistic given \\\delta\\ is
noncentral *t* with \\\mathit{df} = n_1 + n_2 - 2\\ and noncentrality
\\\delta \sqrt{n_1 n_2 / (n_1 + n_2)}\\; equal variances are assumed, as
in the standard JZS development. The raw-scale rows transform the
\\\delta\\ summaries through the pooled standard deviation.

## References

Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G.
(2009). Bayesian t tests for accepting and rejecting the null
hypothesis. *Psychonomic Bulletin & Review, 16*(2), 225–237.
[doi:10.3758/PBR.16.2.225](https://doi.org/10.3758/PBR.16.2.225)

Jeffreys, H. (1961). *Theory of probability* (3rd ed.). Oxford
University Press.

Zellner, A., & Siow, A. (1980). Posterior odds ratios for selected
regression hypotheses. In J. M. Bernardo, M. H. DeGroot, D. V. Lindley,
& A. F. M. Smith (Eds.), *Bayesian statistics: Proceedings of the First
International Meeting* (pp. 585–603). University of Valencia Press.

## See also

[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md)
and
[`bayes_paired_t`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md)
for the other designs;
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) and
[`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md) for
the frequentist analyses of the same comparison.

Other Bayesian t analyses:
[`bayes_one_sample_t()`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
[`bayes_paired_t()`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
g1 <- rnorm(35, 105, 15)
g2 <- rnorm(35, 100, 15)
bayes_independent_t(g1, g2)
#>  term                   value
#>  delta_posterior_median 0.243
#>  delta_posterior_mean   0.246
#>  delta_lower            -0.19
#>  delta_upper            0.693
#>  p_delta_positive       0.864
#>  raw_posterior_median   4    
#>  raw_lower              -3.13
#>  raw_upper              11.4 
#>  bf_10                  0.454
#>  bf_01                  2.2  
#>  t                      1.2  
#>  df                     68   
#>  prior_scale            0.707
#>  n_1                    35   
#>  n_2                    35   
#> 
#> Confidence level: 95%

# The probability that the effect is positive is read directly off the
# p_delta_positive row: a probability statement, not a p-value.
```
