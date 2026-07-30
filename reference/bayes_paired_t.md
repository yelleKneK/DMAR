# Bayesian Paired-Samples *t* Analysis

The Bayesian counterpart of the paired *t* test, the analysis of
[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md)
applied to the within-pair differences. It reports the posterior of the
standardized difference \\\delta = \mu_D / \sigma_D\\ under the default
Jeffreys-Zellner-Siow (JZS) prior, a Cauchy prior on \\\delta\\ with
Jeffreys priors on the nuisance parameters (the variance), summarized by
its median, mean, a credible interval, and the probability statement
\\P(\delta \> 0 \mid \mathrm{data})\\. The JZS default Bayes factor
(Rouder, Speckman, Sun, Morey, & Iverson, 2009) is also reported. See
[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md)
for the model, the package's interpretive stance, and the computational
details (exact quadrature, no Monte Carlo error).

## Usage

``` r
bayes_paired_t(x, y, prior_scale = sqrt(2)/2, conf_level = 0.95)
```

## Arguments

- x, y:

  Numeric vectors of paired observations, the same length, in matching
  order. The analysis is of `x - y`.

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
where \\\delta\\ is the standardized within-pair difference and the raw
rows are on the difference scale; `n` is the number of pairs.

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
for the model and stance;
[`bayes_independent_t`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md)
for unpaired groups;
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)
and
[`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md)
for other paired analyses.

Other Bayesian t analyses:
[`bayes_independent_t()`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md),
[`bayes_one_sample_t()`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
before <- rnorm(30, 100, 12)
after  <- before + rnorm(30, 3, 6)
bayes_paired_t(after, before)
#>  term                   value 
#>  delta_posterior_median 0.623 
#>  delta_posterior_mean   0.626 
#>  delta_lower            0.237 
#>  delta_upper            1.02  
#>  p_delta_positive       0.999 
#>  raw_posterior_median   3.85  
#>  raw_lower              1.46  
#>  raw_upper              6.31  
#>  bf_10                  35.8  
#>  bf_01                  0.0279
#>  t                      3.69  
#>  df                     29    
#>  prior_scale            0.707 
#>  n                      30    
#> 
#> Confidence level: 95%
```
