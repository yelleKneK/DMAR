# Bayesian One-Sample *t* Analysis

The Bayesian counterpart of the one-sample *t* test. It reports the
posterior distribution of the standardized effect \\\delta = (\mu -
\mu_0)/\sigma\\ under the default Jeffreys-Zellner-Siow (JZS) prior, a
Cauchy prior on \\\delta\\ with Jeffreys priors on the nuisance
parameters (the variance), summarized by its median, mean, a credible
interval, and the direct probability statement \\P(\delta \> 0 \mid
\mathrm{data})\\. The JZS default Bayes factor (Rouder, Speckman, Sun,
Morey, & Iverson, 2009) is also reported.

## Usage

``` r
bayes_one_sample_t(x, mu_0 = 0, prior_scale = sqrt(2)/2, conf_level = 0.95)
```

## Arguments

- x:

  Numeric vector of observations.

- mu_0:

  The comparison value for the mean under the point null (and the
  centering value for \\\delta\\). Defaults to 0.

- prior_scale:

  The way to adjust the prior. It is the scale (width) \\r\\ of the
  Cauchy prior on the standardized effect \\\delta\\ (the JZS prior), so
  a user who wants a more or less informative prior sets `prior_scale`:
  larger values say larger effects are plausible a priori, smaller
  values concentrate the prior near zero. The default \\\sqrt{2}/2
  \approx 0.707\\ is the JZS “medium” prior. Fully custom or subjective
  priors beyond the Cauchy family are not supported by the BayesFactor
  engine.

- conf_level:

  Probability mass of the (central) credible interval. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_tbl`) with the posterior summaries of
\\\delta\\ (`delta_posterior_median`, `delta_posterior_mean`,
`delta_lower`, `delta_upper`, `p_delta_positive`), the same summaries
mapped to the raw mean difference scale (`raw_*`), the Bayes factors
(`bf_10`, `bf_01`), the observed `t` and `df`, the `prior_scale`, and
`n`.

## Details

The posterior is a *probability statement* about the parameter given the
model, the prior, and the data, and that is how these functions are
meant to be read. A Bayes factor is a different kind of claim, a
comparison of how well two models predicted the data, and it leans
harder on the prior; it is reported because it may be helpful for some
questions. Neither replaces the estimation-first habits of the rest of
the package;
[`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md) and
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) remain
the frequentist complements.

With a Jeffreys prior on \\(\mu_0, \sigma^2)\\ and \\\delta \sim
\mathrm{Cauchy}(0, r)\\, all inference flows through the observed *t*
statistic, whose likelihood given \\\delta\\ is noncentral *t* with
noncentrality \\\delta \sqrt{n}\\. The posterior of \\\delta\\ is
computed by quadrature (no Monte Carlo error) and the Bayes factor by
the one-dimensional integral of that likelihood against the Cauchy
prior, the exact JZS form. The raw-scale rows transform the \\\delta\\
summaries through the sample standard deviation (a plug-in, as is
conventional for reporting).

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

[`bayes_paired_t`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md)
and
[`bayes_independent_t`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md)
for the two-sample designs;
[`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md) for the
frequentist standardized mean.

Other Bayesian t analyses:
[`bayes_independent_t()`](https://yelleknek.github.io/DMAR/reference/bayes_independent_t.md),
[`bayes_paired_t()`](https://yelleknek.github.io/DMAR/reference/bayes_paired_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
x <- rnorm(40, mean = 0.4, sd = 1)
bayes_one_sample_t(x)
#>  term                   value
#>  delta_posterior_median 0.502
#>  delta_posterior_mean   0.503
#>  delta_lower            0.178
#>  delta_upper            0.832
#>  p_delta_positive       0.999
#>  raw_posterior_median   0.518
#>  raw_lower              0.184
#>  raw_upper              0.858
#>  bf_10                  20   
#>  bf_01                  0.05 
#>  t                      3.39 
#>  df                     39   
#>  prior_scale            0.707
#>  n                      40   
#> 
#> Confidence level: 95%

# Against a nonzero comparison value, with a wider prior.
bayes_one_sample_t(x, mu_0 = 0.1, prior_scale = 1)
#>  term                   value
#>  delta_posterior_median 0.421
#>  delta_posterior_mean   0.423
#>  delta_lower            0.103
#>  delta_upper            0.743
#>  p_delta_positive       0.995
#>  raw_posterior_median   0.435
#>  raw_lower              0.106
#>  raw_upper              0.767
#>  bf_10                  3.87 
#>  bf_01                  0.259
#>  t                      2.78 
#>  df                     39   
#>  prior_scale            1    
#>  n                      40   
#> 
#> Confidence level: 95%
```
