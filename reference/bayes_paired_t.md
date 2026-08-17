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
bayes_paired_t(
  x = NULL,
  y = NULL,
  mean_diff = NULL,
  sd_diff = NULL,
  n = NULL,
  prior_location = 0,
  prior_scale = sqrt(2)/2,
  prior_mean = NULL,
  prior_sd = NULL,
  conf_level = 0.95
)
```

## Arguments

- x, y:

  Numeric vectors of paired observations, the same length, in matching
  order. The analysis is of `x - y`. Omit both to supply summary
  statistics of the differences instead.

- mean_diff, sd_diff, n:

  Summary statistics of the paired differences: their mean, their
  standard deviation, and the number of pairs. These are the quantities
  a paper's paired *t* test reports. Supply either raw data or all three
  summary values, never both.

- prior_location:

  Location of the Cauchy prior on \\\delta\\. Defaults to 0, the JZS
  prior; a nonzero value centers the prior on an expected effect
  (Gronau, Ly, & Wagenmakers, 2020).

- prior_scale:

  The way to adjust the prior. It is the scale (width) of the Cauchy
  prior on the standardized effect \\\delta\\ (the JZS prior), so a user
  who wants a more or less informative prior sets `prior_scale`: larger
  values say larger effects are plausible a priori, smaller values
  concentrate the prior near zero. The default \\\sqrt{2}/2 \approx
  0.707\\ is the JZS “medium” prior. Fully custom or subjective priors
  beyond the Cauchy family are not supported by the BayesFactor engine.

- prior_mean, prior_sd:

  Mean and standard deviation of a normal prior on \\\delta\\, for prior
  beliefs stated as moments. Supplying them selects the normal prior;
  they cannot be combined with the Cauchy arguments. See the prior
  section of Details.

- conf_level:

  Probability mass of the credible interval.

## Value

A `data.frame` (class `dmar_tbl`) with the same rows as
[`bayes_one_sample_t`](https://yelleknek.github.io/DMAR/reference/bayes_one_sample_t.md),
where \\\delta\\ is the standardized within-pair difference and the raw
rows are on the difference scale; `n` is the number of pairs.

**Specifying the prior.** The default prior on the standardized effect
\\\delta\\ is the JZS Cauchy centered at zero. Its `prior_scale` \\r\\
is not a standard deviation: a Cauchy has no mean and no variance (those
integrals diverge), so beliefs stated as prior moments cannot be
expressed through it. What the scale does fix is the quartiles: half the
prior mass lies within \\\pm r\\ of the location, so the default \\r =
\sqrt{2}/2\\ says a 50 percent prior bet that \\\|\delta\| \< 0.71\\. A
directional prior keeps the Cauchy and moves `prior_location` (Gronau,
Ly, & Wagenmakers, 2020). A researcher who thinks in prior moments
instead sets `prior_mean` and `prior_sd`, which use a normal prior with
exactly those moments; the two families are exclusive.

The families are also linked by an exact identity: a Cauchy with
location \\\mu\\ and scale \\r\\ is a normal prior \\N(\mu, r^2/z^2)\\
whose \\z\\ is standard normal, that is, a normal prior whose variance
you are not sure of. Choosing the Cauchy is therefore choosing a normal
prior with built-in doubt about its own width, which is why its tails
are heavier and its Bayes factors more conservative. A normal matched to
the Cauchy's interquartile range has `prior_sd = 1.4826 * prior_scale`.
The full posterior of \\\delta\\ is returned in the `"posterior"`
attribute as a data frame of `delta` and `density`, so any posterior
probability, not only the reported ones, can be computed from it.

## References

Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2020). Informed Bayesian
t-tests. *The American Statistician, 74*(2), 137–143.
[doi:10.1080/00031305.2018.1562983](https://doi.org/10.1080/00031305.2018.1562983)

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
#>  delta_posterior_median 0.624 
#>  delta_posterior_mean   0.626 
#>  delta_lower            0.237 
#>  delta_upper            1.02  
#>  p_delta_positive       0.999 
#>  raw_posterior_median   3.85  
#>  raw_lower              1.47  
#>  raw_upper              6.31  
#>  bf_10                  35.8  
#>  bf_01                  0.0279
#>  t                      3.69  
#>  df                     29    
#>  prior_location         0     
#>  prior_scale            0.707 
#>  n                      30    
#> 
#> Confidence level: 95%
```
