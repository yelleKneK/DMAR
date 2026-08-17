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
bayes_independent_t(
  x = NULL,
  y = NULL,
  mean_1 = NULL,
  sd_1 = NULL,
  n_1 = NULL,
  mean_2 = NULL,
  sd_2 = NULL,
  n_2 = NULL,
  prior_location = 0,
  prior_scale = sqrt(2)/2,
  prior_mean = NULL,
  prior_sd = NULL,
  conf_level = 0.95
)
```

## Arguments

- x, y:

  Numeric vectors: the observations of the two independent groups
  (\\\delta\\ is positive when `x` runs higher). Omit both to supply
  summary statistics instead.

- mean_1, sd_1, n_1:

  Summary statistics of the first group: mean, standard deviation, and
  sample size. The Bayes factor depends on the data only through the *t*
  statistic and the sample sizes, so the summary form is exact, not an
  approximation. Supply either raw data or all six summary values, never
  both.

- mean_2, sd_2, n_2:

  Summary statistics of the second group.

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
plus `n_1` and `n_2`.

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

A standardized effect size enters through the summary form directly: an
observed Cohen's *d* with group sizes \\n_1\\ and \\n_2\\ is
`mean_1 = d, mean_2 = 0, sd_1 = 1, sd_2 = 1`, since *d* is the mean
difference in pooled standard deviation units.

## Details

The likelihood of the pooled-variance *t* statistic given \\\delta\\ is
noncentral *t* with \\\mathit{df} = n_1 + n_2 - 2\\ and noncentrality
\\\delta \sqrt{n_1 n_2 / (n_1 + n_2)}\\; equal variances are assumed, as
in the standard JZS development. The raw-scale rows transform the
\\\delta\\ summaries through the pooled standard deviation.

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
#>  delta_posterior_median 0.244 
#>  delta_posterior_mean   0.246 
#>  delta_lower            -0.189
#>  delta_upper            0.693 
#>  p_delta_positive       0.864 
#>  raw_posterior_median   4.01  
#>  raw_lower              -3.11 
#>  raw_upper              11.4  
#>  bf_10                  0.454 
#>  bf_01                  2.2   
#>  t                      1.2   
#>  df                     68    
#>  prior_location         0     
#>  prior_scale            0.707 
#>  n_1                    35    
#>  n_2                    35    
#> 
#> Confidence level: 95%

# The probability that the effect is positive is read directly off the
# p_delta_positive row: a probability statement, not a p-value.
```
