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
bayes_one_sample_t(
  x = NULL,
  mu_0 = 0,
  mean = NULL,
  sd = NULL,
  n = NULL,
  prior_location = 0,
  prior_scale = sqrt(2)/2,
  prior_mean = NULL,
  prior_sd = NULL,
  conf_level = 0.95
)
```

## Arguments

- x:

  Numeric vector of observations. Omit to supply summary statistics
  instead.

- mu_0:

  The comparison value for the mean under the point null (and the
  centering value for \\\delta\\). Defaults to 0.

- mean, sd, n:

  Summary statistics: the sample mean, standard deviation, and sample
  size. The Bayes factor depends on the data only through the *t*
  statistic and \\n\\, so the summary form is exact, not an
  approximation. Supply either `x` or all three summary values, never
  both.

- prior_location:

  Location of the Cauchy prior on \\\delta\\. Defaults to 0, the JZS
  prior; a nonzero value centers the prior on an expected effect
  (Gronau, Ly, & Wagenmakers, 2020).

- prior_scale:

  The way to adjust the prior. It is the scale (width) \\r\\ of the
  Cauchy prior on the standardized effect \\\delta\\ (the JZS prior), so
  a user who wants a more or less informative prior sets `prior_scale`:
  larger values say larger effects are plausible a priori, smaller
  values concentrate the prior near zero. The default \\\sqrt{2}/2
  \approx 0.707\\ is the JZS “medium” prior. Fully custom or subjective
  priors beyond the Cauchy family are not supported by the BayesFactor
  engine.

- prior_mean, prior_sd:

  Mean and standard deviation of a normal prior on \\\delta\\, for prior
  beliefs stated as moments. Supplying them selects the normal prior;
  they cannot be combined with the Cauchy arguments. See the prior
  section of Details.

- conf_level:

  Probability mass of the (central) credible interval. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_tbl`) with the posterior summaries of
\\\delta\\ (`delta_posterior_median`, `delta_posterior_mean`,
`delta_lower`, `delta_upper`, `p_delta_positive`), the same summaries
mapped to the raw mean difference scale (`raw_*`), the Bayes factors
(`bf_10`, `bf_01`), the observed `t` and `df`, the `prior_scale`, and
`n`.

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
observed *d* relative to `mu_0 = 0` is `mean = d, sd = 1`.

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
#>  delta_posterior_median 0.503
#>  delta_posterior_mean   0.503
#>  delta_lower            0.179
#>  delta_upper            0.832
#>  p_delta_positive       0.999
#>  raw_posterior_median   0.519
#>  raw_lower              0.185
#>  raw_upper              0.859
#>  bf_10                  20   
#>  bf_01                  0.05 
#>  t                      3.39 
#>  df                     39   
#>  prior_location         0    
#>  prior_scale            0.707
#>  n                      40   
#> 
#> Confidence level: 95%

# Against a nonzero comparison value, with a wider prior.
bayes_one_sample_t(x, mu_0 = 0.1, prior_scale = 1)
#>  term                   value
#>  delta_posterior_median 0.422
#>  delta_posterior_mean   0.423
#>  delta_lower            0.105
#>  delta_upper            0.744
#>  p_delta_positive       0.995
#>  raw_posterior_median   0.436
#>  raw_lower              0.108
#>  raw_upper              0.768
#>  bf_10                  3.87 
#>  bf_01                  0.259
#>  t                      2.78 
#>  df                     39   
#>  prior_location         0    
#>  prior_scale            1    
#>  n                      40   
#> 
#> Confidence level: 95%
```
