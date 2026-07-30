# Exact Expected Value of Cohen's *d* (and Hedges' *g* Bias Correction)

Computes \\\mathrm{E}\[\hat d \mid \delta, n_1, n_2\]\\, the exact
expected value of the sample standardized mean difference (Cohen's *d*,
with pooled variance) under bivariate normality and the noncentral *t*
sampling distribution (Hedges, 1981). The sample *d* is upward-biased as
an estimator of the population \\\delta\\: \\\mathrm{E}\[\hat d\] =
\delta / J(\mathit{df})\\, where \$\$J(\mathit{df}) \\=\\
\frac{\Gamma(\mathit{df}/2)} {\sqrt{\mathit{df}/2}\\
\Gamma((\mathit{df}-1)/2)}\$\$ is Hedges' (1981) bias-correction factor
(\\\< 1\\ for finite \\\mathit{df}\\, tending to 1 as \\n \to \infty\\).
Hedges' *g*, the unbiased estimator of \\\delta\\, is then \\g = J \cdot
\hat d\\. The same \\J(\mathit{df})\\ is the workhorse of
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) when
`unbiased = TRUE`.

## Usage

``` r
expected_smd(delta, n_1, n_2 = NULL)
```

## Arguments

- delta:

  Population standardized mean difference. A numeric scalar or vector.

- n_1:

  Sample size in the first group. Scalar or vector.

- n_2:

  Sample size in the second group. Scalar or vector. If omitted,
  defaults to `n_1` (balanced design).

## Value

A `data.frame` with one row per (`delta`, `n_1`, `n_2`) input and the
columns

- `delta`, the population SMD,

- `n_1`, `n_2`, group sample sizes,

- `expected_smd`, \\\mathrm{E}\[\hat d \mid \delta\]\\,

- `bias`, \\\mathrm{E}\[\hat d\] - \delta\\ (the amount by which \\\hat
  d\\ overestimates \\\delta\\ on average),

- `j_correction`, the Hedges (1981) correction factor \\J(\mathit{df})\\
  used to compute the unbiased *g*, where \\\mathit{df} = n_1 + n_2 -
  2\\.

## Details

`expected_smd()` is especially useful at the *design stage* of a study,
where sample size planning typically proceeds from an assumed population
standardized mean difference \\\delta\\. Because the value of \\\hat d\\
a researcher should expect to observe on average is *larger in absolute
value* than \\\delta\\, plugging \\\delta\\ directly into
sampling-distribution machinery (precision of \\\hat d\\, width of a
confidence interval, power of a test of \\H_0\\: \delta = 0\\)
over-promises on the realized precision when the precision is expressed
on the \\\hat d\\ scale. Substituting `expected_smd(delta, n_1, n_2)`
for the bare \\\delta\\ corrects this leading-order bias.

**Derivation.** Under bivariate normality with equal variances, the
observed *t*-statistic \\t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}\\
follows a noncentral *t* distribution with \\\mathit{df} = n_1 + n_2 -
2\\ degrees of freedom and noncentrality parameter \\\lambda = \delta
\sqrt{n_1 n_2 / (n_1 + n_2)}\\. The expected value of a noncentral *t*
variate equals \\\lambda / J(\mathit{df})\\ (Johnson, Kotz, &
Balakrishnan, 1995, Section 31.3), so \\\mathrm{E}\[\hat d\] =
\mathrm{E}\[t\] / \sqrt{n_1 n_2 / (n_1 + n_2)} = \delta /
J(\mathit{df})\\. The bias \\\mathrm{E}\[\hat d\] - \delta =
\delta\\(1 - J)/J\\ is positive when \\\delta \> 0\\.

**Magnitude of the correction.** \\J(\mathit{df}) \approx 1 -
3/(4\\\mathit{df} - 1)\\ to leading order (Hedges & Olkin, 1985, p.\\
81). For \\\delta = 0.5\\, \\n_1 = n_2 = 10\\ (\\\mathit{df} = 18\\),
\\J \approx 0.957\\, so \\\mathrm{E}\[\hat d\] \approx 0.522\\ and the
upward bias is about 4%. For \\n_1 = n_2 = 50\\ the bias is under 1%;
for \\n_1 = n_2 = 5\\ (very small samples) it exceeds 10%.

**Connection to Hedges' *g*.** The natural inverse of this function is
Hedges' *g*: given an observed \\\hat d\\, the unbiased estimator of
\\\delta\\ is \\g = J(\mathit{df}) \hat d\\, which satisfies
\\\mathrm{E}\[g \mid \delta\] = \delta\\ exactly under the same
noncentral *t* model.
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) with
`unbiased = TRUE` returns *g*.

## References

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press. (See Section 5, equations 6 and 9.)

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions, volume 2* (2nd ed.), Section 31.3. Wiley.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced design: delta = 0.5, n = 10 per group.
expected_smd(delta = 0.5, n_1 = 10)
#>  delta n_1 n_2 expected_smd bias   j_correction
#>  0.5   10  10  0.522        0.0221 0.958       

# 2. Bias as a function of n at fixed delta.
expected_smd(delta = 0.5, n_1 = c(5, 10, 20, 50, 100, 500))
#>  delta n_1 n_2 expected_smd bias     j_correction
#>  0.5   5   5   0.554        0.0539   0.903       
#>  0.5   10  10  0.522        0.0221   0.958       
#>  0.5   20  20  0.51         0.0101   0.98        
#>  0.5   50  50  0.504        0.00387  0.992       
#>  0.5   100 100 0.502        0.0019   0.996       
#>  0.5   500 500 0.5          0.000376 0.999       

# 3. Unbalanced design.
expected_smd(delta = 0.5, n_1 = 30, n_2 = 60)
#>  delta n_1 n_2 expected_smd bias    j_correction
#>  0.5   30  60  0.504        0.00431 0.991       

# 4. Design-stage use: an a-priori delta of 0.4, planned n of 40/group.
#        The d we should expect to observe on average is slightly larger.
expected_smd(delta = 0.4, n_1 = 40)
#>  delta n_1 n_2 expected_smd bias   j_correction
#>  0.4   40  40  0.404        0.0039 0.99        
```
