# Minimum Risk Point Estimation of the Population Coefficient of Variation

Sequential estimation of the coefficient of variation with minimum risk.
The function implements the procedure of Chattopadhyay and Kelley
(2016), which considers study cost and the accuracy of the estimated
coefficient of variation together. For the variance of the sample
coefficient of variation it uses an improved analytic solution, Equation
32 of Kelley, Darku, and Chattopadhyay (2018), in place of Equation 22
of the 2016 article, which carries an error in its kurtosis term (see
the part of Details on Equation 22). It returns the pilot sample size,
or it evaluates the risk and the stopping rule on the observations
collected so far and says whether sampling can stop.

## Usage

``` r
mr_cv(
  data,
  A,
  structural_cost,
  epsilon,
  sampling_cost,
  pilot = FALSE,
  m0 = 4,
  gamma = 0.49,
  verbose = FALSE
)
```

## Arguments

- data:

  A numeric vector of the observations collected so far (at least 4),
  from a quantity with a positive mean. It is needed when
  `pilot = FALSE`.

- A:

  The structural cost divided by the square of `epsilon` (see Details),
  a single positive number. Give either `A` or both `structural_cost`
  and `epsilon`.

- structural_cost:

  The structural cost: what one is willing to pay for the study apart
  from the cost of sampling (see Details), a single positive number.

- epsilon:

  The maximum desired difference between the estimated coefficient of
  variation and the population value, a single positive number.

- sampling_cost:

  The cost of collecting one more observation, a single positive number.
  For example, if each survey costs 10 dollars to distribute and score,
  `sampling_cost` would be 10 dollars per additional observation.

- pilot:

  `TRUE` to plan the pilot sample size, or `FALSE` (the default) to
  evaluate the risk function and the stopping rule at the current sample
  size.

- m0:

  The smallest pilot sample size allowed: a whole number of at least 4
  (the default), since the variance estimate needs four observations.

- gamma:

  The exponent \\\gamma\\ of the correction term \\n^{-2\gamma}\\ in the
  stopping rule, which keeps sampling from stopping too early: a single
  number strictly between 0 and 1/2. Chattopadhyay and Kelley (2016)
  take \\0 \< \gamma \< 1/2\\ and use .49 for practical purposes, the
  default.

- verbose:

  If `TRUE`, rows with the quantities behind the stopping rule are added
  to the table; defaults to `FALSE`.

## Value

A `dmar_tbl` with columns `term` and `value`. With `pilot = FALSE` the
rows are

- risk:

  The estimated risk at the current sample size, \\A V_n^2 / n + cn\\,
  with \\\max(V_n^2, n^{-3})\\ for \\V_n^2\\.

- N:

  The current sample size (the number of observations in `data`).

- cv:

  The current coefficient of variation.

With `verbose = TRUE` the rows `v2_n`, `mu_3_n` and `mu_4_n` (the
unbiased estimates of the third and fourth central moments it is built
from), and `criterion` follow. The row `v2_n` always holds the estimate
\\V_n^2\\ itself, also when \\n^{-3}\\ was used in its place and when
the estimate is zero or negative. The attribute `is_satisfied` is `TRUE`
when the stopping rule is satisfied, so that sampling can stop; the
attribute `criterion` holds \\(A / c)(\max(V_n^2, n^{-3}) +
n^{-2\gamma})\\, the value that \\N^2\\ is compared with; and the
attribute `variance_unusable` is `TRUE` when the estimate is zero,
negative, or not finite. In that case `is_satisfied` is `FALSE`, and the
row `risk`, the row `criterion`, and the attribute `criterion` are `NA`.
The decision is printed under the table: whether the stopping rule is
met, the two numbers compared, and whether to stop or to collect another
observation. With `pilot = TRUE` the single row is `pilot_N`, the pilot
sample size.

## Details

The value of `epsilon` is context specific; the smaller the value the
closer the estimated value will tend to be to the population value.

When the aim of a study is to estimate a parameter accurately, such as
the coefficient of variation, the structural cost and the maximum
probable error of the estimate, \\\epsilon\\, are combined to form
\\A\\, the structural cost divided by \\\epsilon^2\\: the price paid per
squared unit of maximum probable error. The structural cost is the
investment in the study apart from the cost of sampling; in a grant
application it is the funds requested less those for sampling
observations (Chattopadhyay & Kelley, 2016). A larger structural cost or
a smaller \\\epsilon\\ gives a larger \\A\\, and a larger \\A\\ leads to
a more expensive study, other things being equal. \\A\\ is fixed before
any data are collected, whether it is given directly or through its two
components.

The sampling cost is not fixed in advance, because the sample size
needed to estimate the coefficient of variation accurately is unknown.
The risk function of Chattopadhyay and Kelley (2016) adds the mean
square error of the estimator, weighted by \\A\\, to the sampling cost
\\cn\\, where \\c\\ is `sampling_cost` and \\n\\ the sample size. Their
sequential procedure collects a pilot sample, then adds one observation
at a time, and stops at the first \\n\\ for which \\n^2 \ge (A /
c)(V_n^2 + n^{-2\gamma})\\ (their Equation 23), where \\V_n^2\\
estimates \\\xi^2\\, the variance of the limiting normal distribution of
\\\sqrt{n}(k_n - \kappa)\\, with \\k_n\\ the sample coefficient of
variation and \\\kappa\\ the population value. Both sides of the rule
are on the scale of \\n^2\\: the criterion the function reports is the
right-hand side, and it is compared with the square of the sample size,
not with the sample size. The row `N` of the table holds this \\n\\;
\\N^2\\, here and in the line printed under the table, is its square.
With `pilot = FALSE` the function evaluates the risk and the stopping
rule at the current sample size; with `pilot = TRUE` it returns the
pilot sample size \\\max(m_0, \lceil (A / c)^{1/(2 + 2\gamma)}
\rceil)\\. A sample size is a whole number, so the formula is rounded
up, and a formula value that is itself a whole number is that number:
with \\A = 100{,}000\\, \\c = 1\\, and \\\gamma = .25\\ it is
\\100{,}000^{0.4} = 100\\, which R computes as 100.00000000000003, and
the pilot sample size is 100, not 101.
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md) is the
corresponding procedure for the standardized mean difference
(Chattopadhyay & Kelley, 2017). For a confidence interval for the
coefficient of variation of normally distributed scores, and for the
sample size that makes it narrow, see
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) and
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
(Kelley, 2007).

### The Variance Estimate

For scores with mean \\\mu\\, variance \\\sigma^2\\, and third and
fourth central moments \\\mu_3\\ and \\\mu_4\\, \$\$\xi^2 =
\frac{\mu_4}{4\mu^2\sigma^2} - \frac{\sigma^2}{4\mu^2} -
\frac{\mu_3}{\mu^3} + \frac{\sigma^4}{\mu^4}\$\$ (Kelley, Darku, &
Chattopadhyay, 2018, Equation 30). With the skewness \\\gamma_1\\ and
the excess kurtosis \\\gamma_2\\ of the scores (neither is the exponent
\\\gamma\\ of the stopping rule), this is \$\$\xi^2 =
\frac{\kappa^2}{2} + \kappa^4 - \kappa^3\gamma_1 +
\frac{\kappa^2\gamma_2}{4}.\$\$ \\V_n^2\\ is the expression in moments
with the sample mean, the sample variance, and the unbiased estimates of
the third and fourth central moments in place of the population values
(their Equation 32), each moment computed from the scores centered at
their mean. It is consistent for \\\xi^2\\ whenever the fourth moment is
finite, so the procedure does not assume a specific distribution for the
scores. The scores must have a positive mean: the coefficient of
variation divides by the mean, and the procedure is for a quantity that
cannot be negative (Chattopadhyay & Kelley, 2016, p. 627).

### A Zero or Negative Estimate

\\V_n^2\\ can be zero or negative in a small sample. A positive estimate
enters the stopping rule and the risk as \\\max(V_n^2, n^{-3})\\, the
estimator of Equation 31 of Kelley, Darku, and Chattopadhyay (2018). The
bound is small beside the correction term of the rule: \\n^{-3}\\ is
\\n^{2\gamma - 3}\\ times \\n^{-2\gamma}\\, which is \\n^{-2.02}\\ times
it at the default \\\gamma\\ and less than \\n^{-2}\\ times it for any
\\\gamma\\. It replaces a positive estimate below \\n^{-3}\\. With a
small coefficient of variation that is the usual case in a small sample:
for normal scores \\\xi^2\\ is close to \\\kappa^2/2\\, which is below
\\n^{-3}\\ while \\n\\ is below \\(2/\kappa^2)^{1/3}\\, about 27
observations for \\\kappa = .01\\ and 9 for \\\kappa = .05\\. From the
pilot sample size on, the bound adds less than 1 to the criterion, so it
can change the decision only at the pilot sample size, and only when
\\N^2\\ there is within 1 of the criterion; the stop then usually comes
one observation later.

When the estimate is zero or negative (it is taken as 0 for scores that
do not vary, for scores that differ only by rounding error, with a range
of at most \\10^{-10}\\ times the largest absolute score, such as equal
gains computed as posttest minus pretest, and when it is within
\\10^{-10}\\ of the size of its four terms), the function reports the
stopping rule as not met, with the risk and the criterion missing, and
says so under the table; the next observation is collected and the rule
is evaluated again. This is a deliberate departure from Equation 23 read
literally, under which a negative \\V_n^2\\ lowers the criterion and can
end sampling at the pilot sample, and from Equation 31, which would put
\\n^{-3}\\ in place of such an estimate: the rule can then be met by a
sample that says nothing usable about the variance.

### Large-Sample Properties

The properties of the procedure are asymptotic (Chattopadhyay & Kelley,
2016). As \\A / c\\ grows, the final sample size approaches the
large-sample value \\\sqrt{A / c}\\\xi\\, at which the approximate risk
\\A \xi^2 / n + cn\\ is smallest. Short of that limit the procedure can
stop early, on average, when the scores are skewed. In 20,000 simulated
runs with \\A = 5{,}000{,}000\\ and \\c = 10\\, the mean final sample
size was 480 for gamma scores with \\\kappa = .8\\, which is 0.94 of the
large-sample value of 512, and 440 for lognormal scores with \\\kappa =
.6\\, which is 0.81 of the large-sample value of 545.

### Equation 22 of Chattopadhyay and Kelley (2016)

The variance estimate of this function improves on the one the 2016
article uses, and its results therefore differ from the article's.
Chattopadhyay and Kelley (2016) estimate \\\xi^2\\ by their Equation 22,
which the article prints in moments, \$\$\frac{s_n^4}{4\bar{X}\_n^4} +
\frac{\hat{\mu}\_{4n}}{4\bar{X}\_n^4} + \frac{s_n^2}{2\bar{X}\_n^2} -
\frac{\hat{\mu}\_{3n}}{\bar{X}\_n^3},\$\$ and of whose four terms only
the last is a term of Equation 32 of Kelley, Darku, and Chattopadhyay
(2018). Written with the skewness and the excess kurtosis, Equation 22
estimates \\\kappa^2/2 + \kappa^4 - \kappa^3\gamma_1 +
\kappa^4\gamma_2/4\\, the first-order term of Bao (2009, Equation 3) as
printed. That term is in error. Bao posted an erratum to the article on
his web page at Purdue University (Bao, n.d.; the file is dated February
20, 2023), which says that the corrections were missed in proofreading
and makes three: in the first line of the mean squared error of his
Equation 3 (p. 293) the term \\\alpha^4\gamma_2/(4n)\\ should be
\\\alpha^2\gamma_2/(4n)\\, where \\\alpha\\ is the coefficient of
variation and \\\gamma_2\\ the excess kurtosis; in the last line of his
Remark 2 (p. 296) \\\alpha^4\gamma_2\\ should be \\\alpha^2\gamma_2\\;
and in the second line of the expression for \\E(\lambda\_{22})\\ (p.
295) the term \\(3/n^2 + 1/n)\alpha^2\gamma_2\\ should be \\(3/n^2 +
1/(n\alpha^2))\alpha^2\gamma_2\\. It credits Dennis Boos of North
Carolina State University with finding the errors. The kurtosis term of
the first-order variance is therefore \\\kappa^2\gamma_2/4\\, as
Equation 30 of the 2018 article has, and Equation 32, the estimate used
here, is the sample version of that corrected variance. The two
population values are equal only when the excess kurtosis is 0 or
\\\kappa = 1\\; with positive excess kurtosis the value of Equation 22
is too small for \\\kappa\\ below 1 and too large above it. For gamma
scores with \\\kappa = .8\\, \\\xi^2\\ is 0.525 and the value of
Equation 22 is 0.304. In a sample the two estimates differ by
\\(k_n^2/4)(1 - k_n^2)(\hat{\mu}\_{4n}/s_n^4 - 3)\\, so the results of
this function differ from those of the article for nearly any data,
normal data included. For the article's example, the first of the
Examples, the stopping rule is met at 35 observations with Equation 22,
as the article reports, and is not met there with the variance estimate
used here: the criterion is 1261.332 and \\N^2\\ is 1225.

## References

Bao, Y. (n.d.). Erratum: “Finite sample moments of the coefficient of
variation” \[Econometric Theory 25(1) (2009) 291–297\]. Department of
Economics, Purdue University. Retrieved November 15, 2026, from
<https://web.ics.purdue.edu/~ybao/Research/erratum_BaoET.pdf>

Bao, Y. (2009). Finite-sample moments of the coefficient of variation.
*Econometric Theory, 25*(1), 291–297.
[doi:10.1017/S0266466608090555](https://doi.org/10.1017/S0266466608090555)

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods, 22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md),
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The example of Chattopadhyay and Kelley (2016): a structural cost of
# 1000 and a maximum probable error of .05, so that A is 400000, with a
# sampling cost of 75 per observation. The pilot sample size, from A or
# from its two components:
mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75)
#>  term    value
#>  pilot_N 18   
mr_cv(pilot = TRUE, structural_cost = 1000, epsilon = .05, sampling_cost = 75)
#>  term    value
#>  pilot_N 18   

# The 18 pilot observations of the article. The line under the table
# gives the decision.
books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41, 5, 3)
mr_cv(data = books, A = 400000, sampling_cost = 75)
#>  term value   
#>  risk 5354.513
#>  N    18      
#>  cv   0.739   
#> 
#> Stopping rule: not met (N^2 = 324 is below the criterion, 1275.012); collect another observation.

# The rule is not met, so one more observation is collected and the rule
# is checked again, and so on. The article stops at 35 observations; with
# the variance estimate used here the rule is not yet met there (see the
# part of Details on Equation 22). The extra rows show what the decision
# rests on.
books <- c(books, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23, 40, 18)
mr_cv(data = books, A = 400000, sampling_cost = 75, verbose = TRUE)
#>  term      value     
#>  risk      4977.259  
#>  N         35        
#>  cv        0.701     
#>  v2_n      0.206     
#>  mu_3_n    3286.967  
#>  mu_4_n    186275.301
#>  criterion 1261.332  
#> 
#> Stopping rule: not met (N^2 = 1225 is below the criterion, 1261.332); collect another observation.

# A complete run. The simulated scores stand in for observations that a
# study would collect one at a time. The loop starts at the pilot sample
# size, adds an observation whenever the rule is not met, and cannot go
# past the scores on hand.
set.seed(113)
scores <- rgamma(1000, shape = 4, scale = 5)
n <- mr_cv(pilot = TRUE, A = 100000, sampling_cost = 1)$value
repeat {
  check <- mr_cv(data = scores[1:n], A = 100000, sampling_cost = 1)
  if (attr(check, "is_satisfied") || n == length(scores)) break
  n <- n + 1
}
check
#>  term value  
#>  risk 194.843
#>  N    103    
#>  cv   0.506  
#> 
#> Stopping rule: met (N^2 = 10609 is at least the criterion, 10525.008); sampling can stop.
```
