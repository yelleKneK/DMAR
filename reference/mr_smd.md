# Minimum Risk Point Estimation of the Population Standardized Mean Difference

A function for the sequential estimation of the standardized mean
difference with minimum risk. The function implements the procedure of
Chattopadhyay and Kelley (2017), which considers study cost and the
accuracy of the estimated standardized mean difference together.
`mr_smd` assumes normally distributed data, equal sample sizes in the
two groups, and an equal cost of sampling per observation in each group.

## Usage

``` r
mr_smd(
  A,
  structural_cost,
  epsilon,
  d,
  n_per_group,
  sampling_cost,
  pilot = FALSE,
  m0 = 4,
  gamma = 0.49
)
```

## Arguments

- A:

  The structural cost divided by the square of `epsilon`: the price one
  is willing to pay for a maximum allowable difference of `epsilon`
  between the estimate of the standardized mean difference and its
  parameter, a single positive number. Give either `A` or both
  `structural_cost` and `epsilon`.

- structural_cost:

  The structural cost of what one is willing to pay in a study, a single
  positive number.

- epsilon:

  The maximum desired difference between the estimated standardized mean
  difference and the population value, a single positive number.

- d:

  The current estimate of the standardized mean difference, a single
  finite number. It is needed when `pilot = FALSE`.

- n_per_group:

  Current sample size *per group* (thus total sample size is \\2n\\), a
  whole number of at least 2; requires equal sample size *per group*. It
  is needed when `pilot = FALSE`.

- sampling_cost:

  The cost of collecting one more observation, a single positive number.
  For example, if each survey costs 10 dollars to distribute and score,
  `sampling_cost` would be 10 dollars per additional observation.

- pilot:

  `TRUE` to plan the pilot sample size per group, or `FALSE` (the
  default) to evaluate the risk function and the stopping rule at the
  current sample size.

- m0:

  The smallest pilot sample size per group allowed: a whole number of at
  least 4 (the default).

- gamma:

  The exponent \\\gamma\\ of the correction term \\n^{-\gamma}\\ in the
  stopping rule, which keeps sampling from stopping too early: a single
  number strictly between 0 and 1/2. Chattopadhyay and Kelley (2017)
  take \\0 \< \gamma \< 1/2\\ and suggest .49 for practical purposes,
  the default.

## Value

A `dmar_tbl` with columns `term` and `value`. With `pilot = FALSE` the
rows are

- risk:

  The estimated risk at the current sample size, \\A V_n^2 / n + 2cn\\.

- n_1:

  Sample size for group 1 (the input `n_per_group`).

- n_2:

  Sample size for group 2 (the input `n_per_group`).

- d:

  The input estimate of the standardized mean difference.

The attribute `is_satisfied` is `TRUE` when the stopping rule is
satisfied, so that sampling can stop, and the attribute `criterion`
holds \\\sqrt{A / (2c)}\\(V_n + n^{-\gamma})\\, the value that the
sample size per group is compared with. The decision is printed under
the table: whether the stopping rule is met, the two numbers compared,
and whether to stop or to collect another observation in each group.
With `pilot = TRUE` the single row is `pilot_n_per_group` (see Note).

## Details

Chattopadhyay and Kelley (2017) estimate the population standardized
mean difference sequentially, minimizing a risk function that adds the
mean square error of the estimator, weighted by `A`, to the total
sampling cost. With \\n\\ observations per group and \\V_n^2 = 2 +
d^2/4\\, sampling stops at the first \\n\\ for which \\n \ge \sqrt{A /
(2c)}\\(V_n + n^{-\gamma})\\, where \\c\\ is `sampling_cost`. With
`pilot = FALSE` the function evaluates that risk and the stopping rule
at the current sample size; with `pilot = TRUE` it returns the pilot
sample size per group, \\\max(m_0, \lceil (A / (2c))^{1/(2 + 2\gamma)}
\rceil)\\. A formula value that is itself a whole number is that number,
although R can compute it a little above (with \\A = 1{,}000{,}000\\,
\\c = 5\\, and \\\gamma = .25\\ it is 100, not 101; see
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md)). See
Kelley and Rausch (2006) for additional information on the standardized
mean difference.

## Note

When `pilot=TRUE` the function returns one row, `pilot_n_per_group`: the
size of the pilot sample *per group* that should be used (thus, the
total sample size is twice the pilot sample size).

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods, 22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# To obtain pilot sample size in a situation in which A=10000. Note that 'A' is
# 'structural_cost' divided by the square of 'epsilon'.

# From Chattopadhyay and Kelley (2017)
mr_smd(pilot = TRUE, A = 10000, sampling_cost = 2.4, gamma = .49)
#>  term              value
#>  pilot_n_per_group 13   

High.SLS <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
Low.SLS <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)

mr_smd(d = 1.021484, n_per_group = 13, A = 10000, sampling_cost = 2.40, gamma = .49)
#>  term value   
#>  risk 1801.521
#>  n_1  13      
#>  n_2  13      
#>  d    1.021   
#> 
#> Stopping rule: not met (the sample size per group, 13, is below the criterion, 81.618); collect another observation in each group.

# Or, using the smd() function:
mr_smd(d = smd(group_1 = High.SLS, group_2 = Low.SLS)[1,2], n_per_group = 13, A = 10000,
       sampling_cost = 2.40, gamma = .49)
#>  term value   
#>  risk 1801.521
#>  n_1  13      
#>  n_2  13      
#>  d    1.021   
#> 
#> Stopping rule: not met (the sample size per group, 13, is below the criterion, 81.618); collect another observation in each group.

# Here, for this situation, the stopping rule is satisfied:
mr_smd(d = 1.00, n_per_group = 75, A = 10000, sampling_cost = 2.40, gamma = .49)
#>  term value
#>  risk 660  
#>  n_1  75   
#>  n_2  75   
#>  d    1    
#> 
#> Stopping rule: met (the sample size per group, 75, is at least the criterion, 73.968); sampling can stop.
```
