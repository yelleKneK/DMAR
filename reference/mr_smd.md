# Minimum Risk Point Estimation of the Population Standardized Mean Difference

A function for the sequential estimation of the standardized mean
difference with minimum risk. The function implements the ideas of
Chattopadhyay and Kelley (2017), which considers study cost and accuracy
of the estimated standardized mean difference simultaneously. This is
important to specify that `mr_smd` was developed under the assumption of
normally distributed data with equal sample size and equal cost of
sampling per observation for each group.

## Usage

``` r
mr_smd(
  A,
  structural_cost,
  epsilon,
  d,
  n,
  sampling_cost,
  pilot = FALSE,
  m0 = 4,
  gamma = 0.49
)
```

## Arguments

- A:

  The price one is willing to pay in order to have a maximum allowable
  difference of \\epsilon^2\\ between the estimate of the standardized
  mean difference and its corresponding parameter

- structural_cost:

  The structural cost of what one is willing to pay in a study

- epsilon:

  The maximum desired difference between the estimated standardized mean
  difference and the population value

- d:

  The current estimate of the standardized mean difference

- n:

  Current sample size *per group* (thus total sample size is \\2n\\);
  requires equal sample size *per group*

- sampling_cost:

  The sampling cost to collect an additional observation. For example,
  if each survey costs 10 dollars to distribute and score,
  `sampling_cost` would be 10 dollars per additional observation

- pilot:

  `TRUE` or `FALSE` based on whether the users is using the function to
  plan a pilot sample size (TRUE) or if it is being used to assess if
  the optimization criterion has been satisfied (FALSE)

- m0:

  The minimum bound on the initial pilot sample size

- gamma:

  A correction factor in which we suggest .49; see the two Chattopadhyay
  & Kelley articles for more details (ignorable for most users)

## Value

- risk:

  The value of the risk function.

- n1:

  Sample size for group 1 (echos the input value)

- n2:

  Sample size for group 2 (echos the input value)

- d:

  Observed value of the standardized mean difference (i.e., *d*; echos
  the input value)

- is_satisfied:

  A `TRUE` or `FALSE` statement of that evaluates a stopping rule using
  the risk function to determine if the optimization criterion has been
  satisfied (based on the goals of the researcher and current
  information available)

## Details

The standardized mean difference is a widely used measure effect size.
In this article, we developed a general theory for estimating the
population standardized mean difference by minimizing both the mean
square error of the estimator and the total sampling cost. This function
implements our ideas discussed in Chattopadhyay and Kelley (2017). See
also Kelley and Rausch (2006) for additional information on the
standardized mean difference.

## Note

When `pilot=TRUE` the function returns the size of the pilot sample
size, *per group*, that should be used (thus, the total sample size is
twice the pilot sample size).

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
approach. *Psychological Methods, 23*, 226–243.
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
#>  term     value
#>  pilot_ss 13   

High.SLS <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
Low.SLS <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)

mr_smd(d = 1.021484, n = 13, A = 10000, sampling_cost = 2.40, gamma = .49)
#>  term value
#>  risk 1800 
#>  n1   13   
#>  n2   13   
#>  d    1.02 

# Or, using the smd() function:
mr_smd(d = smd(group_1 = High.SLS, group_2 = Low.SLS)[1,2], n = 13, A = 10000,
       sampling_cost = 2.40, gamma = .49)
#>  term value
#>  risk 1800 
#>  n1   13   
#>  n2   13   
#>  d    1.02 

# Here, for this situation, the stopping rule is satisfied:
mr_smd(d = 1.00, n = 75, A = 10000, sampling_cost = 2.40, gamma = .49)
#>  term value
#>  risk 660  
#>  n1   75   
#>  n2   75   
#>  d    1    
```
