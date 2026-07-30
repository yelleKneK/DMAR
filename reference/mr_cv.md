# Minimum Risk Point Estimation of the Population Coefficient of Variation

A function for the sequential estimation of the coefficient of
variations with minimum risk. The function implements the ideas of
Chattopadhyay and Kelley (2016), which considers study cost and accuracy
of the estimated coefficient of variation simultaneously.

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

  the data for which to evaluate the function

- A:

  \\structural_cost/(epsilon^2)\\; this is the structural cost that one
  is willing to pay in a study to estimate the coefficient of variation
  divided by the square of the desired difference (between the estimate
  and the parameter)

- structural_cost:

  The structural cost of what one is willing to pay in a study (see note
  below)

- epsilon:

  The maximum desired difference between the estimated coefficient of
  variation and the population value)

- sampling_cost:

  The sampling cost to collect an additional observation. For example,
  if each survey costs 10 dollars to distribute and score,
  `sampling_cost` would be 10 dollars per additional observation

- pilot:

  `TRUE` or `FALSE` based on whether the users is using the function to
  plan a pilot sample size (`TRUE`) or if it is being used to assess if
  the optimization criterion has been satisfied (`FALSE`)

- m0:

  The minimum bound on the initial pilot sample size

- gamma:

  A correction factor in which we suggest .49; see the two Chattopadhyay
  & Kelley articles for more details (ignorable for most users)

- verbose:

  If `TRUE`, extra information is printed; defaults to `FALSE`

## Value

- risk:

  The value of the risk function

- n:

  The current sample size

- cv:

  The current coefficient of variation

- is_satisfied:

  A TRUE/FALSE statement of whether or not the risk function has been
  satisfied. If TRUE then sampling can stop as the stopping rule has
  been satisfied

## Details

The value of `epsilon` is context specific; the smaller the value the
closer the estimated value will tend to be to the population value.

## Note

When a study's aim is to estimate a parameter accurately, such as the
coefficient of variation, the structural costs and the maximum probable
error of the estimate (i.e., \\\epsilon\\) are combined to form \\A\\.
When we say "what the researcher is willing to pay", we literally mean
the structural cost (\\c\\) the researcher is willing to invest in a
study in order to estimate the parameter of interest with the desired
degree of accuracy. This value is implicitly included (along with
anticipated sampling cost) in grant applications for empirical studies
when a certain amount of money is requested to conduct a study. If a
researcher is willing to pay more and/or desire a smaller value of
\\\epsilon\\, \\A\\ is larger than it would have been. A larger \\A\\
value will translate into a more expensive study, holding everything
else constant. Notice that \\A\\ is a fixed value in any investigation,
as the researcher specifies \\A\\ directly or by specifying its two
components (structural cost and \\\epsilon\\) individually. However,
what is not fixed but rather evaluated in multiple steps throughout the
process is the sampling cost, as it is unknown the necessary sample size
in order to accomplish the study's goal of achieving a sufficiently
accurate estimate of the coefficient of variation. This is the core of
our contributions: minimizing sampling cost, and thereby study cost, by
using a sequential procedure that evaluates a stopping rule using the
risk function to determine if the optimization criterion has been
satisfied (based on the goals of the researcher and current information
available). This function implements the ideas of sampling error and the
study costs are considered simultaneously, so that the cost is not
higher than necessary for the tolerable sampling error.

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

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*, 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md),
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Determine pilot sample size:
mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75, gamma = .49)
#>  term     value
#>  pilot_ss 18   

# Collect data (the size of which is the pilot sample size)
Data <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41, 5, 3)

# Use mr_cv() to assess if the criterion for stopping the sequential study has been satisfied:
mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#>  term value
#>  risk 5960 
#>  n    18   
#>  cv   0.739

# Collect another data (m=1 here) and perform another check:
Data <- c(Data, 44)
mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#>  term value
#>  risk 6220 
#>  n    19   
#>  cv   0.711

# Continue adding obervations, checking each time if m=1, until the minimum risk criteria
# are satisfied:
Data <- c(Data, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23, 40, 18)
mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#>  term value
#>  risk 4890 
#>  n    35   
#>  cv   0.701
```
