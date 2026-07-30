# Community Reinforcement Approach Drinking Trial With Homeless Alcohol-Dependent Individuals (Smith, Meyers, & Delaney, 1998)

Nine-month follow-up drinking outcomes for the *N* = 88 homeless
alcohol-dependent participants in Smith, Meyers, and Delaney's (1998)
randomized clinical trial of the Community Reinforcement Approach (CRA),
published in the *Journal of Consulting and Clinical Psychology*.
Participants were recruited from the Salvation Army Adult Rehabilitation
Center in Albuquerque, New Mexico across two consecutive cohorts and
were randomized to CRA, to CRA augmented with disulfiram, or to standard
care. The outcome is the participant's average number of standard drinks
per week at the nine-month follow-up, reported in both raw form
(markedly right-skewed) and after a base-ten log transformation that
approximately normalizes the distribution used in the original published
analyses. The data are reproduced in Maxwell, Delaney, and Kelley (2027,
*Designing Experiments and Analyzing Data: A Model Comparison
Perspective*, 4th ed., Routledge), Chapter 3, Section 3.10.4, as the
textbook's worked example of a between-subjects analysis of variance
with a heavily skewed outcome.

## Usage

``` r
drinks_trial
```

## Format

A data frame with 88 observations on 5 variables.

- `id`:

  Sequential participant identifier, 1 to 88.

- `cohort`:

  Factor with levels `1` and `2`. The study enrolled two consecutive
  cohorts. Cohort 1 compared three conditions (Standard, CRA, and CRA +
  Disulfiram); Cohort 2 dropped the disulfiram cell on the basis of
  Cohort 1 results and compared Standard against CRA only.

- `treatment`:

  Factor with levels `Standard`, `CRA`, and `CRA + Disulfiram`, the
  randomly assigned treatment condition.

- `drinks_per_week`:

  Average number of standard drinks per week at the nine-month
  follow-up. Bounded below at zero and markedly right-skewed (range 0 to
  624.6, mean 36.9, median 3.8).

- `log_drinks`:

  Common-log transformation \\\log\_{10}(\code{drinks\\per\\week} +
  1)\\, the scale on which Smith, Meyers, and Delaney (1998) ran their
  primary between-groups analyses to obtain approximate normality. The
  plus-one inside the logarithm keeps the zero values finite (and mapped
  to zero).

## Source

Smith, J. E., Meyers, R. J., and Delaney, H. D. (1998). The community
reinforcement approach with homeless alcohol-dependent individuals.
*Journal of Consulting and Clinical Psychology, 66*(3), 541–548.
[doi:10.1037/0022-006X.66.3.541](https://doi.org/10.1037/0022-006X.66.3.541)

Also reproduced in Maxwell, Delaney, and Kelley (2027), Chapter 3.

The data are distributed openly with the companion materials of Maxwell,
Delaney, and Kelley (2027), *Designing Experiments and Analyzing Data: A
Model Comparison Perspective*, and are redistributed here on that basis.

## Details

**Per-cell sample sizes.** The (`cohort`, `treatment`) crosstab is
incomplete by design:

|          |          |     |                  |
|----------|----------|-----|------------------|
|          | Standard | CRA | CRA + Disulfiram |
| Cohort 1 | 17       | 15  | 19               |
| Cohort 2 | 20       | 17  | (not run)        |

Treatment marginals sum to 37 Standard, 32 CRA, and 19 CRA + Disulfiram;
cohort marginals sum to 51 in Cohort 1 and 37 in Cohort 2.

**The authors, the published article, and the textbook.** The trial was
conducted and reported by Jane Ellen Smith, Robert J. Meyers, and Harold
D. Delaney, all then in the Department of Psychology at the University
of New Mexico. Smith and Meyers were the substantive PIs of an extensive
program of CRA research; Meyers (with N. H. Azrin) is widely associated
with the dissemination of CRA and is the developer of the related
Community Reinforcement and Family Training (CRAFT) intervention.
Delaney is a quantitative psychologist who served as the trial's
methodologist. In DMAR the data are included as a worked-example
benchmark because they appear in Maxwell, Delaney, and Kelley (2027),
Chapter 3, Section 3.10.4, as the running example for the consequences
of skipping a normalizing transformation when fitting an analysis of
variance to a markedly skewed outcome.

The original published article is

- Smith, J. E., Meyers, R. J., and Delaney, H. D. (1998). The community
  reinforcement approach with homeless alcohol-dependent individuals.
  *Journal of Consulting and Clinical Psychology, 66*(3), 541–548.
  [doi:10.1037/0022-006X.66.3.541](https://doi.org/10.1037/0022-006X.66.3.541)

and the textbook reproduction (with Section 3.10.4 of MDK 2027 devoted
to the worked analysis) is

- Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
  experiments and analyzing data: A model comparison perspective* (4th
  ed.). Routledge.

**The Community Reinforcement Approach (CRA).** CRA is a behavioral
treatment for alcohol-use disorder developed by Nathan Azrin and
colleagues in the 1970s. It uses operant conditioning principles to
align vocational, marital, social, and recreational reinforcers in the
patient's natural environment with abstinence rather than drinking. CRA
had previously been evaluated in housed populations; Smith, Meyers, and
Delaney (1998) extended its evidence base to homeless alcohol-dependent
individuals, the population represented in these data.

**Disulfiram (Antabuse).** The CRA + Disulfiram condition added
disulfiram, a long-standing pharmacological adjunct that produces an
aversive reaction on alcohol consumption. The Cohort 1 finding that
adding disulfiram offered little incremental benefit over CRA alone
motivated dropping the disulfiram cell in Cohort 2.

**Design and analysis.** The published analyses fit a one-way analysis
of variance to `log_drinks` (the raw `drinks_per_week` variable violates
the normality assumption badly enough that the omnibus inference is
misleading without transformation; see Maxwell, Delaney, and Kelley,
2027, Section 3.10.4). Useful contrasts include CRA versus Standard
within each cohort, CRA versus CRA + Disulfiram within Cohort 1 to
estimate the incremental benefit of disulfiram, and pooling across
cohorts to estimate an overall CRA-versus-Standard effect.

**Use as a DMAR benchmark.** This data set is the canonical DMAR example
for one-way between-subjects analysis with a heavily right-skewed
continuous outcome. Methods that pair naturally with the data set
include the one-way analysis of variance, planned contrasts
([`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)),
the standardized mean difference
([`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)),
Cliff's delta
([`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)),
the probability of superiority, and accuracy in parameter estimation
sample size planning
([`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md)).
The contrast between an ANOVA fit to `drinks_per_week` and one fit to
`log_drinks` is a clean teaching example for the consequences of
skipping a normalizing transformation.

**What is and is not here.** The published paper reports drinking
outcomes at multiple follow-up timepoints (2, 6, 9, and 12 months) along
with several demographic and clinical covariates. The values distributed
here cover the nine-month follow-up only and contain no demographic or
covariate information. Anyone wanting the full longitudinal trajectories
or the covariate set should consult Smith, Meyers, and Delaney (1998)
directly.

## References

Smith, J. E., Meyers, R. J., and Delaney, H. D. (1998). The community
reinforcement approach with homeless alcohol-dependent individuals.
*Journal of Consulting and Clinical Psychology, 66*(3), 541–548.
[doi:10.1037/0022-006X.66.3.541](https://doi.org/10.1037/0022-006X.66.3.541)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on between-subjects analysis of variance
and Section 3.10 on transformations.)

Kelley, K. (2026). *DMAR: Methods for the design, measurement, and
analysis of human-centered outcomes in R* \[R package\].
<https://github.com/yelleKneK/DMAR>

## Author

Ken Kelley

## Examples

``` r
data(drinks_trial)
str(drinks_trial)
#> 'data.frame':    88 obs. of  5 variables:
#>  $ id             : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ cohort         : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ treatment      : Factor w/ 3 levels "Standard","CRA",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ drinks_per_week: num  11 21.7 98.3 6.9 19.1 ...
#>  $ log_drinks     : num  1.078 1.357 1.997 0.898 1.303 ...

# Per-cell sample sizes by cohort and treatment.
with(drinks_trial, table(cohort, treatment))
#>       treatment
#> cohort Standard CRA CRA + Disulfiram
#>      1       17  15               19
#>      2       20  17                0

# Right skew of the raw outcome versus the log-transformed scale
# used in the published analyses.
summary(drinks_trial$drinks_per_week)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.000   0.000   3.846  36.879  36.434 624.615 
summary(drinks_trial$log_drinks)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.0000  0.6834  0.8408  1.5724  2.7963 

# One-way analysis of variance on the log scale, treating the
# five filled (cohort, treatment) cells as the design.
fit <- aov(log_drinks ~ cohort:treatment, data = drinks_trial)
summary(fit)
#>                  Df Sum Sq Mean Sq F value Pr(>F)  
#> cohort:treatment  4   9.04  2.2597    3.48 0.0112 *
#> Residuals        83  53.90  0.6494                 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Standardized mean difference (Cohen's d) for every pairwise
# (two-way) comparison of the three treatments, on the normalizing
# log scale. The loop assembles a table of d with its 95%
# confidence interval for each of the three possible pairs.
groups <- split(drinks_trial$log_drinks, drinks_trial$treatment)
cmp    <- combn(names(groups), 2)
smd_table <- do.call(rbind, lapply(seq_len(ncol(cmp)), function(j) {
  g1 <- groups[[cmp[1, j]]]
  g2 <- groups[[cmp[2, j]]]
  d  <- smd(group_1 = g1, group_2 = g2)$value
  ci <- ci_smd(smd = d, n_1 = length(g1), n_2 = length(g2))
  data.frame(
    comparison = paste(cmp[1, j], "vs", cmp[2, j]),
    d          = round(d, 3),
    ci_lower   = round(ci$value[ci$term == "lower_limit"], 3),
    ci_upper   = round(ci$value[ci$term == "upper_limit"], 3)
  )
}))
smd_table
#>                     comparison     d ci_lower ci_upper
#> 1              Standard vs CRA 0.555    0.071    1.036
#> 2 Standard vs CRA + Disulfiram 0.725    0.152    1.292
#> 3      CRA vs CRA + Disulfiram 0.170   -0.400    0.738
```
