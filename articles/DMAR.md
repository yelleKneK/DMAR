# Introduction to DMAR

## What Is DMAR?

DMAR is a more general and greatly expanded reimagining of the **MBESS**
package (Kelley, 2007a, *Journal of Statistical Software*; 2007b,
*Behavior Research Methods*). MBESS itself was originally framed as
*Methods for the Behavioral, Educational, and Social Sciences*, but the
methods it implements became useful far beyond that scope; DMAR ships
under a new name to make that scope expansion explicit. Its emphasis is
on confidence intervals for effect sizes built from noncentral sampling
distributions, sample size planning for both power and accuracy in
parameter estimation (AIPE), reliability, factor analysis, and the
critical values used in classical inference. DMAR is intended to serve
as the computational companion to Maxwell, Delaney, and Kelley’s
*Designing Experiments and Analyzing Data* (4th ed., 2027; henceforth
**MDK**), for use in both instruction and substantive research.

Although the package’s origins are in psychology, the questions it
answers are not field-specific. Researchers in biostatistics
(clinical-trial effect size CIs and AIPE-based trial sizing),
organizational behavior and management science (intervention, training,
and program evaluation), information systems (user-experience A/B
comparisons; technology-adoption studies), education (instructional
effectiveness, longitudinal achievement gains), sociology (program
evaluation, attitudinal surveys), and the methodological literature
itself share the same three tasks: estimating effect sizes, constructing
confidence intervals for them, and planning the sample size of future
studies, for power or for the precision of an estimate.

What DMAR keeps from MBESS is the philosophy, that is, the same
insistence on *effect sizes with uncertainty quantified*, and on
planning studies for *accurate estimation* rather than only for
sufficient power. What DMAR changes is the implementation: a uniform
`snake_case` interface, tidy `data.frame` returns with consistently
named columns, native `ggplot2` visualization, heavier dependencies kept
in Suggests and loaded only by the functions that use them, and a test
suite that checks results against published values and reference
implementations. The function families have also been broadened,
particularly the `ss_power_*` family, which now covers between-subjects,
within-subjects, mixed, and multilevel designs across the chapters of
MDK.

## Why Effect Size Confidence Intervals?

A *p*-value tells you whether an effect is plausibly nonzero in the
direction you tested. It does not tell you how *large* the effect is,
and it does not tell you how *precisely* you have estimated it. Both of
those questions matter scientifically: a study can produce a
“significant” but trivially small effect, or a “nonsignificant” effect
whose confidence interval is so wide that it is consistent with both a
small and a large true value. Reporting an effect size *and* a
confidence interval for it addresses both shortcomings. The *Publication
Manual* of the American Psychological Association calls for effect sizes
and confidence intervals, and the reporting standards of the American
Educational Research Association call for effect sizes and an indication
of their uncertainty, such as a confidence interval (see the review by
Kelley & Preacher, 2012).

For most standardized effect sizes, the sampling distribution is
*noncentral* (e.g., the noncentral *t* underlying Cohen’s $`d`$, the
noncentral *F* underlying $`R^2`$, the noncentral $`\chi^2`$ underlying
RMSEA), which means the exact interval, under the model’s assumptions,
is asymmetric about the point estimate. The symmetric “estimate $`\pm`$
standard error” approximation ignores that skew. DMAR computes the
intervals for these effect sizes from the noncentral distributions
directly, through the routines that
[`ci_nc_t()`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ci_nc_F()`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
and
[`ci_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)
expose.

## Two Philosophies of Sample Size Planning

Researchers who plan studies in DMAR can do so from either of two
complementary perspectives:

- **Power.** Functions in the `ss_power_*` family answer the question
  *“How many subjects do I need so that, if the effect in the population
  is what I expect, I have at least probability $`\pi`$ of rejecting the
  null?”* This is the familiar Cohen-style framing.
- **Accuracy in parameter estimation (AIPE).** Functions in the
  `ss_aipe_*` family answer a different question: *“How many subjects do
  I need so that the confidence interval around my effect size estimate
  is no wider than $`w`$?”* AIPE is appropriate when estimation is the
  primary goal, for example, when an effect is already known to be
  nonzero and the next study should pin down its magnitude rather than
  re-establish significance.

*Accuracy* has a specific meaning in this framework, and the term is
sometimes used loosely. Accuracy considers precision and bias together:
an accurate estimator is one whose sampling distribution is concentrated
*and* centered on the target parameter. Precision alone, in contrast,
can always be improved by trading bias for precision. A trivial
estimator that returns the same constant value for every sample has zero
variance and is therefore arbitrarily “precise,” but it is biased unless
the constant happens to equal the population parameter. The AIPE framing
asks for sample sizes that buy narrow confidence intervals (a precision
target) *without* trading bias to get there. The confidence interval
width target $`w`$ in `ss_aipe_*` is shorthand for that joint accuracy
goal: a CI that is informative because it is both narrow and centered on
the parameter (Kelley & Rausch, 2006; Kelley & Maxwell, 2003; Maxwell,
Kelley, & Rausch, 2008).

The two perspectives often imply different sample sizes for the same
study, and which is larger depends on the width chosen and on the size
of the effect. DMAR provides both side by side so that the choice can be
made deliberately, in light of the scientific question.

## When Do You Reach for DMAR?

The package is organized around six families of tasks. Pick the one that
matches the question you are trying to answer.

### 1. *I Just Ran an ANOVA, t-Test, or Regression. What’s the Effect Size and a CI for It?*

Use the `ci_*` family.

``` r

# Confidence interval for omega squared (proportion of variance accounted
# for) from a one-way ANOVA with F = 11.221 on (4, 50) df, total N = 55.
# F determines the partial value, which in a one-way design is the total.
ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                 estimand = "partial")
```

| effect | omega_squared_partial | omega_squared_partial_ci_lower | omega_squared_partial_ci_upper | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|
| overall | 0.426 | 0.226 | 0.587 | 11.221 | 4 | 50 | 55 |

Confidence level: 95%

``` r


# More naturally, from a fitted model:
fit <- aov(log_drinks ~ treatment, data = drinks_trial)
ci_omega_squared(fit)
```

| effect | omega_squared_total | omega_squared_total_ci_lower | omega_squared_total_ci_upper | omega_squared_partial | omega_squared_partial_ci_lower | omega_squared_partial_ci_upper | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| treatment | 0.069 | 0.002 | 0.207 | 0.069 | 0.002 | 0.207 | 4.261 | 2 | 85 | 88 |

Confidence level: 95%

The fitted model is a one-way ANOVA of the `drinks_trial` data, from the
randomized clinical trial of the Community Reinforcement Approach (CRA)
with 88 homeless alcohol-dependent individuals reported by Smith,
Meyers, and Delaney (1998) and reproduced in MDK (Chapter 3).
Participants were randomized to CRA, to CRA augmented with disulfiram,
or to standard care. The outcome, `log_drinks`, is the common log of one
plus the average number of standard drinks per week at the nine-month
follow-up; the raw averages are markedly right-skewed, and the primary
published analyses were run on this scale. The trial enrolled two
consecutive cohorts and ran the disulfiram arm in the first only. To
keep the demonstration a one-way design, the analysis here and in
section 3 pools the two cohorts. That is a simplification: the contrasts
that
[`?drinks_trial`](https://yelleknek.github.io/DMAR/reference/drinks_trial.md)
suggests compare CRA with CRA plus disulfiram within Cohort 1, the only
cohort in which all three arms ran, and pool the cohorts only for CRA
against standard care.

The output is a tidy `data.frame`: an effect label, the point estimate
(`omega_squared`), the lower and upper limits of the CI, and the inputs
that produced it (`F_value`, `df_effect`, `df_error`, and *N*). That
format can be passed directly to a table or to a forest plot
([`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)
takes it as input).

CI variants exist for many other effect sizes:
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) for
the standardized mean difference (Cohen’s $`d`$);
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md)
for the standardized mean difference using the control-group SD as the
divisor;
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) for
the proportion of variance in the outcome accounted for by the levels of
a factor;
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md) and
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md)
for signal-to-noise ratios;
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
and
[`ci_R()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
for the correlation and the multiple correlation;
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) for the
squared multiple correlation;
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) and
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md) for
unstandardized and standardized regression coefficients;
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)
for the RMSEA model-fit index; and others.

### 2. *I’m Planning a Study. How Many Subjects Do I Need?*

Use `ss_power_*` for sufficient power, or `ss_aipe_*` for sufficient
precision. The same effect can be approached from either angle.

``` r

# Power-based: per-group n needed for 90% power to detect a contrast among
# 4 groups with population means (90, 92, 88, 81) and within-group
# variance 144. The contrast compares the average of the first three
# groups against the fourth.
ss_power_contrast(
  c_weights           = c(1/3, 1/3, 1/3, -1),
  true_means          = c(90, 92, 88, 81),
  true_var_error = 144,
  desired_power       = 0.90
)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 26    |
| total_N               | 104   |
| actual_power          | 0.907 |
| noncentral_t_parm     | 3.312 |
| true_f                | 0.325 |

``` r

# AIPE-based: n needed so a 95% CI around a population standardized mean
# difference of 0.50 is no wider than 0.30 units (a tight interval).
ss_aipe_smd(true_smd = 0.50, conf_level = 0.95, width = 0.30)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 353   |
| true_smd              | 0.500 |
| width                 | 0.300 |

Confidence level: 95%

The `ss_power_*` family spans most of the standard designs covered in
MDK: one-way ANOVA (`ss_power_one_way_anova`), factorial ANOVA
(`ss_power_factorial_anova`), repeated measures ANOVA with sphericity
adjustment (`ss_power_rm_anova`), split-plot mixed ANOVA
(`ss_power_split_plot_anova`), Pearson correlation (`ss_power_r`),
two-independent-groups standardized mean difference (`ss_power_smd`),
contrasts standardized and unstandardized (`ss_power_c`, `ss_power_sc`,
`ss_power_c_ancova`, `ss_power_contrast`), regression coefficients
(`ss_power_reg_coef`), and two-level random-intercept models for
cluster-randomized trials (`ss_power_mixed_effects`). For SEM-style
designs, use
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md).
AIPE counterparts exist for the correspondingly named effects.

### 3. *I Want to Test a Specific Contrast or Set of Pairwise Comparisons.*

Use
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
for arbitrary contrasts in a one-way ANOVA, with optional
multiple-comparison adjustment.

``` r

fit <- aov(log_drinks ~ treatment, data = drinks_trial)
# All pairwise comparisons, Tukey-protected family-wise error rate:
contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
```

| contrast | estimate | se | t_value | df | p_value | p_adjusted | ci_lower | ci_upper |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| CRA - Standard | -0.454 | 0.198 | -2.294 | 85 | 0.0242 | 0.0620 | -0.927 | 0.018 |
| CRA + Disulfiram - Standard | -0.594 | 0.232 | -2.567 | 85 | 0.0120 | 0.0319 | -1.147 | -0.042 |
| CRA + Disulfiram - CRA | -0.140 | 0.238 | -0.589 | 85 | 0.5575 | 0.8264 | -0.707 | 0.427 |

Confidence level: 95%, simultaneous interval (Tukey)

For the bare critical value of a particular procedure, see the `cv_*`
family: [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) (the
studentized maximum modulus),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
and [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md).
These return the threshold value alone, suitable for inclusion in user
code or hand-checking textbook examples.

### 4. *I Have Within-Subjects (Repeated Measures) Data.*

Start with
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
which combines the univariate within-subjects *F* test with a sphericity
diagnostic (Mauchly’s *W*) and the three standard
$`\varepsilon`$-corrected *p*-values (Greenhouse-Geisser, Huynh-Feldt,
lower-bound). The motivation is that an uncorrected within-subjects *F*
test can be substantially anticonservative when the sphericity
assumption fails, and
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md)
puts the diagnostic and the correction in one place.

The `pygmalion` data are the teacher-expectancy data from Rosenthal and
Jacobson’s (1968) *Pygmalion in the Classroom*: intelligence test scores
for 310 elementary school children in grades 1 through 6, of whom 64
were randomly designated to their teachers as likely “intellectual
bloomers” while the remaining 246 served as controls. Each child’s total
IQ was measured three times: at the pretest (`iq_pre`), at an
intermediate follow-up assessment (`iq_4`), and at the end-of-study
follow-up assessment (`iq_8`). The three occasions are the levels of a
within-subjects factor.
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md)
takes the scores in long form, one row per child per occasion, with a
column that identifies the child:

``` r

n_children <- nrow(pygmalion)
occasions  <- c("pretest", "intermediate", "end_of_study")
iq_long <- data.frame(
  id       = factor(rep(seq_len(n_children), times = 3)),
  occasion = factor(rep(occasions, each = n_children), levels = occasions),
  iq       = c(pygmalion$iq_pre, pygmalion$iq_4, pygmalion$iq_8)
)
res <- anova_within(iq_long, id = "id", time = "occasion", outcome = "iq")
res
```

| adjustment         | F_value | df_1  | df_2    | p_value   | epsilon |
|:-------------------|:--------|:------|:--------|:----------|:--------|
| none               | 83.113  | 2     | 618     | \< 0.0001 | NA      |
| Greenhouse-Geisser | 83.113  | 1.928 | 595.813 | \< 0.0001 | 0.964   |
| Huynh-Feldt        | 83.113  | 1.940 | 599.482 | \< 0.0001 | 0.970   |
| lower_bound        | 83.113  | 1     | 309     | \< 0.0001 | 0.500   |

``` r

attr(res, "mauchly")
```

| W     | statistic | df  | p_value | n_subjects | n_levels | method                       |
|:------|:----------|:----|:--------|:-----------|:---------|:-----------------------------|
| 0.963 | 11.689    | 2   | 0.0029  | 310        | 3        | Mauchly’s test of sphericity |

With 310 children, Mauchly’s test is sensitive to small departures from
sphericity (*W* = 0.963, *p* = 0.0029), and the departure here is small:
the Greenhouse-Geisser and Huynh-Feldt estimates of $`\varepsilon`$ are
0.964 and 0.970, near their upper bound of 1, so the corrected tests
differ little from the uncorrected one.

[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md)
keeps the partial $`\eta^2`$ for occasion, a point estimate, in its
`eta_squared_partial` attribute. Its confidence interval comes from
`ci_eta_squared(estimand = "partial")`, which takes the same analysis
fitted by [`aov()`](https://rdrr.io/r/stats/aov.html) with an `Error()`
term for the children and returns one row per effect, naming the error
stratum that the row’s *F* was tested against:

``` r

fit_rm <- aov(iq ~ occasion + Error(id / occasion), data = iq_long)
ci_occasion <- ci_eta_squared(fit_rm, estimand = "partial")
ci_occasion
```

| effect | eta_squared_partial | eta_squared_partial_ci_lower | eta_squared_partial_ci_upper | stratum | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| occasion | 0.212 | 0.157 | 0.265 | id:occasion | 83.113 | 2 | 618 | 930 |

Confidence level: 95%

The interval is built from the uncorrected *F* on 2 and 618 degrees of
freedom (the `none` row of the
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md)
table), so like the uncorrected test it takes sphericity as given; with
$`\varepsilon`$ this near 1, the corrected degrees of freedom are close
to these. The partial $`\eta^2`$ is the proportion of the within-child
variance in IQ that the occasions account for, since its denominator
leaves out the variation due to stable differences among the children.

Mean IQ rose from 98.5 at the pretest to 101.7 at the intermediate
follow-up and 107.9 at the end of the study, Huynh-Feldt corrected
*F*(1.94, 599.48) = 83.11, *p* \< 0.0001, partial $`\eta^2`$ = 0.212,
95% CI \[0.157, 0.265\].

The analysis pools the bloomers and the controls, so it asks whether the
children’s mean IQ changed across the three occasions, not whether the
expectancy manipulation changed it. That question compares the two
conditions, and the analysis of covariance these data are best known for
takes it up in
[`vignette("pygmalion")`](https://yelleknek.github.io/DMAR/articles/pygmalion.md).

For the underlying components separately, see
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md)
and
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md).
For visualizing individual subject trajectories, see
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
and
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).

### 5. *I Want Descriptive Statistics, a Correlation Table, or Distribution-Shape Diagnostics.*

[`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md)
returns a tidy summary suitable for screening. Here it screens six tests
from `holzinger_swineford`, the scores of 301 seventh and eighth grade
pupils at two Chicago elementary schools on 26 ability tests (Holzinger
& Swineford, 1939): three spatial tests and three verbal ones, the six
that the measurement model of the next section takes up. Short names for
the tests keep the correlation table below narrow.

``` r

ability_tests <- setNames(
  holzinger_swineford[, c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                          "t6_paragraph_comprehension", "t7_sentence",
                          "t9_word_meaning")],
  c("visual", "cubes", "lozenges", "paragraph", "sentence", "word_meaning")
)
descriptives(ability_tests)$descriptives
#>       variable    type   n n_missing prop_missing      mean median       sd min
#> 1       visual integer 301         0            0 29.614618     30 7.004593   4
#> 2        cubes integer 301         0            0 24.352159     24 4.709802   9
#> 3     lozenges integer 301         0            0 18.003322     17 9.047835   2
#> 4    paragraph integer 301         0            0  9.182724      9 3.492349   0
#> 5     sentence integer 301         0            0 17.362126     18 5.161889   4
#> 6 word_meaning integer 301         0            0 15.299003     14 7.669222   1
#>   max q25 q75   skewness   kurtosis
#> 1  51  25  34 -0.2569003  0.3553640
#> 2  37  21  27  0.4747983  0.3808109
#> 3  36  11  25  0.3872808 -0.8883744
#> 4  19   7  11  0.2701735  0.1225896
#> 5  28  14  21 -0.3533097 -0.5250063
#> 6  43  10  19  0.8665663  0.8764015
```

Adding `correlations = TRUE` appends a Pearson correlation matrix. For a
publication-ready correlation table with *p*-values, confidence
intervals, and significance stars (and, optionally, HTML or LaTeX
export), use
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md):

``` r

correlations_test(ability_tests, stars = TRUE)
#> Correlations (Pearson, 95% CI)
#> 
#>               visual          cubes           lozenges        paragraph       sentence        word_meaning    
#> --------------------------------------------------------------------------------------------------------------
#> visual        -                                                                                               
#>                                                                                                               
#>                                                                                                               
#>                                                                                                               
#> 
#> cubes         .30***          -                                                                               
#>               p < .0001                                                                                       
#>               [.19, .40]                                                                                      
#>               N = 301                                                                                         
#> 
#> lozenges      .44***          .34***          -                                                               
#>               p < .0001       p < .0001                                                                       
#>               [.34, .53]      [.24, .44]                                                                      
#>               N = 301         N = 301                                                                         
#> 
#> paragraph     .37***          .15**           .16**           -                                               
#>               p < .0001       p = .0079       p = .0058                                                       
#>               [.27, .47]      [.04, .26]      [.05, .27]                                                      
#>               N = 301         N = 301         N = 301                                                         
#> 
#> sentence      .29***          .14*            .08             .73***          -                               
#>               p < .0001       p = .0155       p = .1816       p < .0001                                       
#>               [.19, .39]      [.03, .25]      [-.04, .19]     [.68, .78]                                      
#>               N = 301         N = 301         N = 301         N = 301                                         
#> 
#> word_meaning  .36***          .19***          .20***          .70***          .72***          -               
#>               p < .0001       p = .0008       p = .0006       p < .0001       p < .0001                       
#>               [.25, .45]      [.08, .30]      [.09, .30]      [.64, .76]      [.66, .77]                      
#>               N = 301         N = 301         N = 301         N = 301         N = 301                         
#> 
#> Note. * p < .05, ** p < .01, *** p < .001.
```

The three verbal tests correlate more strongly with one another than any
of them does with a spatial test. The spatial tests hang together less
tightly, so the clustering that the two-factor model of the next section
describes is clearest on the verbal side.

For interrater agreement, use
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
(two raters),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
(multiple raters), or
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md) (intraclass
correlation, all six classical forms). For distribution-shape statistics
in isolation, use
[`skewness()`](https://yelleknek.github.io/DMAR/reference/skewness.md)
and
[`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md).

### 6. *I Am Developing a Measurement Model.*

Use a common-factor model for the shared variation among items.
[`fa_exploratory()`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
returns rotated loadings, factor correlations, and communalities with
confidence intervals, and prints them as a solution: the loadings matrix
with intervals, each item’s complexity, the variance each factor
accounts for, the factor correlations, each factor’s determinacy and
coefficient H with intervals, and the fit, whose RMSEA interval is at 90
percent, the conventional level for that index (`rmsea_conf_level`),
where the other intervals are at 95 percent. With continuous incomplete
observations it fits the raw data by FIML, rather than treating an
estimated correlation matrix as a complete-data sample covariance.

``` r

measurement_items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
exploratory <- fa_exploratory(holzinger_swineford, nfactors = 2,
                              items = measurement_items, seed = 113)
exploratory
```

Exploratory factor analysis: 2 factors, 6 items, N = 301

Estimator: ML; missing data: FIML; rotation: oblimin (oblique); SE:
standard

| item | f1 | f2 | communality | complexity |
|:---|:---|:---|:---|:---|
| t1_visual_perception | 0.546 \[0.391, 0.701\] | 0.229 \[0.088, 0.370\] | 0.426 \[0.286, 0.566\] | 1.341 |
| t2_cubes | 0.443 \[0.304, 0.582\] | 0.052 \[-0.081, 0.185\] | 0.213 \[0.103, 0.323\] | 1.028 |
| t4_lozenges | 0.766 \[0.629, 0.903\] | -0.072 \[-0.125, -0.019\] | 0.559 \[0.333, 0.784\] | 1.018 |
| t6_paragraph_comprehension | 0.040 \[-0.032, 0.112\] | 0.833 \[0.779, 0.888\] | 0.716 \[0.641, 0.791\] | 1.005 |
| t7_sentence | -0.088 \[-0.139, -0.038\] | 0.899 \[0.850, 0.949\] | 0.768 \[0.690, 0.846\] | 1.019 |
| t9_word_meaning | 0.085 \[0.009, 0.161\] | 0.805 \[0.747, 0.863\] | 0.696 \[0.621, 0.771\] | 1.022 |

(standardized pattern loadings and communalities with 95% intervals, and
Hofmann’s complexity; tidy(x) has every parameter in both metrics)

Variance accounted for (sums of squared standardized loadings, the
diagonal of Phi Lambda’ Lambda):

     factor sum_sq_loadings proportion_variance proportion_common
         f1           1.133               0.189             0.335
         f2           2.246               0.374             0.665
      total           3.379               0.563             1.000

Factor correlations with 95% intervals:

                         f1                   f2
    f1                    1 0.302 [0.145, 0.459]
    f2 0.302 [0.145, 0.459]                    1

Factor determinacy and coefficient H with 95% intervals:

     factor          determinacy                min_r                    H
         f1 0.833 [0.775, 0.878] 0.386 [0.201, 0.541] 0.678 [0.563, 0.775]
         f2 0.944 [0.932, 0.955] 0.784 [0.736, 0.824] 0.894 [0.865, 0.917]

(determinacy: the correlation of the regression factor score with its
factor, whose square is the factor’s maximal reliability; min_r: the
lowest correlation two sets of scores equally consistent with the model
can have, 2 x determinacy^2 - 1 (Guttman, 1955); H: Hancock and
Mueller’s reliability of the optimally weighted composite from the
factor’s loadings alone; limits from the logit scale)

Fit: chi square 2.109, df 4, p = 0.7157; RMSEA with 90% interval 0.000
\[0.000, 0.064\]; CFI 1.000; TLI 1.011; SRMR 0.008

Version: robust (auto: the fit uses FIML (missing = “ml”))

Chi square: standard; robust measures: FIML corrected indices of Zhang
and Savalei (2023), version 6 (normal theory)

The fit_measures attribute holds every measure with its version lines.

Confidence level: 95%

Intervals condition on the selected factor count and rotation.

These intervals condition on the chosen factor count and rotation. The
count itself is a judgment, and
[`factor_count()`](https://yelleknek.github.io/DMAR/reference/factor_count.md)
sets what bears on it side by side: for one factor, two, and so on, the
fit of that model and the test of the previous count against it, its
information criteria, the reduced eigenvalue at that rank against its
parallel analysis reference, and Velicer’s minimum average partial. The
parallel analysis reference percentile is a null threshold, not a
sampling confidence interval for an observed eigenvalue or for the
count, and the table applies no cutoff.

``` r

factor_count(holzinger_swineford, items = measurement_items, B = 200,
             seed = 113)
```

Factor counts 1 to 2 for 6 items, N = 301

Estimator: ML; missing data: FIML; unrotated solutions (the fit of a
count does not depend on the rotation)

| factors | solution | chi_square | df | p_value | rmsea | cfi | tli | srmr | bic | delta | weight | p_delta | eigenvalue | reference | map |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 1 | proper | 103.230 | 9 | \< 0.0001 | 0.187 \[0.155, 0.220\] | 0.856 | 0.760 | 0.114 | 11054.5 | 72.6 | 0.000 | NA | 2.444 | 1.032 | 0.110 |
| 2 | proper | 2.109 | 4 | 0.7157 | 0.000 \[0.000, 0.064\] | 1.000 | 1.011 | 0.008 | 10981.9 | 0.0 | 1.000 | \< 0.0001 | 0.568 | 0.170 | 0.102 |

(each count’s model: chi square test, RMSEA with 90% interval, CFI, TLI,
SRMR, BIC with delta (its gap from the lowest) and weight (its relative
weight among the proper solutions); p_delta: the test of the previous
count against this one; eigenvalue: of the reduced correlation matrix,
with its parallel analysis 95th percentile; map: Velicer’s average
squared partial correlation after this many components)

Information criteria over the proper solutions: lowest AIC at 2 factors,
lowest BIC at 2 factors, lowest SABIC at 2 factors.

BIC comparison over the proper solutions (compare_sem()): lowest at 2
factors, weight 1.000; the nearest other count, 1 factor, is 72.6 higher
(weight 0.000, log Bayes factor 36.3 for 2 factors over it).

Parallel analysis (normal reference, B = 200, 95th percentile, reference
N = 301): the reduced eigenvalues exceed their reference through rank 2.

Velicer’s MAP: the average squared partial correlation is smallest after
2 components (0.102; before any component is removed it is 0.164).

Confidence level: 90%

No cutoff is applied and no count is chosen: the count is a judgment
that combines the fit of each model, the difference tests, the
information criteria, the eigenvalues against their reference, and the
MAP with theory and the interpretability of the rotated solutions, which
fa_exploratory() gives at a chosen count. Each model’s interval and
p-values condition on its own count, and the criterion weights are
conditional on the set of proper counts compared.

Once a measurement structure is supported, a latent SEM estimates
relations among constructs without treating factor scores as error-free.
Local structural-after-measurement estimation is another option.
[`sem_sam()`](https://yelleknek.github.io/DMAR/reference/sem_sam.md)
fits measurement blocks first and propagates their uncertainty through
two-step SEs for the structural parameters (Rosseel & Loh, 2024). Its
current scope is complete continuous data; incomplete-data uncertainty
and arbitrary cross-stage derived effects are refused rather than
silently approximated.

``` r

model <- "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
          verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
          verbal ~ visual"
sem_sam(model, holzinger_swineford, estimator = "MLR")
```

| term | lhs | op | rhs | estimate | se | statistic | p_value | ci_lower | ci_upper | conf_level | stage |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| visual =~ t1_visual_perception | visual | =~ | t1_visual_perception | 1.000 | 0.000 | NA | NA | 1.000 | 1.000 | 0.950 | measurement |
| visual =~ t2_cubes | visual | =~ | t2_cubes | 0.519 | 0.103 | 5.045 | \< 0.0001 | 0.317 | 0.720 | 0.950 | measurement |
| visual =~ t4_lozenges | visual | =~ | t4_lozenges | 1.476 | 0.305 | 4.843 | \< 0.0001 | 0.879 | 2.074 | 0.950 | measurement |
| verbal =~ t6_paragraph_comprehension | verbal | =~ | t6_paragraph_comprehension | 1.000 | 0.000 | NA | NA | 1.000 | 1.000 | 0.950 | measurement |
| verbal =~ t7_sentence | verbal | =~ | t7_sentence | 1.511 | 0.093 | 16.218 | \< 0.0001 | 1.328 | 1.693 | 0.950 | measurement |
| verbal =~ t9_word_meaning | verbal | =~ | t9_word_meaning | 2.156 | 0.138 | 15.614 | \< 0.0001 | 1.886 | 2.427 | 0.950 | measurement |
| verbal ~ visual | verbal | ~ | visual | 0.264 | 0.062 | 4.251 | \< 0.0001 | 0.142 | 0.386 | 0.950 | structural |
| t1_visual_perception  t1_visual_perception | t1_visual_perception |  | t1_visual_perception | 30.047 | 4.833 | 6.217 | \< 0.0001 | 20.574 | 39.520 | 0.950 | measurement |
| t2_cubes  t2_cubes | t2_cubes |  | t2_cubes | 17.039 | 1.742 | 9.778 | \< 0.0001 | 13.624 | 20.454 | 0.950 | measurement |
| t4_lozenges  t4_lozenges | t4_lozenges |  | t4_lozenges | 40.497 | 8.063 | 5.023 | \< 0.0001 | 24.694 | 56.300 | 0.950 | measurement |
| t6_paragraph_comprehension  t6_paragraph_comprehension | t6_paragraph_comprehension |  | t6_paragraph_comprehension | 3.435 | 0.462 | 7.434 | \< 0.0001 | 2.530 | 4.341 | 0.950 | measurement |
| t7_sentence  t7_sentence | t7_sentence |  | t7_sentence | 6.658 | 0.945 | 7.047 | \< 0.0001 | 4.807 | 8.510 | 0.950 | measurement |
| t9_word_meaning  t9_word_meaning | t9_word_meaning |  | t9_word_meaning | 18.068 | 2.359 | 7.659 | \< 0.0001 | 13.444 | 22.692 | 0.950 | measurement |
| visual  visual | visual |  | visual | 18.854 | 5.278 | 3.572 | 0.0004 | 8.509 | 29.199 | 0.950 | structural |
| verbal  verbal | verbal |  | verbal | 7.405 | 0.977 | 7.583 | \< 0.0001 | 5.491 | 9.319 | 0.950 | structural |

Confidence level: 95%

Local SAM: MLR; SE: twostep.robust; complete continuous data;
gamma_unbiased = TRUE. Intervals concern individual original-unit
parameters. Cross-stage covariance cells are incomplete; do not use the
attached fit for arbitrary derived or standardized structural intervals.

ML and MLR use the same likelihood objective here. Their SEs differ
because normal-theory and sandwich inference make different covariance
assumptions, not because one point estimate was computed iteratively. A
bootstrap is a separate uncertainty procedure and should be named
explicitly when used.

## What DMAR Is *Not*

DMAR deliberately overlaps with, but does not duplicate, several
specialized packages.

- For full-featured factor analysis and SEM, use `lavaan`.
- For coefficient $`\alpha`$ / $`\omega`$ with extensive options, use
  `psych::omega()`.
- For mixed-effects models, use `lme4` or `nlme`.
- For robust ANOVA, use `WRS2`.

DMAR layers tidy effect size CIs, sample size planning, and
publication-ready output over those tools, rather than replacing them.

## A Note on Argument Naming

The package uses `snake_case` throughout. Capital letters in identifiers
are reserved for two uses. The first is statistically meaningful
symbols, for example `R2` for the squared multiple correlation, `N` for
total sample size, `S` for a covariance matrix, `F_value` for an *F*
statistic, `B` and `G` for the numbers of bootstrap and Monte Carlo
replications, and `M` for item means. The second is the labels of
variables, factors, and groups, for example `X`, `Y`, and `Z` (as in
`true_mean_Y` and `true_cor_YX`), `A` and `B` for the factors of a
two-factor design (as in `factor_A`), and `C` for a control group (as in
`mean_C`). A population value given to a planner carries the prefix
`true_` and an English name for the quantity, and the planning value a
sensitivity function chooses its design under carries the prefix
`supposed_`, so examples read consistently regardless of discipline:

``` r

ci_R2(R2 = 0.25, N = 100, p = 5, random_predictors = TRUE)
ss_power_reg_coef(true_R2_Y_X = 0.78, true_R2_Y_X_without_j = 0.74,
                  p = 5, desired_power = 0.85)
ss_aipe_cv_sensitivity(true_cv = 0.25, supposed_cv = 0.25,
                       width = 0.10, conf_level = 0.95, G = 200)
```

DMAR does not accept the MBESS dot.case style (`conf.level`,
`Random.Predictors`, `Specified.N`). Scripts written against MBESS will
need to convert dots to underscores in argument names and lowercase any
prefix words that are not themselves statistical notation. Many function
names carry over with underscores in place of dots (`ci.smd()` is
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
`ss.aipe.reg.coef()` is
[`ss_aipe_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)),
and results are returned as data frames rather than lists (see
[`vignette("mbess_to_dmar")`](https://yelleknek.github.io/DMAR/articles/mbess_to_dmar.md)).

If you have scripts written against MBESS and you see errors of the form

    Error in foo(...) : unused argument (Group.1 = ...)

the first step is to replace dots with underscores in the argument name
(`Group_1`, `n_1`, `Mean_1`, `Specified_N`, etc.) and then lowercase the
capitalized prefix words that are not statistical notation (`group_1`,
`mean_1`, `specified_N`, …). Some arguments are also named for the role
the value plays: the population standardized mean difference that
MBESS’s `ss.aipe.smd()` takes as `delta` is `true_smd` in
[`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md).

The single convention exists for one reason: a mix of dot.case and
snake_case arguments is confusing in side-by-side examples, especially
in teaching.

## Where to Next

- For a worked walkthrough of the visualization functions and effect
  size confidence intervals on the `pygmalion` teacher expectancy data,
  see
  [`vignette("effect-size-visualization", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/effect-size-visualization.md).
- For when full information maximum likelihood regression actually helps
  over [`lm()`](https://rdrr.io/r/stats/lm.html), with
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md),
  see
  [`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md).
- For the function reference, `?DMAR-package` lists the exported
  functions grouped by family.
- Citation: `citation("DMAR")`.
- Author web site: <https://kenkelley.org>; related publications:
  <https://kenkelley.org/publications/>.

## References

Holzinger, K. J., & Swineford, F. (1939). *A study in factor analysis:
The stability of a bi-factor solution* (Supplementary Educational
Monographs, No. 48). University of Chicago Press.

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. <https://doi.org/10.18637/jss.v020.i08>

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
<https://doi.org/10.3758/BF03192993>

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*(2), 137–152.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
<https://doi.org/10.1146/annurev.psych.59.103006.093735>

Rosenthal, R., & Jacobson, L. (1968). *Pygmalion in the classroom:
Teacher expectation and pupils’ intellectual development*. Holt,
Rinehart and Winston.

Rosseel, Y., & Loh, W. W. (2024). A structural after measurement
approach to structural equation modeling. *Psychological Methods,
29*(3), 561–588. <https://doi.org/10.1037/met0000503>

Smith, J. E., Meyers, R. J., & Delaney, H. D. (1998). The community
reinforcement approach with homeless alcohol-dependent individuals.
*Journal of Consulting and Clinical Psychology, 66*(3), 541–548.
<https://doi.org/10.1037/0022-006X.66.3.541>
