# Visualizing Effect Sizes and Distributions With DMAR

## Introduction

DMAR (*Design, Measurement, and Analysis in R*) is a modern and greatly
expanded reimagining of the widely used MBESS package (Kelley, 2007a,
2007b). Both packages share the same statistical foundations (confidence
intervals for standardized effect sizes, sample size planning from the
accuracy in parameter estimation (AIPE) perspective, and measurement),
but DMAR adopts a tidy, data-frame-first output style that plays well
with modern R workflows.

Although the worked example below draws from educational psychology, the
same tools apply equally to organizational behavior research (comparing
intervention groups on engagement, commitment, or performance),
biostatistics (group comparisons of clinical outcomes), information
systems (user-experience or adoption studies), management science
(process or training-program evaluations), sociology (program
evaluations on attitudinal scales), and any other discipline where the
inferential question is “how large is the effect, and how precisely have
we estimated it?”

MBESS (Kelley, 2007a, 2007b) continues to be available on CRAN and
remains a reliable choice, especially for users with existing scripts.
DMAR is designed for new projects and teaching, where a consistent, tidy
interface makes it easier to integrate effect size estimation with data
visualization and reporting.

This vignette shows how to:

1.  Screen data with
    [`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md).
2.  Visualize group distributions with raincloud plots (via the
    **ggrain** package).
3.  Compute and visualize standardized mean differences, $`R^2`$, and
    ANOVA effect sizes with
    [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
    [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
    and
    [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md).

Every plot produced by DMAR includes a confidence interval and the
sample size on which the estimate is based, two details that are
essential for transparent reporting (Cumming, 2012; Kelley & Preacher,
2012). Both annotations are on by default and can be turned off with
`show_ci = FALSE` and `show_n = FALSE`.

The implicit standard DMAR enforces is that an effect size, on its own,
is not enough information to interpret. Without an interval, the reader
cannot judge how precisely the effect has been estimated; without a
sample size, the reader cannot judge how seriously the interval should
be taken. The default annotations make these three pieces of information
travel together. When you do choose to suppress them (typically for a
slide or a teaching figure that needs to be visually minimal), the
suppression is a deliberate editorial decision rather than an oversight.

## Setup

``` r

library(DMAR)
```

``` r

library(ggplot2)
```

## The Data

We simulate data from a hypothetical study on the effect of a
growth-mindset intervention on academic motivation among college
students. One hundred and twenty participants were randomly assigned to
a mindfulness-based growth mindset intervention ($`n = 60`$) or a
waitlist control ($`n = 60`$). Four measures were collected at
post-test:

| Variable        | Description                    | Scale |
|:----------------|:-------------------------------|:------|
| `motivation`    | Academic Motivation Scale      | 1–7   |
| `self_efficacy` | General Self-Efficacy Scale    | 1–5   |
| `test_anxiety`  | Test Anxiety Inventory         | 20–80 |
| `gpa`           | Cumulative grade-point average | 0–4   |

``` r

set.seed(113)
n_per_group <- 60

study_data <- data.frame(
  participant   = seq_len(2 * n_per_group),
  group         = factor(rep(c("Control", "Intervention"), each = n_per_group)),
  motivation    = c(rnorm(n_per_group, mean = 4.20, sd = 0.90),
                    rnorm(n_per_group, mean = 4.70, sd = 0.90)),
  self_efficacy = c(rnorm(n_per_group, mean = 3.40, sd = 0.65),
                    rnorm(n_per_group, mean = 3.60, sd = 0.65)),
  test_anxiety  = c(rnorm(n_per_group, mean = 48, sd = 10),
                    rnorm(n_per_group, mean = 44, sd = 10))
)

# GPA as a function of motivation and test anxiety.
study_data$gpa <- with(study_data,
  2.0 + 0.25 * motivation - 0.015 * test_anxiety + rnorm(2 * n_per_group, sd = 0.30)
)
study_data$gpa <- pmin(pmax(study_data$gpa, 0), 4.0)

knitr::kable(head(study_data), digits = 2)
```

| participant | group   | motivation | self_efficacy | test_anxiety |  gpa |
|------------:|:--------|-----------:|--------------:|-------------:|-----:|
|           1 | Control |       4.32 |          2.90 |        50.25 | 2.15 |
|           2 | Control |       5.44 |          3.59 |        47.34 | 2.74 |
|           3 | Control |       4.87 |          3.88 |        49.61 | 2.12 |
|           4 | Control |       3.04 |          2.44 |        58.87 | 2.43 |
|           5 | Control |       3.70 |          3.33 |        46.00 | 2.34 |
|           6 | Control |       2.64 |          2.40 |        32.41 | 1.77 |

## Descriptive Statistics

The
[`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md)
function provides a compact summary useful for data screening and
psychometric work. It reports per-variable sample size, missingness,
central tendency, spread, skewness, and excess kurtosis, all in a single
tidy data frame.

``` r

desc <- descriptives(
  study_data[, c("motivation", "self_efficacy", "test_anxiety", "gpa")]
)
knitr::kable(desc$descriptives, digits = 3)
```

| variable | type | n | n_missing | prop_missing | mean | median | sd | min | max | q25 | q75 | skewness | kurtosis |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| motivation | numeric | 120 | 0 | 0 | 4.465 | 4.503 | 0.932 | 2.168 | 6.450 | 3.893 | 5.056 | -0.204 | -0.336 |
| self_efficacy | numeric | 120 | 0 | 0 | 3.441 | 3.430 | 0.637 | 1.875 | 5.413 | 3.011 | 3.789 | 0.204 | 0.141 |
| test_anxiety | numeric | 120 | 0 | 0 | 45.878 | 45.533 | 10.012 | 21.558 | 69.855 | 38.652 | 53.067 | 0.080 | -0.372 |
| gpa | numeric | 120 | 0 | 0 | 2.460 | 2.461 | 0.413 | 1.369 | 3.501 | 2.182 | 2.743 | -0.222 | -0.267 |

Skewness and kurtosis values close to zero suggest approximate
normality. As a rough guideline, $`|\text{skewness}| > 2`$ or
$`|\text{kurtosis}| > 7`$ may warrant concern for normal-theory methods.

Adding `correlations = TRUE` appends a correlation matrix, which is
handy during scale development or when checking multicollinearity.

``` r

desc_cor <- descriptives(
  study_data[, c("motivation", "self_efficacy", "test_anxiety", "gpa")],
  correlations = TRUE
)
knitr::kable(desc_cor$correlations, digits = 3)
```

|               | motivation | self_efficacy | test_anxiety |    gpa |
|:--------------|-----------:|--------------:|-------------:|-------:|
| motivation    |      1.000 |         0.215 |       -0.018 |  0.646 |
| self_efficacy |      0.215 |         1.000 |        0.001 |  0.101 |
| test_anxiety  |     -0.018 |         0.001 |        1.000 | -0.327 |
| gpa           |      0.646 |         0.101 |       -0.327 |  1.000 |

## Visualizing Distributions: Raincloud Plots

Raincloud plots (Allen et al., 2019) combine a half-violin, jittered raw
data points, and a boxplot into a single display. They show the full
distributional shape, individual observations, and summary statistics at
once. The **ggrain** package (Patil, 2023) makes them easy to produce
with a single
[`geom_rain()`](https://rdrr.io/pkg/ggrain/man/geom_rain.html) call.

``` r

library(ggrain)

# Source group colors from DMAR's own colorblind-safe palette engine.
group_fills <- setNames(unname(grDevices::palette.colors(2)), c("Control", "Intervention"))

ggplot(study_data, aes(x = group, y = motivation, fill = group)) +
  geom_rain(alpha = 0.5) +
  scale_fill_manual(values = group_fills) +
  labs(
    title = "Academic Motivation by Group",
    x     = NULL,
    y     = "Academic Motivation Scale (1\u20137)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))
```

![Raincloud plot of academic motivation by group, combining a
half-violin density, jittered raw scores, and a boxplot summary (Allen
et al.,
2019).](effect-size-visualization_files/figure-html/raincloud-motivation-1.png)

Raincloud plot of academic motivation by group, combining a half-violin
density, jittered raw scores, and a boxplot summary (Allen et al.,
2019).

``` r

ggplot(study_data, aes(x = group, y = test_anxiety, fill = group)) +
  geom_rain(alpha = 0.5) +
  scale_fill_manual(values = group_fills) +
  labs(
    title = "Test Anxiety by Group",
    x     = NULL,
    y     = "Test Anxiety Inventory (20\u201380)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))
```

![Raincloud plot of test anxiety by group. Reading across the two
distributions, the intervention group's anxiety distribution is shifted
somewhat lower and the spread is
comparable.](effect-size-visualization_files/figure-html/raincloud-anxiety-1.png)

Raincloud plot of test anxiety by group. Reading across the two
distributions, the intervention group’s anxiety distribution is shifted
somewhat lower and the spread is comparable.

The raincloud plots give an immediate qualitative impression of group
differences. The next sections quantify those differences with
standardized effect sizes and confidence intervals.

## Standardized Mean Difference

The standardized mean difference (Cohen’s $`d`$) expresses a group
difference in standard-deviation units, making it comparable across
studies that use different measurement scales. DMAR’s
[`smd()`](https://yelleknek.github.io/DMAR/reference/smd.md) computes
the biased (Cohen) or unbiased (Hedges) estimate;
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) wraps
it in a noncentral $`t`$-based confidence interval.

``` r

# Split the data by group.
ctrl <- study_data$motivation[study_data$group == "Control"]
intv <- study_data$motivation[study_data$group == "Intervention"]

# Point estimate (Hedges' unbiased g).
smd(group_1 = intv, group_2 = ctrl, unbiased = TRUE)
```

| term | value |
|:-----|:------|
| smd  | 0.207 |

``` r


# 95% confidence interval.
d_hat <- smd(group_1 = intv, group_2 = ctrl)$value
d_ci  <- ci_smd(smd = d_hat, n_1 = length(intv), n_2 = length(ctrl))
d_ci
```

| term        | value  |
|:------------|:-------|
| lower_limit | -0.151 |
| smd         | 0.208  |
| upper_limit | 0.567  |

Confidence level: 95%

``` r


# The results-section sentence, written by the package so the numbers
# can never drift from the table they came from.
results_sentence(d_ci, label = "the standardized mean difference")
#> [1] "the standardized mean difference = 0.21, 95% CI [-0.15, 0.57]"
```

Reported as one would in a results section, the intervention group did
not differ reliably from control on academic motivation: the
standardized mean difference = 0.21, 95% CI \[-0.15, 0.57\]. Because the
interval includes zero, the data are consistent with no effect, which
matches the omnibus test below ($`p = 0.256`$). The point estimate alone
would overstate what the study established; the interval is what makes
the uncertainty visible.

### Visualizing With `plot_smd()`

[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)
draws two unit-variance normal distributions separated by $`d`$. The
overlap makes the practical meaning of the effect size immediately
visible. A confidence interval and the per-group sample sizes are
displayed by default.

``` r

plot_smd(
  group_1      = intv,
  group_2      = ctrl,
  group_labels = c("Intervention", "Control"),
  title        = "Motivation: Intervention vs. Control"
)
```

![Cohen's \\d\\ with 95\\ confidence interval and per-group sample
sizes, displayed as two unit-variance normal curves separated by the
observed effect
size.](effect-size-visualization_files/figure-html/plot-smd-1.png)

Cohen’s $`d`$ with 95% confidence interval and per-group sample sizes,
displayed as two unit-variance normal curves separated by the observed
effect size.

``` r

# You can also pass a known d and sample sizes directly.
plot_smd(smd = 0.80, n_1 = 30, n_2 = 30,
         title = "What Does d = 0.80 Look Like?")
```

![An effect of d = 0.80, near the top of the values Cohen (1988)
tabulated, still leaves substantial overlap between the two
unit-variance normal
distributions.](effect-size-visualization_files/figure-html/plot-smd-known-1.png)

An effect of d = 0.80, near the top of the values Cohen (1988)
tabulated, still leaves substantial overlap between the two
unit-variance normal distributions.

The visual makes it clear that even an effect of $`d = 0.80`$, near the
top of the values Cohen (1988) tabulated, leaves substantial overlap
between the two distributions, a useful corrective to the common
misconception that a significant $`p`$-value implies complete
separation.

Numeric reference values for $`d`$ such as $`0.20`$, $`0.50`$, and
$`0.80`$ are sometimes treated as defaults when no field-specific
calibration is available, but they are not a substitute for substantive
interpretation. A $`d`$ of $`0.50`$ in a randomized educational
intervention may be remarkable; the same $`d`$ in a basic perception
experiment may be uninteresting. Pairing the plot with a confidence
interval prevents the second mistake that often accompanies the first,
that is, reporting a point estimate as if it were known with certainty.
Pairing it with the per-group sample sizes makes the precision visible
at a glance.

## ANOVA Effect Sizes

### Omega Squared ($`\omega^2`$)

Omega squared estimates the proportion of variance in the population
accounted for by the fixed effect, correcting for the upward bias of
$`\eta^2`$ (Olejnik & Algina, 2003). DMAR’s
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
accepts either raw ANOVA summary values or a fitted
[`aov()`](https://rdrr.io/r/stats/aov.html) object.

``` r

fit <- aov(motivation ~ group, data = study_data)
print_anova(summary(fit)[[1]])
#>              Df     Sum Sq   Mean Sq  F value Pr(>F)
#> group         1   1.130034 1.1300340 1.303211 0.2559
#> Residuals   118 102.319620 0.8671154       NA   <NA>
```

``` r

omega_result <- ci_omega_squared(fit)
#> Warning: The observed F_value is below the alpha_lower critical value of the
#> central F-distribution; the lower noncentrality limit has been clamped to 0 and
#> the reported 'prob_greater' on the lower_limit row reflects the actual
#> upper-tail probability at lambda = 0.
omega_result
```

| effect | omega_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N   |
|:-------|:--------------|:------------|:------------|:--------|:----------|:---------|:----|
| group  | 0.00252       | 0           | 0.0743      | 1.3     | 1         | 118      | 120 |

### Visualizing With `plot_ci()`

[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)
creates a forest-plot-style display that works with any DMAR confidence
interval output. When it receives output from
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
it automatically extracts effect names, point estimates, confidence
bounds, and sample sizes.

``` r

plot_ci(omega_result,
        reference_line = 0,
        xlab  = expression(omega^2),
        title = "ANOVA Effect Size: Motivation ~ Group")
```

![Forest-plot display of partial \\\omega^2\\ with 95\\ noncentral \\F\\
confidence interval (Steiger, 2004) for the group effect on
motivation.](effect-size-visualization_files/figure-html/plot-ci-omega-1.png)

Forest-plot display of partial $`\omega^2`$ with 95% noncentral $`F`$
confidence interval (Steiger, 2004) for the group effect on motivation.

For factorial designs,
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)
produces a multi-row forest plot with one row per effect:

``` r

# Using the built-in warpbreaks data (2 x 3 factorial).
fit_factorial <- aov(breaks ~ wool * tension, data = warpbreaks)
omega_factorial <- ci_omega_squared(fit_factorial)
#> Warning: The observed F_value is below the alpha_lower critical value of the
#> central F-distribution; the lower noncentrality limit has been clamped to 0 and
#> the reported 'prob_greater' on the lower_limit row reflects the actual
#> upper-tail probability at lambda = 0.
omega_factorial
```

| effect | omega_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N |
|:---|:---|:---|:---|:---|:---|:---|:---|
| wool | 0.0487 | 0 | 0.222 | 3.77 | 1 | 48 | 54 |
| tension | 0.217 | 0.0558 | 0.411 | 8.5 | 2 | 48 | 54 |
| wool:tension | 0.106 | 0.00191 | 0.298 | 4.19 | 2 | 48 | 54 |

``` r


plot_ci(omega_factorial,
        reference_line = 0,
        xlab  = expression(omega^2),
        title = "Warpbreaks: Partial Omega Squared per Effect")
```

![Multi-row forest plot of partial \\\omega^2\\ for each effect in a 2x3
factorial ANOVA on the warpbreaks data. Reading down the rows is the
recommended diagnostic for ANOVA model summaries (Maxwell, Delaney, &
Kelley,
2027).](effect-size-visualization_files/figure-html/plot-ci-factorial-1.png)

Multi-row forest plot of partial $`\omega^2`$ for each effect in a 2x3
factorial ANOVA on the warpbreaks data. Reading down the rows is the
recommended diagnostic for ANOVA model summaries (Maxwell, Delaney, &
Kelley, 2027).

A multi-row forest plot of this kind is read row by row. Each row gives
the partial $`\omega^2`$ for one effect, with the horizontal bar
covering the confidence interval and the dot marking the point estimate.
The `reference_line = 0` argument draws a vertical guide at zero
variance explained: any interval that crosses that line is consistent
with no effect of that term. Reading down the rows lets you compare the
relative magnitudes of the effects at a glance, and the per-row
confidence intervals communicate which of those comparisons are reliable
and which are not.

## Proportion of Variance Explained: $`R^2`$

The squared multiple correlation $`R^2`$ answers: “What proportion of
the variance in the outcome is linearly predictable from the set of
predictors?”
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)
displays this as a horizontal proportion bar, making the magnitude
immediately interpretable.

``` r

reg_fit <- lm(gpa ~ motivation + test_anxiety, data = study_data)
print_summary(reg_fit)
#> Coefficients:
#>                 Estimate  Std. Error   t value Pr(>|t|)
#> (Intercept)   1.79056606 0.179526638  9.973818 < 0.0001
#> motivation    0.28374794 0.028478156  9.963705 < 0.0001
#> test_anxiety -0.01302498 0.002652155 -4.911094 < 0.0001
#> 
#> Residual standard error: 0.2896 on 117 degrees of freedom
#> Multiple R-squared: 0.5169,  Adjusted R-squared: 0.5087
#> F-statistic:  62.6 on 2 and 117 DF, p-value: < 0.0001
```

The number of predictors is a single quantity here (there are two:
motivation and test anxiety), and every sibling function names it the
same way: `p`, whether in
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md), or
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md). We store
it once and pass it along.

``` r

R2_obs <- summary(reg_fit)$r.squared
N      <- nrow(study_data)
p      <- 2  # number of predictors

ci_R2(R2 = R2_obs, N = N, p = p, random_predictors = TRUE)
```

| term        | value | prob_less | prob_greater |
|:------------|:------|:----------|:-------------|
| lower_limit | 0.377 | 0.025     | 0.975        |
| R2          | 0.517 | NA        | NA           |
| upper_limit | 0.628 | 0.975     | 0.025        |

Confidence level: 95%

``` r

plot_R2(R2 = R2_obs, N = N, p = p,
        title = "GPA ~ Motivation + Test Anxiety")
```

![Squared multiple correlation \\R^2\\ with 95% confidence interval for
the regression of GPA on motivation and test
anxiety.](effect-size-visualization_files/figure-html/plot-r2-1.png)

Squared multiple correlation $`R^2`$ with 95% confidence interval for
the regression of GPA on motivation and test anxiety.

The bar makes the “glass half-full, glass half-empty” nature of $`R^2`$
vivid: even an $`R^2`$ that is statistically significant may leave
substantial unexplained variance. The accompanying confidence interval
makes a second, equally important point. With moderate sample sizes, the
CI on $`R^2`$ is typically wide. It is not unusual for an observed
$`R^2 = 0.25`$ from $`N = 100`$ to be consistent with population values
ranging from roughly $`0.10`$ to $`0.40`$. Treating the point estimate
as the answer, while ignoring that range, overstates what the data
actually show.

### Correlation Confidence Interval

For the multiple correlation $`R`$ itself,
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md) provides
a confidence interval that can be displayed with
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md):

``` r

R_obs <- sqrt(R2_obs)
r_ci  <- ci_r(R = R_obs, N = N, p = p)
r_ci
```

| term        | value | prob_less | prob_greater |
|:------------|:------|:----------|:-------------|
| lower_limit | 0.614 | 0.025     | 0.975        |
| R           | 0.719 | NA        | NA           |
| upper_limit | 0.792 | 0.975     | 0.025        |

Confidence level: 95%

``` r


plot_ci(r_ci,
        estimate = R_obs,
        n = N,
        reference_line = 0,
        xlab  = "Multiple R",
        title = "Confidence Interval for the Multiple Correlation")
```

![Confidence interval for the multiple correlation \\R\\, constructed
from the random-predictor sampling distribution of \\R^2\\ (Lee, 1971),
the default for
\`ci_r()\`.](effect-size-visualization_files/figure-html/ci-r-1.png)

Confidence interval for the multiple correlation $`R`$, constructed from
the random-predictor sampling distribution of $`R^2`$ (Lee, 1971), the
default for
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md).

## Combining Multiple Effect Sizes

[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)
also accepts explicit vectors, so you can build a side-by-side
comparison of different effects:

``` r

# Effect sizes from our study.
d_motivation  <- smd(group_1 = intv, group_2 = ctrl)$value
d_self_eff    <- smd(
  group_1 = study_data$self_efficacy[study_data$group == "Intervention"],
  group_2 = study_data$self_efficacy[study_data$group == "Control"]
)$value
d_anxiety     <- smd(
  group_1 = study_data$test_anxiety[study_data$group == "Intervention"],
  group_2 = study_data$test_anxiety[study_data$group == "Control"]
)$value

# CIs.
ci_m <- ci_smd(smd = d_motivation, n_1 = 60, n_2 = 60)
ci_s <- ci_smd(smd = d_self_eff,   n_1 = 60, n_2 = 60)
ci_a <- ci_smd(smd = d_anxiety,    n_1 = 60, n_2 = 60)

plot_ci(
  estimate = c(d_motivation, d_self_eff, d_anxiety),
  lower    = c(ci_m$value[1], ci_s$value[1], ci_a$value[1]),
  upper    = c(ci_m$value[3], ci_s$value[3], ci_a$value[3]),
  names    = c("Motivation", "Self-Efficacy", "Test Anxiety"),
  n        = 120,
  reference_line = 0,
  xlab  = "Standardized Mean Difference (d)",
  title = "Intervention Effects Across Outcomes"
)
```

![Cohen's \\d\\ with 95\\ confidence intervals for the intervention
effect on three distinct outcomes, displayed as a combined forest plot.
The \\N\\ annotations remind the reader that the precision of every
comparison depends on its sample
size.](effect-size-visualization_files/figure-html/combined-forest-1.png)

Cohen’s $`d`$ with 95% confidence intervals for the intervention effect
on three distinct outcomes, displayed as a combined forest plot. The
$`N`$ annotations remind the reader that the precision of every
comparison depends on its sample size.

This combined display makes it easy to see which outcomes show the
strongest intervention effects and where the confidence intervals
include zero (suggesting a non-significant difference at the chosen
$`\alpha`$ level).

## Turning Off Annotations

Every DMAR plot function shows the confidence interval and sample size
by default, because these are essential for transparent scientific
reporting. However, for presentations or simplified displays, both can
be suppressed:

``` r

plot_smd(smd = 0.50, show_ci = FALSE, show_n = FALSE,
         title = "Minimal Display: d = 0.50")
```

![Minimal plot suitable for a slide or teaching figure, with confidence
interval and sample size suppressed via \`show_ci = FALSE\` and \`show_n
=
FALSE\`.](effect-size-visualization_files/figure-html/plot-smd-minimal-1.png)

Minimal plot suitable for a slide or teaching figure, with confidence
interval and sample size suppressed via `show_ci = FALSE` and
`show_n = FALSE`.

In a manuscript or report, the recommendation is to leave both
annotations on. In a slide deck where the same information appears in
the surrounding text, suppression is reasonable.

## Customizing the Plots

Every DMAR plot function returns a `ggplot2` object, which means the
output is a fully editable plot rather than a fixed image. Any layer,
scale, theme, or annotation supported by `ggplot2` can be added on top.

``` r

plot_smd(smd = 0.65, n_1 = 80, n_2 = 80,
         group_labels = c("Treatment", "Control"),
         title = "Reading Intervention") +
  ggplot2::theme(legend.position = "top",
                 plot.title = ggplot2::element_text(size = 14,
                                                    face = "bold")) +
  ggplot2::scale_fill_manual(
    values = c("Treatment" = "#1B7837", "Control" = "#762A83"), name = NULL
  )
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.
```

![Every DMAR plot returns a \`ggplot2\` object, so additional
\`ggplot2\` layers can be added with \`+\`. Here the legend is moved to
the top, the title is bold, and the group fills follow a
colorblind-friendly
palette.](effect-size-visualization_files/figure-html/customize-1.png)

Every DMAR plot returns a `ggplot2` object, so additional `ggplot2`
layers can be added with `+`. Here the legend is moved to the top, the
title is bold, and the group fills follow a colorblind-friendly palette.

The same idiom (add `+` followed by another `ggplot2` layer) works for
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md) and
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md) as
well. Because the underlying object is just a `ggplot`, it can be saved
with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
embedded in an R Markdown document, or further composed into multi-panel
figures with packages like `patchwork` or `cowplot`.

## A Note on MBESS

DMAR is a more modern, more general, and greatly expanded reimagining of
the MBESS package (Kelley, 2007a, 2007b), which remains available on
CRAN. MBESS provides the same core capabilities (confidence intervals
for standardized effect sizes, sample size planning, and measurement)
implemented in a conventional S3 style. DMAR reimagines that work with:

- **Tidy output**: all functions return data frames with `term` and
  `value` columns (or similarly structured output), making results easy
  to pipe into downstream analysis and visualization.
- **Integrated visualization**:
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  and
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)
  produce publication-quality `ggplot2` graphics with confidence
  intervals and sample sizes shown by default.
- **Modern R conventions**: `TRUE`/`FALSE` (never `T`/`F`),
  [`warning()`](https://rdrr.io/r/base/warning.html) instead of
  [`print()`](https://rdrr.io/r/base/print.html) for diagnostics, and
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) for
  optional dependencies.

For existing scripts and reproducibility, MBESS continues to work
exactly as it always has. For new projects and teaching, DMAR offers a
cleaner, more consistent interface.

## References

Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit, R. A.
(2019). Raincloud plots: A multi-platform tool for robust data
visualization. *Wellcome Open Research, 4*, 63.

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Lawrence Erlbaum.

Cumming, G. (2012). *Understanding the new statistics: Effect sizes,
confidence intervals, and meta-analysis*. Routledge.

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press.

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*(2), 137–152.

Lee, Y.-S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society. Series B, 33*(1), 117–130.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.

Patil, I. (2023). *ggrain: A ‘ggplot2’ extension for raincloud plots*. R
package. <https://CRAN.R-project.org/package=ggrain>

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
