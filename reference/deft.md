# Kish's Design Effect (DEFF), DEFT, and the Effective Sample Size

Computes the design effect (DEFF) and its square root (DEFT) for a
clustered or multistage sample, given a vector of per-cluster sample
sizes and a value of the intraclass correlation. Returns Kish's (1965)
classic formula together with the effective sample size and a
description of the clustering, including the number of empty clusters
(no observations) and the number of singleton clusters (one
observation), which carry different amounts of within-cluster
information in a mixed-effects context. The function is also available
under the spelled-out alias `design_effect()`.

## Usage

``` r
deft(cluster_sizes, icc)

design_effect(cluster_sizes, icc)
```

## Arguments

- cluster_sizes:

  Numeric vector of per-cluster sample sizes. Each element is the number
  of observations in one cluster (so the length of the vector is the
  number of clusters). Zero values are allowed and counted as empty
  clusters; they do not affect the computation of DEFF.

- icc:

  Intraclass correlation coefficient, in \\\[0, 1)\\. Use the population
  value when planning prospectively, or the sample estimate (e.g., from
  [`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) or
  [`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md))
  when describing an observed clustered sample.

## Value

A `data.frame` with rows for the design effect, its square root, the
effective sample size, the total observation and cluster counts
(including separate counts of empty and singleton clusters), the mean
cluster size, Kish's design-weighted mean cluster size, and the input
`icc`.

## Details

**Definition (Kish, 1965).** For a clustered sample with per- cluster
sizes \\m_1, m_2, \ldots, m_K\\ and intraclass correlation \\\rho\\, the
design effect on the variance of the mean is \$\$\mathrm{DEFF} \\=\\ 1 +
(\bar m^{\*} - 1) \rho,\$\$ where \\\bar m^{\*} = \sum_k m_k^2 / \sum_k
m_k\\ is the design- weighted average cluster size (Kish, 1965, eq. 5.4;
sometimes called the "Kish weighted average" or "effective cluster
size"). For equal cluster sizes \\m_k = m\\, this reduces to the
classroom form \\1 + (m - 1)\rho\\. The DEFT is the square root of DEFF
and is the inflation factor on the *standard error* of the mean (whereas
DEFF inflates the variance).

**Effective sample size.** The number of observations from a simple
random sample that would yield the same standard error as the clustered
sample is \$\$N\_{\mathrm{eff}} \\=\\ N / \mathrm{DEFF} \\=\\ N /
\mathrm{DEFT}^2,\$\$ where \\N = \sum_k m_k\\ is the total observations.

**Why empty and singleton clusters are reported separately.** In
mixed-effects / multilevel modeling, clusters with zero observations
carry no information (they should be dropped before fitting), and
clusters with one observation contribute to \\N\\ and to the
fixed-effect estimate but contribute nothing to the estimation of the
random-effect variance or to the within-cluster residual. Hox et al.
(2017) note that singleton-heavy designs have a design effect close to 1
even at moderate \\\rho\\ because the weighted cluster size is small.
The output reports the counts of empty and singleton clusters so the
user can see at a glance how much of the nominal sample size carries
clustering information.

If a fitted `lmerMod` object is available, the typical workflow is to
extract `icc` via
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md) and
the per-cluster sample sizes via `table(cluster_id)`, then pass both to
`deft()`; see the second example.

## References

Hox, J. J., Moerbeek, M., & van de Schoot, R. (2017). *Multilevel
analysis: Techniques and applications* (3rd ed.). Routledge.

Kish, L. (1965). *Survey sampling*. Wiley.

Kish, L. (1992). Weighting for unequal Pi. *Journal of Official
Statistics, 8*(2), 183–200.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 15 and 16 on mixed-effects models and
nested designs.)

Snijders, T. A. B., & Bosker, R. J. (2012). *Multilevel analysis: An
introduction to basic and advanced multilevel modeling* (2nd ed.). Sage.

## See also

[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)

Other design utilities:
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
[`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
[`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced design: K = 30 clusters of size 20, ICC = 0.10.
deft(cluster_sizes = rep(20, 30), icc = 0.10)
#>  term                   value
#>  design_effect          2.9  
#>  deft                   1.7  
#>  effective_n            207  
#>  n_total                600  
#>  n_clusters_total       30   
#>  n_clusters_with_data   30   
#>  n_clusters_empty       0    
#>  n_clusters_singletons  0    
#>  n_clusters_informative 30   
#>  m_bar                  20   
#>  m_kish                 20   
#>  icc                    0.1  
# DEFF = 1 + (20 - 1) * 0.10 = 2.9; DEFT = 1.7; effective N = 600 / 2.9 = 207.

# 2. Unbalanced design with some empty and some singleton clusters.
# Suppose K = 25 schools with attendance ranging from 0 (closed)
# through 1 (single student showed up) to 30 (full class).
set.seed(113)
sizes <- c(0, 0, 1, 1, 1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 22,
           22, 25, 26, 28, 28, 30, 30, 30, 30, 30)
deft(cluster_sizes = sizes, icc = 0.10)
#>  term                   value
#>  design_effect          3.33 
#>  deft                   1.82 
#>  effective_n            119  
#>  n_total                397  
#>  n_clusters_total       25   
#>  n_clusters_with_data   23   
#>  n_clusters_empty       2    
#>  n_clusters_singletons  3    
#>  n_clusters_informative 20   
#>  m_bar                  17.3 
#>  m_kish                 24.3 
#>  icc                    0.1  

# 3. From a cluster-id vector: tabulate, then call deft().
cluster_id <- rep(1:8, times = c(20, 15, 22, 1, 0, 30, 18, 25))
deft(cluster_sizes = as.numeric(table(factor(cluster_id, levels = 1:8))),
     icc = 0.15)
#>  term                   value
#>  design_effect          4.24 
#>  deft                   2.06 
#>  effective_n            30.9 
#>  n_total                131  
#>  n_clusters_total       8    
#>  n_clusters_with_data   7    
#>  n_clusters_empty       1    
#>  n_clusters_singletons  1    
#>  n_clusters_informative 6    
#>  m_bar                  18.7 
#>  m_kish                 22.6 
#>  icc                    0.15 
```
