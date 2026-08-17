#' @rdname ss_aipe_crd_es
#' @name ss_aipe_crd_es
#' @aliases ss_aipe_crd_es_n_clusters_fixed_width
#' @aliases ss_aipe_crd_es_n_individuals_fixed_width
#' @aliases ss_aipe_crd_es_n_clusters_fixed_budget
#' @aliases ss_aipe_crd_es_n_individuals_fixed_budget
#' @aliases ss_aipe_crd_es_both_fixed_budget
#' @aliases ss_aipe_crd_es_both_fixed_width
#'
#' @title Find Target Sample Sizes for the Accuracy in Standardized Conditions Means Estimation in CRD
#'
#' @description
#' Find target sample sizes (the number of clusters, cluster size, or both) for the accuracy in standardized conditions
#' means estimation in CRD. If users wish to seek for both types of sample sizes simultaneously, an additional
#' constraint is required, such as a desired width or a desired budget. This function uses the likelihood-based
#' confidence interval (Cheung, 2009) by the \code{OpenMx} package (Boker et al., 2011). See further details at
#' Pornprasertmanit and Schneider (2014).
#'
#' @param width The desired width of the confidence interval of the unstandardized means difference
#' @param budget The desired amount of budget
#' @param n_clusters The desired number of clusters
#' @param n_individuals The number of individuals in each cluster (cluster size)
#' @param pr_treat The proportion of treatment clusters
#' @param clus_cost The cost of collecting a new cluster regardless of the number of individuals collected in each cluster
#' @param indiv_cost The cost of collecting a new individual
#' @param icc_Y The intraclass correlation of the dependent variable
#' @param icc_Z The intraclass correlation of the covariate (used when \code{covariate = TRUE}). If \code{icc_Z = 0}, the within-level covariate will be only used. If \code{icc_Z = 1}, the between-level covariate will be only used
#' @param es The amount of effect size
#' @param es_type The type of effect size. There are only three possible options: 0 = the effect size using total standard deviation, 1 = the effect size using the individual-level standard deviation (level 1), 2 = the effect size using the cluster-level standard deviation (level 2)
#' @param R2_within The proportion of variance explained in the within level (used when \code{covariate = TRUE})
#' @param R2_between The proportion of variance explained in the between level (used when \code{covariate = TRUE})
#' @param num_predictors The number of predictors used in the between level
#' @param assurance The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)
#' @param conf_level The desired level of confidence for the confidence interval
#' @param diff_size Difference cluster size specification. The differences in cluster sizes can be specified in two ways, and the specified vector is recycled across the clusters. First, users may specify differences as integers, which can be negative or positive; the resulting cluster sizes add the specified values to the estimated cluster size. For example, if the cluster size is 25, the number of clusters is 10, and \code{diff_size = c(-1, 0, 1)}, the cluster sizes will be 24, 25, 26, 24, 25, 26, 24, 25, 26, and 24. Second, users may specify multipliers of the cluster size as positive decimals; at least one value must be non-integer, which is what selects the multiplicative form. The resulting cluster sizes multiply the estimated cluster size by the specified values and round to the nearest integer. For example, if the cluster size is 25, the number of clusters is 10, and \code{diff_size = c(0.8, 1, 1.2)}, the cluster sizes will be 20, 25, 30, 20, 25, 30, 20, 25, 30, and 20. In either form a resulting cluster size below 1 is set to 1. If \code{NULL}, the cluster size is equal across clusters
#' @param nrep The number of replications used in a priori Monte Carlo simulation
#' @param seed An optional integer seed for the a priori Monte Carlo simulation. The default \code{NULL} uses the current state of the random number generator and leaves it unchanged, so repeated calls reflect the genuine sampling variability of the simulation. Supply an integer for reproducible results, in which case the generator state is restored on exit
#' @param multicore Use multiple processors within a computer. Specify as \code{TRUE} to use it
#' @param num_proc The number of processors to be used when \code{multicore = TRUE}. If it is not specified, the package will use the maximum number of processors in a machine
#'
#' @details
#' Here are the functions' descriptions:
#' \describe{
#'   \item{\code{ss_aipe_crd_es_n_clusters_fixed_width}}{Find the number of clusters given a specified width of the confidence interval and the cluster size}
#'   \item{\code{ss_aipe_crd_es_n_individuals_fixed_width}}{Find the cluster size given a specified width of the confidence interval and the number of clusters}
#'   \item{\code{ss_aipe_crd_es_n_clusters_fixed_budget}}{Find the number of clusters given a budget and the cluster size}
#'   \item{\code{ss_aipe_crd_es_n_individuals_fixed_budget}}{Find the cluster size given a budget and the number of clusters}
#'   \item{\code{ss_aipe_crd_es_both_fixed_budget}}{Find the sample size combinations (the number of clusters and that cluster size) providing the narrowest confidence interval given the fixed budget}
#'   \item{\code{ss_aipe_crd_es_both_fixed_width}}{Find the sample size combinations (the number of clusters and that cluster size) providing the lowest cost given the specified width of the confidence interval}
#' }
#'
#' @return
#' The \code{ss_aipe_crd_es_n_clusters_fixed_width} and \code{ss_aipe_crd_es_n_clusters_fixed_budget} functions provide the number of clusters.
#' The \code{ss_aipe_crd_es_n_individuals_fixed_width} and \code{ss_aipe_crd_es_n_individuals_fixed_budget} functions provide the cluster size.
#' The \code{ss_aipe_crd_es_both_fixed_budget} and \code{ss_aipe_crd_es_both_fixed_width} provide the number of clusters and the
#' cluster size, respectively.
#'
#' @references
#' Boker, S. M., Neale, M. C., Maes, H. H., Wilde, M., Spiegel, M.,
#'   Brick, T. R., ... Fox, J. (2011). OpenMx: An open source extended
#'   structural equation modeling framework. \emph{Psychometrika, 76}(2),
#'   306--317. \doi{10.1007/s11336-010-9200-6}
#'
#' Cheung, M. W.-L. (2009). Constructing approximate confidence intervals
#'   for parameters with structural equation models. \emph{Structural
#'   Equation Modeling, 16}(2), 267--294. \doi{10.1080/10705510902751291}
#'
#' Pornprasertmanit, S., & Schneider, W. J. (2010). \emph{Efficient sample size for power and desired accuracy in Cohen's d
#' estimation in two-group cluster randomized design} (Master Thesis). Illinois State University, Normal, IL.
#'
#' Pornprasertmanit, S., & Schneider, W. J. (2014). Accuracy in parameter
#'   estimation in cluster randomized designs. \emph{Psychological
#'   Methods, 19}(3), 356--379. \doi{10.1037/a0037036}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Two of these planners answer a question the budget alone settles. Given
#' # what it costs to open a cluster and what it costs to collect one more
#' # individual, the first reports how many clusters a budget buys at a fixed
#' # cluster size and the second reports how large each cluster can be at a
#' # fixed number of clusters. Clusters cost nothing to open here and each
#' # individual costs 1, so the budget buys 1000 individuals and the only
#' # question is how to arrange them.
#' ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, n_individuals = 20,
#'   clus_cost = 0, indiv_cost = 1)
#'
#' ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 200,
#'   clus_cost = 0, indiv_cost = 1)
#'
#' # The interval width these planners can report, and every answer that
#' # targets a width, rests on an a priori Monte Carlo simulation: a candidate
#' # is evaluated by generating nrep data sets and reading the
#' # likelihood-based confidence interval on the standardized effect size from
#' # OpenMx, and the planners that search over the number of clusters evaluate
#' # many candidates in turn. Those calls run for seconds to minutes apiece,
#' # so they are shown below but not run; the package's tests exercise them.
#' # Each one describes a population standardized effect size of 0.5, with
#' # es_type = 1 putting that effect size in individual-level standard
#' # deviation units and a quarter of the outcome variance lying between
#' # clusters.
#' #
#' # Supplying nrep and the population values to a budget planner adds the
#' # expected width of the interval the affordable design buys:
#' #   ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, n_individuals = 20,
#' #     clus_cost = 0, indiv_cost = 1, es = 0.5, es_type = 1, icc_Y = 0.25,
#' #     pr_treat = 0.5, nrep = 1000, seed = 113)
#' #
#' # Cluster size needed for a target width, given the number of clusters.
#' # With 250 clusters the planner settles on the smallest cluster size it
#' # will consider, two individuals per cluster, and the expected width still
#' # comes in well under the target: for a contrast between conditions that
#' # are assigned at the cluster level, precision is bought with clusters
#' # rather than with what happens inside them.
#' #   ss_aipe_crd_es_n_individuals_fixed_width(width = 0.5, n_clusters = 250,
#' #     es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 1000,
#' #     seed = 113)
#' #
#' # Once recruiting a cluster costs 5, the number of clusters and the cluster
#' # size trade off against each other, and this planner searches the
#' # combinations the budget allows for the narrowest expected interval:
#' #   ss_aipe_crd_es_both_fixed_budget(budget = 1000, clus_cost = 5,
#' #     indiv_cost = 1, es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5,
#' #     nrep = 1000, seed = 113)
#' #
#' # Number of clusters needed for a target width, given the cluster size:
#' #   ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20,
#' #     es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 1000,
#' #     seed = 113)
#' #
#' # Both quantities under a target width, taking the least costly combination
#' # that reaches it:
#' #   ss_aipe_crd_es_both_fixed_width(width = 0.5, clus_cost = 5,
#' #     indiv_cost = 1, es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5,
#' #     nrep = 1000, seed = 113)
#' #
#' # Unequal cluster sizes: diff_size gives each cluster's deviation from
#' # n_individuals (additive) or its multiplicative factor.
#' #   ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20,
#' #     es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 1000,
#' #     seed = 113, diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))
#' #
#' #   ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20,
#' #     es = 0.5, es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 1000,
#' #     seed = 113, diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export
ss_aipe_crd_es_n_clusters_fixed_width <- function(width, n_individuals, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, clus_cost = NULL, indiv_cost = NULL, diff_size = NULL) {
    suppressWarnings(result <- .find_n_clus_crd_es(width = width, n_individuals = n_individuals, es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size))
    calculated_cost <- NULL
    if (!is.null(clus_cost) && !is.null(indiv_cost)) calculated_cost <- .cost_crd(result[1], n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(result[1], n_individuals, result[2], cost = calculated_cost, es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "cluster_size"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd_es
#' @export
ss_aipe_crd_es_n_individuals_fixed_width <- function(width, n_clusters, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, clus_cost = NULL, indiv_cost = NULL, diff_size = NULL) {
    suppressWarnings(result <- .find_n_indiv_crd_es(width = width, n_clusters = n_clusters, es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size))
    calculated_cost <- NULL
    if (!is.null(clus_cost) && !is.null(indiv_cost)) calculated_cost <- .cost_crd(n_clusters, result[1], clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, result[1], result[2], cost = calculated_cost, es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "necessary_n_clusters"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd_es
#' @export
ss_aipe_crd_es_n_clusters_fixed_budget <- function(budget, n_individuals, clus_cost, indiv_cost, nrep = NULL, pr_treat = NULL, icc_Y = NULL, es = NULL, es_type = 1, num_predictors = 0, icc_Z = NULL, R2_within = NULL, R2_between = NULL, assurance = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- .find_n_clus_crd_budget(budget = budget, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    calculated_width <- NULL
    if (!is.null(nrep) && !is.null(pr_treat) && !is.null(n_individuals) && !is.null(icc_Y)) {
        suppressWarnings(calculated_width <- .find_width_crd_es(nrep = nrep, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = 1, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, assurance = assurance, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size))
    }
    calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "cluster_size"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd_es
#' @export
ss_aipe_crd_es_n_individuals_fixed_budget <- function(budget, n_clusters, clus_cost, indiv_cost, nrep = NULL, pr_treat = NULL, icc_Y = NULL, es = NULL, es_type = 1, num_predictors = 0, icc_Z = NULL, R2_within = NULL, R2_between = NULL, assurance = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, conf_level = 0.95, diff_size = NULL) {
    n_individuals <- .find_n_indiv_crd_budget(budget = budget, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    calculated_width <- NULL
    if (!is.null(nrep) && !is.null(pr_treat) && !is.null(n_individuals) && !is.null(icc_Y)) {
        suppressWarnings(calculated_width <- .find_width_crd_es(nrep = nrep, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = 1, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, assurance = assurance, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size))
    }
    calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "necessary_n_clusters"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd_es
#' @export
ss_aipe_crd_es_both_fixed_budget <- function(budget, clus_cost = 0, indiv_cost = 1, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    suppressWarnings(result <- .find_min_width_crd_es(budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size))
    calculated_cost <- .cost_crd(result[1], result[2], clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    .as_dmar_tbl(.report_crd(result[1], result[2], result[3], cost = calculated_cost, es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size), conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd_es
#' @export
ss_aipe_crd_es_both_fixed_width <- function(width, clus_cost = 0, indiv_cost = 1, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    suppressWarnings(result <- .find_min_cost_crd_es(width = width, clus_cost = clus_cost, indiv_cost = indiv_cost, es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size))
    .as_dmar_tbl(.report_crd(result[1], result[2], result[4], cost = result[3], es = TRUE, es_type = es_type, assurance = assurance, diff_size = diff_size), conf_level = conf_level, subclass = "dmar_ss_aipe")
}
