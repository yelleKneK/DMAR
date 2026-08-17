#' @rdname ss_aipe_crd
#' @name ss_aipe_crd
#' @aliases ss_aipe_crd_n_clusters_fixed_width
#' @aliases ss_aipe_crd_n_individuals_fixed_width
#' @aliases ss_aipe_crd_n_clusters_fixed_budget
#' @aliases ss_aipe_crd_n_individuals_fixed_budget
#' @aliases ss_aipe_crd_both_fixed_budget
#' @aliases ss_aipe_crd_both_fixed_width
#'
#' @title Find Target Sample Sizes for the Accuracy in Unstandardized Conditions Means Estimation in CRD
#'
#' @description
#' Find target sample sizes (the number of clusters, cluster size, or both) for the accuracy in unstandardized conditions
#' means estimation in CRD. If users wish to seek for both types of sample sizes simultaneously, an additional constraint
#' is required, such as a desired width or a desired budget.
#'
#' @param width The desired width of the confidence interval of the unstandardized means difference
#' @param budget The desired amount of budget
#' @param n_clusters The desired number of clusters
#' @param n_individuals The number of individuals in each cluster (cluster size)
#' @param pr_treat The proportion of treatment clusters
#' @param clus_cost The cost of collecting a new cluster regardless of the number of individuals collected in each cluster
#' @param indiv_cost The cost of collecting a new individual
#' @param tau_Y The residual variance in the between level before accounting for the covariate
#' @param sigma2_Y The residual variance in the within level before accounting for the covariate
#' @param total_var The total residual variance before accounting for the covariate
#' @param icc_Y The intraclass correlation of the dependent variable
#' @param R2_within The proportion of variance explained in the within level (used when \code{covariate = TRUE})
#' @param R2_between The proportion of variance explained in the between level (used when \code{covariate = TRUE})
#' @param num_predictors The number of predictors used in the between level
#' @param assurance The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)
#' @param conf_level The desired level of confidence for the confidence interval
#' @param diff_size Difference cluster size specification. The differences in cluster sizes can be specified in two ways, and the specified vector is recycled across the clusters. First, users may specify differences as integers, which can be negative or positive; the resulting cluster sizes add the specified values to the estimated cluster size. For example, if the cluster size is 25, the number of clusters is 10, and \code{diff_size = c(-1, 0, 1)}, the cluster sizes will be 24, 25, 26, 24, 25, 26, 24, 25, 26, and 24. Second, users may specify multipliers of the cluster size as positive decimals; at least one value must be non-integer, which is what selects the multiplicative form. The resulting cluster sizes multiply the estimated cluster size by the specified values and round to the nearest integer. For example, if the cluster size is 25, the number of clusters is 10, and \code{diff_size = c(0.8, 1, 1.2)}, the cluster sizes will be 20, 25, 30, 20, 25, 30, 20, 25, 30, and 20. In either form a resulting cluster size below 1 is set to 1. If \code{NULL}, the cluster size is equal across clusters
#'
#' @details
#' Here are the functions' descriptions:
#' \describe{
#'   \item{\code{ss_aipe_crd_n_clusters_fixed_width}}{Find the number of clusters given a specified width of the confidence interval and the cluster size}
#'   \item{\code{ss_aipe_crd_n_individuals_fixed_width}}{Find the cluster size given a specified width of the confidence interval and the number of clusters}
#'   \item{\code{ss_aipe_crd_n_clusters_fixed_budget}}{Find the number of clusters given a budget and the cluster size}
#'   \item{\code{ss_aipe_crd_n_individuals_fixed_budget}}{Find the cluster size given a budget and the number of clusters}
#'   \item{\code{ss_aipe_crd_both_fixed_budget}}{Find the sample size combinations (the number of clusters and that cluster size) providing the narrowest confidence interval given the fixed budget}
#'   \item{\code{ss_aipe_crd_both_fixed_width}}{Find the sample size combinations (the number of clusters and that cluster size) providing the lowest cost given the specified width of the confidence interval}
#' }
#'
#' @return
#' The \code{ss_aipe_crd_n_clusters_fixed_width} and \code{ss_aipe_crd_n_clusters_fixed_budget} functions provide the number of clusters.
#' The \code{ss_aipe_crd_n_individuals_fixed_width} and \code{ss_aipe_crd_n_individuals_fixed_budget} functions provide the cluster size.
#' The \code{ss_aipe_crd_both_fixed_budget} and \code{ss_aipe_crd_both_fixed_width} provide the number of clusters and the
#' cluster size, respectively.
#'
#' @references
#' Pornprasertmanit, S., & Schneider, W. J. (2014). Accuracy in parameter
#'   estimation in cluster randomized designs. \emph{Psychological
#'   Methods, 19}(3), 356--379. \doi{10.1037/a0037036}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Examples for each function
#' ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 30,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#'
#' ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#'
#' ss_aipe_crd_n_clusters_fixed_budget(budget = 10000, n_individuals = 20,
#'   clus_cost = 20, indiv_cost = 1)
#'
#' ss_aipe_crd_n_individuals_fixed_budget(budget = 10000, n_clusters = 30,
#'   clus_cost = 20, indiv_cost = 1,
#'   pr_treat = 0.5, tau_Y = 0.05, sigma2_Y = 0.95, assurance = 0.8)
#'
#' ss_aipe_crd_both_fixed_budget(budget = 10000, clus_cost = 30, indiv_cost = 1,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#'
#' ss_aipe_crd_both_fixed_width(width = 0.3, clus_cost = 0, indiv_cost = 1,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#'
#' # Examples for different cluster size
#' set.seed(113)
#' ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 30,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75,
#'   diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))
#'
#' # Examples for different number of clusters
#' ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
#'   pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75,
#'   diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export
ss_aipe_crd_n_clusters_fixed_width <- function(width, n_individuals, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, clus_cost = NULL, indiv_cost = NULL, diff_size = NULL) {
    n_clusters <- .find_n_clus_crd_diff(width = width, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    calculated_width <- .find_width_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    calculated_cost <- NULL
    if (!is.null(clus_cost) && !is.null(indiv_cost)) calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost, indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "cluster_size"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd
#' @export
ss_aipe_crd_n_individuals_fixed_width <- function(width, n_clusters, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, clus_cost = NULL, indiv_cost = NULL, diff_size = NULL) {
    n_individuals <- .find_n_indiv_crd_diff(width = width, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (n_individuals == "> 100000") stop("With the current number of clusters, it is impossible to achieve the target width. Please increase the number of clusters.")
    calculated_width <- .find_width_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    calculated_cost <- NULL
    if (!is.null(clus_cost) && !is.null(indiv_cost)) calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost, indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "necessary_n_clusters"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd
#' @export
ss_aipe_crd_n_clusters_fixed_budget <- function(budget, n_individuals, clus_cost = 0, indiv_cost = 1, pr_treat = NULL, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- .find_n_clus_crd_budget(budget = budget, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    calculated_width <- NULL
    if (!is.null(pr_treat)) {
        calculated_width <- .find_width_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    }
    calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "cluster_size"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd
#' @export
ss_aipe_crd_n_individuals_fixed_budget <- function(budget, n_clusters, clus_cost = 0, indiv_cost = 1, pr_treat = NULL, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_individuals <- .find_n_indiv_crd_budget(budget = budget, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    calculated_width <- NULL
    if (!is.null(pr_treat)) {
        calculated_width <- .find_width_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    }
    calculated_cost <- .cost_crd(n_clusters, n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    re <- .report_crd(n_clusters, n_individuals, calculated_width, cost = calculated_cost, es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size)
    .as_dmar_tbl(re[which(!re$term == "necessary_n_clusters"), ], conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd
#' @export
ss_aipe_crd_both_fixed_budget <- function(budget, clus_cost = 0, indiv_cost = 1, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    result <- .find_min_width_crd_diff(budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    calculated_cost <- .cost_crd(result[1], result[2], clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    .as_dmar_tbl(.report_crd(result[1], result[2], result[3], cost = calculated_cost, es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size), conf_level = conf_level, subclass = "dmar_ss_aipe")
}

#' @rdname ss_aipe_crd
#' @export
ss_aipe_crd_both_fixed_width <- function(width, clus_cost = 0, indiv_cost = 1, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    result <- .find_min_cost_crd_diff(width = width, clus_cost = clus_cost, indiv_cost = indiv_cost, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    calculated_width <- .find_width_crd_diff(n_clusters = result[1], n_individuals = result[2], pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    .as_dmar_tbl(.report_crd(result[1], result[2], calculated_width, cost = result[3], es = FALSE, es_type = 0, assurance = assurance, diff_size = diff_size), conf_level = conf_level, subclass = "dmar_ss_aipe")
}
