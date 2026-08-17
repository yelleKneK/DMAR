# Internal helpers shared by R/ss_aipe_crd.R and R/ss_aipe_crd_es.R.
# Not exported. The "*_diff" helpers handle the unstandardized
# conditions difference; the "*_es" helpers handle the standardized
# effect size variants and use a priori Monte Carlo simulation via
# OpenMx for likelihood-based confidence intervals.

# Get a nice report for a public function
.report_crd <- function(n_clusters, n_individuals, width = NULL, cost = NULL, es = FALSE, es_type = 0, assurance = NULL, diff_size = NULL) {
    term <- "necessary_n_clusters"
    value <- n_clusters
    if (is.null(diff_size)) {
        term <- c(term, "cluster_size")
        value <- c(value, n_individuals)
    } else {
        out <- data.frame(table(.find_n_indiv_vec(as.numeric(n_individuals), diff_size, n_clusters)))
        colnames(out) <- c("cluster_size", "freq")
        print(out)
    }
    if (!is.null(width)) {
        if (es) {
            es_label <- NULL
            if (es_type == 0) {
                es_label <- "total"
            } else if (es_type == 1) {
                es_label <- "individual-level"
            } else {
                es_label <- "cluster-level"
            }
            if (is.null(assurance)) {
                term <- c(term, paste("exp_width_of_", es_label, "_effect_size", sep = ""))
                value <- c(value, width)
            } else {
                term <- c(term, paste("width_of_", es_label, "_effect_size_with_", round(assurance, 2), "_assurance", sep = ""))
                value <- c(value, width)
            }
        } else {
            if (is.null(assurance)) {
                term <- c(term, paste("exp_width_of_unstd_conditions_diff", sep = ""))
                value <- c(value, width)
            } else {
                term <- c(term, paste("width_of_unstd_conditions_diff_with_", round(assurance, 2), "_assurance", sep = ""))
                value <- c(value, width)
            }
        }
    }
    if (!is.null(cost)) {
        term <- c(term, "budget")
        value <- c(value, cost)
    }
    return(data.frame(term, value))
}

.find_n_indiv_vec <- function(n_individuals, diff_size, n_clusters) {
    is_integer <- all(round(diff_size) == diff_size)
    result <- NULL
    if (is_integer) {
        result <- rep(n_individuals + diff_size, length.out = n_clusters)
    } else {
        result <- rep(round(n_individuals * diff_size), length.out = n_clusters)
    }
    result[result < 1] <- 1
    result
}

# Find the number of clusters given the specified width of ES and the cluster size
.find_n_clus_crd_es <- function(width, n_individuals, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    if (num_predictors > 0 && is.null(icc_Z)) icc_Z <- icc_Y
    if (num_predictors > 1) stop("Only one predictor is allowed.")
    total_var <- 1
    if (es_type == 0) {
        total_var <- 1
    } else if (es_type == 1) {
        total_var <- 1 / (1 - icc_Y)
    } else if (es_type == 2) {
        total_var <- 1 / icc_Y
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    startval <- .find_n_clus_crd_diff(width = width, n_individuals = n_individuals, pr_treat = pr_treat, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    startval <- as.numeric(startval)
    startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = startval, n_treat_clus = round(startval * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
    if (startwidth < width) {
        repeat {
            startval <- startval - 1
            if (round(startval * pr_treat) == 1 || (startval - round(startval * pr_treat)) == 1) return(c(startval + 1, startwidth))
            savedwidth <- startwidth
            startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = startval, n_treat_clus = round(startval * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
            if (startwidth > width) return(c(startval + 1, savedwidth))
        }
    } else if (startwidth > width) {
        repeat {
            startval <- startval + 1
            startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = startval, n_treat_clus = round(startval * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
            if (startwidth < width) return(c(startval, startwidth))
        }
    } else {
        return(c(startval, startwidth))
    }
}

# Simulate multiple datasets and find the average of the width of CI of ES (with or without the degree of assurance)
.find_width_crd_es <- function(nrep, n_clusters, n_treat_clus, n_individuals, icc_Y, es, es_type = 1, total_var = 1, covariate = FALSE, icc_Z = NULL, R2_within = NULL, R2_between = NULL, total_var_Z = 1, assurance = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, conf_level = 0.95, diff_size = NULL) {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        .old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
      } else {
        on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }

    seed_list <- as.list(sample(1:999999, nrep))
    result_l <- NULL

    if (multicore) {
        if (!requireNamespace("parallel", quietly = TRUE)) stop("The package 'parallel' is needed; please install the package and try again.")
        sys <- .Platform$OS.type
        if (is.null(num_proc))
            num_proc <- parallel::detectCores()
        if (sys == "windows") {
            cl <- parallel::makeCluster(num_proc, type = "PSOCK")
            result_l <- parallel::clusterApplyLB(cl, seed_list, .runrep_width_es_crd, n_clusters = n_clusters, n_treat_clus = n_treat_clus, n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, conf_level = conf_level, diff_size = diff_size)
            parallel::stopCluster(cl)
        } else {
            result_l <- parallel::mclapply(seed_list, .runrep_width_es_crd, n_clusters = n_clusters, n_treat_clus = n_treat_clus, n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, conf_level = conf_level, diff_size = diff_size)
        }
    } else {
        result_l <- lapply(seed_list, .runrep_width_es_crd, n_clusters = n_clusters, n_treat_clus = n_treat_clus, n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, conf_level = conf_level, diff_size = diff_size)
    }
    result <- do.call(c, result_l)
    if (is.null(assurance)) {
        return(mean(result, na.rm = TRUE))
    } else {
        return(quantile(result, assurance, na.rm = TRUE))
    }
}

# Create data and find the width of the likelihood-based CI of ES
.runrep_width_es_crd <- function(seed, n_clusters, n_treat_clus, n_individuals, icc_Y, es, es_type = 1, total_var = 1, covariate = FALSE, icc_Z = NULL, R2_within = NULL, R2_between = NULL, total_var_Z = 1, conf_level = 0.95, diff_size = NULL) {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        .old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
      } else {
        on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }
    datawide <- .create_data_crd_wide(n_clusters = n_clusters, n_treat_clus = n_treat_clus, n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, diff_size = diff_size)
    y_label <- NULL
    if (!is.null(diff_size)) {
        size <- as.numeric(names(datawide))
        y_label <- lapply(size, function(x) paste("y", 1:x, sep = ""))
    } else {
        y_label <- paste("y", 1:n_individuals, sep = "")
    }

    x_label <- "x"
    z_within_label <- NULL
    if (covariate && icc_Z != 1) {
        if (!is.null(diff_size)) {
            size <- as.numeric(names(datawide))
            z_within_label <- lapply(size, function(x) paste("zw", 1:x, sep = ""))
        } else {
            z_within_label <- paste("zw", 1:n_individuals, sep = "")
        }
    }
    z_between_label <- NULL
    if (covariate && icc_Z != 0) z_between_label <- "zb"
    if (!is.null(diff_size)) {
        screencapture <- utils::capture.output(
            result <- .lik_ci_es_crd_unequal(datawide = datawide, y_label = y_label, x_label = x_label, z_within_label = z_within_label, z_between_label = z_between_label, es_type = es_type, icc_Y = icc_Y, es = es, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, conf_level = conf_level))
    } else {
        screencapture <- utils::capture.output(
            result <- .lik_ci_es_crd(datawide = datawide, y_label = y_label, x_label = x_label, z_within_label = z_within_label, z_between_label = z_between_label, es_type = es_type, icc_Y = icc_Y, es = es, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, conf_level = conf_level)
        )
    }
    return(result[2] - result[1])
}

# Calculate a likelihood-based CI of ES
.lik_ci_es_crd <- function(datawide, y_label, x_label, z_within_label = NULL, z_between_label = NULL, es_type = 1, icc_Y = 0.25, es = 0.5, total_var = 1, covariate = FALSE, icc_Z = 0.25, R2_within = 0.5, R2_between = 0.5, total_var_Z = 1, conf_level = 0.95) {
    tau <- icc_Y * total_var
    sigma <- (1 - icc_Y) * total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    if (covariate) {
        if (icc_Z == 0 && R2_between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
        if (icc_Z == 1 && R2_within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
        tauz <- total_var_Z * icc_Z
        sigmaz <- total_var_Z * (1 - icc_Z)
        gammazb <- sqrt(R2_between * tau / tauz)
        gammazw <- sqrt(R2_within * sigma / sigmaz)
        tau <- (1 - R2_between) * tau
        sigma <- (1 - R2_within) * sigma
    }
    probx <- sum(datawide[, x_label]) / nrow(datawide)

    if (!requireNamespace("OpenMx", quietly = TRUE)) stop("The package 'OpenMx' is needed; please install the package and try again.")

    latent_label <- c("intcept", "slope")
    if (is.null(z_within_label)) latent_label <- "intcept"

    frowlab <- c(y_label, x_label, z_between_label)
    fcollab <- c(frowlab, latent_label)
    lenrow <- length(frowlab)
    lencol <- length(fcollab)

    Alab <- matrix(NA, lencol, lencol)
    Aval <- matrix(0, lencol, lencol)
    Afree <- matrix(FALSE, lencol, lencol)
    colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab

    Alab["intcept", x_label] <- "groupdiff"
    if (!is.null(z_between_label)) Alab["intcept", z_between_label] <- "zbeffect"
    if (!is.null(z_within_label)) Alab[y_label, "slope"] <- paste("data.", z_within_label, sep = "")

    Aval["intcept", x_label] <- gamma1
    if (!is.null(z_between_label)) Aval["intcept", z_between_label] <- gammazb
    Aval[y_label, latent_label] <- 1

    Afree["intcept", x_label] <- TRUE
    if (!is.null(z_between_label)) Afree["intcept", z_between_label] <- TRUE

    Slab <- matrix(NA, lencol, lencol)
    Sval <- matrix(0, lencol, lencol)
    Sfree <- matrix(FALSE, lencol, lencol)
    colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
    diag(Slab)[1:length(y_label)] <- "l1error"
    diag(Sval)[1:length(y_label)] <- sigma
    diag(Sfree)[1:length(y_label)] <- TRUE
    Slab["intcept", "intcept"] <- "l2error"
    Sval["intcept", "intcept"] <- tau
    Sfree["intcept", "intcept"] <- TRUE
    Slab[x_label, x_label] <- "varx"
    Sval[x_label, x_label] <- probx * (1 - probx)
    Sfree[x_label, x_label] <- TRUE
    if (!is.null(z_between_label)) {
        Slab[c(x_label, z_between_label), c(x_label, z_between_label)] <- "covxzb"
        Slab[x_label, x_label] <- "varx"
        Slab[z_between_label, z_between_label] <- "varzb"
        Sval[c(x_label, z_between_label), c(x_label, z_between_label)] <- 0
        Sval[x_label, x_label] <- probx * (1 - probx)
        Sval[z_between_label, z_between_label] <- tauz
        Sfree[c(x_label, z_between_label), c(x_label, z_between_label)] <- TRUE
    }

    Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latent_label)))
    Flab <- matrix(NA, lenrow, lencol)
    Ffree <- matrix(FALSE, lenrow, lencol)
    colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
    rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

    Mlab <- c(rep(NA, length(y_label)), "meanX")
    Mval <- c(rep(0, length(y_label)), probx)
    Mfree <- c(rep(FALSE, length(y_label)), TRUE)

    if (!is.null(z_between_label)) {
        Mlab <- c(Mlab, "meanzb")
        Mval <- c(Mval, 0)
        Mfree <- c(Mfree, TRUE)
    }
    Mlab <- c(Mlab, "meanctrl")
    Mval <- c(Mval, 0)
    Mfree <- c(Mfree, TRUE)
    if (!is.null(z_within_label)) {
        Mlab <- c(Mlab, "zweffect")
        Mval <- c(Mval, gammazw)
        Mfree <- c(Mfree, TRUE)
    }

    Mlab <- matrix(Mlab, 1, lencol)
    Mval <- matrix(Mval, 1, lencol)
    Mfree <- matrix(Mfree, 1, lencol)
    colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

    varzw <- 0
    if (!is.null(z_within_label)) varzw <- var(as.vector(as.matrix(datawide[, z_within_label])))

    # These lines prevent "no visible binding" notes when OpenMx evaluates the mxAlgebra expressions.
    groupdiff <- NULL
    l1error <- NULL
    l2error <- NULL
    zbeffect <- NULL
    varzb <- NULL
    zweffect <- NULL

    constraint <- NULL
    if (es_type == 0) {
        if (is.null(z_within_label)) {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        } else {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        }
    } else if (es_type == 1) {
        if (is.null(z_within_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw)), name = "es")
        }
    } else if (es_type == 2) {
        if (is.null(z_between_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
        }
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    onecov <- OpenMx::mxModel("effect size CRD", type = "RAM",
                              OpenMx::mxData(datawide, type = "raw"),
                              OpenMx::mxMatrix(type = "Full", nrow = lencol, ncol = lencol, values = Aval, free = Afree, labels = Alab, name = "A"),
                              OpenMx::mxMatrix(type = "Symm", nrow = lencol, ncol = lencol, values = Sval, free = Sfree, labels = Slab, name = "S"),
                              OpenMx::mxMatrix(type = "Full", nrow = lenrow, ncol = lencol, values = Fval, free = Ffree, labels = Flab, name = "F"),
                              OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = lencol, values = Mval, free = Mfree, labels = Mlab, name = "M"),
                              OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"),
                              OpenMx::mxRAMObjective("A", "S", "F", "M", dimnames = fcollab), constraint, OpenMx::mxCI(c("es"), interval = conf_level)
    )

    onecovfit <- OpenMx::mxRun(onecov, intervals = TRUE, silent = TRUE)
    return(onecovfit@output$confidenceIntervals)
}

# Calculate a likelihood-based CI of ES (unequal cluster sizes)
.lik_ci_es_crd_unequal <- function(datawide, y_label, x_label, z_within_label = NULL, z_between_label = NULL, es_type = 1, icc_Y = 0.25, es = 0.5, total_var = 1, covariate = FALSE, icc_Z = 0.25, R2_within = 0.5, R2_between = 0.5, total_var_Z = 1, conf_level = 0.95) {
    if (!requireNamespace("OpenMx", quietly = TRUE)) stop("The package 'OpenMx' is needed; please install the package and try again.")

    tau <- icc_Y * total_var
    sigma <- (1 - icc_Y) * total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }

    if (covariate) {
        if (icc_Z == 0 && R2_between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
        if (icc_Z == 1 && R2_within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
        tauz <- total_var_Z * icc_Z
        sigmaz <- total_var_Z * (1 - icc_Z)
        gammazb <- sqrt(R2_between * tau / tauz)
        gammazw <- sqrt(R2_within * sigma / sigmaz)
        tau <- (1 - R2_between) * tau
        sigma <- (1 - R2_within) * sigma
    }

    ntreat <- sum(sapply(datawide, function(x) sum(x[, x_label])))
    totaln <- sum(sapply(datawide, nrow))
    probx <- ntreat / totaln

    FUNgroupsize <- function(dat, y, zw = NULL) {
        latent_label <- c("intcept", "slope")
        if (is.null(zw)) latent_label <- "intcept"

        frowlab <- c(y, x_label, z_between_label)
        fcollab <- c(frowlab, latent_label)
        lenrow <- length(frowlab)
        lencol <- length(fcollab)

        Alab <- matrix(NA, lencol, lencol)
        Aval <- matrix(0, lencol, lencol)
        Afree <- matrix(FALSE, lencol, lencol)
        colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab

        Alab["intcept", x_label] <- "groupdiff"
        if (!is.null(z_between_label)) Alab["intcept", z_between_label] <- "zbeffect"
        if (!is.null(zw)) Alab[y, "slope"] <- paste("data.", zw, sep = "")

        Aval["intcept", x_label] <- gamma1
        if (!is.null(z_between_label)) Aval["intcept", z_between_label] <- gammazb
        Aval[y, latent_label] <- 1

        Afree["intcept", x_label] <- TRUE
        if (!is.null(z_between_label)) Afree["intcept", z_between_label] <- TRUE

        Slab <- matrix(NA, lencol, lencol)
        Sval <- matrix(0, lencol, lencol)
        Sfree <- matrix(FALSE, lencol, lencol)
        colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
        diag(Slab)[1:length(y)] <- "l1error"
        diag(Sval)[1:length(y)] <- sigma
        diag(Sfree)[1:length(y)] <- TRUE
        Slab["intcept", "intcept"] <- "l2error"
        Sval["intcept", "intcept"] <- tau
        Sfree["intcept", "intcept"] <- TRUE
        Slab[x_label, x_label] <- "varx"
        Sval[x_label, x_label] <- probx * (1 - probx)
        Sfree[x_label, x_label] <- TRUE
        if (!is.null(z_between_label)) {
            Slab[c(x_label, z_between_label), c(x_label, z_between_label)] <- "covxzb"
            Slab[x_label, x_label] <- "varx"
            Slab[z_between_label, z_between_label] <- "varzb"
            Sval[c(x_label, z_between_label), c(x_label, z_between_label)] <- 0
            Sval[x_label, x_label] <- probx * (1 - probx)
            Sval[z_between_label, z_between_label] <- tauz
            Sfree[c(x_label, z_between_label), c(x_label, z_between_label)] <- TRUE
        }

        Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latent_label)))
        Flab <- matrix(NA, lenrow, lencol)
        Ffree <- matrix(FALSE, lenrow, lencol)
        colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
        rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

        Mlab <- c(rep(NA, length(y)), "meanX")
        Mval <- c(rep(0, length(y)), probx)
        Mfree <- c(rep(FALSE, length(y)), TRUE)

        if (!is.null(z_between_label)) {
            Mlab <- c(Mlab, "meanzb")
            Mval <- c(Mval, 0)
            Mfree <- c(Mfree, TRUE)
        }
        Mlab <- c(Mlab, "meanctrl")
        Mval <- c(Mval, 0)
        Mfree <- c(Mfree, TRUE)
        if (!is.null(zw)) {
            Mlab <- c(Mlab, "zweffect")
            Mval <- c(Mval, gammazw)
            Mfree <- c(Mfree, TRUE)
        }

        Mlab <- matrix(Mlab, 1, lencol)
        Mval <- matrix(Mval, 1, lencol)
        Mfree <- matrix(Mfree, 1, lencol)
        colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

        onecov <- OpenMx::mxModel(paste0("group", length(y)), type = "RAM",
                                  OpenMx::mxData(dat, type = "raw"),
                                  OpenMx::mxMatrix(type = "Full", nrow = lencol, ncol = lencol, values = Aval, free = Afree, labels = Alab, name = "A"),
                                  OpenMx::mxMatrix(type = "Symm", nrow = lencol, ncol = lencol, values = Sval, free = Sfree, labels = Slab, name = "S"),
                                  OpenMx::mxMatrix(type = "Full", nrow = lenrow, ncol = lencol, values = Fval, free = Ffree, labels = Flab, name = "F"),
                                  OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = lencol, values = Mval, free = Mfree, labels = Mlab, name = "M"),
                                  OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"),
                                  OpenMx::mxRAMObjective("A", "S", "F", "M", dimnames = fcollab)
        )
        return(onecov)
    }
    varzw <- 0
    if (!is.null(z_within_label)) varzw <- weighted.mean(do.call(c, mapply(function(x, y) var(as.vector(as.matrix(x[, y]))), x = datawide, y = z_within_label, SIMPLIFY = FALSE)), as.numeric(names(datawide)))

    # These lines prevent "no visible binding" notes when OpenMx evaluates the mxAlgebra expressions.
    groupdiff <- NULL
    l1error <- NULL
    l2error <- NULL
    zbeffect <- NULL
    varzb <- NULL
    zweffect <- NULL

    constraint <- NULL
    if (es_type == 0) {
        if (is.null(z_within_label)) {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        } else {
            if (is.null(z_between_label)) {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
            } else {
                constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
            }
        }
    } else if (es_type == 1) {
        if (is.null(z_within_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l1error + (zweffect^2 * varzw)), name = "es")
        }
    } else if (es_type == 2) {
        if (is.null(z_between_label)) {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error), name = "es")
        } else {
            constraint <- OpenMx::mxAlgebra(expression = groupdiff / sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
        }
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    list_model <- NULL
    if (!is.null(z_within_label)) {
        list_model <- mapply(FUNgroupsize, dat = datawide, y = y_label, zw = z_within_label)
    } else {
        list_model <- mapply(FUNgroupsize, dat = datawide, y = y_label)
    }
    title <- "Effect Size CRD"
    # Build the summed objective from the per-group objectives as one
    # parsed expression. (The old idiom of constructing an empty
    # mxAlgebra("") and patching its @formula slot fails on current
    # OpenMx, which validates the expression eagerly.)
    groupnames <- paste0("group", names(datawide), ".objective")
    algebra <- OpenMx::mxAlgebraFromString(
        paste0("sum(", paste(groupnames, collapse = ", "), ")"),
        name = "allobjective")
    objective <- OpenMx::mxAlgebraObjective("allobjective")
    finalmodel <- OpenMx::mxModel(title, OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, values = varzw, free = FALSE, labels = "varzw", name = "J"), unlist(list_model), constraint, algebra, objective, OpenMx::mxCI(c("es"), interval = conf_level))
    finalmodelfit <- OpenMx::mxRun(finalmodel, intervals = TRUE, silent = TRUE)
    return(finalmodelfit@output$confidenceIntervals)
}

# Create dataset from the CRD model and transform it to the wide format
.create_data_crd_wide <- function(n_clusters, n_treat_clus, n_individuals, icc_Y, es, es_type = 1, total_var = 1, covariate = FALSE, icc_Z = NULL, R2_within = NULL, R2_between = NULL, total_var_Z = 1, diff_size = NULL) {
    dat <- .create_data_crd(n_clusters = n_clusters, n_treat_clus = n_treat_clus, n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = covariate, icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = total_var_Z, diff_size = diff_size)
    datawide <- NULL
    if (!is.null(diff_size)) {
        if (covariate) {
            if (icc_Z == 0) {
                datawide <- .wide_format_unequal(dat, 3, c(2, 4), 1)
            } else if (icc_Z == 1) {
                datawide <- .wide_format_unequal(dat, c(3, 5), 2, 1)
            } else {
                datawide <- .wide_format_unequal(dat, c(3, 5), c(2, 4), 1)
            }
        } else {
            datawide <- .wide_format_unequal(dat, 3, 2, 1)
        }
    } else {
        if (covariate) {
            if (icc_Z == 0) {
                datawide <- .wide_format(dat, 3, c(2, 4), 1)
            } else if (icc_Z == 1) {
                datawide <- .wide_format(dat, c(3, 5), 2, 1)
            } else {
                datawide <- .wide_format(dat, c(3, 5), c(2, 4), 1)
            }
        } else {
            datawide <- .wide_format(dat, 3, 2, 1)
        }
    }
    return(datawide)
}

# Reshape long-format data to wide format
.wide_format <- function(data, between_col, within_col, id_col) {
    temp <- split(data[, within_col], data[, id_col])
    temp <- lapply(temp, function(x) as.vector(as.matrix(x)))
    dataw <- do.call(rbind, temp)
    datab <- as.matrix(data[match(unique(data[, id_col]), data[, id_col]), between_col])
    colnames(datab) <- colnames(data)[between_col]
    n_individuals <- nrow(data) / nrow(datab)
    colnames(dataw) <- paste(rep(colnames(data)[within_col], each = n_individuals), rep(1:n_individuals, length(within_col)), sep = "")
    return(data.frame(datab, dataw))
}

# Reshape long-format data to wide format with unequal cluster sizes
.wide_format_unequal <- function(data, between_col, within_col, id_col) {
    temp <- split(data[, within_col], data[, id_col])
    temp <- lapply(temp, function(x) as.vector(as.matrix(x)))

    size <- sapply(temp, length) / length(within_col)
    dataw <- lapply(split(temp, size), function(x) do.call(rbind, x))
    datab <- split(data[match(unique(data[, id_col]), data[, id_col]), between_col], size)
    resultdat <- mapply(data.frame, datab, dataw)
    varnamesw <- lapply(sapply(dataw, ncol) / length(within_col), function(x) paste(rep(colnames(data)[within_col], each = x), rep(1:x, length(within_col)), sep = ""))
    varnames <- lapply(varnamesw, function(x) c(colnames(data)[between_col], x))
    resultdat <- mapply(function(x, y) { colnames(x) <- y; x }, x = resultdat, y = varnames)
    return(resultdat)
}

# Create dataset from the CRD model
.create_data_crd <- function(n_clusters, n_treat_clus, n_individuals, icc_Y, es, es_type = 1, total_var = 1, covariate = FALSE, icc_Z = NULL, R2_within = NULL, R2_between = NULL, total_var_Z = 1, diff_size = NULL) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

    if (!is.null(diff_size)) {
        n_individuals <- .find_n_indiv_vec(n_individuals, diff_size, n_clusters)
    } else {
        n_individuals <- rep(n_individuals, each = n_clusters)
    }
    id <- do.call(c, mapply(rep, 1:n_clusters, n_individuals, SIMPLIFY = FALSE))
    x <- c(rep(1, n_treat_clus), rep(0, n_clusters - n_treat_clus))
    tau <- icc_Y * total_var
    sigma <- (1 - icc_Y) * total_var
    gamma1 <- NULL
    if (es_type == 0) {
        gamma1 <- es * sqrt(total_var)
    } else if (es_type == 1) {
        gamma1 <- es * sqrt(sigma)
    } else if (es_type == 2) {
        gamma1 <- es * sqrt(tau)
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    if (covariate) {
        if (icc_Z == 0 && R2_between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
        if (icc_Z == 1 && R2_within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
        tauz <- total_var_Z * icc_Z
        sigmaz <- total_var_Z * (1 - icc_Z)
        gammazb <- sqrt(R2_between * tau / tauz)
        gammazw <- sqrt(R2_within * sigma / sigmaz)
        tau <- (1 - R2_between) * tau
        sigma <- (1 - R2_within) * sigma
    }
    ybetween <- (gamma1 * x) + rnorm(n_clusters, 0, sqrt(tau))
    if (covariate && icc_Z != 0) {
        zbetween <- rnorm(n_clusters, 0, sqrt(tauz))
        ybetween <- ybetween + gammazb * zbetween
    }
    ybetween <- do.call(c, mapply(rep, ybetween, n_individuals, SIMPLIFY = FALSE))
    y <- ybetween + rnorm(sum(n_individuals), 0, sqrt(sigma))
    if (covariate && icc_Z != 1) {
        zwithin <- rnorm(sum(n_individuals), 0, sqrt(sigmaz))
        y <- y + gammazw * zwithin
    }
    x <- do.call(c, mapply(rep, x, n_individuals, SIMPLIFY = FALSE))
    dat <- data.frame(id, y, x)
    if (covariate) {
        z <- NULL
        if (icc_Z == 0) {
            z <- zwithin
        } else if (icc_Z == 1) {
            z <- do.call(c, mapply(rep, zbetween, n_individuals, SIMPLIFY = FALSE))
        } else {
            z <- do.call(c, mapply(rep, zbetween, n_individuals, SIMPLIFY = FALSE)) + zwithin
        }
        if (icc_Z == 0) {
            zb <- 0
            zw <- z
        } else {
            zlist <- split(z, id)
            zb <- do.call(c, mapply(rep, sapply(zlist, mean), n_individuals, SIMPLIFY = FALSE))
            zw <- do.call(c, lapply(zlist, scale, scale = FALSE))
        }
        z <- data.frame(zw = zw, zb = zb)
        dat <- data.frame(dat, z)
    }
    rownames(dat) <- NULL
    return(dat)
}

# Find the number of clusters given the specified width and the cluster size
.find_n_clus_crd_diff <- function(width, n_individuals, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- seq(100, 1000, 100)
    result <- sapply(n_clusters, .find_width_crd_diff, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)

    if (all(width > result)) {
        n_clusters <- seq(ceiling(2 * (1 / min(pr_treat, 1 - pr_treat))) + num_predictors, 100, 1)
        result <- sapply(n_clusters, .find_width_crd_diff, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width > result)) {
            return(3 + num_predictors)
        } else if (all(width < result)) {
            return(100)
        } else {
            return(n_clusters[which(width > result)[1]])
        }
    } else if (all(width < result)) {
        start <- 1000
        repeat {
            n_clusters <- seq(start, start + 1000, 100)
            result <- sapply(n_clusters, .find_width_crd_diff, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
            if (all(width > result)) {
                return(start)
            } else if (all(width < result)) {
                start <- start + 1000
            } else {
                minval <- n_clusters[which(width > result)[1] - 1]
                maxval <- n_clusters[which(width > result)[1]]
                n_clusters <- seq(minval, maxval, 1)
                result <- sapply(n_clusters, .find_width_crd_diff, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
                if (all(width > result)) {
                    return(minval)
                } else if (all(width < result)) {
                    return(maxval)
                } else {
                    return(n_clusters[which(width > result)[1]])
                }
            }
        }
    } else {
        minval <- n_clusters[which(width > result)[1] - 1]
        maxval <- n_clusters[which(width > result)[1]]
        n_clusters <- seq(minval, maxval, 1)
        result <- sapply(n_clusters, .find_width_crd_diff, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width > result)) {
            return(minval)
        } else if (all(width < result)) {
            return(maxval)
        } else {
            return(n_clusters[which(width > result)[1]])
        }
    }
}

# Find width of the CI of the unstandardized conditions difference
.find_width_crd_diff <- function(n_clusters, n_individuals, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- as.numeric(n_clusters)
    pr_treat <- round(n_clusters * pr_treat) / n_clusters
    n_individuals <- as.numeric(n_individuals)
    if (!is.null(diff_size)) n_individuals <- .find_n_indiv_vec(n_individuals, diff_size, n_clusters)
    v <- .var_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance)
    alpha <- 1 - conf_level
    width <- 2 * sqrt(v) * qt(1 - alpha / 2, n_clusters - 2 - num_predictors)
    return(width)
}

# Find the variance of the conditions difference with or without a covariate
.var_crd_diff <- function(n_clusters, n_individuals, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL) {
    if (is.null(tau_Y)) {
        if (is.null(icc_Y)) {
            tau_Y <- total_var - sigma2_Y
        } else {
            tau_Y <- icc_Y * total_var
        }
    }
    if (is.null(sigma2_Y)) {
        if (is.null(icc_Y)) {
            sigma2_Y <- total_var - tau_Y
        } else {
            sigma2_Y <- (1 - icc_Y) * total_var
        }
    }
    if (length(n_individuals) > 1) {
        sigma2_Y <- sigma2_Y * (1 - R2_within)
        tau_Y <- tau_Y * (1 - R2_between)
        D <- 1 / sigma2_Y
        F <- 0 - (tau_Y) / (sigma2_Y * (sigma2_Y + n_individuals * tau_Y))
        ntreat <- round(pr_treat * n_clusters)
        pr_treat <- ntreat / n_clusters
        terms <- n_individuals * (D + n_individuals * F)
        termstotal <- sum(terms)
        termstreatment <- sum(terms[1:ntreat])
        vargamma1 <- termstotal / (termstotal * termstreatment - termstreatment^2)
        varclusmean <- vargamma1 * (n_clusters * pr_treat * (1 - pr_treat))
    } else {
        ntreat <- round(pr_treat * n_clusters)
        pr_treat <- ntreat / n_clusters
        varclusmean <- ((sigma2_Y * (1 - R2_within)) / n_individuals) + (tau_Y * (1 - R2_between))
    }
    if (!is.null(assurance)) {
        df <- n_clusters - 2 - num_predictors
        varclusmean <- varclusmean * qchisq(assurance, df) / df
    }
    denominator <- n_clusters * pr_treat * (1 - pr_treat)
    return(varclusmean / denominator)
}

# Find the total cost of the whole cluster randomized design
.cost_crd <- function(n_clusters, n_individuals, clus_cost, indiv_cost, diff_size = NULL) {
    n_individuals[n_individuals == "> 100000"] <- Inf
    n_individuals <- as.numeric(n_individuals)
    if (!is.null(diff_size)) n_individuals <- .find_n_indiv_vec(n_individuals, diff_size, n_clusters)
    result <- NULL
    if (length(n_individuals) > 1) {
        result <- (n_clusters * clus_cost) + sum(n_individuals * indiv_cost)
    } else {
        result <- n_clusters * (clus_cost + (n_individuals * indiv_cost))
    }
    return(result)
}

# Find the cluster size given the specified width of ES and the number of clusters
.find_n_indiv_crd_es <- function(width, n_clusters, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    if (num_predictors > 0 && is.null(icc_Z)) icc_Z <- icc_Y
    if (num_predictors > 1) stop("Only one predictor is allowed.")
    total_var <- 1
    if (es_type == 0) {
        total_var <- 1
    } else if (es_type == 1) {
        total_var <- 1 / (1 - icc_Y)
    } else if (es_type == 2) {
        total_var <- 1 / icc_Y
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    startval <- .find_n_indiv_crd_diff(width = width, n_clusters = n_clusters, pr_treat = pr_treat, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (startval == "> 100000") stop("The starting number of individuals is > 100,000. With the specified number of clusters, it seems impossible to get the specified width.")
    startval <- as.numeric(startval)
    startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = startval, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
    if (startwidth < width) {
        repeat {
            startval <- startval - 1
            if (startval == 1) return(c(startval + 1, startwidth))
            savedwidth <- startwidth
            startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = startval, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
            if (startwidth > width) return(c(startval + 1, savedwidth))
        }
    } else if (startwidth > width) {
        repeat {
            startval <- startval + 1
            startwidth <- .find_width_crd_es(nrep, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = startval, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
            if (startwidth < width) return(c(startval, startwidth))
        }
    } else {
        return(c(startval, startwidth))
    }
}

# Find the cluster size given the specified width and the number of clusters
.find_n_indiv_crd_diff <- function(width, n_clusters, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    baremaximum <- .find_width_crd_diff(n_clusters = n_clusters, n_individuals = 100000, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (width < baremaximum) return("> 100000")
    n_individuals <- seq(100, 1000, 100)
    result <- sapply(n_individuals, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    if (all(width > result)) {
        n_individuals <- seq(2, 100, 1)
        result <- sapply(n_individuals, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width > result)) {
            return(2)
        } else if (all(width < result)) {
            return(100)
        } else {
            return(n_individuals[which(width > result)[1]])
        }
    } else if (all(width < result)) {
        start <- 1000
        repeat {
            n_individuals <- seq(start, start + 1000, 100)
            result <- sapply(n_individuals, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
            if (all(width > result)) {
                return(start)
            } else if (all(width < result)) {
                start <- start + 1000
            } else {
                minval <- n_individuals[which(width > result)[1] - 1]
                maxval <- n_individuals[which(width > result)[1]]
                n_individuals <- seq(minval, maxval, 1)
                result <- sapply(n_individuals, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
                if (all(width > result)) {
                    return(minval)
                } else if (all(width < result)) {
                    return(maxval)
                } else {
                    return(n_individuals[which(width > result)[1]])
                }
            }
        }
    } else {
        minval <- n_individuals[which(width > result)[1] - 1]
        maxval <- n_individuals[which(width > result)[1]]
        n_individuals <- seq(minval, maxval, 1)
        result <- sapply(n_individuals, .find_width_crd_diff, n_clusters = n_clusters, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        if (all(width > result)) {
            return(minval)
        } else if (all(width < result)) {
            return(maxval)
        } else {
            return(n_individuals[which(width > result)[1]])
        }
    }
}

# Find the number of clusters given the budget and the cluster size
.find_n_clus_crd_budget <- function(budget, n_individuals, clus_cost, indiv_cost, diff_size = NULL) {
    n_clusters <- seq(100, 1000, 100)
    result <- sapply(n_clusters, .cost_crd, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    if (all(budget < result)) {
        n_clusters <- seq(4, 100, 1)
        result <- sapply(n_clusters, .cost_crd, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(NA)
        } else if (all(budget > result)) {
            return(100)
        } else {
            index <- which(budget >= result)
            return(n_clusters[index[length(index)]])
        }
    } else if (all(budget > result)) {
        start <- 1000
        repeat {
            n_clusters <- seq(start, start + 1000, 100)
            result <- sapply(n_clusters, .cost_crd, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            if (all(budget < result)) {
                return(start)
            } else if (all(budget > result)) {
                start <- start + 1000
            } else {
                minval <- n_clusters[which(budget < result)[1] - 1]
                maxval <- n_clusters[which(budget < result)[1]]
                n_clusters <- seq(minval, maxval, 1)
                result <- sapply(n_clusters, .cost_crd, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
                if (all(budget < result)) {
                    return(minval)
                } else if (all(budget > result)) {
                    return(maxval)
                } else {
                    index <- which(budget >= result)
                    return(n_clusters[index[length(index)]])
                }
            }
        }
    } else {
        minval <- n_clusters[which(budget < result)[1] - 1]
        maxval <- n_clusters[which(budget < result)[1]]
        n_clusters <- seq(minval, maxval, 1)
        result <- sapply(n_clusters, .cost_crd, n_individuals = n_individuals, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(minval)
        } else if (all(budget > result)) {
            return(maxval)
        } else {
            index <- which(budget >= result)
            return(n_clusters[index[length(index)]])
        }
    }
}

# Find the cluster size given the budget and the number of clusters
.find_n_indiv_crd_budget <- function(budget, n_clusters, clus_cost, indiv_cost, diff_size = NULL) {
    n_individuals <- seq(100, 1000, 100)
    result <- sapply(n_individuals, .cost_crd, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    if (all(budget < result)) {
        n_individuals <- seq(2, 100, 1)
        result <- sapply(n_individuals, .cost_crd, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(NA)
        } else if (all(budget > result)) {
            return(100)
        } else {
            index <- which(budget >= result)
            return(n_individuals[index[length(index)]])
        }
    } else if (all(budget > result)) {
        start <- 1000
        repeat {
            n_individuals <- seq(start, start + 1000, 100)
            result <- sapply(n_individuals, .cost_crd, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            if (all(budget < result)) {
                return(start)
            } else if (all(budget > result)) {
                start <- start + 1000
            } else {
                minval <- n_individuals[which(budget < result)[1] - 1]
                maxval <- n_individuals[which(budget < result)[1]]
                n_individuals <- seq(minval, maxval, 1)
                result <- sapply(n_individuals, .cost_crd, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
                if (all(budget < result)) {
                    return(minval)
                } else if (all(budget > result)) {
                    return(maxval)
                } else {
                    index <- which(budget >= result)
                    return(n_individuals[index[length(index)]])
                }
            }
        }
    } else {
        minval <- n_individuals[which(budget < result)[1] - 1]
        maxval <- n_individuals[which(budget < result)[1]]
        n_individuals <- seq(minval, maxval, 1)
        result <- sapply(n_individuals, .cost_crd, n_clusters = n_clusters, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        if (all(budget < result)) {
            return(minval)
        } else if (all(budget > result)) {
            return(maxval)
        } else {
            index <- which(budget >= result)
            return(n_individuals[index[length(index)]])
        }
    }
}

# Find the least width of ES combination of the number of clusters and cluster size given the budget
.find_min_width_crd_es <- function(budget, clus_cost = 0, indiv_cost = 1, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    if (num_predictors > 0 && is.null(icc_Z)) icc_Z <- icc_Y
    if (num_predictors > 1) stop("Only one predictor is allowed.")
    total_var <- 1
    if (es_type == 0) {
        total_var <- 1
    } else if (es_type == 1) {
        total_var <- 1 / (1 - icc_Y)
    } else if (es_type == 2) {
        total_var <- 1 / icc_Y
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    FUN <- function(n_clusters, n_individuals) {
        .find_width_crd_es(nrep = nrep, assurance = assurance, n_clusters = n_clusters, n_treat_clus = round(n_clusters * pr_treat), n_individuals = n_individuals, icc_Y = icc_Y, es = es, es_type = es_type, total_var = total_var, covariate = as.logical(num_predictors), icc_Z = icc_Z, R2_within = R2_within, R2_between = R2_between, total_var_Z = 1, seed = seed, multicore = multicore, num_proc = num_proc, conf_level = conf_level, diff_size = diff_size)
    }
    startval <- .find_min_width_crd_diff(budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, pr_treat = pr_treat, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)

    start_n_clusters <- c(startval[1] - 1, startval[1], startval[1] + 1)
    result_n_individuals <- sapply(start_n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    startwidth <- mapply(FUN, n_clusters = start_n_clusters, n_individuals = result_n_individuals)

    if (which(startwidth == min(startwidth)) == 1) {
        repeat {
            start_n_clusters <- start_n_clusters - 1
            result_n_individuals[2:3] <- result_n_individuals[1:2]
            startwidth[2:3] <- startwidth[1:2]
            result_n_individuals[1] <- .find_n_indiv_crd_budget(start_n_clusters[1], budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            startwidth[1] <- FUN(n_clusters = start_n_clusters[1], n_individuals = result_n_individuals[1])
            if (which(startwidth == min(startwidth)) != 1) return(c(start_n_clusters[2], result_n_individuals[2], startwidth[2]))
        }
    } else if (which(startwidth == min(startwidth)) == 3) {
        repeat {
            start_n_clusters <- start_n_clusters + 1
            result_n_individuals[1:2] <- result_n_individuals[2:3]
            startwidth[1:2] <- startwidth[2:3]
            result_n_individuals[3] <- .find_n_indiv_crd_budget(start_n_clusters[3], budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            startwidth[3] <- FUN(n_clusters = start_n_clusters[3], n_individuals = result_n_individuals[3])
            if (which(startwidth == min(startwidth)) != 3) return(c(start_n_clusters[2], result_n_individuals[2], startwidth[2]))
        }
    } else {
        return(c(start_n_clusters[2], result_n_individuals[2], startwidth[2]))
    }
}

# Find the least width combination of the number of clusters and cluster size given the budget
.find_min_width_crd_diff <- function(budget, clus_cost = 0, indiv_cost = 1, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    FUN <- function(n_clusters, n_individuals) {
        .find_width_crd_diff(n_clusters = n_clusters, n_individuals = n_individuals, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    }
    n_clusters <- seq(100, 1100, 100)
    n_individuals <- sapply(n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
    if (all(is.na(n_individuals))) {
        posmin <- 1
    } else {
        resultwidth <- rep(NA, length(n_clusters))
        resultwidth[!is.na(n_individuals)] <- mapply(FUN, n_clusters = n_clusters[!is.na(n_individuals)], n_individuals = n_individuals[!is.na(n_individuals)], SIMPLIFY = TRUE)
        posmin <- which(resultwidth == min(resultwidth, na.rm = TRUE))
    }
    if (posmin == 1) {
        n_clusters <- seq(ceiling(2 * (1 / min(pr_treat, 1 - pr_treat))) + num_predictors, 200, 1)
        n_individuals <- sapply(n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        resultwidth <- rep(NA, length(n_clusters))
        resultwidth[!is.na(n_individuals)] <- mapply(FUN, n_clusters = n_clusters[!is.na(n_individuals)], n_individuals = n_individuals[!is.na(n_individuals)])
        posmin <- which(resultwidth == min(resultwidth, na.rm = TRUE))
        return(c(n_clusters[posmin], n_individuals[posmin], resultwidth[posmin]))
    } else if (posmin == length(resultwidth)) {
        start <- n_clusters[length(n_clusters) - 1]
        repeat {
            n_clusters <- seq(start, start + 1100, 100)
            n_individuals <- sapply(n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            resultwidth <- rep(NA, length(n_clusters))
            resultwidth[!is.na(n_individuals)] <- mapply(FUN, n_clusters = n_clusters[!is.na(n_individuals)], n_individuals = n_individuals[!is.na(n_individuals)])
            posmin <- which(resultwidth == min(resultwidth, na.rm = TRUE))
            if (posmin == 1) {
                return(c(n_clusters[posmin], n_individuals[posmin], resultwidth[posmin]))
            } else if (posmin == length(resultwidth)) {
                start <- start + 1000
            } else {
                minval <- n_clusters[posmin - 1]
                maxval <- n_clusters[posmin + 1]
                n_clusters <- seq(minval, maxval, 1)
                n_individuals <- sapply(n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
                resultwidth <- rep(NA, length(n_clusters))
                resultwidth[!is.na(n_individuals)] <- mapply(FUN, n_clusters = n_clusters[!is.na(n_individuals)], n_individuals = n_individuals[!is.na(n_individuals)])
                posmin <- which(resultwidth == min(resultwidth, na.rm = TRUE))
                return(c(n_clusters[posmin], n_individuals[posmin], resultwidth[posmin]))
            }
        }
    } else {
        minval <- n_clusters[posmin - 1]
        maxval <- n_clusters[posmin + 1]
        n_clusters <- seq(minval, maxval, 1)
        n_individuals <- sapply(n_clusters, .find_n_indiv_crd_budget, budget = budget, clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
        resultwidth <- rep(NA, length(n_clusters))
        resultwidth[!is.na(n_individuals)] <- mapply(FUN, n_clusters = n_clusters[!is.na(n_individuals)], n_individuals = n_individuals[!is.na(n_individuals)])
        posmin <- which(resultwidth == min(resultwidth, na.rm = TRUE))
        return(c(n_clusters[posmin], n_individuals[posmin], resultwidth[posmin]))
    }
}

# Find the least expensive combination of the number of clusters and cluster size given the specified width of ES
.find_min_cost_crd_es <- function(width, clus_cost = 0, indiv_cost = 1, es, es_type = 1, icc_Y, pr_treat, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, nrep = 1000, icc_Z = NULL, seed = NULL, multicore = FALSE, num_proc = NULL, diff_size = NULL) {
    if (num_predictors > 0 && is.null(icc_Z)) icc_Z <- icc_Y
    if (num_predictors > 1) stop("Only one predictor is allowed.")
    total_var <- 1
    if (es_type == 0) {
        total_var <- 1
    } else if (es_type == 1) {
        total_var <- 1 / (1 - icc_Y)
    } else if (es_type == 2) {
        total_var <- 1 / icc_Y
    } else {
        stop("'es_type' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
    }
    startval <- .find_min_cost_crd_diff(width = width, clus_cost = clus_cost, indiv_cost = indiv_cost, pr_treat = pr_treat, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
    startval <- as.numeric(startval)
    start_n_individuals <- c(startval[2] - 1, startval[2], startval[2] + 1)
    result <- sapply(start_n_individuals, .find_n_clus_crd_es, width = width, es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size)
    result_n_clusters <- result[1, ]
    resultwidth <- result[2, ]
    startbudget <- mapply(.cost_crd, n_clusters = result_n_clusters, n_individuals = start_n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size))
    if (which(startbudget == min(startbudget)) == 1) {
        repeat {
            start_n_individuals <- start_n_individuals - 1
            result_n_clusters[2:3] <- result_n_clusters[1:2]
            startbudget[2:3] <- startbudget[1:2]
            resultwidth[2:3] <- resultwidth[1:2]
            result <- .find_n_clus_crd_es(width = width, n_individuals = start_n_individuals[1], es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size)
            result_n_clusters[1] <- result[1]
            resultwidth[1] <- result[2]
            startbudget[1] <- .cost_crd(n_clusters = result_n_clusters[1], n_individuals = start_n_individuals[1], clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            if (which(startbudget == min(startbudget)) != 1) return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
        }
    } else if (which(startbudget == min(startbudget)) == 3) {
        repeat {
            start_n_individuals <- start_n_individuals + 1
            result_n_clusters[1:2] <- result_n_clusters[2:3]
            startbudget[1:2] <- startbudget[2:3]
            resultwidth[1:2] <- resultwidth[2:3]
            result <- .find_n_clus_crd_es(width = width, n_individuals = start_n_individuals[3], es = es, es_type = es_type, icc_Y = icc_Y, pr_treat = pr_treat, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, nrep = nrep, icc_Z = icc_Z, seed = seed, multicore = multicore, num_proc = num_proc, diff_size = diff_size)
            result_n_clusters[3] <- result[1]
            resultwidth[3] <- result[2]
            startbudget[3] <- .cost_crd(n_clusters = result_n_clusters[3], n_individuals = start_n_individuals[3], clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size)
            if (which(startbudget == min(startbudget)) != 3) return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
        }
    } else {
        return(c(result_n_clusters[2], start_n_individuals[2], startbudget[2], resultwidth[2]))
    }
}

# Find the least expensive combination of the number of clusters and cluster size given the specified width
.find_min_cost_crd_diff <- function(width, clus_cost = 0, indiv_cost = 1, pr_treat, tau_Y = NULL, sigma2_Y = NULL, total_var = NULL, icc_Y = NULL, R2_between = 0, R2_within = 0, num_predictors = 0, assurance = NULL, conf_level = 0.95, diff_size = NULL) {
    n_clusters <- seq(100, 1100, 100)
    repeat {
        n_individuals <- sapply(n_clusters, .find_n_indiv_crd_diff, width = width, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        cost <- mapply(.cost_crd, n_clusters = n_clusters, n_individuals = n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size), SIMPLIFY = TRUE)
        if (!all(cost == Inf)) break
        n_clusters <- n_clusters + 1000
    }
    posmin <- which(cost == min(cost))
    if (length(posmin) > 1) posmin <- posmin[1]
    if (posmin == 1) {
        n_clusters <- seq(ceiling(2 * (1 / min(pr_treat, 1 - pr_treat))) + num_predictors, 200, 1)
        n_individuals <- sapply(n_clusters, .find_n_indiv_crd_diff, width = width, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        cost <- mapply(.cost_crd, n_clusters = n_clusters, n_individuals = n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size), SIMPLIFY = TRUE)
        posmin <- which(cost == min(cost))
        if (length(posmin) > 1) posmin <- posmin[1]
        return(c(n_clusters[posmin], n_individuals[posmin], cost[posmin]))
    } else if (posmin == length(cost)) {
        start <- n_clusters[length(n_clusters) - 1]
        repeat {
            n_clusters <- seq(start, start + 1100, 100)
            n_individuals <- sapply(n_clusters, .find_n_indiv_crd_diff, width = width, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
            cost <- mapply(.cost_crd, n_clusters = n_clusters, n_individuals = n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size), SIMPLIFY = TRUE)
            posmin <- which(cost == min(cost))
            if (length(posmin) > 1) posmin <- posmin[1]
            if (posmin == 1) {
                return(c(n_clusters[posmin], n_individuals[posmin], cost[posmin]))
            } else if (posmin == length(cost)) {
                start <- start + 1000
            } else {
                minval <- n_clusters[posmin - 1]
                maxval <- n_clusters[posmin + 1]
                n_clusters <- seq(minval, maxval, 1)
                n_individuals <- sapply(n_clusters, .find_n_indiv_crd_diff, width = width, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
                cost <- mapply(.cost_crd, n_clusters = n_clusters, n_individuals = n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size), SIMPLIFY = TRUE)
                posmin <- which(cost == min(cost))
                if (length(posmin) > 1) posmin <- posmin[1]
                return(c(n_clusters[posmin], n_individuals[posmin], cost[posmin]))
            }
        }
    } else {
        minval <- n_clusters[posmin - 1]
        maxval <- n_clusters[posmin + 1]
        n_clusters <- seq(minval, maxval, 1)
        n_individuals <- sapply(n_clusters, .find_n_indiv_crd_diff, width = width, pr_treat = pr_treat, tau_Y = tau_Y, sigma2_Y = sigma2_Y, total_var = total_var, icc_Y = icc_Y, R2_between = R2_between, R2_within = R2_within, num_predictors = num_predictors, assurance = assurance, conf_level = conf_level, diff_size = diff_size)
        cost <- mapply(.cost_crd, n_clusters = n_clusters, n_individuals = n_individuals, MoreArgs = list(clus_cost = clus_cost, indiv_cost = indiv_cost, diff_size = diff_size), SIMPLIFY = TRUE)
        posmin <- which(cost == min(cost))
        if (length(posmin) > 1) posmin <- posmin[1]
        return(c(n_clusters[posmin], n_individuals[posmin], cost[posmin]))
    }
}
