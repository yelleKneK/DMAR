#' Heterotrait-Monotrait Ratio of Correlations (HTMT)
#'
#' Computes the HTMT discriminant-validity index of Henseler, Ringle, and
#' Sarstedt (2015) for every pair of constructs: the average correlation
#' between items of \emph{different} constructs, divided by the geometric
#' mean of the average correlations among items \emph{within} each
#' construct. Two constructs whose HTMT approaches 1 are empirically
#' indistinguishable however cleanly the model draws them; the customary
#' red flags are 0.85 (strict) or 0.90 (liberal). An optional bootstrap
#' gives a one-sided upper confidence bound, the quantity actually
#' compared against the cutoff in the validity literature.
#'
#' @param data A \code{data.frame} of item responses.
#' @param blocks Named list of character vectors: each element names a
#'   construct and gives its item columns (two or more per construct, two
#'   or more constructs). The blocks are disjoint: each item is listed
#'   once, in one block, and each construct has a name of its own (see
#'   \emph{Disjoint blocks} under \emph{Details}).
#' @param B Number of bootstrap resamples for the upper confidence bound,
#'   a single whole number; \code{0} (default) skips the bootstrap and
#'   reports the point estimates only. When positive it must be at least
#'   200 at the default \code{conf_level} (100 at 0.90, 1000 at 0.99): the
#'   bound is one quantile of the resamples, reported only when at least 10
#'   of them lie beyond it and never from fewer than 100. A smaller
#'   positive \code{B} is refused before anything is computed (class
#'   \code{dmar_count_invalid}). See \emph{Bootstrap replications} under
#'   \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param conf_level Confidence level for the one-sided upper bound.
#'   Defaults to 0.95.
#' @param seed Optional integer seed for the bootstrap, used locally (the
#'   caller's random number generator state is restored on exit).
#'
#' @details
#' All correlations are Pearson, computed on pairwise-complete
#' observations, and every one of them enters its average as an absolute
#' value, the monotrait correlations among the items of a block as well
#' as the heterotrait correlations between blocks. This is the ratio of
#' arithmetic means that \pkg{semTools}' \code{htmt()} computes with
#' \code{htmt2 = FALSE} and its default \code{absolute = TRUE}; on
#' incomplete data the two agree only when \pkg{semTools} is given the
#' same pairwise-complete correlations (through its \code{sample.cov}
#' argument), since by default it deletes incomplete rows listwise, and
#' its \code{missing = "pairwise"} option estimates the correlations
#' through \pkg{lavaan}, whose values differ from
#' \code{stats::cor(use = "pairwise.complete.obs")}. Absolute
#' values make the ratio indifferent to the direction in which an item
#' is keyed, since reversing an item changes the signs of its
#' correlations and not their sizes. With signed correlations, a reverse
#' keyed item that was not recoded lowers its block's average
#' within-construct correlation, which raises the ratio, while its
#' negative correlations with the items of the other block partly cancel
#' the positive ones in the numerator, so the ratio would depend on how
#' the items happen to be scored. Values near or above 1 indicate that
#' the two item sets correlate across constructs about as strongly as
#' within them. HTMT is a correlation-based screen, deliberately
#' model-free; the confirmatory companion is the latent correlation
#' between the two factors (see
#' \code{\link{correction_for_attenuation}} and its factor-model
#' discussion, and the \code{phi} rows of \code{\link{fa_confirmatory}} with
#' \code{output = "measurement"}).
#'
#' \strong{Disjoint blocks.} The ratio compares two sets of items, so the
#' blocks are disjoint sets of distinct items under distinct names. An
#' item listed twice in one block would enter that block's average with
#' its correlation of 1 with itself, and its other correlations would be
#' counted twice; an item placed in two blocks would enter the average
#' correlation between them with its correlation of 1 with itself; and
#' two blocks under one name could not be told apart. Each of these, and
#' an item whose name belongs to two or more columns of \code{data}, is
#' refused with an error of class \code{dmar_htmt_blocks_invalid} (see
#' \code{\link{dmar_metadata}}) before any correlation is computed or any
#' resample drawn, whatever \code{B}. The condition's \code{reason}
#' element names the check that failed (\code{"construct_names"},
#' \code{"repeated"}, \code{"shared"}, or \code{"data_names"}), and its
#' \code{blocks} and \code{items} elements the names at issue. A ratio
#' for item sets that overlap would be a different quantity, which this
#' function does not compute.
#'
#' A pair's ratio is undefined when a correlation it averages cannot be
#' computed: an item with no variance, or two items never observed
#' together, as in a planned missing data design whose forms give the
#' items of two blocks to different respondents. That pair's
#' \code{htmt} (and \code{ci_upper}) is \code{NA}, with a warning that
#' names it, and the other pairs are reported; the call stops only when
#' no pair is defined.
#'
#' The bootstrap, when requested (\code{B > 0}), resamples the rows of
#' \code{data} with replacement \code{B} times and recomputes every
#' pairwise HTMT on each resample; the reported \code{ci_upper} is
#' the \code{conf_level} empirical quantile of each pair's bootstrap
#' distribution, a one-sided upper percentile bound (Efron &
#' Tibshirani, 1993). That bound is the only interval offered, matching
#' how the validity literature uses HTMT (the question is whether the
#' ratio credibly exceeds the cutoff); no two-sided or bias-corrected
#' and accelerated (BCa) variant is provided. The ratio of each resample is
#' the one reported, with every correlation an absolute value, so a
#' resample that draws a negative correlation within a block enters with
#' its size. Bootstrap results
#' vary from run to run; supply \code{seed} for reproducibility.
#'
#' \strong{Bootstrap replications.} The bound of each pair is the
#' \code{conf_level} quantile of its \code{B} resampled ratios, so it is
#' run only when \code{B} puts at least 10 resamples beyond it (200 at the
#' 95 percent level). A replication in which the HTMT of some pair defined
#' on the data is undefined (an item that is constant in the resample, or
#' two items never observed together) is dropped and counted, never
#' replaced, and the bound is computed from the rest. The result records
#' \code{B} (requested), \code{B_used} (the number the bound rests on),
#' \code{n_dropped}, and \code{B_required} (the fewest replications a
#' request may ask for; that count is held to \code{B}, not to
#' \code{B_used}) as attributes, with \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}). One
#' warning of class \code{dmar_bootstrap_replications_dropped} states the
#' counts when any replication was dropped, with a caution when the dropped
#' replications are at least the 5 percent of \code{B} that lies beyond one
#' limit at the 95 percent level, because they could then decide by
#' themselves where that limit falls. When more than half were dropped, the
#' bound of every pair is
#' \code{NA}, the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row
#'   per construct pair and the same four columns whatever \code{B}:
#'   \code{construct_1}, \code{construct_2}, \code{htmt}, and
#'   \code{ci_upper}, the one-sided \code{conf_level} bootstrap
#'   percentile bound, which is \code{NA} when \code{B = 0} (no bootstrap
#'   was run). Both are \code{NA} for a pair whose ratio is undefined on
#'   the data (see Details). With \code{B > 0} the table also carries the
#'   attributes \code{conf_level}, \code{B} (requested), \code{B_used}
#'   (the number of resamples the bound was computed from),
#'   \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#'   (see \emph{Bootstrap replications} under Details).
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
#'   assessing discriminant validity in variance-based structural equation
#'   modeling. \emph{Journal of the Academy of Marketing Science, 43}(1),
#'   115--135. \doi{10.1007/s11747-014-0403-8}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{average_variance_extracted}} for the convergent
#'   side of the validity ledger; \code{\link{reliability_omega}} for the
#'   composite reliability of each block;
#'   \code{\link{correction_for_attenuation}} for the latent correlation
#'   route.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Two clean constructs and their six items.
#' set.seed(113)
#' n <- 300
#' f1 <- rnorm(n); f2 <- 0.3 * f1 + sqrt(1 - 0.09) * rnorm(n)
#' d <- data.frame(
#'   a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
#'   a3 = .6 * f1 + rnorm(n, 0, .8),
#'   b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
#'   b3 = .6 * f2 + rnorm(n, 0, .8))
#' h <- htmt(d, blocks = list(A = c("a1", "a2", "a3"),
#'                            B = c("b1", "b2", "b3")))
#' h
#'
#' # The broom verbs: one row per construct pair.
#' generics::tidy(h)
#' generics::glance(h)
#'
#' # The upper confidence bound, which is the quantity the validity
#' # literature compares against 0.85 or 0.90, comes from a bootstrap
#' # that recomputes every pairwise ratio on each of B resamples and
#' # fills the ci_upper column, NA above. B = 500 keeps the example
#' # quick; a claim about discriminant validity deserves the bound at
#' # B = 10000 rather than the point estimate alone.
#' htmt(d, blocks = list(A = c("a1", "a2", "a3"),
#'                       B = c("b1", "b2", "b3")),
#'      B = 500, seed = 113)
#'
#' @export
#' @importFrom stats cor quantile
htmt <- function(data, blocks, B = 0, conf_level = 0.95, seed = NULL) {
  if (!is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  if (!is.list(blocks) || length(blocks) < 2L ||
      is.null(names(blocks)) || anyNA(names(blocks)) ||
      any(names(blocks) == "")) {
    .dmar_refuse_argument(
      "'blocks' must be a named list of two or more constructs.",
      arg = "blocks")
  }
  # The blocks are disjoint sets of distinct items under distinct names,
  # checked before any correlation is computed or any resample drawn
  # (.htmt_check_blocks()).
  .htmt_check_blocks(data, blocks)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # B = 0 means no bootstrap. Otherwise the rule of ?dmar_bootstrap for a
  # one-sided bound at the conf_level quantile of the resamples: at least
  # 10 of them beyond it (200 at the 95 percent level).
  boot_tail <- min(conf_level, 1 - conf_level)
  B <- .dmar_check_B(B, tail = boot_tail, conf_level = conf_level,
                     zero_ok = TRUE, what = "bound", position = conf_level)

  # Every correlation enters as an absolute value, within the blocks as
  # well as between them, so the ratio does not depend on how an item is
  # keyed. Through 1.0.1 only the heterotrait correlations were absolute:
  # a reverse keyed item left unrecoded then shrank its block's average
  # within-construct correlation (or turned it negative, and the call
  # stopped) while the numerator kept its size, which inflated the ratio.
  # A correlation that cannot be computed (a constant item, or two items
  # never observed together) leaves the ratio of a pair undefined, and so
  # do items of a block uncorrelated with one another. cor()'s warning
  # about a zero standard deviation is muffled because the warning below
  # says what went wrong; htmt_from() returns NA for such a pair, with the
  # reason in its "undefined" attribute.
  htmt_from <- function(d) {
    R <- withCallingHandlers(
      abs(stats::cor(d[, unlist(blocks), drop = FALSE],
                     use = "pairwise.complete.obs")),
      warning = function(w) {
        if (grepl("standard deviation is zero", conditionMessage(w),
                  fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })
    pairs <- utils::combn(names(blocks), 2, simplify = FALSE)
    val <- rep(NA_real_, length(pairs))
    why <- character(length(pairs))
    for (k in seq_along(pairs)) {
      i <- blocks[[pairs[[k]][1]]]; j <- blocks[[pairs[[k]][2]]]
      hetero <- mean(R[i, j])
      mono_i <- mean(R[i, i][lower.tri(R[i, i])])
      mono_j <- mean(R[j, j][lower.tri(R[j, j])])
      if (anyNA(c(hetero, mono_i, mono_j))) {
        why[k] <- "cannot_compute"
      } else if (mono_i == 0 || mono_j == 0) {
        why[k] <- "uncorrelated"
      } else {
        val[k] <- hetero / sqrt(mono_i * mono_j)
      }
    }
    attr(val, "undefined") <- why
    val
  }

  est <- htmt_from(data)
  pairs <- utils::combn(names(blocks), 2)
  # A pair undefined on the data is NA, with one warning that names it;
  # the other pairs are reported. Through 1.0.1 a pair whose items were
  # never observed together across its two blocks was NA without a word,
  # and an undefined block (a constant item, or two of its items never
  # observed together) stopped the call with R's "missing value where
  # TRUE/FALSE needed"; the first versions of 1.1.0 stopped the call for
  # any undefined pair, which cost fa_confirmatory(output = "measurement") its
  # whole table. The call stops only when no pair is defined.
  why <- attr(est, "undefined")
  est <- as.numeric(est)
  undefined <- nzchar(why)
  if (any(undefined)) {
    msg <- .htmt_undefined_message(pairs[, undefined, drop = FALSE],
                                   why[undefined])
    if (all(undefined)) {
      .dmar_refuse(
        msg,
        "htmt_undefined")
    }
    .dmar_warn(c(
      msg, " ", if (sum(undefined) == 1L) "Its row is" else
"Their rows are", " NA."), "result_qualified")
  }

  # One shape whatever B: ci_upper is NA unless the bootstrap runs.
  out <- data.frame(
    construct_1 = pairs[1, ],
    construct_2 = pairs[2, ],
    htmt        = unname(est),
    ci_upper = NA_real_,
    stringsAsFactors = FALSE
  )

  if (B > 0) {
    .dmar_local_seed(seed)
    N <- nrow(data)
    # Only the pairs defined on the data are bootstrapped; an undefined
    # pair keeps an NA bound. A resample can draw a case mix in which some
    # correlation of a defined pair cannot be computed (an item constant
    # in the resample), which leaves HTMT undefined there; that resample
    # is dropped rather than the whole call aborted.
    boots <- vapply(seq_len(B), function(i) {
      tryCatch(
        as.numeric(htmt_from(data[sample.int(N, N, replace = TRUE), ,
                                  drop = FALSE])),
        error = function(e) rep(NA_real_, length(est)))
    }, numeric(length(est)))
    boots <- matrix(boots, nrow = length(est))[!undefined, , drop = FALSE]
    ok <- colSums(!is.finite(boots)) == 0L
    boots <- boots[, ok, drop = FALSE]
    # The verdict of the rule (R/dmar_bootstrap.R): more than half of
    # the resamples dropped leave the bound NA under a classed warning
    # (through 1.0.x the call stopped when fewer than 100 remained).
    verdict <- .dmar_boot_settle(
      B, sum(ok), tail = boot_tail, conf_level = conf_level, fun = "htmt",
      quantity = "the HTMT ratio of some pair",
      dropped_because = paste0("left HTMT undefined for some pair (a ",
                               "correlation that could not be computed, ",
                               "as when an item is constant in the ",
                               "resample)"),
      what = "bound", position = conf_level,
      replications = "bootstrap resamples")
    if (verdict$formed) {
      out$ci_upper[!undefined] <- apply(boots, 1, stats::quantile,
                                           probs = conf_level, names = FALSE)
    }
  }
  res <- .as_dmar_tbl(out)
  if (B > 0) {
    attr(res, "conf_level") <- conf_level
    res <- .dmar_boot_attrs(res, verdict)
  }
  res
}

# The message for the pairs whose HTMT is undefined on the data, one
# sentence per reason: a correlation among their items that cannot be
# computed (an item with no variance, or two items never observed
# together), or a block whose items are uncorrelated with one another.
# `pairs` is a two-row matrix of construct names, one column per pair;
# `why` holds "cannot_compute" or "uncorrelated" for each.
.htmt_undefined_message <- function(pairs, why) {
  each <- sprintf("'%s' and '%s'", pairs[1L, ], pairs[2L, ])
  reasons <- c(
    cannot_compute = paste0("a correlation among their items cannot be ",
                            "computed (an item with no variance, or two ",
                            "items never observed together)"),
    uncorrelated = "the items of a block are uncorrelated with one another")
  parts <- vapply(unique(why), function(w) {
    paste0("HTMT is undefined for ", paste(each[why == w], collapse = ", for "),
           ": ", reasons[[w]], ".")
  }, character(1L))
  paste(parts, collapse = " ")
}

# .htmt_check_blocks(data, blocks)
#
# The blocks of htmt() are disjoint sets of distinct items under distinct
# names. `blocks` is already a named list of two or more elements with no
# empty or NA name. Refuses, before any correlation is computed or any
# resample drawn:
#   - a block that is not two or more columns of `data` (class
#     dmar_argument_invalid, as through 1.0.x);
# and, with class dmar_htmt_blocks_invalid, whose `reason` says which
# check failed and whose `blocks` and `items` hold the names at issue:
#   - "construct_names": a construct name given to two or more blocks.
#     blocks[["A"]] reads the first block of that name, so both would
#     have been read as the first block, and their pair would have
#     compared that block with itself (1.533 on the same fixture);
#   - "repeated": an item listed twice in one block. Its correlation of 1
#     with itself would enter the block's average and its other
#     correlations would be counted twice (on the seed-113 fixture of the
#     deep review, x1 listed twice moved HTMT from .615 to .544);
#   - "shared": an item in two or more blocks. Its correlation of 1 with
#     itself would enter the heterotrait average (x3 in both blocks moved
#     the same ratio to .944);
#   - "data_names": an item that names two or more columns of `data`, so
#     the column meant is ambiguous (the first would have been read).
# Through 1.0.x none of the four was checked, and each gave a number
# without a word.
.htmt_check_blocks <- function(data, blocks) {
  constructs <- names(blocks)
  quoted <- function(x) paste0("'", x, "'")
  and <- .measurement_invariance_and
  refuse <- function(message, reason, blocks, items = character(0)) {
    .dmar_refuse(message, "htmt_blocks_invalid", reason = reason,
                 blocks = unique(blocks), items = unique(items),
                 arg = "blocks")
  }
  twice <- unique(constructs[duplicated(constructs)])
  if (length(twice)) {
    refuse(c(
      "Each block of 'blocks' needs a name of its own, but ",
      and(quoted(twice)), if (length(twice) == 1L) " names" else " each name",
      " more than one block. Rename the blocks so that each construct has ",
      "one name."), "construct_names", blocks = twice)
  }
  for (k in seq_along(blocks)) {
    bl <- blocks[[k]]
    if (!is.character(bl) || length(bl) < 2L || !all(bl %in% names(data))) {
      .dmar_refuse_argument(
        sprintf("Block '%s' must name two or more columns of 'data'.",
                constructs[k]), arg = "data")
    }
  }
  repeated <- lapply(blocks, function(v) unique(v[duplicated(v)]))
  at <- lengths(repeated) > 0L
  if (any(at)) {
    refuse(c(
      "Each item is listed once in its block, but ",
      paste0(vapply(repeated[at], function(v) and(v), character(1L)),
             ifelse(lengths(repeated[at]) == 1L, " is", " are"),
             " listed more than once in block ", quoted(constructs[at]),
             collapse = "; "),
      ". A repeated item would enter the block's average with its ",
      "correlation of 1 with itself and with its other correlations ",
      "counted twice. List each item once."),
      "repeated", blocks = constructs[at],
      items = unlist(repeated[at], use.names = FALSE))
  }
  items <- unlist(blocks, use.names = FALSE)
  owner <- rep(constructs, lengths(blocks))
  shared <- unique(items[duplicated(items)])
  if (length(shared)) {
    where <- vapply(shared, function(i) and(quoted(owner[items == i])),
                    character(1L))
    refuse(c(
      "HTMT compares blocks that share no item, but ",
      paste0(shared, " is in blocks ", where, collapse = "; "),
      ". An item in two blocks would enter the average correlation ",
      "between them with its correlation of 1 with itself. Give each item ",
      "to one block."),
      "shared", blocks = owner[items %in% shared], items = shared)
  }
  dup_cols <- unique(names(data)[duplicated(names(data))])
  ambiguous <- intersect(items, dup_cols)
  if (length(ambiguous)) {
    where <- unique(owner[items %in% ambiguous])
    one <- length(ambiguous) == 1L
    refuse(c(
      if (one) "The name " else "The names ", and(ambiguous),
      if (one) " belongs" else " each belong",
      " to more than one column of 'data', and ",
      if (length(where) == 1L) "block " else "blocks ", and(quoted(where)),
      if (length(where) == 1L) " lists " else " list ",
      if (one) "it" else "them", ", so the column meant is ambiguous. Give ",
      "the columns of 'data' distinct names."),
      "data_names", blocks = where, items = ambiguous)
  }
  invisible(NULL)
}
