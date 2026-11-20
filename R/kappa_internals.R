# Internal helpers shared by cohen_kappa() and fleiss_kappa().
#
#   .kappa_boot_limits()   the percentile or BCa limits of a kappa from its
#                          bootstrap replicates, under the rule of
#                          R/dmar_bootstrap.R


# .kappa_boot_limits(boots, B, kappa, kappa_at, N, conf_level, ci_method,
#                    fun)
#
# The limits of a bootstrap interval for a kappa coefficient.
#   boots     : the B replicated values, NA where chance agreement was 1
#               in the resample and kappa is undefined.
#   kappa     : the sample estimate.
#   kappa_at  : function(idx), kappa on the subjects `idx` (NA when
#               undefined), for the jackknife of the BCa acceleration.
#   N         : the number of subjects.
#   ci_method : "percentile" or "bca".
#   fun       : the calling function's name, for the messages.
# Returns list(limits, verdict). A replication with an undefined kappa is
# dropped and counted. The verdict (.dmar_boot_settle()) withholds the
# interval (limits NA) when more than half were dropped, and for BCa when
# the B requested puts fewer than 10 replications beyond a limit where the
# bias correction and the acceleration place it (the count rule is held to
# the request, which the caller checked).
# A BCa interval whose bias correction or acceleration is undefined is
# refused (class dmar_bca_undefined), as mediate() refuses it; through the
# first versions of 1.1.0 an undefined leave-one-out kappa was removed
# from the jackknife without a word.
.kappa_boot_limits <- function(boots, B, kappa, kappa_at, N, conf_level,
                               ci_method, fun) {
  alpha_2 <- (1 - conf_level) / 2
  kept <- boots[is.finite(boots)]
  verdict <- .dmar_boot_settle(
    B, length(kept), kind = if (ci_method == "bca") "bca" else "quantile",
    tail = alpha_2, conf_level = conf_level, fun = fun, quantity = "kappa",
    dropped_because = paste0("left kappa undefined (chance agreement of 1 ",
                             "in the resample)"),
    signal = FALSE)
  limits <- c(NA_real_, NA_real_)
  if (verdict$formed && ci_method == "percentile") {
    limits <- stats::quantile(kept, c(alpha_2, 1 - alpha_2), names = FALSE)
  } else if (verdict$formed) {
    # BCa: median bias from the bootstrap distribution, acceleration
    # from the jackknife (Efron & Tibshirani, 1993).
    z0 <- stats::qnorm(mean(kept < kappa))
    jack <- vapply(seq_len(N), function(i) kappa_at(seq_len(N)[-i]),
                   numeric(1))
    jm <- mean(jack)
    acc <- sum((jm - jack)^3) / (6 * (sum((jm - jack)^2))^1.5)
    zq <- stats::qnorm(c(alpha_2, 1 - alpha_2))
    adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
    if (!all(is.finite(jack)) || !is.finite(z0) || !is.finite(acc) ||
        !all(is.finite(adj))) {
      .dmar_refuse(c(
        fun, "(): the BCa interval is not available for these data: ",
        if (!all(is.finite(jack))) {
          paste0("its acceleration needs kappa with each case left out in ",
                 "turn, and kappa is not defined with case ",
                 paste(utils::head(which(!is.finite(jack)), 5L),
                       collapse = ", "),
                 " left out (chance agreement 1).")
        } else if (!is.finite(z0)) {
          paste0("every bootstrap estimate falls on one side of the ",
                 "sample estimate, so the bias correction is infinite.")
        } else {
          "its acceleration is not finite."
        },
        " ci_method = \"percentile\" does not need these quantities."),
        "bca_undefined")
    }
    verdict <- .dmar_boot_bca_positions(verdict, adj)
    if (verdict$formed) limits <- stats::quantile(kept, adj, names = FALSE)
  }
  .dmar_boot_signal(verdict)
  list(limits = limits, verdict = verdict)
}
