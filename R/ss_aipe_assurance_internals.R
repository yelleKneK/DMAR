# The argument assurance_method of the planners with an assurance step
# (ss_aipe_smd(), ss_aipe_sc(), ss_aipe_sc_ancova(), ss_aipe_reg_coef(),
# ss_aipe_rc(), ss_aipe_src(), and their sensitivity functions), and the
# record an assurance plan carries of how its probability was found. None
# is exported.
#
# Two methods.
#   "analytic", the default everywhere, computes without simulation. Where
#     the probability that the realized interval is no wider than `width`
#     has an exact form (ss_aipe_smd(), ss_aipe_sc(), ss_aipe_rc(), and
#     ss_aipe_reg_coef() with noncentral = FALSE), the plan is the smallest
#     sample size, at or above the plan without assurance, whose exact
#     probability is at least `assurance`, and that probability is
#     assurance_attained. Where it has none (ss_aipe_sc_ancova(),
#     ss_aipe_reg_coef() with noncentral = TRUE, and both intervals of
#     ss_aipe_src()), the plan is the published approximation, which does
#     not compute the probability it attains, and assurance_attained is NA.
#   "simulation", offered only where the probability has no exact form,
#     estimates it from G simulated samples and returns the smallest sample
#     size, at or above the plan without assurance, whose estimate reaches
#     `assurance`; assurance_attained is that estimate.

# .ss_aipe_assurance_method(assurance_method, fn, simulation, where = "",
#                           both = TRUE)
#   assurance_method : the value given, or the formal's default.
#   fn               : the exported function the user called, for the
#                      message.
#   simulation       : TRUE when this call offers "simulation".
#   where            : the words that follow "fn()" when the offer depends
#                      on another argument (" with noncentral = FALSE").
#   both             : TRUE when the formal's default is c("analytic",
#                      "simulation"); FALSE where it is "analytic", the only
#                      method the function offers.
# Returns "analytic" or "simulation". The default of a formal that offers
# both, c("analytic", "simulation"), is "analytic"; where the formal is
# "analytic", that vector is refused as any other value is. A value is
# matched exactly, not partially, and anything other than one of the two
# strings is refused with a condition of class dmar_argument_invalid (arg =
# "assurance_method") that carries the value given. "simulation" where it is
# not offered is refused the same way, and the message says why and where
# it is offered.
.ss_aipe_assurance_method <- function(assurance_method, fn, simulation,
                                      where = "", both = TRUE) {
  choices <- c("analytic", "simulation")
  if (both && identical(assurance_method, choices)) return("analytic")
  if (!(is.character(assurance_method) && length(assurance_method) == 1L &&
        !is.na(assurance_method) && assurance_method %in% choices))
    .dmar_refuse_argument(c(
      if (simulation)
        "'assurance_method' must be \"analytic\" (the default) or \"simulation\"; "
      else c("'assurance_method' must be \"analytic\", the only method ", fn,
             "()", where, " offers; "),
      .dmar_given(assurance_method), "."),
      arg = "assurance_method", value = assurance_method)
  if (assurance_method == "simulation" && !simulation) {
    # A sensitivity function whose offer depends on noncentral offers it
    # itself with noncentral = TRUE.
    offered <- paste0("by ", c(
      if (nzchar(where) && fn != "ss_aipe_reg_coef")
        paste0(fn, "() with noncentral = TRUE"),
      "ss_aipe_sc_ancova()", "ss_aipe_reg_coef() with noncentral = TRUE",
      "ss_aipe_src()"))
    last <- length(offered)
    .dmar_refuse_argument(c(
      "assurance_method = \"simulation\" is not offered by ", fn, "()",
      where, ": the analytic probability that the interval is no wider ",
      "than 'width' is exact for this interval, so a simulated estimate ",
      "would only add Monte Carlo error. Leave 'assurance_method' at its ",
      "default, \"analytic\". Simulation is offered ",
      paste(offered[-last], collapse = ", "), ", and ", offered[last],
      ", whose probabilities have no exact form."),
      arg = "assurance_method", value = assurance_method)
  }
  assurance_method
}

# .ss_aipe_assurance_record(out, method, probability)
#   out         : an assurance plan, a dmar_tbl whose rows include
#                 assurance and assurance_attained.
#   method      : "analytic" or "simulation", the assurance_method used.
#   probability : how assurance_attained was found: "exact"; "simulation",
#                 an estimate from G simulated samples (the row G), whose
#                 Monte Carlo standard error is the attribute
#                 assurance_mc_se; or "not computed", for the published
#                 approximation, whose row assurance_attained is NA.
# Sets the attributes assurance_method and assurance_probability, which the
# printed footer reads (.ss_aipe_assurance_line()), and returns `out`.
.ss_aipe_assurance_record <- function(out, method, probability) {
  attr(out, "assurance_method") <- method
  attr(out, "assurance_probability") <- probability
  out
}

# .ss_aipe_assurance_line(x): the footer line of an assurance plan
# (.dmar_tbl_footer.default() in R/dmar_output_helpers.R), or NULL when x
# records no assurance method.
.ss_aipe_assurance_line <- function(x) {
  method <- attr(x, "assurance_method", exact = TRUE)
  probability <- attr(x, "assurance_probability", exact = TRUE)
  if (!is.character(method) || length(method) != 1L ||
      !is.character(probability) || length(probability) != 1L) return(NULL)
  switch(probability,
    exact = paste0(
      "Assurance method: analytic; assurance_attained is the exact ",
      "probability that the interval is no wider than the width"),
    simulation = {
      G <- x$value[x$term == "G"]
      se <- attr(x, "assurance_mc_se", exact = TRUE)
      paste0(
        "Assurance method: simulation; assurance_attained is estimated from ",
        if (length(G) == 1L)
          paste0("G = ", format(G, scientific = FALSE), " "),
        "simulated samples",
        if (is.numeric(se) && length(se) == 1L)
          paste0(" (Monte Carlo SE ",
                 formatC(se, format = "f",
                         digits = getOption("dmar.digits", 3L)), ")"))
    },
    `not computed` = paste0(
      "Assurance method: analytic; the plan is the published approximation, ",
      "which does not compute the probability that the interval is no wider ",
      "than the width (assurance_attained is NA) and can attain less than ",
      "the assurance. assurance_method = \"simulation\" computes it and ",
      "returns the smallest plan that attains it"),
    NULL)
}
