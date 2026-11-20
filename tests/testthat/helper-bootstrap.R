# local_small_B()
#
# Lets a test run a bootstrap with a small B to exercise its mechanics
# (counts, attributes, progress, cancellation, agreement with a reference
# computed on the same draws). It replaces .dmar_boot_required(), the
# fewest replications a request may ask for (R/dmar_bootstrap.R),
# by a function that returns 2, for the calling frame only: a test_that()
# block, or a whole test file when called at file scope, where testthat
# undoes it when the file ends. It also sets to 0 the count of
# replications required beyond each limit of a BCa interval where the
# bias correction and the acceleration place it (.dmar_boot_beyond()),
# which no small B meets. The share rule (no interval when more than half
# of the replications were dropped) is not affected. Tests of the rule
# itself run without this.
local_small_B <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    .dmar_boot_required = function(kind = "quantile", tail = NULL,
                                   floor = 100L) 2L,
    .dmar_boot_beyond = function() 0,
    .package = "DMAR", .env = env)
}


# local_count_unmet()
#
# Lets a test see what a function returns and prints when the count rule
# leaves no interval. The rule is held to the number requested, which
# every exported function checks before anything is computed, so a
# request that fails it never reaches the replications in a call a user
# can make; the verdict "insufficient_count" then remains for a BCa
# interval withheld at its adjusted positions and for an internal caller.
# This replaces .dmar_boot_count_met() (R/dmar_bootstrap.R) by a function
# that returns FALSE, for the calling frame only, so that a run whose
# request was accepted is settled as one whose request was below the
# requirement.
local_count_unmet <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    .dmar_boot_count_met = function(B, required) FALSE,
    .package = "DMAR", .env = env)
}
