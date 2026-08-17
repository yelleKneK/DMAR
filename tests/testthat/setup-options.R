# Fail-loudly options for the whole test run.
#
# Partial matching of argument names, list elements ($), and attributes is a
# silent correctness hazard: code that works today breaks when a function later
# gains an argument sharing the abbreviated prefix. Enabling the warnings here
# means any internal reliance on partial matching surfaces as a test warning
# the day it is introduced rather than as a user-facing bug years later. The
# options are process-local to the testthat run and do not affect users.
options(
  warnPartialMatchArgs   = TRUE,
  warnPartialMatchDollar = TRUE,
  warnPartialMatchAttr   = TRUE
)
