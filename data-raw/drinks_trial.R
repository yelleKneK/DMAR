## data-raw/drinks_trial.R
##
## Construction script for the drinks_trial data set.
##
## Source: Smith, Meyers, and Delaney (1998), "The community
## reinforcement approach with homeless alcohol-dependent
## individuals," Journal of Consulting and Clinical Psychology,
## 66(3), 541-548. The raw nine-month follow-up values used here
## are also distributed via the AMCP companion package for
## Maxwell, Delaney, and Kelley (2027) as chapter_3_table_7_raw.
##
## What this script does:
##   1. Pulls the 88 raw nine-month drinks-per-week values from
##      AMCP's chapter_3_table_7_raw.
##   2. Decodes the AMCP integer group code into a substantive
##      (cohort, treatment) factor pair, with treatment labels
##      rewritten for clarity ("Standard", "CRA",
##      "CRA + Disulfiram") rather than the AMCP CRA-D / CRA+D
##      shorthand that can be misread.
##   3. Drops the AMCP integer Group column. The (cohort,
##      treatment) factor pair carries all the design information
##      and is easier to use in formulas.
##   4. Sanity-checks per-cell counts against the published paper.
##
## AMCP is only needed at build time and is not a runtime
## dependency of DMAR; the .rda is self-contained.
##
## To rebuild from scratch, run from the package root:
##   source("data-raw/drinks_trial.R")

stopifnot(requireNamespace("AMCP", quietly = TRUE))

data("chapter_3_table_7_raw", package = "AMCP", envir = environment())
raw <- chapter_3_table_7_raw

stopifnot(nrow(raw) == 88L, ncol(raw) == 3L,
          identical(names(raw), c("Group", "Drinks", "LgDrinks")))

# AMCP group coding (from the AMCP help page):
#   1 = Cohort 1 CRA-D    (CRA without disulfiram)
#   2 = Cohort 1 CRA+D    (CRA with disulfiram)
#   3 = Cohort 1 Standard
#   4 = Cohort 2 CRA-D
#   5 = Cohort 2 Standard
# Cohort 2 dropped the CRA + Disulfiram cell after Cohort 1
# results suggested disulfiram added little to CRA alone.
cohort_lookup    <- c("1", "1", "1", "2", "2")
treatment_lookup <- c("CRA", "CRA + Disulfiram", "Standard",
                      "CRA", "Standard")
cohort    <- factor(cohort_lookup[raw$Group],
                    levels = c("1", "2"))
treatment <- factor(treatment_lookup[raw$Group],
                    levels = c("Standard", "CRA", "CRA + Disulfiram"))

drinks_trial <- data.frame(
  id              = seq_len(nrow(raw)),
  cohort          = cohort,
  treatment       = treatment,
  drinks_per_week = raw$Drinks,
  log_drinks      = raw$LgDrinks,
  stringsAsFactors = FALSE
)

# Reconciliation checks against the published study.
stopifnot(
  identical(dim(drinks_trial), c(88L, 5L)),
  identical(sum(drinks_trial$cohort == "1"), 51L),
  identical(sum(drinks_trial$cohort == "2"), 37L),
  identical(sum(drinks_trial$treatment == "Standard"),         37L),
  identical(sum(drinks_trial$treatment == "CRA"),              32L),
  identical(sum(drinks_trial$treatment == "CRA + Disulfiram"), 19L),
  # Per-cell sizes, cohort by treatment.
  identical(sum(drinks_trial$cohort == "1" &
                drinks_trial$treatment == "CRA"),              15L),
  identical(sum(drinks_trial$cohort == "1" &
                drinks_trial$treatment == "CRA + Disulfiram"), 19L),
  identical(sum(drinks_trial$cohort == "1" &
                drinks_trial$treatment == "Standard"),         17L),
  identical(sum(drinks_trial$cohort == "2" &
                drinks_trial$treatment == "CRA"),              17L),
  identical(sum(drinks_trial$cohort == "2" &
                drinks_trial$treatment == "Standard"),         20L),
  # Outcome range.
  min(drinks_trial$drinks_per_week)  == 0,
  abs(max(drinks_trial$drinks_per_week) - 624.615) < 0.001,
  # log_drinks is log10(drinks + 1); the zeros map to zeros.
  abs(min(drinks_trial$log_drinks)) < 1e-9,
  abs(max(drinks_trial$log_drinks) - 2.7963) < 0.001
)

save(drinks_trial,
     file = "data/drinks_trial.rda",
     version = 2L, compress = "bzip2")
