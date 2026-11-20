## data-raw/bessel_errors.R
##
## Construction script for the bessel_errors data set.
##
## Source: Friedrich Wilhelm Bessel (1818), "Fundamenta
## astronomiae pro anno MDCCLV deducta ex observationibus viri
## incomparabilis James Bradley in specula astronomica
## Grenovicensi per annos 1750--1762 institutis." Friedrich
## Nicolovius. The 9-bin grouped frequency distribution of the
## absolute errors that Bessel computed from 300 of British
## Astronomer Royal James Bradley's stellar observations.
## Reproduced in Maxwell, Delaney, and Kelley (2027, "Designing
## Experiments and Analyzing Data: A Model Comparison
## Perspective," 4th ed., Routledge), Table 1.4.
##
## The numbers below are taken from MDK4, Table 1.4. Both the
## observed and expected columns sum to 300, matching the
## published table.
##
## To rebuild from scratch, run from the package root:
##   source("data-raw/bessel_errors.R")

bessel_errors <- data.frame(
  bin       = 1:9,
  lower     = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
  upper     = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  observed  = c(114L, 84L, 53L, 24L, 14L, 6L, 3L, 1L, 1L),
  expected  = c(107L, 87L, 57L, 30L, 13L, 5L, 1L, 0L, 0L),
  stringsAsFactors = FALSE
)
# Compute midpoint from the bin edges to avoid IEEE 754 drift from
# literal 0.05 / 0.15 / ... values that look identical in print but
# do not satisfy bit-exact equality with (lower + upper) / 2.
bessel_errors$midpoint <- (bessel_errors$lower + bessel_errors$upper) / 2

# Reorder columns to keep the midpoint between upper and observed,
# matching the documented layout in the Rd \format{} block.
bessel_errors <- bessel_errors[, c("bin", "lower", "upper",
                                   "midpoint", "observed",
                                   "expected")]

# Reconciliation checks against MDK4 Table 1.4.
stopifnot(
  identical(dim(bessel_errors), c(9L, 6L)),
  identical(sum(bessel_errors$observed), 300L),
  identical(sum(bessel_errors$expected), 300L),
  all.equal(bessel_errors$midpoint,
            (bessel_errors$lower + bessel_errors$upper) / 2)
)

save(bessel_errors,
     file = "data/bessel_errors.rda",
     version = 2L, compress = "bzip2")
