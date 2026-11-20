# car and lme4 both register an S3 method for na.action on merMod objects,
# so whichever namespace loads second emits a "Registered S3 method
# overwritten" message. Load both here, car first, with the registration
# message muffled so it does not leak into the output of whichever test file
# happens to load lme4 first. The curated subset that runs on CRAN (NOT_CRAN
# unset; cran_files in tests/testthat.R) never loads lme4, so there is no
# message to muffle, and the half second that loading lme4 takes is not
# spent on CRAN's clock.
if (isTRUE(as.logical(Sys.getenv("NOT_CRAN")))) {
  suppressMessages({
    requireNamespace("car", quietly = TRUE)
    requireNamespace("lme4", quietly = TRUE)
  })
}
