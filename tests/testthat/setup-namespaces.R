# car and lme4 both register an S3 method for na.action on merMod objects,
# so whichever namespace loads second emits a "Registered S3 method
# overwritten" message. Load both here, car first, with the registration
# message muffled so it does not leak into the output of whichever test file
# happens to load lme4 first.
suppressMessages({
  requireNamespace("car", quietly = TRUE)
  requireNamespace("lme4", quietly = TRUE)
})
