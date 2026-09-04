# Seed discipline for every function that takes a `seed` argument.
#
# A supplied seed is set for the duration of the calling function only. The
# caller's random number generator state, including the generator kind, is
# restored when that function exits, so a call made reproducible by `seed`
# does not perturb the random stream of the script around it. The saving and
# restoring is delegated to withr::local_seed(), the accepted way to do this
# on CRAN; the package itself never reads or writes `.Random.seed`. A NULL
# seed is a no-op: the draws come from the user's current generator state,
# exactly as if the function had no `seed` argument. Call it from the body of
# the function whose exit should restore the state (the default `envir` is
# that function's frame), never from a helper called by that function unless
# the helper is where the random draws begin and end.
.dmar_local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be a single number or NULL.", call. = FALSE)
  }
  withr::local_seed(as.integer(seed), .local_envir = envir)
  invisible(seed)
}
