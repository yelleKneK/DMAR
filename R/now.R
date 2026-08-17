#' A Friendly Date / Time Stamp and a Simple Stopwatch
#'
#' \code{now()} is a small utility for two related tasks: (1)
#' printing the current date and time in a form that reads
#' naturally in a console log or report ("September 21, 2026
#' (3:42 PM)"), and (2) measuring how long a piece of work takes
#' by capturing the time before and after the work and subtracting
#' the two stamps.
#'
#' The function returns an object of class \code{"dmar_now"} that
#' carries the underlying \code{\link{Sys.time}} value and prints
#' in the human-readable form. Two \code{dmar_now} objects can be
#' subtracted with the ordinary \code{-} operator; the result is a
#' \code{\link{difftime}} object with units chosen automatically
#' from the magnitude of the elapsed interval. The pattern is the
#' R analog of a stopwatch: capture before, capture after, take the
#' difference.
#'
#' For analyses that should record \emph{when} a long-running
#' computation was performed (a bootstrap confidence interval, a
#' Monte Carlo simulation, an \code{ss_aipe_*_sensitivity} run),
#' inserting a \code{now()} call at the start and end of the work
#' creates a self-documenting log of when the computation ran and
#' how long it took. The print method's natural-language form is
#' designed to be copy-pasted into a methods section or an analysis
#' note without further formatting.
#'
#' The name coincides with \code{lubridate::now()}, and the masking
#' is deliberate: both functions return the current time as a
#' \code{POSIXct} object, so a script written for either remains
#' correct with the other attached; only the printed form differs.
#'
#' @param time Logical; if \code{TRUE} (default) the hour, minute,
#'   and AM/PM are appended to the printed date. \code{FALSE}
#'   returns the date only.
#' @param tidy Logical; if \code{TRUE}, the timestamp is returned
#'   as a \code{data.frame} with columns \code{term} and
#'   \code{value} (numeric components only: day, year, hour,
#'   minute) and the non-numeric components (month name, AM/PM)
#'   attached as attributes. Defaults to \code{FALSE}, which
#'   returns the printable \code{dmar_now} stamp. Use \code{tidy =
#'   TRUE} only when the components are needed individually for
#'   programmatic processing; the \code{dmar_now} default supports
#'   subtraction and prettily prints in any context.
#'
#' @return When \code{tidy = FALSE} (default), an object of class
#'   \code{c("dmar_now", "POSIXct", "POSIXt")} whose print method
#'   yields the natural-language form ("September 21, 2026 (3:42
#'   PM)" or, with \code{time = FALSE}, "September 21, 2026"). The
#'   underlying numeric value is the \code{Sys.time()} stamp at the
#'   moment of the call, so the \code{-} operator gives the elapsed
#'   time between two captures as a \code{\link{difftime}} object.
#'
#'   When \code{tidy = TRUE}, a \code{data.frame} with columns
#'   \code{term} and \code{value}, where \code{value} is numeric
#'   (day, year, hour, minute) and the month name and AM/PM marker
#'   are attached as the \code{"month"} and \code{"am_pm"}
#'   attributes.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{Sys.time}} for the underlying timestamp;
#'   \code{\link[base]{difftime}} for the elapsed-time class
#'   returned by the \code{-} operator.
#'
#' @examples
#' # Print the current date and time.
#' now()
#'
#' # Time how long a piece of work takes. The pattern is the same
#' # whether the work is a bootstrap, a simulation, or a numeric
#' # search: capture a stamp before, capture one after, subtract.
#' # The "-" operator returns the elapsed time as a difftime, with
#' # units chosen automatically. Here the work is a descriptive
#' # summary of three Holzinger and Swineford cognitive tests, small
#' # enough that the elapsed time is a fraction of a second.
#' start <- now()
#' d <- descriptives(holzinger_swineford[, c("t1_visual_perception",
#'                                           "t2_cubes", "t4_lozenges")])
#' end <- now()
#' end - start
#'
#' # A deliberate wait shows the same pattern on a longer interval.
#' # Not run here because the only way to demonstrate a wait is to
#' # make the example wait; the calls are:
#' # start <- now()
#' # Sys.sleep(0.5)
#' # end <- now()
#' # end - start          # elapsed time, here about 0.5 seconds
#'
#' # Date only.
#' now(time = FALSE)
#'
#' # Tidy data.frame form, when the components are needed
#' # individually for programmatic processing (e.g., embedding in a
#' # report's metadata block).
#' now(tidy = TRUE)
#'
#' @export
now <- function(time = TRUE, tidy = FALSE) {
  present <- Sys.time()
  if (isTRUE(tidy)) {
    return(.now_tidy_frame(present, time = time))
  }
  out <- present
  attr(out, "print_time") <- isTRUE(time)
  class(out) <- c("dmar_now", class(present))
  out
}


#' @export
print.dmar_now <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}


#' @export
format.dmar_now <- function(x, ...) {
  month  <- format(unclass_now(x), "%B")
  day    <- as.integer(format(unclass_now(x), "%d"))
  year   <- as.integer(format(unclass_now(x), "%Y"))
  hour   <- as.integer(format(unclass_now(x), "%I"))
  minute <- format(unclass_now(x), "%M")
  am_pm  <- toupper(paste0(
    substr(format(unclass_now(x), "%P"), 1, 1), "M"
  ))
  mdy  <- paste0(month, " ", day, ", ", year)
  if (isTRUE(attr(x, "print_time"))) {
    paste0(mdy, " (", hour, ":", minute, " ", am_pm, ")")
  } else {
    mdy
  }
}


#' @export
"-.dmar_now" <- function(e1, e2) {
  # Both operands are dmar_now objects (the supported stopwatch
  # pattern: capture before, capture after, subtract). Return the
  # elapsed time as a difftime with units chosen by R's default
  # heuristic (seconds, minutes, hours, days, weeks).
  t1 <- unclass_now(e1)
  t2 <- unclass_now(e2)
  out <- difftime(t1, t2)
  out
}


# Internal: drop the dmar_now class so base R arithmetic sees the
# POSIXct value rather than dispatching back to our method.
unclass_now <- function(x) {
  cls <- class(x)
  class(x) <- setdiff(cls, "dmar_now")
  x
}


# Internal: build the data.frame return for tidy = TRUE.
.now_tidy_frame <- function(present, time) {
  month  <- format(present, "%B")
  day    <- as.numeric(format(present, "%d"))
  year   <- as.numeric(format(present, "%Y"))
  hour   <- as.numeric(format(present, "%I"))
  minute <- as.numeric(format(present, "%M"))
  am_pm  <- toupper(paste0(
    substr(format(present, "%P"), 1, 1), "M"
  ))
  if (isTRUE(time)) {
    out <- data.frame(
      term  = c("day", "year", "hour", "minute"),
      value = c(day, year, hour, minute)
    )
    attr(out, "month") <- month
    attr(out, "am_pm") <- am_pm
  } else {
    out <- data.frame(
      term  = c("day", "year"),
      value = c(day, year)
    )
    attr(out, "month") <- month
  }
  out
}
