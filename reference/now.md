# A Friendly Date / Time Stamp and a Simple Stopwatch

`now()` is a small utility for two related tasks: (1) printing the
current date and time in a form that reads naturally in a console log or
report ("September 21, 2026 (3:42 PM)"), and (2) measuring how long a
piece of work takes by capturing the time before and after the work and
subtracting the two stamps.

## Usage

``` r
now(time = TRUE, tidy = FALSE)
```

## Arguments

- time:

  Logical; if `TRUE` (default) the hour, minute, and AM/PM are appended
  to the printed date. `FALSE` returns the date only.

- tidy:

  Logical; if `TRUE`, the timestamp is returned as a `data.frame` with
  columns `term` and `value` (numeric components only: day, year, hour,
  minute) and the non-numeric components (month name, AM/PM) attached as
  attributes. Defaults to `FALSE`, which returns the printable
  `dmar_now` stamp. Use `tidy = TRUE` only when the components are
  needed individually for programmatic processing; the `dmar_now`
  default supports subtraction and prettily prints in any context.

## Value

When `tidy = FALSE` (default), an object of class
`c("dmar_now", "POSIXct", "POSIXt")` whose print method yields the
natural-language form ("September 21, 2026 (3:42 PM)" or, with
`time = FALSE`, "September 21, 2026"). The underlying numeric value is
the [`Sys.time()`](https://rdrr.io/r/base/Sys.time.html) stamp at the
moment of the call, so the `-` operator gives the elapsed time between
two captures as a [`difftime`](https://rdrr.io/r/base/difftime.html)
object.

When `tidy = TRUE`, a `data.frame` with columns `term` and `value`,
where `value` is numeric (day, year, hour, minute) and the month name
and AM/PM marker are attached as the `"month"` and `"am_pm"` attributes.

## Details

The function returns an object of class `"dmar_now"` that carries the
underlying [`Sys.time`](https://rdrr.io/r/base/Sys.time.html) value and
prints in the human-readable form. Two `dmar_now` objects can be
subtracted with the ordinary `-` operator; the result is a
[`difftime`](https://rdrr.io/r/base/difftime.html) object with units
chosen automatically from the magnitude of the elapsed interval. The
pattern is the R analog of a stopwatch: capture before, capture after,
take the difference.

For analyses that should record *when* a long-running computation was
performed (a bootstrap confidence interval, a Monte Carlo simulation, an
`ss_aipe_*_sensitivity` run), inserting a `now()` call at the start and
end of the work creates a self-documenting log of when the computation
ran and how long it took. The print method's natural-language form is
designed to be copy-pasted into a methods section or an analysis note
without further formatting.

The name coincides with
[`lubridate::now()`](https://lubridate.tidyverse.org/reference/now.html),
and the masking is deliberate: both functions return the current time as
a `POSIXct` object, so a script written for either remains correct with
the other attached; only the printed form differs.

## See also

[`Sys.time`](https://rdrr.io/r/base/Sys.time.html) for the underlying
timestamp; [`difftime`](https://rdrr.io/r/base/difftime.html) for the
elapsed-time class returned by the `-` operator.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Print the current date and time.
now()
#> August 17, 2026 (12:51 PM)

# Time how long a piece of work takes. The pattern is the same
# whether the work is a bootstrap, a simulation, or a numeric
# search: capture a stamp before, capture one after, subtract.
# The "-" operator returns the elapsed time as a difftime, with
# units chosen automatically. Here the work is a descriptive
# summary of three Holzinger and Swineford cognitive tests, small
# enough that the elapsed time is a fraction of a second.
start <- now()
d <- descriptives(holzinger_swineford[, c("t1_visual_perception",
                                          "t2_cubes", "t4_lozenges")])
end <- now()
end - start
#> Time difference of 0.003771305 secs

# A deliberate wait shows the same pattern on a longer interval.
# Not run here because the only way to demonstrate a wait is to
# make the example wait; the calls are:
# start <- now()
# Sys.sleep(0.5)
# end <- now()
# end - start          # elapsed time, here about 0.5 seconds

# Date only.
now(time = FALSE)
#> August 17, 2026

# Tidy data.frame form, when the components are needed
# individually for programmatic processing (e.g., embedding in a
# report's metadata block).
now(tidy = TRUE)
#>     term value
#> 1    day    17
#> 2   year  2026
#> 3   hour    12
#> 4 minute    51
```
