# .match_which_width(which_width)
#
# The `which_width` of ss_aipe_R2() and ss_aipe_reg_coef() (and through it
# ss_aipe_rc() and ss_aipe_src()): "Full", "Lower", or "Upper", or an
# abbreviation of one ("F", "L", "Up"), returned in full, as match.arg()
# returns it for the planners whose which_width lists the three values.
# Anything else (another string, "", NA, a value that is not a character
# string, or more than one value) is refused with a condition of class
# dmar_argument_invalid (arg = "which_width"). DMAR 1.0.x checked the
# value with char.expand() but kept it as given, so ss_aipe_reg_coef()
# planned "L" and "Up" for the full width and ss_aipe_R2() returned NULL
# for "F", "L", and "Up".
.match_which_width <- function(which_width) {
  choices <- c("Full", "Lower", "Upper")
  if (is.character(which_width) && length(which_width) == 1L &&
      !is.na(which_width) && nzchar(which_width)) {
    m <- charmatch(which_width, choices)
    if (!is.na(m) && m > 0L) return(choices[m])
  }
  given <- if (is.character(which_width) && length(which_width) == 1L) {
    if (is.na(which_width)) "NA was given" else
      paste0("\"", which_width, "\" was given")
  } else {
    paste0("an object of class ", class(which_width)[1L], " and length ",
           length(which_width), " was given")
  }
  .dmar_refuse_argument(c(
    "'which_width' must be \"Full\", \"Lower\", or \"Upper\", or an ",
    "abbreviation of one such as \"L\"; ", given, "."),
    arg = "which_width", value = which_width)
}
