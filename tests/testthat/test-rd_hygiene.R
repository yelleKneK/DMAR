# Documentation hygiene that R CMD check only catches with the manuals built.
#
# The 2026-08-19 CRAN pretest returned the package for a single help page
# whose equation carried Rd markup inside \eqn{}: \code{} in LaTeX math
# becomes \texttt, the PDF manual stops on it, and the HTML manual reports
# a math rendering problem. Nothing else in the toolchain objects (roxygen,
# checkRd, the text help renderer, examples, and tests all pass), so the
# defect class gets its own permanent detector here, run on every test
# pass locally and on CRAN.

rd_sources <- function() {
  man <- file.path("..", "..", "man")
  if (dir.exists(man)) {
    tools::Rd_db(dir = file.path("..", ".."))
  } else {
    tools::Rd_db("DMAR")
  }
}

# Text of every \eqn{} / \deqn{} node (all arguments), keyed by Rd file.
math_text <- function(rd) {
  out <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% c("\\eqn", "\\deqn")) {
      out <<- c(out, paste(unlist(x), collapse = ""))
    }
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  out
}

# Tags present anywhere in the Rd object.
rd_tags <- function(rd) {
  tags <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag)) tags <<- c(tags, tag)
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  unique(tags)
}

# Text of the \examples section.
examples_text <- function(rd) {
  ex <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\examples"), rd)
  if (!length(ex)) return("")
  paste(unlist(ex), collapse = "")
}

markup_in_math <- "\\\\(code|emph|link|strong|pkg|bold|verb|var|samp|file|env|dQuote|sQuote|href|url)\\{"

test_that("the math-markup detector flags the construct CRAN rejected", {
  bad <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\code{delta_lower}, \\code{delta_upper})}.}")))
  expect_true(any(grepl(markup_in_math, math_text(bad))))
  good <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\delta_L, +\\delta_U)}, with \\code{delta_lower} as ",
    "\\eqn{\\delta_L}.}")))
  expect_false(any(grepl(markup_in_math, math_text(good))))
})

test_that("no help page carries Rd markup inside \\eqn{} or \\deqn{}", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) any(grepl(markup_in_math, math_text(rd))), db))
  expect_identical(offenders, character(0))
})

test_that("no help page uses \\donttest{} or \\dontrun{}", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) any(c("\\donttest", "\\dontrun") %in% rd_tags(rd)), db))
  expect_identical(offenders, character(0))
})

test_that("no example block gates on requireNamespace()", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) grepl("requireNamespace", examples_text(rd), fixed = TRUE), db))
  expect_identical(offenders, character(0))
})
