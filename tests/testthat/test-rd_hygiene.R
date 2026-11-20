# Documentation hygiene that R CMD check only catches with the manuals built.
#
# The 2026-08-19 CRAN pretest returned the package for a single help page
# whose equation carried Rd markup inside \eqn{}: \code{} in LaTeX math
# becomes \texttt, the PDF manual stops on it, and the HTML manual reports
# a math rendering problem. Nothing else in the toolchain objects (roxygen,
# checkRd, the text help renderer, examples, and tests all pass), so the
# defect class gets its own permanent detector here, run on every test
# pass locally and on CRAN.

rd_source_cache <- new.env(parent = emptyenv())
rd_sources <- function() {
  if (exists("db", envir = rd_source_cache, inherits = FALSE)) {
    return(rd_source_cache$db)
  }
  man <- file.path("..", "..", "man")
  if (dir.exists(man)) {
    rd_source_cache$db <- tools::Rd_db(dir = file.path("..", ".."))
  } else {
    rd_source_cache$db <- tools::Rd_db("DMAR")
  }
  rd_source_cache$db
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

# Tags present anywhere in the Rd object, each once, in the order first met.
# A tag is kept only the first time it is seen. Appending every node's tag
# and taking unique() at the end gives the same vector, but each append
# copies all the tags before it, and over the whole manual that cost about
# 0.4 s of the tests CRAN runs.
rd_tags <- function(rd) {
  tags <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && !(tag %in% tags)) tags <<- c(tags, tag)
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  tags
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

# Comment lines inside \examples whose text is R code. The 2026-09-04 CRAN
# review returned the package for exactly this ("Some code lines in examples
# are commented out. Please never do that."), after the previous round had
# adopted a plain comment as the idiom for a call shown but not run. Every
# contiguous window of comment lines is parsed with the comment markers
# stripped; a window that parses to a call or an assignment is commented-out
# code. Prose does not parse, so a comment that reads as a sentence is never
# flagged, while a bare call, an assignment, or an equation chain is.
commented_code <- function(rd) {
  ex <- examples_text(rd)
  if (!nzchar(ex)) return(character())
  lines <- strsplit(ex, "\n", fixed = TRUE)[[1]]
  cm <- grep("^\\s*#", lines)
  if (!length(cm)) return(character())
  is_code <- function(txt) {
    parsed <- tryCatch(parse(text = txt, keep.source = FALSE),
                       error = function(e) NULL)
    if (is.null(parsed) || !length(parsed)) return(FALSE)
    is_call <- vapply(parsed, is.call, logical(1))
    is_formula <- vapply(parsed, function(e)
      is.call(e) && identical(e[[1L]], as.name("~")), logical(1))
    any(is_call & !is_formula)
  }
  hits <- integer()
  for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
    txt <- sub("^\\s*#+\\s?", "", lines[run])
    n <- length(run)
    for (a in seq_len(n)) for (b in a:n) {
      window <- paste(txt[a:b], collapse = "\n")
      if (nzchar(trimws(window)) && is_code(window)) hits <- c(hits, run[a:b])
    }
  }
  lines[sort(unique(hits))]
}

test_that("the commented-code detector separates code from prose", {
  rd <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# The interval is reported for comparison.\n",
    "# ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)\n",
    "# start <- now()\n",
    "# Sys.sleep(0.2)\n",
    "ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)\n}")))
  hits <- commented_code(rd)
  expect_length(hits, 3L)
  expect_false(any(grepl("reported for comparison", hits)))
})

test_that("no help page carries commented-out code in its examples", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- Filter(length, lapply(db, commented_code))
  expect_identical(names(offenders), character(0),
                   info = paste(unlist(lapply(names(offenders), function(f)
                     paste(f, offenders[[f]], sep = ": "))), collapse = "\n"))
})

# The same CRAN review named three run-time policies: no writing to the
# global environment (which its checker reads as any assign() or rm() aimed
# at .GlobalEnv, and any <<-), no options(warn = -1), and no default path in
# a function that writes a file. The package's own functions are the
# evidence, so the namespace is deparsed and inspected here. The option is
# `warn` itself, followed by "=": `warning.length`, which the refusal of a
# fit with no maximum raises while it prints its lavaan code
# (.no_max_refuse()), is another option.
namespace_functions <- function() {
  ns <- asNamespace("DMAR")
  nms <- ls(ns, all.names = TRUE)
  fns <- Filter(is.function, mget(nms, envir = ns, inherits = FALSE))
  fns[!vapply(fns, is.primitive, logical(1))]
}

test_that("no function touches the global environment, uses <<-, or sets options(warn)", {
  fns <- namespace_functions()
  text <- vapply(fns, function(f) paste(deparse(f, width.cutoff = 500L), collapse = "\n"),
                 character(1))
  offenders <- names(text)[grepl("\\.GlobalEnv|globalenv\\(\\)|\\.Random\\.seed|<<-|options\\(\\s*warn\\s*=", text)]
  expect_identical(offenders, character(0))
})

test_that("no function that writes a file has a default path", {
  fns <- namespace_functions()
  bad <- character()
  for (nm in names(fns)) {
    fm <- formals(fns[[nm]])
    for (arg in intersect(names(fm), c("file", "filename", "path", "dir", "directory"))) {
      if (!is.null(fm[[arg]]) && !identical(fm[[arg]], quote(expr = ))) {
        bad <- c(bad, sprintf("%s(%s = %s)", nm, arg, deparse(fm[[arg]])))
      }
    }
    if ("save" %in% names(fm) && grepl("_sensitivity$", nm)) {
      bad <- c(bad, sprintf("%s still has a 'save' argument", nm))
    }
  }
  expect_identical(bad, character(0))
})

# A spaced double hyphen in Rd prose renders as an en dash used as a
# dash, the em dash style the package does not use; in an example comment
# it reads the same way. A double hyphen with no space beside it, between
# numbers (123--124, 1997--1999, \code{1}--\code{6}) or names
# (Bryant--Paulson), is a range or a joined pair and is left alone, as are
# the reference lists, whose titles are reproduced as published, and code
# and math.
prose_dashes <- function(rd) {
  dash <- "\\s--(\\s|$)|(^|\\s)--\\s"
  skip <- c("\\references", "\\examples", "\\eqn", "\\deqn", "\\code",
            "\\preformatted", "\\verb", "\\href", "\\url")
  out <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% skip) return(invisible())
    if (is.list(x)) {
      for (el in x) walk(el)
    } else if (is.character(x) && any(grepl(dash, x))) {
      out <<- c(out, trimws(x[grepl(dash, x)]))
    }
  }
  # One vectorized pass over every text node of the page first; only a
  # page with a spaced double hyphen somewhere (most have none) is walked
  # node by node to set aside the sections where one is allowed.
  if (any(grepl(dash, unlist(rd, use.names = FALSE)))) walk(rd)
  ex <- strsplit(examples_text(rd), "\n", fixed = TRUE)[[1]]
  comments <- ex[grepl("^\\s*#", ex)]
  c(out, trimws(comments[grepl(dash, comments)]))
}

test_that("the prose dash detector separates a dash from a range", {
  rd <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}",
    "\\description{One of three tests -- the main effect -- is chosen.\n",
    "Pages 123--124 cover the Bryant--Paulson range, levels ",
    "\\code{1}--\\code{6}.}",
    "\\references{Communications in Statistics -- Simulation.}",
    "\\examples{\n# A null of 0.2 -- looking for more\n# rows 1--2\nx <- 1\n}")))
  hits <- prose_dashes(rd)
  expect_length(hits, 2L)
  expect_match(hits[1L], "^One of three tests")
  expect_match(hits[2L], "^# A null of 0.2")
})

test_that("no help page uses a double hyphen as a dash", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- Filter(length, lapply(db, prose_dashes))
  expect_identical(names(offenders), character(0),
                   info = paste(unlist(lapply(names(offenders), function(f)
                     paste(f, offenders[[f]], sep = ": "))), collapse = "\n"))
})
