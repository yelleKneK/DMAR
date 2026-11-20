# Check dependent-package contracts against the recorded snapshot.
#
# This reads dependent R sources and NAMESPACE entries, compares referenced
# DMAR objects and formals, and executes the snapshot's canonical calls to
# check result columns, terms, classes and consumed attributes. Dynamically
# constructed names must be declared in the snapshot.
#
# Usage, from the package root:
#   Rscript tools/revdep_contract.R
#   Rscript tools/revdep_contract.R --lib=<isolated-library>
#   Rscript tools/revdep_contract.R --update
# Options:
#   --siever=PATH     dependent checkout (default: ../SIEVER)
#   --smsd=PATH       dependent checkout (default: ../SMSD)
#   --snapshot=FILE   contract snapshot (default: tools/revdep_contract.json)
#   --pkg=DIR         DMAR source checkout (default: .)
#
# Exit status: 0 for a complete match; 1 for a contract difference; 2 for
# usage, load or snapshot errors; 3 when a dependent could not be read.
# Missing dependents do not count as a pass. --update refuses incomplete or
# inconsistent inventories and records portable API data, not checkout paths.
# Review contract changes and coordinate them with consumers before release.
# The full release gate checks this against its installed artifact and a
# planted incompatible formal. tools/ is excluded from the source package.

args <- commandArgs(trailingOnly = TRUE)
known <- c("--update", "--lib", "--pkg", "--snapshot", "--siever", "--smsd")
usage <- paste("usage: Rscript tools/revdep_contract.R [--update] [--lib=DIR] [--pkg=DIR]",
               "[--snapshot=FILE] [--siever=PATH] [--smsd=PATH]")
bad <- args[!sub("=.*", "", args) %in% known]
if (length(bad)) { cat("unknown argument(s):", bad, "\n", usage, "\n"); quit(status = 2L) }
flag <- function(name, default) {
  v <- sub(paste0("^", name, "="), "", grep(paste0("^", name, "="), args, value = TRUE))
  if (length(v)) v[length(v)] else default
}
update <- "--update" %in% args
lib <- flag("--lib", "")
pkg <- flag("--pkg", ".")
snapshot_file <- flag("--snapshot", if (pkg == ".") file.path("tools", "revdep_contract.json")
                                    else file.path(pkg, "tools", "revdep_contract.json"))
dependents <- c(SIEVER = flag("--siever", file.path(normalizePath("..", mustWork = TRUE), "SIEVER")),
                SMSD = flag("--smsd", file.path(normalizePath("..", mustWork = TRUE), "SMSD")))
if (!requireNamespace("jsonlite", quietly = TRUE)) { cat("jsonlite is needed to read the snapshot\n"); quit(status = 2L) }

## ---- The canonical calls and declared names, used when there is no snapshot ----
# Each call is evaluated with DMAR loaded, in a fresh environment, after
# set.seed(113). "reads" lists what a dependent reads from the result: rows
# of the term column ("terms"), columns ("columns"), and attributes.
default_calls <- list(
  list(label = "mr_cv, verbose", used_by = "SMSD",
       call = paste("DMAR::mr_cv(data = c(10, 12, 9, 15, 11, 14, 8, 13, 12, 16, 10, 11, 14, 9, 12, 13, 15, 11),",
                    "A = 400000, sampling_cost = 75, gamma = 0.25, verbose = TRUE)"),
       reads = list(terms = c("cv", "risk", "v2_n", "mu_3_n", "mu_4_n"),
                    attributes = c("criterion", "is_satisfied", "variance_unusable"))),
  list(label = "mr_cv, pilot", used_by = "SMSD",
       call = "DMAR::mr_cv(pilot = TRUE, A = 1e5, sampling_cost = 1, m0 = 4, gamma = 0.25)",
       reads = list(terms = "pilot_N")),
  list(label = "smd, unbiased = FALSE", used_by = "SMSD",
       call = paste("DMAR::smd(group_1 = c(48.1, 52.3, 45.9, 50.2, 55.0, 47.7, 49.8, 51.4),",
                    "group_2 = c(53.2, 57.9, 50.1, 56.4, 54.8, 59.3, 52.0, 55.6), unbiased = FALSE)"),
       reads = list(terms = "smd")),
  list(label = "mr_smd, at a sample size", used_by = "SMSD",
       call = "DMAR::mr_smd(A = 10000, d = 0.5, n_per_group = 50, sampling_cost = 2, gamma = 0.25)",
       reads = list(terms = "risk", attributes = c("criterion", "is_satisfied"))),
  list(label = "mr_smd, pilot", used_by = "SMSD",
       call = "DMAR::mr_smd(pilot = TRUE, A = 1e6, sampling_cost = 5, m0 = 4, gamma = 0.25)",
       reads = list(terms = "pilot_n_per_group")),
  list(label = "ci_crossed", used_by = "SMSD",
       call = paste("DMAR::ci_crossed(matrix(c(7, 9, 6, 8, 10, 6, 8, 5, 9, 9, 8, 9, 7, 8, 11, 5, 7, 6, 7, 9,",
                    "6, 9, 7, 9, 10, 7, 8, 6, 8, 10, 8, 10, 8, 9, 12, 6, 8, 7, 8, 10, 7, 8, 7, 9, 11, 6, 9,",
                    "6, 8, 9, 9, 10, 8, 9, 11, 7, 8, 7, 9, 10), nrow = 12, ncol = 5), conf_level = 0.95)"),
       reads = list(terms = c("ci_lower", "grand_mean", "ci_upper", "se", "df", "var_unit", "var_rater",
                              "var_residual"))),
  list(label = "ci_proportion", used_by = "SIEVER",
       call = "DMAR::ci_proportion(successes = 50L, N = 100L)",
       reads = list(terms = c("ci_lower", "proportion", "ci_upper"))),
  list(label = "irt_information", used_by = "SIEVER",
       call = "DMAR::irt_information(discrimination = 1, threshold = 0, trait_value = c(-1, 0, 1))",
       reads = list(columns = c("trait_value", "test_information", "se"),
                    attributes = c("discrimination", "trait_value_max_information", "item_information"))),
  list(label = "reliability_alpha", used_by = "SIEVER",
       call = paste("DMAR::reliability_alpha(data = DMAR::holzinger_swineford[, c(\"t1_visual_perception\",",
                    "\"t2_cubes\", \"t3_paper_form_board\", \"t4_lozenges\")])"),
       reads = list(terms = c("estimate", "se", "ci_lower", "ci_upper"))),
  list(label = "sem_fit_measures, from summary statistics", used_by = "SIEVER",
       call = "DMAR::sem_fit_measures(chisq = 85.3, df = 24, N = 301, baseline_chisq = 918.9, baseline_df = 36)",
       reads = list())
)
default_declared <- list(
  SIEVER = list(
    missing_patterns = "named at run time: paste0(\"missing_\", kind) in R/fct_missing_data.R, after a test against getNamespaceExports(\"DMAR\")",
    missing_coverage = "named at run time: paste0(\"missing_\", kind) in R/fct_missing_data.R, after a test against getNamespaceExports(\"DMAR\")"),
  SMSD = list()
)

## ---- Load DMAR -------------------------------------------------------------------
loaded <- tryCatch({
  if (nzchar(lib)) {
    lib <- normalizePath(lib, mustWork = TRUE)
    ns <- loadNamespace("DMAR", lib.loc = lib)
    path <- normalizePath(getNamespaceInfo(ns, "path"))
    if (!startsWith(path, lib)) stop("DMAR was loaded from ", path, ", not from ", lib)
    list(ns = ns, from = sprintf("DMAR %s installed in %s", getNamespaceVersion(ns), lib), commit = NA_character_,
         origin = sprintf("DMAR %s, installed", getNamespaceVersion(ns)))
  } else {
    pkg <- normalizePath(pkg, mustWork = TRUE)
    suppressMessages(pkgload::load_all(pkg, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE,
                                       quiet = TRUE))
    ns <- asNamespace("DMAR")
    commit <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(pkg), "rev-parse", "--short=8", "HEAD"),
                                                stdout = TRUE, stderr = FALSE)), error = function(e) character())
    dirty <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(pkg), "status", "--porcelain",
                                                        "--untracked-files=no", "--", "R", "NAMESPACE", "data"),
                                               stdout = TRUE, stderr = FALSE)), error = function(e) character())
    commit <- if (length(commit) == 1L) paste0(commit, if (length(dirty)) " with uncommitted changes in R/, NAMESPACE, or data/" else "") else NA_character_
    list(ns = ns, from = sprintf("DMAR %s from the source tree %s (pkgload)", getNamespaceVersion(ns), pkg),
         commit = commit, origin = sprintf("DMAR %s, source tree loaded with pkgload", getNamespaceVersion(ns)))
  }
}, error = function(e) e)
if (inherits(loaded, "error")) {
  cat("DMAR could not be loaded:", conditionMessage(loaded), "\n")
  cat("revdep contract: 0 of", length(dependents), "dependents read (DMAR not loaded); not run\n")
  quit(status = 2L)
}
ns <- loaded$ns
dmar_exports <- sort(getNamespaceExports(ns))
dmar_datasets <- sort(ls(getNamespaceInfo(ns, "lazydata"), all.names = TRUE))
dmar_all <- ls(ns, all.names = TRUE)

## ---- What a dependent reaches ------------------------------------------------------------
# One row per reference: the object's name, how it is reached, and where.
is_empty_arg <- function(x) is.symbol(x) && !nzchar(as.character(x))
lit <- function(x) if (is.character(x) && length(x) == 1L && !is.na(x)) x else NA_character_
read_dependent <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(file.path(path, "DESCRIPTION")) || !dir.exists(file.path(path, "R")))
    return(list(read = FALSE, path = path, why = sprintf("no package checkout at %s (no DESCRIPTION or R/)", path)))
  desc <- read.dcf(file.path(path, "DESCRIPTION"))
  field <- function(f) if (f %in% colnames(desc)) desc[1L, f] else ""
  depends_dmar <- grepl("(^|[,[:space:]])DMAR([[:space:](,]|$)", field("Depends"))
  nsf <- file.path(path, "NAMESPACE")
  ns_calls <- if (file.exists(nsf)) as.list(parse(nsf, keep.source = FALSE)) else list()
  head_is <- function(e, f) is.call(e) && is.name(e[[1L]]) && as.character(e[[1L]]) == f
  ns_args <- function(e) vapply(as.list(e)[-1L], function(a) as.character(a), "")
  import_all <- any(vapply(ns_calls, function(e) head_is(e, "import") && "DMAR" %in% ns_args(e), logical(1)))
  refs <- list()
  add <- function(name, how, where) refs[[length(refs) + 1L]] <<- data.frame(name = name, how = how, where = where)
  other_imports <- character()
  for (e in ns_calls) if (head_is(e, "importFrom")) {
    a <- ns_args(e)
    if (a[1L] == "DMAR") for (n in a[-1L]) add(n, "importFrom(DMAR)", "NAMESPACE") else other_imports <- c(other_imports, a[-1L])
  }
  files <- list.files(file.path(path, "R"), pattern = "[.][RrSsq]$", full.names = TRUE)
  dynamic <- FALSE
  assigned <- character()
  fun_strings <- list(); bare_calls <- list()   # per file
  for (f in files) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) e)
    if (inherits(exprs, "error")) return(list(read = FALSE, path = path,
                                             why = sprintf("%s does not parse: %s", f, conditionMessage(exprs))))
    where <- file.path("R", basename(f))
    names_dmar_ns <- FALSE
    get_literals <- character()
    heads <- character()
    walk <- function(x) {
      if (is.pairlist(x) && length(x)) { for (k in seq_along(x)) walk(x[[k]]); return(invisible()) }
      if (is_empty_arg(x) || !is.call(x)) return(invisible())
      h <- x[[1L]]
      if (is.name(h)) {
        hn <- as.character(h)
        if (hn %in% c("::", ":::") && length(x) == 3L && as.character(x[[2L]]) == "DMAR") {
          add(as.character(x[[3L]]), paste0("DMAR", hn), where)
          return(invisible())
        }
        if (hn %in% c("<-", "=", "<<-") && is.name(x[[2L]])) assigned <<- c(assigned, as.character(x[[2L]]))
        if (hn == "function" && !is.null(x[[2L]])) assigned <<- c(assigned, names(x[[2L]]))
        if (hn == "asNamespace" && length(x) >= 2L && identical(lit(x[[2L]]), "DMAR")) names_dmar_ns <<- TRUE
        if (hn == "getExportedValue" && length(x) >= 3L) {
          m <- tryCatch(match.call(function(ns, name) NULL, x), error = function(e) NULL)
          if (!is.null(m) && identical(lit(m$ns), "DMAR")) {
            if (!is.na(lit(m$name))) add(lit(m$name), "getExportedValue(\"DMAR\", \"name\")", where) else dynamic <<- TRUE
          }
        }
        if (hn == "getFromNamespace" && length(x) >= 3L) {
          m <- tryCatch(match.call(function(x, ns, pos, envir) NULL, x), error = function(e) NULL)
          if (!is.null(m) && identical(lit(m$ns), "DMAR") && !is.na(lit(m$x)))
            add(lit(m$x), "getFromNamespace(\"name\", \"DMAR\")", where)
        }
        if (hn %in% c("get", "get0", "exists") && length(x) >= 2L && !is.na(lit(x[[2L]])))
          get_literals <<- c(get_literals, lit(x[[2L]]))
        heads <<- c(heads, hn)
      }
      nm <- names(x)
      if (!is.null(nm)) for (k in which(nm == "fun")) {
        v <- lit(x[[k]])
        if (!is.na(v)) fun_strings[[length(fun_strings) + 1L]] <<- data.frame(name = v, where = where)
      }
      for (k in seq_along(x)) walk(x[[k]])
      invisible()
    }
    for (e in exprs) walk(e)
    heads <- unique(heads)
    if (length(heads)) bare_calls[[length(bare_calls) + 1L]] <- data.frame(name = heads, where = where)
    if (names_dmar_ns) for (g in intersect(get_literals, c(dmar_exports, dmar_datasets)))
      add(g, "get(\"name\", <DMAR namespace>)", where)
  }
  fs <- do.call(rbind, fun_strings)
  if (dynamic && !is.null(fs)) {
    fs <- fs[fs$name %in% dmar_exports, , drop = FALSE]
    for (i in seq_len(nrow(fs))) add(fs$name[i], "fun = \"name\" (called through getExportedValue(\"DMAR\", fun))", fs$where[i])
  }
  bc <- do.call(rbind, bare_calls)
  if ((import_all || depends_dmar) && !is.null(bc)) {
    bc <- bc[bc$name %in% dmar_exports & !bc$name %in% c(assigned, other_imports, ls(baseenv(), all.names = TRUE)), , drop = FALSE]
    for (i in seq_len(nrow(bc))) add(bc$name[i], if (import_all) "bare call (import(DMAR))" else "bare call (Depends: DMAR)", bc$where[i])
  }
  refs <- do.call(rbind, refs)
  if (is.null(refs)) refs <- data.frame(name = character(), how = character(), where = character())
  commit <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(path), "rev-parse", "--short", "HEAD"),
                                              stdout = TRUE, stderr = FALSE)), error = function(e) character())
  dirty <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(path), "status", "--porcelain", "--untracked-files=no"),
                                             stdout = TRUE, stderr = FALSE)), error = function(e) character())
  list(read = TRUE, path = path, refs = unique(refs),
       commit = if (length(commit) == 1L) commit else NA_character_, n_dirty = length(dirty))
}

## ---- The current state of an object, and of a call's result ----------------------------
deparse1_all <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
object_state <- function(name) {
  if (name %in% dmar_datasets) {
    d <- get(name, envir = getNamespaceInfo(ns, "lazydata"))
    return(list(kind = "dataset", class = class(d), names = if (is.null(names(d))) list() else names(d)))
  }
  if (!name %in% dmar_all) return(list(kind = "absent"))
  obj <- get(name, envir = ns, inherits = FALSE)
  if (!is.function(obj)) return(list(kind = class(obj)[1L], exported = name %in% dmar_exports))
  fm <- formals(obj)
  list(kind = "function", exported = name %in% dmar_exports,
       formals = stats::setNames(as.list(vapply(seq_along(fm), function(i) deparse1_all(fm[[i]]), "")), names(fm)))
}
shape <- function(a) list(class = class(a), length = length(a))
call_state <- function(cl) {
  set.seed(113)
  env <- new.env(parent = globalenv())
  res <- tryCatch(suppressWarnings(suppressMessages(eval(parse(text = cl$call, keep.source = FALSE)[[1L]], env))),
                  error = function(e) e)
  if (inherits(res, "error")) return(list(error = conditionMessage(res)))
  rn <- if (is.data.frame(res)) attr(res, "row.names") else rownames(res)
  auto_rn <- is.null(rn) || (is.numeric(rn) && identical(as.integer(rn), seq_len(NROW(res)))) ||
    (is.integer(rn) && length(rn) == 2L && is.na(rn[1L]))
  att <- list()
  for (a in unlist(cl$reads$attributes)) att[[a]] <- if (is.null(attr(res, a, exact = TRUE))) list(present = FALSE)
                                                    else c(list(present = TRUE), shape(attr(res, a, exact = TRUE)))
  list(class = class(res),
       names = if (is.null(names(res))) list() else names(res),
       term = if (is.data.frame(res) && "term" %in% names(res)) as.character(res$term) else list(),
       row_names = if (auto_rn) list() else as.character(rn),
       attributes = att)
}

## ---- Read the snapshot and the dependents ------------------------------------------------
snap <- if (file.exists(snapshot_file)) tryCatch(jsonlite::read_json(snapshot_file, simplifyVector = FALSE),
                                                 error = function(e) e) else NULL
if (inherits(snap, "error")) { cat("the snapshot", snapshot_file, "could not be read:", conditionMessage(snap), "\n"); quit(status = 2L) }
if (is.null(snap) && !update) {
  cat("no snapshot at", snapshot_file, "; create one with --update\n")
  cat("revdep contract: 0 of", length(dependents), "dependents read (no snapshot); not run\n")
  quit(status = 2L)
}
calls <- if (!is.null(snap$calls)) snap$calls else default_calls
declared <- if (!is.null(snap$declared)) snap$declared else default_declared

cat(sprintf("revdep contract: %s%s\n", loaded$from, if (!is.na(loaded$commit)) sprintf(" (commit %s)", loaded$commit) else ""))
if (!is.null(snap)) cat(sprintf("snapshot: %s (taken %s from %s)\n", snapshot_file, snap$taken$date, snap$taken$from))

dep <- lapply(dependents, read_dependent)
problems <- character()     # differences and broken reliances: exit 1
relied <- list()            # name -> dependents that rely on it, now
for (d in names(dep)) {
  x <- dep[[d]]
  if (!x$read) { cat(sprintf("[NOT RUN] %s: %s; its current use of DMAR was not read\n", d, x$why)); next }
  dec <- names(declared[[d]])
  names_now <- sort(unique(c(x$refs$name, dec)))
  hows <- table(factor(x$refs$how[!duplicated(x$refs[c("name", "how")])]))
  cat(sprintf("%s at %s (commit %s%s): relies on %d DMAR objects (%s%s)\n", d, x$path, x$commit,
              if (x$n_dirty) sprintf(", %d uncommitted change(s)", x$n_dirty) else "",
              length(names_now), paste(sprintf("%d by %s", as.integer(hows), names(hows)), collapse = "; "),
              if (length(dec)) sprintf("; %d declared in the snapshot", length(dec)) else ""))
  for (n in names_now) relied[[n]] <- c(relied[[n]], d)
  # A reliance on something this DMAR does not provide is a break whatever the snapshot says.
  ns_refs <- x$refs[x$refs$how == "DMAR::" | grepl("^importFrom|^bare|^getExported|^fun", x$refs$how), , drop = FALSE]
  for (n in unique(ns_refs$name)) if (!n %in% c(dmar_exports, dmar_datasets)) {
    w <- ns_refs[ns_refs$name == n, ]
    problems <- c(problems, sprintf("%s reaches %s (%s, %s), which this DMAR does not export", d, n, w$how[1L], w$where[1L]))
  }
  for (n in unique(x$refs$name[x$refs$how == "DMAR:::"])) {
    if (!n %in% dmar_all) problems <- c(problems, sprintf("%s reaches DMAR:::%s, which this DMAR does not have", d, n))
    else cat(sprintf("  note: %s reaches the internal DMAR:::%s, which no rule keeps stable\n", d, n))
  }
  for (n in dec) if (!n %in% c(dmar_exports, dmar_datasets))
    problems <- c(problems, sprintf("%s relies on %s (declared in the snapshot), which this DMAR does not export", d, n))
}
# What a dependent that was not read relies on is taken from the snapshot, so it is still checked.
for (d in names(dep)[!vapply(dep, `[[`, logical(1), "read")]) {
  for (n in unlist(snap$dependents[[d]]$objects)) relied[[n]] <- c(relied[[n]], d)
}

## ---- --update: rewrite the snapshot -------------------------------------------------------
if (update) {
  missing_dep <- names(dep)[!vapply(dep, `[[`, logical(1), "read")]
  if (length(missing_dep)) {
    cat(sprintf("--update refused: %s not read, so the snapshot would drop what it relies on\n", paste(missing_dep, collapse = ", ")))
    quit(status = 2L)
  }
  if (length(problems)) {
    cat(paste0("[DIFF] ", problems, "\n"), sep = "")
    cat("--update refused: a dependent reaches an object this DMAR does not provide; the snapshot would record a broken contract\n")
    quit(status = 2L)
  }
  objects <- lapply(sort(names(relied)), function(n) c(object_state(n), list(used_by = sort(unique(relied[[n]])))))
  names(objects) <- sort(names(relied))
  new_calls <- lapply(calls, function(cl) {
    st <- call_state(cl)
    if (!is.null(st$error)) { cat(sprintf("--update refused: the call '%s' fails: %s\n", cl$label, st$error)); quit(status = 2L) }
    c(cl[c("label", "used_by", "call", "reads")], list(result = st))
  })
  out <- list(
    about = paste("What DMAR's dependents rely on, as read and run by tools/revdep_contract.R; rewrite only",
                  "with --update after the differences it reports have been reviewed and the dependents told."),
    taken = list(date = format(Sys.Date()),
                 from = paste0(loaded$origin, if (!is.na(loaded$commit)) sprintf(", commit %s", loaded$commit) else "")),
    dependents = lapply(dep, function(x) list(commit = x$commit,
                                               objects = sort(unique(c(x$refs$name))),
                                               references = lapply(split(x$refs, x$refs$name), function(r)
                                                 unique(sprintf("%s in %s", r$how, r$where))))),
    declared = declared, objects = objects, calls = new_calls)
  for (d in names(out$dependents)) out$dependents[[d]]$objects <- sort(unique(c(out$dependents[[d]]$objects, names(declared[[d]]))))
  jsonlite::write_json(out, snapshot_file, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
  cat(sprintf("snapshot written: %s (%d objects, %d calls)\n", snapshot_file, length(objects), length(new_calls)))
  cat(sprintf("revdep contract: %d of %d dependents read (%s); 0 differences; snapshot rewritten\n",
              length(dep), length(dep), paste(names(dep), collapse = ", ")))
  quit(status = 0L)
}

## ---- Compare with the snapshot --------------------------------------------------------------
fmt <- function(x) if (length(x)) paste(unlist(x), collapse = ", ") else "(none)"
same <- function(a, b) identical(unlist(a), unlist(b))
# The set of objects relied on, per dependent read.
for (d in names(dep)) if (dep[[d]]$read) {
  was <- unlist(snap$dependents[[d]]$objects)
  now <- sort(unique(c(dep[[d]]$refs$name, names(declared[[d]]))))
  if (length(x <- setdiff(now, was))) problems <- c(problems, sprintf("%s now relies on %s, which the snapshot does not record", d, fmt(x)))
  if (length(x <- setdiff(was, now))) problems <- c(problems, sprintf("%s no longer relies on %s, which the snapshot records", d, fmt(x)))
}
n_fun <- 0L; n_data <- 0L
for (n in sort(names(relied))) {
  who <- paste(sort(unique(relied[[n]])), collapse = ", ")
  rec <- snap$objects[[n]]
  cur <- object_state(n)
  if (identical(cur$kind, "function")) n_fun <- n_fun + 1L else if (identical(cur$kind, "dataset")) n_data <- n_data + 1L
  if (is.null(rec)) next   # named above as not recorded
  if (!identical(cur$kind, rec$kind)) {
    problems <- c(problems, sprintf("%s (%s): recorded as %s, now %s", n, who, rec$kind, cur$kind)); next
  }
  if (identical(cur$kind, "function")) {
    if (!identical(isTRUE(cur$exported), isTRUE(rec$exported)))
      problems <- c(problems, sprintf("%s (%s): %s", n, who, if (isTRUE(cur$exported)) "now exported" else "no longer exported"))
    rf <- rec$formals; cf <- cur$formals
    if (!identical(names(rf), names(cf))) {
      problems <- c(problems, sprintf("%s (%s): formals differ: recorded (%s), now (%s)%s", n, who,
                                      fmt(names(rf)), fmt(names(cf)),
                                      if (length(x <- setdiff(names(rf), names(cf)))) sprintf("; gone: %s", fmt(x)) else ""))
    }
    for (a in intersect(names(rf), names(cf))) if (!identical(unlist(rf[[a]]), unlist(cf[[a]])))
      problems <- c(problems, sprintf("%s (%s): default of %s was %s, now %s", n, who, a,
                                      if (nzchar(unlist(rf[[a]]))) unlist(rf[[a]]) else "(none)",
                                      if (nzchar(cf[[a]])) cf[[a]] else "(none)"))
  } else if (identical(cur$kind, "dataset")) {
    if (!same(cur$class, rec$class)) problems <- c(problems, sprintf("%s (%s): class was %s, now %s", n, who, fmt(rec$class), fmt(cur$class)))
    if (!same(cur$names, rec$names)) problems <- c(problems, sprintf("%s (%s): columns differ; gone: %s; new: %s", n, who,
                                                                    fmt(setdiff(unlist(rec$names), unlist(cur$names))),
                                                                    fmt(setdiff(unlist(cur$names), unlist(rec$names)))))
  }
}
n_obj_problems <- length(problems)
plural <- function(n, what) sprintf("%d %s%s", n, what, if (n == 1L) "" else "s")
cat(sprintf("[%s] %s and %s relied on: present, exported, and formals or columns as recorded\n",
            if (n_obj_problems) "DIFF" else "OK", plural(n_fun, "function"), plural(n_data, "dataset")))
call_problems <- character()
for (cl in calls) {
  who <- fmt(cl$used_by)
  st <- call_state(cl)
  rec <- cl$result
  tag <- sprintf("call '%s' (%s)", cl$label, who)
  if (!is.null(st$error)) { call_problems <- c(call_problems, sprintf("%s fails: %s", tag, st$error)); next }
  if (is.null(rec)) { call_problems <- c(call_problems, sprintf("%s has no recorded result", tag)); next }
  if (!same(st$class, rec$class)) call_problems <- c(call_problems, sprintf("%s: class was %s, now %s", tag, fmt(rec$class), fmt(st$class)))
  for (part in c("names", "term", "row_names")) if (!same(st[[part]], rec[[part]])) {
    gone <- setdiff(unlist(rec[[part]]), unlist(st[[part]]))
    read_gone <- intersect(gone, unlist(c(cl$reads$terms, cl$reads$columns)))
    call_problems <- c(call_problems, sprintf("%s: %s differ: recorded (%s), now (%s)%s", tag,
                                              c(names = "columns", term = "terms", row_names = "row names")[[part]],
                                              fmt(rec[[part]]), fmt(st[[part]]),
                                              if (length(read_gone)) sprintf("; gone and read by %s: %s", who, fmt(read_gone)) else ""))
  }
  for (a in names(st$attributes)) if (!identical(unlist(st$attributes[[a]]), unlist(rec$attributes[[a]])))
    call_problems <- c(call_problems, sprintf("%s: attribute %s (read by %s) was %s, now %s", tag, a, who,
                                              if (isTRUE(rec$attributes[[a]]$present)) sprintf("%s of length %s", fmt(rec$attributes[[a]]$class), rec$attributes[[a]]$length) else "absent",
                                              if (isTRUE(st$attributes[[a]]$present)) sprintf("%s of length %s", fmt(st$attributes[[a]]$class), st$attributes[[a]]$length) else "absent"))
}
cat(sprintf("[%s] %d canonical calls: class, columns, terms or row names, and the attributes read as recorded\n",
            if (length(call_problems)) "DIFF" else "OK", length(calls)))
problems <- c(problems, call_problems)
if (length(problems)) cat(paste0("[DIFF] ", problems, "\n"), sep = "")
read_ok <- vapply(dep, `[[`, logical(1), "read")
cat(sprintf("revdep contract: %d of %d dependents read (%s); %d difference%s\n", sum(read_ok), length(dep),
            paste(sprintf("%s %s", names(dep), ifelse(read_ok, "read", "not run")), collapse = ", "), length(problems),
            if (length(problems) == 1L) "" else "s"))
if (length(problems)) quit(status = 1L)
if (!all(read_ok)) quit(status = 3L)
quit(status = 0L)
