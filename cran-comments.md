# DMAR 1.0.0 submission (resubmission)

Outcome: accepted. The automated pretest passed, Uwe Ligges moved the
submission to the manual inspection queue for first submissions on
2026-09-08 ("Not so much you can do, hence moved to the manual newbie
inspection queue"), and DMAR 1.0.0 was published on CRAN on 2026-09-21
from the tarball packaged 2026-09-07 23:48 UTC (the gate run of that
evening; CRAN adds Repository and Date/Publication fields to the
DESCRIPTION it publishes, so the published file's checksum differs from
the certified one by design). The farm's first results (2026-09-22) are
OK on every flavor, with total check times of 399 to 403 s on the Linux
flavors and 632 to 681 s on the Intel macOS flavors; `tools/cran_status.R`
carries those values and reports any drift. GitHub carries this release
as tag `v1.0.0.0` with Version 1.0.0.0, which R orders as the same
version as CRAN's 1.0.0: the same code, with the fourth component zero
at a release; the commit whose tree matches the CRAN tarball exactly is
bcc0624 (2026-09-07), and the version, README, and NEWS wording were
brought up to date after publication. The next submission starts a new
letter; this one stays as the record of the four rounds.

## Response to the CRAN review of 2026-09-04

The review of the 2026-08-21 submission named five things. Each is
answered below; the four code-level items each have a mechanical
detector in the test suite (`tests/testthat/test-rd_hygiene.R`, which
runs on CRAN) and in the maintainer's release tooling, so none can
return, and the first (the title) was settled by the email exchange
that granted the exception.

* **"in R" in the title and description.** The description no longer
  says it: it opens "Methods for design, measurement, and analysis,
  with the aim of ...". The title keeps it because the title is the
  expansion of the package name: DMAR stands for Design, Measurement,
  and Analysis in R. An exception was requested by email on 2026-09-04
  and granted on 2026-09-05 (K. Lauseker), with the request that the
  acronym be added in parentheses directly after the title so the
  connection is obvious; the Title field is now "Design, Measurement,
  and Analysis in R (DMAR)".

* **Code lines in examples commented out.** Every example line in the
  package now runs. The comment idiom was adopted in the 2026-07-31
  round in place of `\donttest{}` (which makes `R CMD check --as-cran`
  run the whole example corpus a second time), and it is gone from
  every page: the eighteen pages the review listed and twenty-nine
  others that a parse-based detector found (every window of comment
  lines in every examples block is parsed; a window that parses to a
  call or an assignment is commented-out code). The detector is now a
  permanent test. Calls that were commented out because they were
  slow run with a small replication count, stated in a comment along
  with the count a reported analysis deserves; nothing was replaced by
  smaller or artificial data. Three pages in the list (`cov_sem`,
  `dmar_tbl`, `holzinger_swineford`) carried prose comments only;
  their comments were reworded anyway so that no comment line reads
  like code. The package still contains no `\donttest{}` and no
  `\dontrun{}`: with every line live, the slowest help page takes
  1.3 s locally, and all 312 pages with examples run in 29 s in a
  single pass.

* **Functions writing to the home filespace, and default paths.** The
  twenty-seven sensitivity functions that carried
  `save = FALSE, filename = "<name>.csv"` no longer have a default path:
  the `save` switch is gone and `filename` defaults to `NULL`, so
  nothing is written unless the user supplies a path.
  `correlations_test()` already took `file = NULL`. No other function
  writes a file. The examples and tests that exercise writing use
  `tempfile()` and remove the file afterwards; no example, test, or
  vignette writes anywhere but `tempdir()`.

* **Modifying `.GlobalEnv`.** The package restored the user's random
  number generator state after a seeded call with
  `assign(".Random.seed", ..., envir = .GlobalEnv)` in twenty-five
  functions, and counted warnings inside calling handlers with `<<-`
  (each modifying a local of the enclosing function, but the operator
  is the operator). Both are gone. A supplied seed now goes through
  `withr::local_seed()` in one internal helper, which is why `withr`
  (which has no dependencies of its own) joined Imports, the only
  dependency change; the counters live in local environments. A test
  deparses every function in the namespace and fails on any
  `.GlobalEnv`, `globalenv()`, `.Random.seed`, `<<-`, or `options(warn`.

* **`options(warn = -1)`.** Removed from the six functions that set
  it (each had restored the previous value on exit). Each now muffles
  only the specific warnings its Monte Carlo loop is known to raise
  (lavaan convergence and variance warnings on borderline replicates,
  the noncentral clamp), with `withCallingHandlers()` or
  `suppressWarnings()` on the call that emits them; every other
  warning reaches the user.

The local check for this round, `R CMD check --as-cran` on the tarball
with both manuals built, ran with 0 errors, 0 warnings, and the two
NOTEs described under "R CMD check results" below. The full local test
suite runs 8,479 expectations with no failures and none skipped. On
CRAN the check runs a curated subset of the suite (see "Response to
the incoming pretest of 2026-09-05" below): 81 test files carrying
2,068 expectations, comprising the mechanical detectors promised
above, the published-value and oracle anchor files that pin the
package's numerics across platforms (every closed-form `ss_aipe_*`
planner file among them), and per-family smoke tests. The full suite,
Monte Carlo blocks included, runs locally and in the maintainer's
release gate before any upload; no test was deleted or weakened.

## Response to the incoming pretest of 2026-09-05

The 2026-09-05 pretest returned an overall checktime of 14 minutes on
r-devel-windows against the 10 minute threshold. The check used to
take much longer, and the reason is deliberate: the aim of this
package is a comprehensive and robust treatment of design,
measurement, and analysis, and that breadth carries real check weight
(313 documented help pages with every example line live, eighteen
teaching vignettes, nine data sets, a test suite of 8,479
expectations). Earlier rounds brought the checktime from 18 minutes
to 14 by removing the `\donttest{}` double pass and trimming
replication counts; this round moves the remaining computation off
CRAN's clock without removing anything from the package.

* All eighteen vignettes are now precomputed. Their executable
  sources are maintained in the repository and regenerated by a
  maintainer script; the shipped vignettes are statically knitted
  twins, so the check re-runs no vignette computation. (Four were
  already built this way; the other fourteen now follow the same
  convention.)

* The CRAN check now runs the curated test subset described above
  (81 files, 2,068 expectations, about 19 seconds serially on the
  local machine against 53 for the previous CRAN path), selected so
  that every mechanical detector promised in the previous round, every
  published-value and oracle anchor family, and every closed-form
  sample size planner keeps cross-platform verification on CRAN's
  farm. The full suite continues to run locally and in the release
  gate.

* Examples were already single-pass and fast (slowest page 1.2
  seconds) and are unchanged.


## Response to the incoming pretest of 2026-08-19

The 2026-08-19 pretest (Debian and Windows) returned the package for
a LaTeX error in the PDF manual, traced to a single help page:
`ss_aipe_equivalence_smd_sensitivity.Rd` carried `\code{}` markup
inside an `\eqn{}`, which renders as `\texttt` inside LaTeX math and
fails the manual build ("Missing $ inserted"); the HTML manual
reported the same line as a math rendering problem, and the leftover
`DMAR-manual.tex` NOTE followed from the failed build. The equation
is now written in symbols, with the argument names in the
surrounding prose, and no `\eqn{}` or `\deqn{}` in the package
contains Rd markup. The same pretest asked for a trailing slash on
the package website URL in README.md; added. This round's local
check ran `R CMD check --as-cran` with the PDF and HTML manuals
built (the earlier rounds had used `--no-manual` locally, which is
how the line slipped through): both manuals now build clean.

## Response to the incoming pretest of 2026-07-31

Third pretest (overall checktime 18 minutes on Windows). The previous
round tried to buy time by wrapping the slow examples in
`\donttest{}`, which was the wrong lever, and this round removes that
construct from the package entirely.

Under `R CMD check --as-cran`, a single `\donttest{}` block anywhere
in a package causes the entire example corpus to be run twice: once
as `checking examples`, with the donttest blocks commented out, and
again as `checking examples with --run-donttest`, with everything
included. With 59 of the help pages carrying a wrapper, all 304 pages with
examples were being run a second time. The wrapper was costing time
rather than saving it.

The package now contains no `\donttest{}` and no `\dontrun{}` at all.
The examples that used to be slow are handled in one of three ways.
Replication counts in the cheap demonstrations were lowered, each with
a comment naming what a reported analysis deserves. No example runs a
bootstrap confidence interval; those calls are carried as commented
code so a reader of `?fn` still sees the exact syntax, and the help
page prose explains the interval and when to ask for it. The
randomization tests keep their resampling, since permuting the data is
what those functions do and their examples cost a tenth of a second.
Anything else that was expensive is likewise commented out, with a
sentence saying what it does and why it is not run, so the code
remains available to a reader working through the page. The data
behind every example are unchanged; nothing was replaced by a smaller
or artificial data set. The slowest help page now takes under a
second, and all 313 pages with examples together take 42 seconds
locally (122 seconds on win-builder), in a single pass.

Alongside that, that round skipped the remaining expensive test
blocks on CRAN at the block level: every file it touched kept a fast
published-value anchor on the CRAN path, and its tarball ran 6,881 of
the suite's then 8,251 expectations there. (The current round
replaces that block-level convention with the curated file-level
subset described at the top of this letter.)

The package now ships eighteen vignettes. Three that were previously
included establish how the methods perform rather than how they are
used: a 10,000-replication Monte Carlo study of the AIPE planners, a
simulation study of the Bryant-Paulson intervals, and the reproduction
of the published Bryant-Paulson critical value tables. That is material
for a paper, not package documentation, so all three are maintained
outside the package. How to use those functions stays in the package,
in the `bryant_paulson_ancova` and `critical_values` vignettes.

In that round, four of the vignettes were precomputed, the four whose
code was genuinely expensive: their executable sources are the
`vignettes/*.Rmd.orig` files kept in the maintained repository, and
the shipped `.Rmd` files are statically knitted twins that execute no
computation at check time. As of the current round all eighteen
vignettes follow that convention (see "Response to the incoming
pretest of 2026-09-05" above).

The two composite planning vignettes had their Monte Carlo replication
counts cut, which is what took `composite_sem_planning` from eleven
minutes to twenty-two seconds. Both open with a note saying plainly
that the printed values are an illustration of the workflow and not
planning values, and giving the simulation standard error the count
implies. In `composite_power_ancova`, which runs at 25 replications,
every simulated quantity is additionally reported beside the same
quantity at 10,000 replications, so a reader sees both the method and
the answer; the reference values are regenerated by
`tools/composite_power_reference.R`.

On win-builder (R-release, 2026-08-17, against this exact tarball)
the corresponding steps ran examples in 122 seconds, tests in 115,
and the vignette rebuild in 53, with no step flagged for time; the `--run-donttest` pass no longer
exists because the package contains no `\donttest{}`.

The tarball is 4.8 Mb. The three simulation studies that left the
package in the previous round remain outside it; the growth since
then is the nonlinear growth curve family (five simulators, the
`analysis_of_change()` fitter, and their vignette).

Suggests is down from thirty packages to fifteen. Every package that
appeared only as a test oracle (MBESS, emmeans, semTools, metafor,
mirt, sirt, multcomp, performance, r2mlm, irr, irrCAC, psych,
BayesFactor, gsl) has been replaced by pinned reference values, each
carrying a provenance comment naming the package, its version, and the
capture date; the live comparisons run at release time from
tools/oracle_checks.R, which is not shipped. Three runtime uses were
removed rather than pinned: BiasedUrn by computing Fisher's noncentral
hypergeometric density directly, kableExtra by building the HTML and
LaTeX tables with knitr alone, and AMCP by shipping the Chapter 9
depression data as depression_bdi. What remains in Suggests is what
the features genuinely use, and none of it requires a compiled system
library; lme4 joined for the optional mixed-effects method of
`analysis_of_change()`. Imports is unchanged and remains base R plus
generics.

The word "reimagining" flagged by the spellchecker is standard
American English (Merriam-Webster).

First pretest, all items addressed:

* Test failure (tests, both platforms): the comparison of iterative
  maximum likelihood against closed-form least squares coefficients in
  the no-intercept mlmr() test assumed agreement at 1e-6, which the
  r-devel Windows and Debian optimizer builds do not deliver; the
  tolerance is now 1e-5 with a comment recording why.
* Vignette rebuild failure (both platforms): contrast weight
  validation required the coefficient sum to equal zero exactly, and
  a thirds-based contrast sums to zero on some platforms but to about
  1e-16 under long double accumulation on others. All nine strict
  checks now use a 1e-8 tolerance, matching the package's other
  contrast validations; genuinely wrong weights still stop.
* Examples over 10 seconds (Windows): the slow portions of the four
  flagged pages (bryant_paulson, average_variance_extracted,
  ci_c_ancova_bp, mlmr) are wrapped in donttest; each page's
  remaining examples run in about a second locally.
* Invalid URL (301): the site address carries its trailing slash in
  DESCRIPTION, the package help page, and README.


First submission of DMAR, a greatly expanded reimagining of the
MBESS package (on CRAN since 2004, by the same author).

## Test environments

* local macOS (Apple Silicon), R 4.6.1 (2026-09-07): full
  `R CMD check --as-cran` with both manuals built, from the release
  tarball, built from a clean archive of the repository
* win-builder, R-release 4.6.1 on Windows Server 2022: the 2026-08-17
  check of the previous tarball returned 1 NOTE (the new-submission
  NOTE, with the five domain terms addressed under "Possibly
  misspelled words" below). [Run win-builder on this round's tarball
  before upload and record its date and status here.]

## R CMD check results

0 errors, 0 warnings, 1 NOTE (plus one that is environmental):

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ken Kelley <kkelley@nd.edu>'

New submission
```

The "New submission" note is expected for the first submission.

The two addresses that the previous pretest flagged as possibly
invalid are no longer written as links. Both are archival citations
rather than resources a reader needs to follow from the help page: the
original Indiana House Bill 1166, which authorized the Prime Time
program the `prime_time_achievement` data set documents, and the
Education Week article on the Project STAR class size work. The
Indiana General Assembly archive returns HTTP 403 to non-browser user
agents, so the address is now set with `\verb{}` in
`man/prime_time_achievement.Rd`. It is still displayed in full, and
still resolves in a browser; it is simply not a link the URL checker
follows.

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent
enough HTML Tidy.
```

Environmental, on the local machine only; the CRAN check farm has
a current HTML Tidy and this NOTE does not appear there.

The check may also report an installed-size NOTE; the size is
justified by the 18 vignettes and nine teaching data sets that carry
the package's pedagogical purpose. Win-builder's R-release check
reported no size NOTE.

## Check time

Locally the release-gate check (2026-09-07, on the submitted tarball)
ran examples in 59 seconds across 312 pages in a single pass (the
package contains no `\donttest{}` block), the curated test subset in
about 19 seconds serially, and the vignette rebuild in 14 seconds
(rendering only; every vignette is precomputed), with no step flagged
for time; the slowest help page took 1.3 seconds. The
2026-08-17 win-builder run of the previous tarball took about three
times the local elapsed time at each step (examples 122 seconds,
tests 115, vignette rebuild 53), which keeps every page well under
the 5 second flag there. The four computation-heavy vignettes are
precomputed, and the Monte Carlo, repeated model fitting, and numeric
integration test blocks are skipped on CRAN with their
published-value anchors retained.

## Possibly misspelled words in DESCRIPTION

The CRAN incoming-feasibility spellchecker may flag the following
domain-specific terms; all are correctly spelled:

* **AIPE** is the established acronym for "accuracy in parameter
  estimation," used throughout the methodological literature
  (Kelley & Maxwell, 2003, *Psychological Methods*; Kelley &
  Rausch, 2006, *Psychological Methods*; and others).
* **noncentral** is the standard spelling for the family of
  probability distributions used to build confidence intervals
  from F, t, and chi-square statistics under the alternative
  hypothesis.
* **nonstandard** is standard American English (Merriam-Webster);
  the Description uses it for methods outside the routine toolkit.
* **MBCO** is the literature's acronym for model-based constrained
  optimization (Tofighi & Kelley, 2020, Psychological Methods).
* **reimagining** is standard American English (Merriam-Webster).
* **ANCOVA** and **noninferiority** are standard statistical terms;
  **methodologists** is standard American English.

## Downstream dependencies

This is a new package; no reverse dependencies on CRAN.

## Relationship to MBESS

DMAR is a greatly expanded reimagining of MBESS, with extensive
notational and programming changes that would have broken existing MBESS code,
examples, and downstream packages. MBESS remains stable on CRAN;
DMAR is the recommended path forward for new users. The two
packages are designed to coexist.
