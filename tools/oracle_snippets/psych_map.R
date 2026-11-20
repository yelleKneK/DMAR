# Live oracle comparison against psych, run at release time after DMAR is
# installed (it calls DMAR::factor_count) when psych is installed locally,
# as the other snippets in this folder are.
## Written against psych 2.6.5.
## from tests/testthat/test-factor_count.R
local({
  # Velicer's (1976) minimum average partial: the average squared partial
  # correlation after m principal components. psych's VSS() computes the
  # same average (its `map`) for m = 1, 2, ...; factor_count() keeps it for
  # m = 0 (the average squared correlation) to p - 2 in its `map` attribute.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  d <- DMAR::holzinger_swineford[items]
  fc <- DMAR::factor_count(d, B = 0)
  vss <- suppressWarnings(suppressMessages(
    psych::VSS(stats::cor(d), n = 4, n.obs = nrow(d), rotate = "none",
               fm = "minres", plot = FALSE)))
  stopifnot(isTRUE(all.equal(unname(attr(fc, "map")[2:5]), unname(vss$map[1:4]),
                             tolerance = 1e-8)))
  # The eigenvalue column is the reduced correlation matrix's, as
  # parallel_analysis() forms it, and the chi square of each count is the
  # unrotated model's.
  pa <- DMAR::parallel_analysis(d, B = 200, seed = 1)
  stopifnot(isTRUE(all.equal(fc$eigenvalue, pa$estimate[1:2])))
  two <- DMAR::fa_exploratory(d, 2, rotation = "none")
  fm <- attr(two, "fit_measures")
  stopifnot(isTRUE(all.equal(fc$chi_square[2], fm$value[fm$term == "chi_square"])))
})
