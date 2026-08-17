## data-raw/pygmalion.R
##
## Construction script for the pygmalion data set.
##
## Source: pygmalion_raw.csv, a faithful 310-row transcription of the
## teacher-expectancy ("Pygmalion in the Classroom") data as used in
## Maxwell, Delaney, and Kelley, Designing Experiments and Analyzing
## Data: A Model Comparison Perspective (Routledge), where it is the
## running example for the analysis of covariance with heterogeneity
## of regression (Chapter 9; the same data also appear as a Chapter 3
## exercise). The transcription preserves the book's column names
## (Grade, Treatment, IQPre, IQ4, IQ8, IQGain) and the 0/1 treatment
## coding so that the source file is an exact copy of the analysis
## data; this script does the (purely cosmetic) renaming and factor
## labeling and verifies that every published quantity reproduces.
##
## The same numeric data ship in the AMCP package (companion to the
## book) as chapter_9_exercise_15 / chapter_9_extension_exercise_3
## (n = 310, with IQGain) and chapter_3_exercise_22 (n = 310, without
## IQGain). No value is altered here; only the categorical treatment
## variable is given human-readable labels and the columns are renamed
## to DMAR's descriptive snake_case house style.
##
## To rebuild from scratch, run this script from the package root:
##   source("data-raw/pygmalion.R")

raw_path <- "data-raw/pygmalion_raw.csv"
stopifnot(file.exists(raw_path))

raw <- read.csv(raw_path, stringsAsFactors = FALSE)
stopifnot(
  identical(names(raw), c("Grade", "Treatment", "IQPre", "IQ4", "IQ8", "IQGain")),
  nrow(raw) == 310L
)

# Treatment: in the source, 1 = the (randomly chosen) children whose
# teachers were led to expect unusual intellectual growth -- the
# "bloomers" or experimental group (n = 64); 0 = the remaining
# children -- the control group (n = 246). Control is made the
# reference level so that the model intercept and the IQPre slope
# correspond to the control group, matching the book's parameterization.
treatment <- factor(
  ifelse(raw$Treatment == 1L, "Bloomer", "Control"),
  levels = c("Control", "Bloomer")
)

pygmalion <- data.frame(
  grade     = as.integer(raw$Grade),
  treatment = treatment,
  iq_pre    = as.integer(raw$IQPre),
  iq_4      = as.integer(raw$IQ4),
  iq_8      = as.integer(raw$IQ8),
  iq_gain   = as.integer(raw$IQGain),
  stringsAsFactors = FALSE
)

# ---- Structural checks -----------------------------------------------
stopifnot(
  identical(dim(pygmalion), c(310L, 6L)),
  identical(sum(pygmalion$treatment == "Control"), 246L),
  identical(sum(pygmalion$treatment == "Bloomer"), 64L),
  identical(as.vector(table(pygmalion$grade)),
            c(52L, 58L, 53L, 59L, 34L, 54L)),
  all(pygmalion$grade %in% 1:6),
  # IQGain is, by construction, the total IQ change from pretest to IQ8.
  all(pygmalion$iq_gain == pygmalion$iq_8 - pygmalion$iq_pre)
)

# ---- Content spot checks (first source row) --------------------------
stopifnot(
  pygmalion$grade[1]   == 1L,
  pygmalion$treatment[1] == "Control",
  pygmalion$iq_pre[1]  == 45L,
  pygmalion$iq_4[1]    == 58L,
  pygmalion$iq_8[1]    == 76L,
  pygmalion$iq_gain[1] == 31L,
  sum(pygmalion$iq_pre)  == 30524L,
  sum(pygmalion$iq_8)    == 33448L,
  sum(pygmalion$iq_gain) ==  2924L
)

# ---- Reproduce the published analysis --------------------------------
# Heterogeneity-of-regression ANCOVA: separate IQ8-on-IQpre slopes for
# the two groups. The book reports control slope 0.77799 and bloomer
# slope 0.96895; the pooled within-group residual variance is 175.3251
# and the (sample) variance of the covariate is 348.9099. These are the
# inputs to the "variance of the estimated treatment effect" worked
# example (see MBESS::var.ete / the var_ete development code).
fit_het <- lm(iq_8 ~ iq_pre * treatment, data = pygmalion)
b <- coef(fit_het)
slope_control <- unname(b["iq_pre"])
slope_bloomer <- unname(b["iq_pre"] + b["iq_pre:treatmentBloomer"])
resid_var     <- sum(resid(fit_het)^2) / fit_het$df.residual

# Compare to the published constants with a small tolerance: lm gives a
# bloomer slope of 0.96894, whereas MBESS::var.ete lists 0.96895 (a
# fifth-decimal rounding difference), and var.ete lists the covariate
# variance as 348.9099 vs. 348.9097 from var().
stopifnot(
  abs(slope_control - 0.77799)  < 5e-5,
  abs(slope_bloomer - 0.96895)  < 5e-5,
  abs(resid_var     - 175.3251) < 1e-3,
  abs(var(pygmalion$iq_pre) - 348.9099) < 1e-3,
  fit_het$df.residual == 306L
)

save(pygmalion,
     file = "data/pygmalion.rda",
     version = 2L, compress = "bzip2")

message("pygmalion: built and verified (310 x 6); ",
        "control slope ", round(slope_control, 5),
        ", bloomer slope ", round(slope_bloomer, 5),
        ", residual var ", round(resid_var, 4), ".")
