## data-raw/holzinger_swineford.R
##
## Construction script for the holzinger_swineford data set and its
## HS_Data alias.
##
## Source: hs39all_raw.txt, a 33-column ASCII transcription of the
## complete Holzinger and Swineford (1939) data, one row per pupil,
## with "." in place of NA for tests 25 and 26 in the Pasteur sample.
## The values in this transcription agree, on every test cell, with
## MBESS::HS (since the post-4.6.0 correction) and with
## psychTools::holzinger.raw.
##
## Two cells in the ASCII source have age_months that disagree with
## the obvious formula (ids 225 and 277 list age = 14 years 0 months
## but age_months = 156 rather than 168). This script recomputes
## age_months as 12 * age + month_since_birthday so the column is
## internally consistent for all 301 cases. No test scores are
## altered.
##
## To rebuild from scratch, run this script from the package root:
##   source("data-raw/holzinger_swineford.R")

raw_path <- "data-raw/hs39all_raw.txt"
stopifnot(file.exists(raw_path))

cols <- c(
  "id", "school_code", "grade", "sex_code", "age",
  "month_since_birthday", "age_months_source",
  "t1_visual_perception", "t2_cubes", "t3_paper_form_board",
  "t4_lozenges", "t5_general_information",
  "t6_paragraph_comprehension", "t7_sentence",
  "t8_word_classification", "t9_word_meaning",
  "t10_addition", "t11_code",
  "t12_counting_groups_of_dots",
  "t13_straight_and_curved_capitals",
  "t14_word_recognition", "t15_number_recognition",
  "t16_figure_recognition", "t17_object_number",
  "t18_number_figure", "t19_figure_word",
  "t20_deduction", "t21_numerical_puzzles",
  "t22_problem_reasoning", "t23_series_completion",
  "t24_woody_mccall", "t25_paper_form_board_r", "t26_flags"
)

raw <- read.table(raw_path, header = FALSE, na.strings = ".",
                  col.names = cols, stringsAsFactors = FALSE)
stopifnot(nrow(raw) == 301L, ncol(raw) == 33L)

# School: 0 = Pasteur, 1 = Grant-White in the source.
school <- factor(ifelse(raw$school_code == 0L, "Pasteur", "Grant-White"),
                 levels = c("Grant-White", "Pasteur"))

# Sex: 0 = Male, 1 = Female in the source.
sex <- factor(ifelse(raw$sex_code == 0L, "Male", "Female"),
              levels = c("Female", "Male"))

# Recompute age_months from age and month_since_birthday so the
# column is internally consistent (the ASCII source has two rows
# where age_months disagrees with 12 * age + month_since_birthday).
age_months <- 12L * raw$age + raw$month_since_birthday
age_years  <- raw$age + raw$month_since_birthday / 12

holzinger_swineford <- data.frame(
  id                   = raw$id,
  sex                  = sex,
  grade                = raw$grade,
  age                  = raw$age,
  month_since_birthday = raw$month_since_birthday,
  age_months           = age_months,
  age_years            = age_years,
  school               = school,
  raw[, grep("^t[0-9]+_", names(raw), value = TRUE)],
  stringsAsFactors     = FALSE
)

# Structural and content checks.
stopifnot(
  identical(dim(holzinger_swineford), c(301L, 34L)),
  identical(sum(holzinger_swineford$school == "Pasteur"), 156L),
  identical(sum(holzinger_swineford$school == "Grant-White"), 145L),
  all(is.na(holzinger_swineford$t25_paper_form_board_r[
        holzinger_swineford$school == "Pasteur"])),
  all(is.na(holzinger_swineford$t26_flags[
        holzinger_swineford$school == "Pasteur"])),
  all(!is.na(holzinger_swineford$t25_paper_form_board_r[
        holzinger_swineford$school == "Grant-White"])),
  all(!is.na(holzinger_swineford$t26_flags[
        holzinger_swineford$school == "Grant-White"]))
)

# Spot-check against several cells that distinguish the corrected
# data from the pre-4.6.0 MBESS snapshot still shipped in
# sem::HS.data and OpenMx::HS.ability.data.
stopifnot(
  holzinger_swineford$t20_deduction[holzinger_swineford$id == 2L]        == -3L,
  holzinger_swineford$t22_problem_reasoning[holzinger_swineford$id == 2L] == 21L,
  holzinger_swineford$t24_woody_mccall[holzinger_swineford$id == 2L]      == 12L,
  holzinger_swineford$t15_number_recognition[holzinger_swineford$id == 1L] == 86L
)

save(holzinger_swineford,
     file = "data/holzinger_swineford.rda",
     version = 2L, compress = "bzip2")

# The short MBESS-era alias HS_Data is not saved as a second .rda:
# it is an active binding defined in R/aliases.R, so the bare name
# resolves to the same data frame after library(DMAR). Saving a
# duplicate here would ship two copies of the data.
