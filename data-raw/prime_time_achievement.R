## data-raw/prime_time_achievement.R
##
## Construction script for the prime_time_achievement data set.
##
## Source: the original Indiana Department of Education program
## evaluation file used in Lapsley, Daytner, Kelley, and Maxwell
## (2002, ERIC ED466679). The SPSS system file shipped to the
## evaluation team is preserved in data-raw/ verbatim
## ("Original New Achievement Data, DO NOT EDIT 7.23.01.sav") so
## that anyone who wants to rebuild the data set has the canonical
## source on hand.
##
## What this script does:
##   1. Reads the SPSS .sav file. The 999 student-level and 888
##      Not-Applicable codes are already represented in R as NA
##      after read_sav() because SPSS marks them as system-missing
##      via user-defined missingness ranges; the script verifies
##      this and recodes the 888 codes that survive (typmulti,
##      ptstatus) to NA.
##   2. Lowercases all variable names to match the published
##      multilevel-modeling code base (Finch, Bolin, and Kelley,
##      2019, 2nd ed.). Preserves the Indiana DOE typos calender,
##      hispanc1, hispanc2 for code compatibility with that book.
##   3. Drops the SPSS "Select Cases" filter variable FILTER_$
##      (geog == 4 indicator), which is not real study data.
##   4. Preserves the SPSS variable labels as a `label` attribute
##      on each column so users can recover the IDOE variable
##      description with attr(prime_time_achievement$X, "label").
##   5. Adds three derived unique-cluster identifiers that respect
##      the nesting irregularities documented in the data:
##        - corp_id   = paste(region, corp, sep = "_")
##        - school_id = paste(corp, school, sep = "_")
##        - class_id  = paste(corp, school, class, sep = "_")
##      corp 2400 spans regions 2 and 3 in the source file, so
##      paste(corp) alone would understate the cluster count by 1.
##      With corp_id the corp count matches the paper's "61
##      corporations." school_id and class_id are already unique
##      without region because the school numbers are unique
##      within corp; using corp keeps the IDs short and matches
##      common practice in the multilevel literature.
##
## To rebuild from scratch, run from the package root:
##   source("data-raw/prime_time_achievement.R")

stopifnot(requireNamespace("haven", quietly = TRUE))

sav_path <- "data-raw/prime_time_achievement.sav"
stopifnot(file.exists(sav_path))

sav <- haven::read_sav(sav_path)
stopifnot(nrow(sav) == 10927L)

# Capture variable labels before stripping haven_labelled class.
var_labels <- vapply(sav, function(x) {
  lab <- attr(x, "label")
  if (is.null(lab)) NA_character_ else as.character(lab)
}, character(1))

# Lowercase column names; rename FILTER_$ to a legal R name and
# drop it after the rename.
names(sav) <- tolower(names(sav))
sav[["filter_dol"]] <- sav[["filter_$"]]
sav[["filter_$"]] <- NULL
sav[["filter_dol"]] <- NULL

# Strip haven_labelled but keep variable label as attribute.
df <- as.data.frame(lapply(sav, function(x) {
  out <- if (inherits(x, "haven_labelled")) as.numeric(haven::zap_labels(x)) else x
  lab <- attr(x, "label")
  if (!is.null(lab)) attr(out, "label") <- as.character(lab)
  out
}), stringsAsFactors = FALSE)

# 888 codes in typmulti and ptstatus mark "Not Applicable" and
# survive the SPSS-to-R conversion. Recode to NA for consistency.
df$typmulti[df$typmulti == 888] <- NA
df$ptstatus[df$ptstatus == 888] <- NA

# Sanity: 999 was the IDOE student-level missing code and should
# already be NA after the haven import.
miss_999_vars <- c("age", "gender", "race", "geread", "gevocab",
                   "gereadcm", "gelang", "gelangmc", "gelangcm",
                   "gemath", "gemathcp", "gemathcm", "getotal",
                   "ncread", "ncvocab", "ncreadcm", "nclang",
                   "nclangmc", "nclangcm", "ncmath", "ncmathcp",
                   "ncmathcm", "nctotal",
                   "aaread", "aavocab", "aareadcm", "aalang",
                   "aalangmc", "aalangcm", "aamath", "aamathcp",
                   "aamathcm", "aatotal",
                   "npanverb", "npamem", "npaverb", "npatotal",
                   "csi")
for (v in miss_999_vars) {
  stopifnot(!any(df[[v]] == 999, na.rm = TRUE))
}

# Derived unique-cluster IDs.
df$corp_id   <- paste(df$region, df$corp, sep = "_")
df$school_id <- paste(df$corp,   df$school, sep = "_")
df$class_id  <- paste(df$corp,   df$school, df$class, sep = "_")

# Place identifiers up front; keep the rest of the column order
# (which mirrors the IDOE source file layout).
id_block <- c("id", "region", "corp", "school", "class",
              "corp_id", "school_id", "class_id")
rest     <- setdiff(names(df), id_block)
prime_time_achievement <- df[, c(id_block, rest)]

# Reattach the SPSS variable labels. Add labels for the three
# derived identifiers.
for (v in setdiff(names(prime_time_achievement), c("corp_id","school_id","class_id"))) {
  src_name <- v
  if (src_name %in% names(var_labels)) {
    lab <- var_labels[[src_name]]
    if (!is.na(lab) && nzchar(lab)) {
      attr(prime_time_achievement[[v]], "label") <- lab
    }
  }
}
attr(prime_time_achievement$corp_id,   "label") <-
  "Unique corporation identifier (region_corp), 61 distinct values"
attr(prime_time_achievement$school_id, "label") <-
  "Unique school identifier (corp_school), 163 distinct values"
attr(prime_time_achievement$class_id,  "label") <-
  "Unique classroom identifier (corp_school_class), 586 distinct values"

# Structural / reconciliation checks against ED466679.
stopifnot(
  nrow(prime_time_achievement) == 10927L,
  ncol(prime_time_achievement) == 113L,                          # 110 SPSS vars + 3 IDs
  length(unique(prime_time_achievement$corp_id))   == 61L,       # paper: 61 corporations
  length(unique(prime_time_achievement$school_id)) == 163L,      # paper: 163 schools
  length(unique(prime_time_achievement$class_id))  == 586L,      # paper: 573 (see Details)
  sum(prime_time_achievement$gender == 1L, na.rm = TRUE) == 5425L,
  sum(prime_time_achievement$gender == 2L, na.rm = TRUE) == 5457L,
  sum(prime_time_achievement$race == 5L, na.rm = TRUE)   == 9207L, # Caucasian
  sum(prime_time_achievement$race == 2L, na.rm = TRUE)   == 995L,  # African American
  sum(prime_time_achievement$race == 4L, na.rm = TRUE)   == 348L,  # Hispanic
  sum(prime_time_achievement$race == 1L, na.rm = TRUE)   == 16L,   # American Indian/Alaskan
  sum(prime_time_achievement$race == 3L, na.rm = TRUE)   == 62L,   # Asian American
  sum(prime_time_achievement$race == 6L, na.rm = TRUE)   == 188L,  # Multi-racial
  sum(prime_time_achievement$ptia == 1L) == 4021L,
  sum(prime_time_achievement$ptia == 2L) == 6789L,
  sum(prime_time_achievement$ptia == 3L) == 117L
)

save(prime_time_achievement,
     file = "data/prime_time_achievement.rda",
     version = 2L, compress = "bzip2")

# The short alias Prime_Time is not saved as a second .rda: it is an
# active binding defined in R/aliases.R, so the bare name resolves to
# the same data frame after library(DMAR). Saving a duplicate here
# would ship two copies of the data.
