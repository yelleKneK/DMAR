# Build data/teacher_expectancy.rda
#
# Source of record: Raudenbush (1984), Table 1 (Journal of Educational
# Psychology, 76, 85-97): effect sizes d (treatment effect in IQ points
# divided by the control group posttest standard deviation), one-tailed
# p-values, estimated weeks of teacher-student contact prior to expectancy
# induction, group vs. individual testing, and aware vs. blind test
# administrator, for 18 experiments; Pellegrini and Hicks (1972) contributes
# two conditions (aware and blind), so the table has 19 effect rows. The
# study-level Pellegrini and Hicks values used in Raudenbush's 18-study
# analyses (d = 0.52, one-tailed p = .010) are documented on the help page.
#
# Per-condition sample sizes are not printed in the 1984 table; they are
# taken from the same studies as tabulated in Raudenbush and Bryk (1985)
# (the version distributed as dat.raudenbush1985 in metafor), matched by
# author, year, weeks, testing, and tester. Note the 1985 chapter
# standardizes differently, so its yi values are not equal to the 1984 d
# column; the d and p columns here are the 1984 paper's.

teacher_expectancy <- data.frame(
  study  = 1:19,
  author = c("Rosenthal, Baratz, & Hall", "Conn, Edwards, Rosenthal, & Crowne",
             "Jose & Cody", "Pellegrini & Hicks (tester aware)",
             "Pellegrini & Hicks (tester blind)", "Evans & Rosenthal",
             "Fielder, Cohen, & Feeney", "Fleming & Anttonen", "Claiborn",
             "Kester", "Maxwell", "Carter", "Flowers", "Keshock", "Henrikson",
             "Fine", "Ginsburg", "Grieger", "Rosenthal & Jacobson"),
  year   = c(1974L, 1968L, 1971L, 1972L, 1972L, 1969L, 1971L, 1971L, 1969L,
             1969L, 1970L, 1970L, 1966L, 1970L, 1970L, 1972L, 1970L, 1970L,
             1968L),
  weeks  = c(2L, 21L, 19L, 0L, 0L, 3L, 17L, 2L, 24L, 0L, 1L, 0L, 0L, 1L, 2L,
             17L, 7L, 5L, 1L),
  testing = factor(c("group", "group", "group", "group", "group", "group",
                     "group", "group", "group", "group", "individual",
                     "group", "group", "individual", "individual", "group",
                     "group", "group", "group")),
  tester  = factor(c("aware", "aware", "aware", "aware", "blind", "aware",
                     "blind", "blind", "aware", "aware", "blind", "blind",
                     "blind", "blind", "blind", "aware", "aware", "blind",
                     "aware")),
  n_experimental = c(77L, 60L, 72L, 11L, 11L, 129L, 110L, 233L, 26L, 75L,
                     32L, 22L, 43L, 24L, 19L, 80L, 65L, 72L, 65L),
  n_control      = c(339L, 198L, 72L, 22L, 22L, 348L, 636L, 224L, 99L, 74L,
                     32L, 22L, 38L, 24L, 32L, 79L, 67L, 72L, 255L),
  d            = c(0.02, 0.14, -0.03, 0.85, 0.19, -0.04, -0.02, 0.05, -0.13,
                   0.27, 0.55, 0.30, 0.18, -0.01, 0.16, -0.13, -0.02, -0.06,
                   0.21),
  p_one_tailed = c(.401, .206, .791, .003, .242, .709, .595, .224, .928,
                   .050, .002, .043, .210, .528, .250, .877, .519, .637,
                   .016),
  stringsAsFactors = FALSE
)

stopifnot(nrow(teacher_expectancy) == 19L)
# Study-level reconstruction (Pellegrini & Hicks merged to d = .52, p = .010)
# must reproduce the paper's own summary: M = 0.11, SD = 0.20, range .55 to -.13.
d18 <- append(teacher_expectancy$d[-c(4, 5)], 0.52, after = 3)
stopifnot(abs(mean(d18) - 0.11) < 0.005, abs(sd(d18) - 0.20) < 0.005,
          max(d18) == 0.55, min(d18) == -0.13)

save(teacher_expectancy, file = "data/teacher_expectancy.rda",
     compress = "xz")
