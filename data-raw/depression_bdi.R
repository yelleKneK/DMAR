# Build data/depression_bdi.rda: the hypothetical three-group depression
# study of Maxwell, Delaney, and Kelley (2027), Chapter 9, Table 9.7.
# Thirty depressive individuals randomly assigned to a selective
# serotonin reuptake inhibitor (SSRI), placebo, or wait list control,
# with the Beck Depression Inventory administered before the study and
# again at its end. The values are transcribed from the book table (they
# also ship in the book's data companion, the AMCP package, as
# chapter_9_table_7); they are embedded as literals here so regenerating
# the dataset depends on nothing outside this file.

depression_bdi <- data.frame(
  condition = factor(rep(c("ssri", "placebo", "wait_list"), each = 10),
                     levels = c("ssri", "placebo", "wait_list")),
  bdi_pre  = c(18, 16, 16, 15, 14, 20, 14, 21, 25, 11,
               18, 16, 15, 14, 20, 25, 11, 25, 11, 22,
               15, 19, 10, 29, 24, 15,  9, 18, 22, 13),
  bdi_post = c(12,  0, 10,  9,  0, 11,  2,  4, 15, 10,
               11,  4, 19, 15,  3, 14, 10, 16, 10, 20,
               17, 25, 10, 22, 23, 10,  2, 10, 14,  7)
)

stopifnot(nrow(depression_bdi) == 30,
          all(table(depression_bdi$condition) == 10),
          # The fingerprints the ANCOVA help pages quote.
          isTRUE(all.equal(as.vector(tapply(depression_bdi$bdi_pre,
                                            depression_bdi$condition, mean)),
                           c(17, 17.7, 17.4))),
          sum(tapply(depression_bdi$bdi_pre, depression_bdi$condition,
                     function(x) sum((x - mean(x))^2))) == 752.5)

save(depression_bdi, file = "data/depression_bdi.rda",
     compress = "xz", version = 2)
