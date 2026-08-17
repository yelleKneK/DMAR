# .dmar_palette() is the internal color source for the plot_* family. It
# only forwards to base R's palette.colors(), so these tests pin that it
# stays a thin, neutral, colorblind-safe wrapper with no palette of its own.

test_that(".dmar_palette() is neutral by default (base R Okabe-Ito)", {
  expect_equal(.dmar_palette(1), "#000000")
  expect_equal(.dmar_palette(2), c("#000000", "#E69F00"))
  expect_equal(.dmar_palette(4),
               c("#000000", "#E69F00", "#56B4E9", "#009E73"))
  # base R's palette exactly, unnamed, not a private copy
  expect_equal(.dmar_palette(3),
               unname(grDevices::palette.colors(3, "Okabe-Ito")))
  expect_null(names(.dmar_palette(3)))
})

test_that(".dmar_palette() with NULL returns the canonical anchor sets", {
  expect_length(.dmar_palette(), 9L)                     # nine Okabe-Ito
  expect_length(.dmar_palette(palette = "tableau"), 10L) # ten Tableau
})

test_that(".dmar_palette() interpolates beyond the anchors", {
  cols <- .dmar_palette(20)
  expect_length(cols, 20L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

test_that(".dmar_palette() selects the tableau alternative", {
  expect_equal(.dmar_palette(2, "tableau"),
               unname(grDevices::palette.colors(2, "Tableau 10")))
  # "okabe_ito" is the default, so naming it explicitly agrees.
  expect_equal(.dmar_palette(3, "okabe_ito"), .dmar_palette(3))
})

test_that(".dmar_palette() honors reverse and the okabe-ito alias", {
  expect_equal(.dmar_palette(3, reverse = TRUE), rev(.dmar_palette(3)))
  expect_equal(.dmar_palette(reverse = TRUE), rev(.dmar_palette()))
  expect_equal(.dmar_palette(3, palette = "okabe-ito"), .dmar_palette(3))
})

test_that(".dmar_palette() rejects bad input and unknown palettes", {
  expect_error(.dmar_palette(0), "positive integer")
  expect_error(.dmar_palette(2.5), "positive integer")
  expect_error(.dmar_palette(-1), "positive integer")
  expect_error(.dmar_palette(c(1, 2)), "positive integer")
  expect_error(.dmar_palette(NA), "positive integer")
  expect_error(.dmar_palette(2, reverse = NA), "TRUE or FALSE")
  # the Notre Dame palette is gone: it is now simply an unknown name
  expect_error(.dmar_palette(2, palette = "dmar_ND"), "Unknown 'palette'")
  expect_error(.dmar_palette(2, palette = "viridis"), "Unknown 'palette'")
})
