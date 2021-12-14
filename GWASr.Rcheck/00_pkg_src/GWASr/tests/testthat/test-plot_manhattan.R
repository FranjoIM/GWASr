library(testthat)
library(GWASr)

# Expect error if DF is missing in plot_manhattan.

test_that("test if f function is correct", {
  expect_error(plot_manhattan("Test"))
}
)
