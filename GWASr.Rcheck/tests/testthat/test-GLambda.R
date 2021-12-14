library(testthat)
library(GWASr)

# Test if GLambda works properly.

test_that("test if f function is correct", {
  expect_equal(round(GLambda(c(0.0001,0.0002,0.0003))), 30)
}
)
