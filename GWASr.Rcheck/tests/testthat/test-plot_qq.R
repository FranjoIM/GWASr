library(testthat)
library(GWASr)

# Check if inherits is right for simulated QQ plot
test_that("test if f function is correct", {
  expect_s3_class(plot_qq(data.frame(P = ppoints(1000))), "ggplot")
}
)
