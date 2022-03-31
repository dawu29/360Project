# testmars.R
library(testthat)
load("testmars.RData")
test_that("mars() returns the correct object", {
  expect_equal(mars(y~., marstestdata, testmc), testmars, ignore_attr=TRUE)
})
