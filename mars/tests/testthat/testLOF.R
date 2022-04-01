library(mars)
load("testLOF.RData")
test_that("LOF() returns the correct object", {
  expect_equal(LOF(y~.-1,dat, testmc), testLOF)
})
