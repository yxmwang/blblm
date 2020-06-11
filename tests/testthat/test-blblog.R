test_that("test blblog() fit the logistic regression model with bag of little bootstrap when using one CPU correctly", {
  bi <- sample(iris[1:100,], 1000, replace = TRUE)  #binomial iris
  suppressWarnings(fit <- blblog(Species ~ Sepal.Width + Sepal.Length, data = bi, m = 3, B = 100, parallel = FALSE))
  #suppress the warning of "not meaningful for factors"
  expect_s3_class(fit, "blblog")
})

test_that("test blblog() fit the logistic regression model with bag of little bootstrap when using parallelization correctly", {
  bi <- sample(iris[1:100,], 1000, replace = TRUE)  #binomial iris
  suppressWarnings(fit <- blblog(Species ~ Sepal.Width + Sepal.Length, data = bi, m = 3, B = 100, parallel = TRUE))
  #suppress the warning of "not meaningful for factors"
  expect_s3_class(fit, "blblog")
})