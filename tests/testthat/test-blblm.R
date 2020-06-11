test_that("test blblm() fit the linear regression model with bag of little bootstrap when using one CPU correctly", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
  expect_s3_class(fit, "blblm")

  co <- coef(fit)
  expect_length(co, 4)
})

test_that("test blblm() fit the linear regression model with bag of little bootstrap when using parallelization correctly", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)
  expect_s3_class(fit, "blblm")

  co <- coef(fit)
  expect_length(co, 4)
})