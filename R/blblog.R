#' @import purrr
#' @import furrr
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Logistic Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Compute logistic regression model with Bag of Little Bootstraps
#'
#' Find logistic regression model using subsamples and bootstrap samples with parallelization.
#' The user can specify the number of workers by plan(multisession, workers = 4) if they choose to run in parallezation.
#' For example, when they write parallel = TRUE, they need can change the number of clusters as they needed.
#'
#' @param formula logistic regression formula
#' @param data data frame
#' @param m integer indicating the number of subsets
#' @param B integer indicating the number of bootstrap
#' @param parallel boolean which contains only two value, TRUE or FALSE
#'
#' @return fitted model by using blb and parallezation
#' @export
blblog <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (parallel == TRUE) {
    mpa_func <- future_map
  } else {
    mpa_func <- map
  }
  estimates <- mpa_func(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
  )
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblog"
  invisible(res)
}

#' split data
#'
#' split data into m parts of approximated equal sizes
#'
#' @param data data frame
#' @param m numeric
#'
#' @return splitted dataframe
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates in each subsample
#'
#' Compute the regression estimate in each subsample
#' @param formula logistic regression formula
#' @param data data frame
#' @param n numeric
#' @param B numeric
#'
#' @return the regression estimate in each subsample
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the estimates for each boot
#'
#' Compute the regression estimate in each boot
#'
#' @param formula logistic regression formula
#' @param data data frame
#' @param n integer
#'
#' @return  the regression estimate for each blb dataset
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, n)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula logistic regression formula
#' @param data data frame
#' @param freqs number of repetitions
#' @param n integer
#'
#' @return the regression estimates
glm1 <- function(formula, data, freqs, n) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- suppressWarnings(glm(formula, family = binomial, data, weights = freqs, maxit = 100))
  while(!fit$converged){
    freqs <- rmultinom(1, n, rep(1, nrow(data)))
    fit <- suppressWarnings(glm(formula, family= binomial, data, weights = freqs, maxit=100))
  }
  list(coef = blblogcoef(fit), sigma = blblogsigma(fit))
}

#' compute the coefficients from object
#'
#' @param fit fitted model from blblog
#' @return the coefficients
blblogcoef <- function(fit) {
  coef(fit)
}

#' compute sigma from object
#'
#' @param fit fitted model from blblog
#' @return the standard errors
blblogsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the result
#'
#' @param x input
#' @param ... additional arguements/inputs
#'
#' @export
#' @return the printed result
#' @method print blblog
print.blblog <- function(x, ...) {
  cat("blblog model:", capture.output(x$formula))
  cat("\n")
}


#' find the standard eror for each blb
#'
#' @param object fitted model from blblog
#' @param confidence boolean
#' @param level numeric, the confidence interval level
#' @param ... additional arguements/inputs
#'
#' @return the standard error
#' @export
#' @method sigma blblog
sigma.blblog <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' find the coefficient in blblog
#'
#' @param object fitted model from blblog
#' @param ... additional arguements/inputs
#'
#' @return the estimated coefficient
#' @export
#' @method coef blblog
coef.blblog <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' find the confidence interval in blblog
#'
#' @param object fitted model from blbog
#' @param parm numerical value indicating the parameters you want to construct the ci for
#' @param level numeric, the confidence interval level
#' @param ... additional arguements/inputs
#'
#' @return the confidence interval of the estimates
#' @export
#' @method confint blblog
confint.blblog <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' the prediction interval in blblog
#'
#' @param object fitted model from blblog
#' @param new_data a new data frame
#' @param confidence boolean
#' @param level numeric, the confidence interval level
#' @param ... additional arguements/inputs
#'
#' @return the confidence interval of the estimates
#' @export
#' @method predict blblog
predict.blblog <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>% rowMeans())
  }
}

invlogit <- function(x) {
  1 / (1 + exp(-x))
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(object = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
