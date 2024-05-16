#' Set up a lambda sequence
#'
#' A function that sets up a lambda sequence for a sparse group penalty.
#'
#' @param X The design matrix without intercept with the variables to be selected.
#' @param y The response vector.
#' @param group A vector indicating the group membership of each variable in X.
#' @param Z The design matrix of the variables to be included in the model without penalization.
#' @param type A string indicating the type of regression model (linear or binomial).
#' @param alpha Tuning parameter for the mixture of penalties at group and variable level.
#' A value of 0 results in a selection at group level, a value of 1
#' results in a selection at variable level and everything in between
#' is bi-level selection.
#' @param lambda.min An integer multiplied by the maximum lambda to define the end of the lambda sequence.
#' @param log.lambda A Boolean value that specifies whether the values of the lambda
#' sequence should be on the log scale.
#' @param nlambda An integer that specifies the length of the lambda sequence.
#' @param group.weight A vector specifying weights that are multiplied by the group
#' penalty to account for different group sizes.
#' @param ada_mult An integer that defines the multiplier for adjusting the convergence threshold.
#'
#' @returns A vector with values for lambda.
#'
process.lambda <- function(X, y, group, Z, type, alpha, lambda.min, log.lambda, nlambda, group.weight, ada_mult) {

  # Validation and correction of lambda related input
  if (nlambda < 2){
    warning("nlambda must be at least 2 and was set to its default value.")
    nlambda <- 100
  }
  if (ada_mult < 1){
    warning("ada_mult must be >= 1 and was set to its default value.")
    ada_mult <- 2
  }

  # GLM with unpenalized variables
  n <- length(y)
  if (is.null(Z)) {
    fit <- stats::glm(y ~ 1, family = ifelse(type == "linear", "gaussian", "binomial"))
  } else {
    fit <- stats::glm(y ~ Z, family = ifelse(type == "linear", "gaussian", "binomial"))
  }

  # Determine lambda.max
  if (type == "linear") {
    r <- fit$residuals
  } else {
    w <- fit$weights
    w <- w + 3 * w * (1 - alpha) * alpha
    if (max(w) < 1e-4) stop("Intercept or variables in Z result in a saturated model: no residuals left for selection", call. = FALSE)
    r <- stats::residuals(fit, "working") * w
  }

  lambda.max <- max_cor(X, r, c(0, cumsum(table(group))), as.double(group.weight), alpha)

  # Generate lambda sequence
  if (log.lambda) {
    if (lambda.min == 0) {
      lambdas <- c(exp(seq(log(lambda.max), log(.001 * lambda.max), length = nlambda - 1)), 0)
    } else {
      lambdas <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), length = nlambda))
    }
  } else {
    if (lambda.min == 0) {
      lambdas <- c(seq(lambda.max, 0.001 * lambda.max, length = nlambda - 1), 0)
    } else {
      lambdas <- seq(lambda.max, lambda.min * lambda.max, length = nlambda)
    }
  }
  lambdas
}
