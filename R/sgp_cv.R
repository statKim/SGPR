#' Cross-validation for sparse group penalties
#'
#' A function that performs k-fold cross-validation for sparse group penalties for a lambda sequence.
#'
#' @param X The design matrix without intercept with the variables to be selected.
#' @param y The response vector.
#' @param group A vector indicating the group membership of each variable in X.
#' @param Z The design matrix of the variables to be included in the model without penalization.
#' @param nfolds The number of folds for cross-validation.
#' @param seed A seed provided by the user for the random number generator.
#' @param fold A vector of folds specified by the user (default is a random assignment).
#' @param type A string indicating the type of regression model (linear or binomial).
#' @param returnY A Boolean value indicating whether the fitted values should be returned.
#' @param print.trace A Boolean value that specifies whether the beginning of a fold should be printed.
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns A list containing:
#' \describe{
#' \item{cve}{The average cross-validation error for each value of lambda.}
#' \item{cvse}{The estimated standard error for each value of cve.}
#' \item{lambdas}{The sequence of lambda values.}
#' \item{fit}{The sparse group penalty model fitted to the entire data.}
#' \item{fold}{The fold assignments for each observation for the cross-validation procedure.}
#' \item{min}{The index of lambda corresponding to the minimum cross-validation error.}
#' \item{lambda.min}{The value of lambda with the minimum cross-validation error.}
#' \item{null.dev}{The deviance for the empty model.}
#' \item{pe}{The cross-validation prediction error for each value of lambda (for binomial only).}
#' \item{pred}{The fitted values from the cross-validation folds.}
#' }
#'
#' @examples
#' \donttest{
#' # Generate data
#'  n <- 100
#'  p <- 200
#'  nr <- 10
#'  g <- ceiling(1:p / nr)
#'  X <- matrix(rnorm(n * p), n, p)
#'  b <- c(-3:3)
#'  y_lin <- X[, 1:length(b)] %*% b + 5 * rnorm(n)
#'  y_log <- rbinom(n, 1, exp(y_lin) / (1 + exp(y_lin)))
#'
#' # Linear regression
#'  lin_fit <- sgp.cv(X, y_lin, g, type = "linear", penalty = "sgl")
#'  plot(lin_fit)
#'  predict(lin_fit, extract = "vars")
#'  lin_fit <- sgp.cv(X, y_lin, g, type = "linear", penalty = "sgs")
#'  plot(lin_fit)
#'  predict(lin_fit, extract = "vars")
#'  lin_fit <- sgp.cv(X, y_lin, g, type = "linear", penalty = "sgm")
#'  plot(lin_fit)
#'  predict(lin_fit, extract = "vars")
#'  lin_fit <- sgp.cv(X, y_lin, g, type = "linear", penalty = "sge")
#'  plot(lin_fit)
#'  predict(lin_fit, extract = "vars")
#'
#' # Logistic regression
#'  log_fit <- sgp.cv(X, y_log, g, type = "logit", penalty = "sgl")
#'  plot(log_fit)
#'  predict(log_fit, extract = "vars")
#'  log_fit <- sgp.cv(X, y_log, g, type = "logit", penalty = "sgs")
#'  plot(log_fit)
#'  predict(log_fit, extract = "vars")
#'  log_fit <- sgp.cv(X, y_log, g, type = "logit", penalty = "sgm")
#'  plot(log_fit)
#'  predict(log_fit, extract = "vars")
#'  log_fit <- sgp.cv(X, y_log, g, type = "logit", penalty = "sge")
#'  plot(log_fit)
#'  predict(log_fit, extract = "vars")
#' }
#'
#' @export
#'
sgp.cv <- function(X, y, group = 1:ncol(X), Z = NULL, ..., nfolds = 10, seed, fold, type,
                   returnY = FALSE, print.trace = FALSE) {

  # Fit a SGP
  tune_all <- list(...)
  tune_all$X <- X
  tune_all$y <- y
  tune_all$Z <- Z
  tune_all$group <- group
  tune_all$returnX <- TRUE
  tune_all$type <- type
  fit_all <- do.call("sgp", tune_all)

  # Extract dimensions
  X <- fit_all$X$X
  y <- fit_all$y
  n <- fit_all$n

  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit_all$XG <- NULL

  # Generate folds
  if (!missing(seed)) set.seed(seed)
  if (missing(fold)) {
    if (fit_all$type == "linear") {
      fold <- sample((1:n %% nfolds) + 1)
    } else {
      fold <- integer(n)
      events <- sum(y)
      fails <- n - events
      fold[y == 1] <- sample((1:events %% nfolds) + 1)
      fold[y == 0] <- sample(((events + 1:fails) %% nfolds) + 1)
    }
  } else {
    nfolds <- max(fold)
  }

  # Performe cross-validation
  Loss <- Pred <- matrix(NA, n, length(fit_all$lambda))
  if (fit_all$type == "logit") class <- Pred

  tune_fold <- list(...)
  tune_fold$group <- fit_all$group
  tune_fold$lambdas <- fit_all$lambdas
  tune_fold$warn <- FALSE
  tune_fold$type <- type

  for (i in 1:nfolds) {
    if (print.trace) cat("Fitting SGP in fold: ", i, sep = "", "\n")

    X_out <- X[fold == i, , drop = FALSE]
    y_out <- y[fold == i]
    Z_out <- Z[fold == i]

    tune_fold$X <- X[fold != i, , drop = FALSE]
    tune_fold$y <- y[fold != i]
    tune_fold$Z <- Z[fold != i]

    fit_fold <- suppressWarnings(do.call("sgp", tune_fold))

    intercept <- fit_fold$beta[1, ]
    beta <- fit_fold$beta[-1, , drop = FALSE]
    eta <- sweep(X_out %*% beta, 2, intercept, "+")

    if (fit_fold$type == "logit") {
      pred <- exp(eta)/(1 + exp(eta))
    } else {
      pred <- eta
    }

    loss <- get.loss(y_out, pred, fit_fold$type)
    pe <- if (fit_all$type == "logit") {(pred < 0.5) == y_out} else NULL
    cv_l_l <- length(fit_fold$lambdas)
    if(cv_l_l>ncol(Pred))cv_l_l <- ncol(Pred)

    Pred[fold == i, 1:cv_l_l] <- pred[, 1:cv_l_l]
    Loss[fold == i, 1:cv_l_l] <- loss[, 1:cv_l_l]
    if (fit_all$type == "logit") class[fold == i, 1:cv_l_l] <- pe[, 1:cv_l_l]
  }

  # Remove failed lambdas
  conv <- which(apply(is.finite(Pred), 2, all))
  Loss <- Loss[, conv, drop = FALSE]
  Pred <- Pred[, conv]
  lambdas <- fit_all$lambdas[conv]

  # Fit null model
  if (is.null(fit_all$Z)) {
    null_model <- stats::glm(y ~ 1, family = ifelse(fit_all$type == "linear", "gaussian", "binomial"))
  } else {
    null_model <- stats::glm(y ~ Z, family = ifelse(fit_all$type == "linear", "gaussian", "binomial"))
  }

  # Summarize
  cve <- apply(Loss, 2, mean)
  cvse <- apply(Loss, 2, stats::sd) / sqrt(n)
  min <- which.min(cve)
  null.dev <- mean(get.loss(y, stats::predict(null_model, type = "response"), fit_all$type))

  res <- list(cve = cve, cvse = cvse, lambdas = lambdas, fit = fit_all, fold = fold,
              min = min, lambda.min = lambdas[min], null.dev = null.dev)
  if (fit_all$type == "logit") res$pe <- apply(class, 2, mean)
  if (returnY) {
    if (fit_all$type == "linear") res$Pred <- Pred + attr(y, "mean")
    else res$Pred <- Pred
  }
  structure(res, class = "sgp.cv")
}
