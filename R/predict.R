#' Predictions based on a SGP model
#'
#' A function that extracts information from a SGP object and performs predictions.
#'
#' @param object A object that was generated with sgp.
#' @param X The design matrix for making predictions.
#' @param extract A string indicating the type of information to return.
#' @param lambda The value of lambda at which predictions should be made.
#' @param index The index that indicates the lambda at which predictions should be made (alternative to specifying 'lambda').
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns Different objects depending on the sting indicated by 'extract'.
#'
#' @examples
#' n <- 100
#' p <- 12
#' nr <- 4
#' g <- paste0("Group ",ceiling(1:p / nr))
#' X <- matrix(rnorm(n * p), n, p)
#' b <- c(-3:3)
#' y_lin <- X[, 1:length(b)] %*% b + 5 * rnorm(n)
#' y_log <- rbinom(n, 1, exp(y_lin) / (1 + exp(y_lin)))
#'
#' lin_fit <- sgp(X, y_lin, g, type = "linear")
#' predict(lin_fit, X = X, extract = "nvars")
#'
#' log_fit <- sgp(X, y_log, g, type = "logit")
#' predict(log_fit, X = X, extract = "nvars")
#'
#' @export
#'
predict.sgp <- function(object, X = NULL, extract = c("link", "response", "class", "coef", "vars", "groups", "nvars", "ngroups", "norm"),
                        lambda, index = 1:length(object$lambda), ...) {

  if (extract == "class" & object$type != "logit") stop("'class' is only applicable for a logistic model", call. = FALSE)
  if (extract == "response" & object$type != "logit") stop("'response' is only applicable for a logistic model", call. = FALSE)

  extract <- match.arg(extract)
  coef <- coef.sgp(object, lambda = lambda, index = index, drop = FALSE)

  if (extract == "coef") return(coef)

  intercept <- coef[1,]
  coef <- coef[-1, , drop = FALSE]

  if (extract == "vars") return(drop(apply(coef != 0, 2, FUN = which)))

  if (extract == "groups") return(drop(apply(coef != 0, 2, function(x) unique(object$group[x]))))

  if (extract == "nvars") {
    nonzero <- drop(apply(coef != 0, 2, FUN = which))
    if (is.list(nonzero)) {return(sapply(nonzero, length))}
    else {return(length(nonzero))}
  }

  if (extract == "ngroups") {
    nonzero <- drop(apply(coef != 0, 2, function(x) unique(object$group[x])))
    if (is.list(nonzero)) {return(sapply(nonzero, length))}
    else {return(length(nonzero))}
  }

  if (extract == "norm") return(drop(apply(coef, 2, function(x) tapply(x, object$group, function(x){sqrt(sum(x^2))}))))

  if (missing(X) | is.null(X)) stop("Please provide X", call. = FALSE)

  eta <- sweep(X %*% coef, 2, intercept, "+")
  if (object$type == "linear" | extract == "link") return(drop(eta))

  resp <- exp(eta) / (1+exp(eta))
  if (extract == "response") return(drop(resp))
  if (extract == "class") return(drop(1 * (eta > 0)))
}

#' Coefficients from an SGP model
#'
#' A function that extracts the estimated coefficients from an SGP object.
#'
#' @param object A object that was generated with sgp.
#' @param lambda The value of lambda at which the coefficients are to be extracted.
#' @param index The index that indicates the lambda at which the coefficients are to be extracted (alternative to specifying 'lambda').
#' @param drop A Boolean value that specifies whether empty dimensions should be removed.
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns A vector or matrix with the estimated coefficients.
#'
#' @examples
#' n <- 100
#' p <- 12
#' nr <- 4
#' g <- paste0("Group ",ceiling(1:p / nr))
#' X <- matrix(rnorm(n * p), n, p)
#' b <- c(-3:3)
#' y_lin <- X[, 1:length(b)] %*% b + 5 * rnorm(n)
#' y_log <- rbinom(n, 1, exp(y_lin) / (1 + exp(y_lin)))
#'
#' lin_fit <- sgp(X, y_lin, g, type = "linear")
#' coef(lin_fit, index = 5:7)
#'
#' log_fit <- sgp(X, y_log, g, type = "logit")
#' coef(log_fit, index = 5:7)
#'
#' @export
#'
coef.sgp <- function(object, lambda, index = 1:length(object$lambda), drop = TRUE, ...) {
  if (!missing(lambda)) {
    if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) stop('Please enter a lambda that lies within the range of the coefficient path', call. = FALSE)
    pos <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    w <- pos %% 1
    coef <- object$beta[, ceiling(pos), drop = FALSE] * w + object$beta[, floor(pos), drop = FALSE] * (1 - w)
    colnames(coef) <- round(lambda, 4)
  } else {
    coef <- object$beta[, index, drop = FALSE]
  }
  if (drop) return(drop(coef)) else return(coef)
}

#' Predictions based on a SGP models
#'
#' A function that extracts information from a cross-validated SGP object and performs predictions.
#'
#' @param object A object that was generated with sgp.cv.
#' @param X The design matrix for making predictions.
#' @param extract A string indicating the type of information to return.
#' @param lambda The value of lambda at which predictions should be made.
#' @param index The index that indicates the lambda at which predictions should be made (alternative to specifying 'lambda').
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns Different objects depending on the sting indicated by 'extract'.
#'
#' @examples
#' n <- 100
#' p <- 12
#' nr <- 4
#' g <- paste0("Group ",ceiling(1:p / nr))
#' X <- matrix(rnorm(n * p), n, p)
#' b <- c(-3:3)
#' y_lin <- X[, 1:length(b)] %*% b + 5 * rnorm(n)
#' y_log <- rbinom(n, 1, exp(y_lin) / (1 + exp(y_lin)))
#'
#' lin_fit <- sgp.cv(X, y_lin, g, type = "linear")
#' predict(lin_fit, X = X, extract = "link")
#'
#' log_fit <- sgp.cv(X, y_log, g, type = "logit")
#' predict(log_fit, X = X, extract = "class")
#'
#' @export
#'
predict.sgp.cv <- function(object, X, lambda = object$lambda.min, index = object$min,
                           extract = c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"), ...) {
  extract <- match.arg(extract)
  return(predict.sgp(object$fit, X = X, lambda = lambda, index = index, extract = extract, ...))
}

#' Coefficients from SGP models
#'
#' A function that extracts the estimated coefficients from an cross-validated SGP object.
#'
#' @param object A object that was generated with sgp.cv.
#' @param lambda The value of lambda at which the coefficients are to be extracted.
#' @param index The index that indicates the lambda at which the coefficients are to be extracted (alternative to specifying 'lambda').
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns A vector or matrix with the estimated coefficients.
#'
#' @examples
#' n <- 100
#' p <- 12
#' nr <- 4
#' g <- paste0("Group ",ceiling(1:p / nr))
#' X <- matrix(rnorm(n * p), n, p)
#' b <- c(-3:3)
#' y_lin <- X[, 1:length(b)] %*% b + 5 * rnorm(n)
#' y_log <- rbinom(n, 1, exp(y_lin) / (1 + exp(y_lin)))
#'
#' lin_fit <- sgp.cv(X, y_lin, g, type = "linear")
#' coef(lin_fit)
#'
#' log_fit <- sgp.cv(X, y_log, g, type = "logit")
#' coef(log_fit)
#'
#' @export
#'
coef.sgp.cv <- function(object, lambda = object$lambda.min, index = object$min, ...) {
  coef.sgp(object$fit, lambda = lambda, index = index, ...)
}
