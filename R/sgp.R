#' Fit a sparse group regularized regression path
#'
#' A function that determines the regularization paths for models with
#' sparse group penalties at a grid of values for the regularization parameter lambda.
#'
#' Two options are available for choosing a penalty. With the argument \code{penalty},
#' the methods Sparse Group LASSO, Sparse Group SCAD, Sparse Group MCP and Sparse Group EP
#' can be selected with the abbreviations \code{sgl}, \code{sgs}, \code{sgm} and \code{sge}.
#' Alternatively, penalties can be combined additively with the arguments \code{pvar}
#' and \code{pgr}, where \code{pvar} is the penalty applied at the variable level and
#' \code{pgr} is the penalty applied at the group level. The options are \code{lasso},
#' \code{scad}, \code{mcp} and \code{exp} for Least Absolute Shrinkage and Selection Operator,
#' Smoothly Clipped Absolute Deviation, Minimax Concave Penalty and Exponential Penalty.
#'
#' @param X The design matrix without intercept with the variables to be selected.
#' @param y The response vector.
#' @param group A vector indicating the group membership of each variable in X.
#' @param penalty A string that specifies the sparse group penalty to be used.
#' @param alpha Tuning parameter for the mixture of penalties at group and variable level.
#' A value of 0 results in a selection at group level, a value of 1
#' results in a selection at variable level and everything in between
#' is bi-level selection.
#' @param type A string indicating the type of regression model (linear or binomial).
#' @param Z The design matrix of the variables to be included in the model without penalization.
#' @param nlambda An integer that specifies the length of the lambda sequence.
#' @param lambda.min An integer multiplied by the maximum lambda to define the end of the lambda sequence.
#' @param log.lambda A Boolean value that specifies whether the values of the lambda
#' sequence should be on the log scale.
#' @param lambdas A user supplied vector with values for lambda.
#' @param prec The convergence threshold for the algorithm.
#' @param ada_mult An integer that defines the multiplier for adjusting the convergence threshold.
#' @param max.iter The convergence threshold for the algorithm.
#' @param standardize An integer that defines the multiplier for adjusting the convergence threshold.
#' @param vargamma An integer that defines the value of gamma for the penalty at the variable level.
#' @param grgamma An integer that specifies the value of gamma for the penalty at the group level.
#' @param vartau An integer that defines the value of tau for the penalty at the variable level.
#' @param grtau An integer that specifies the value of tau for the penalty at the group level.
#' @param pvar A string that specifies the penalty used at the variable level.
#' @param pgr A string that specifies the penalty used at the group level.
#' @param group.weight A vector specifying weights that are multiplied by the group
#' penalty to account for different group sizes.
#' @param returnX A Boolean value that specifies whether standardized design matrix should be returned.
#' @param \dots Other parameters of underlying basic functions.
#'
#' @returns A list containing:
#' \describe{
#' \item{beta}{A vector with estimated coefficients.}
#' \item{type}{A string indicating the type of regression model (linear or binomial).}
#' \item{group}{A vector indicating the group membership of the individual variables in X.}
#' \item{lambdas}{The sequence of lambda values.}
#' \item{alpha}{Tuning parameter for the mixture of penalties at group and variable level.}
#' \item{loss}{A vector containing either the residual sum of squares (linear) or the negative log-likelihood (binomial).}
#' \item{prec}{The convergence threshold used for each lambda.}
#' \item{n}{Number of observations.}
#' \item{penalty}{A string indicating the sparse group penalty used.}
#' \item{df}{A vector of pseudo degrees of freedom for each lambda.}
#' \item{iter}{A vector of the number of iterations for each lambda.}
#' \item{group.weight}{A vector of weights multiplied by the group penalty.}
#' \item{y}{The response vector.}
#' \item{X}{The design matrix without intercept.}
#' }
#'
#' @references
#' \itemize{
#' \item Buch, G., Schulz, A., Schmidtmann, I., Strauch, K., and Wild, P. S. (2024)
#' Sparse Group Penalties for bi-level variable selection. Biometrical Journal, 66, 2200334.
#' \doi{https://doi.org/10.1002/bimj.202200334}
#'
#' \item Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011)
#' A Sparse-Group Lasso. Journal of computational and graphical statistics, 22(2), 231-245.
#' \doi{https://doi.org/10.1080/10618600.2012.681250}
#'
#' \item  Breheny, P., and Huang J. (2009)
#' Penalized methods for bi-level variable selection. Statistics and its interface, 2: 369-380.
#' \doi{10.4310/sii.2009.v2.n3.a10}
#' }
#'
#' @examples
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
#'  lin_fit <- sgp(X, y_lin, g, type = "linear", penalty = "sgl")
#'  plot(lin_fit)
#'  lin_fit <- sgp(X, y_lin, g, type = "linear", penalty = "sgs")
#'  plot(lin_fit)
#'  lin_fit <- sgp(X, y_lin, g, type = "linear", penalty = "sgm")
#'  plot(lin_fit)
#'  lin_fit <- sgp(X, y_lin, g, type = "linear", penalty = "sge")
#'  plot(lin_fit)
#'
#' # Logistic regression
#'  log_fit <- sgp(X, y_log, g, type = "logit", penalty = "sgl")
#'  plot(log_fit)
#'  log_fit <- sgp(X, y_log, g, type = "logit", penalty = "sgs")
#'  plot(log_fit)
#'  log_fit <- sgp(X, y_log, g, type = "logit", penalty = "sgm")
#'  plot(log_fit)
#'  log_fit <- sgp(X, y_log, g, type = "logit", penalty = "sge")
#'  plot(log_fit)
#'
#' @export
#'
sgp <- function(X, y, group = 1:ncol(X), penalty = c("sgl", "sgs", "sgm", "sge"), alpha = 1/3, type = c("linear", "logit"), Z = NULL,
                nlambda = 100, lambda.min = {if(nrow(X) > ncol(X)) 1e-4 else .05}, log.lambda = TRUE, lambdas,
                prec = 1e-4, ada_mult = 2, max.iter = 10000, standardize = TRUE, intercept = TRUE,
                vargamma = ifelse(pvar == "scad"|penalty == "sgs", 4, 3), grgamma = ifelse(pgr == "scad"|penalty == "sgs", 4, 3), vartau = 1, grtau = 1,
                pvar = c("lasso", "scad", "mcp", "exp"), pgr = c("lasso", "scad", "mcp", "exp"),
                group.weight = rep(1, length(unique(group))), returnX = FALSE, ...) {

  type <- match.arg(type)
  if(!is.null(penalty)) penalty <- match.arg(penalty)
  if(!is.null(pvar))    pvar <- match.arg(pvar)
  if(!is.null(pgr))     pgr <- match.arg(pgr)

  # Validation and correction
  sgp <- process.penalty(penalty, pvar, pgr, vargamma, grgamma, vartau, grtau, alpha)

  response <- process.y(y, type)
  grouping <- process.group(group, group.weight)
  predictors <- process.X(X, group)
  covariates <- process.Z(Z)

  n <- length(response)
  if (nrow(predictors$X) != n) stop("Dimensions of X is not compatible with y", call. = FALSE)
  if (!is.null(Z) && (nrow(covariates$Z) != n)) stop("Dimensions of Z is not compatible with y", call. = FALSE)

  # Reorder processed dimensions
  to_order <- F
  if (any(order(grouping) != 1:length(grouping))) {
    to_order <- T
    neworder <- order(group)
    grouping <- grouping[neworder]
    predictors <- predictors[, neworder]
    group.weight <- group.weight[neworder]
  }

  # Combine processed dimensions
  if (is.null(Z)) {
    dat <- predictors
    Z_groups <- grouping
  } else {
    dat <- list(X = cbind(covariates$Z, predictors$X),
                vars = c(covariates$vars, predictors$vars),
                center = c(covariates$center, predictors$center),
                scale = c(covariates$scale, predictors$scale))
    Z_groups <- c(rep(0, ifelse(is.null(covariates), 0, ncol(covariates))), grouping)
  }

  if (!standardize) {
    dat <- list(X = X,
                vars = predictors$vars)
  }
                  
  p <- ncol(dat$X)

  # Setup lambdas
  if (missing(lambdas)) {
    if (is.null(covariates)) {
      lambdas <- process.lambda(dat$X, response, grouping, NULL, type, alpha, lambda.min, log.lambda, nlambda, group.weight, ada_mult)
    } else{lambdas <- process.lambda(dat$X, response, grouping, covariates$Z, type, alpha, lambda.min, log.lambda, nlambda, group.weight, ada_mult)}

    #lam.max <- lambdas[1]
    own_l <- FALSE
  } else {
    #lam.max <- -1
    nlambda <- length(lambdas)
    own_l <- TRUE
  }

  # Indices for groupings
  J <- as.integer(table(Z_groups))
  J0 <- as.integer(if (min(Z_groups) == 0) J[1] else 0)
  JG <- as.integer(if (min(Z_groups) == 0) cumsum(J) else c(0, cumsum(J)))
  if (J0) {
    lambdas[1] <- lambdas[1] + 1e-5
    own_l <- TRUE
  }

  # Call C++
  if (type == "linear") {
    fit <- lcdfit_linear(dat$X, response, JG, J0, lambdas, alpha, prec, ada_mult, grgamma, vargamma, grtau, vartau,
                         as.integer(max.iter), group.weight, as.integer(own_l), as.integer(sgp[[1]]), as.integer(sgp[[2]]))
    b <- rbind(mean(y), matrix(fit[[1]], nrow = p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1
    loss <- fit[[4]]
    ada.prec <- fit[[5]]
  } else {
    fit <- lcdfit_logistic(dat$X, response, JG, J0, lambdas, alpha, prec, ada_mult, grgamma, vargamma, grtau, vartau,
                           as.integer(max.iter), group.weight, as.integer(own_l), as.integer(sgp[[1]]), as.integer(sgp[[2]]))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow = p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
    ada.prec <- fit[[6]]
  }

  # Remove failed lambdas
  conv <- !is.na(iter)
  b <- b[, conv, drop = FALSE]
  iter <- iter[conv]
  lambdas <- lambdas[conv]
  df <- df[conv]
  loss <- loss[conv]

  if (iter[1] == max.iter) stop("The algorithm could not converge already at the first lambda. Try a higher value for alpha or another penalty", call. = FALSE)

  # Unstandardize
  if (to_order) b[-1, ] <- b[1 + neworder, ]

  if(standardize){
    beta <- matrix(0, nrow = 1 + p, ncol = ncol(b))
    beta[-1, ] <- b[-1, ] / dat$scale
    beta[1, ] <- b[1, ] - dat$center %*% beta[-1, , drop = FALSE]
  } else {
    beta <- b
  }
                   
  # Labeling
  if (intercept) { 
    varnames <- c("Intercept", dat$vars)
  } else {
    varnames <-  dat$vars
  }
  

  #dimnames(beta) <- list(varnames, round(lambdas, digits = 4))

  res <- structure(list(beta = beta,
                        type = type,
                        group = factor(group),
                        lambdas = lambdas,
                        alpha = alpha,
                        loss = loss,
                        prec = ada.prec,
                        n = n,
                        penalty = penalty,
                        df = df,
                        iter = iter,
                        group.weight = dat$m),
                   class = "sgp")

  if (type == 'linear') {
    res$y <- response + attr(response, 'mean')
  } else {
    res$y <- response
  }
  if (returnX) res$X <- dat

  return(res)
}
