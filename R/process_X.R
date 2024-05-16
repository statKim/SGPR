#' Process X for a sparse group penalty
#'
#' A function that checks the design matrix X for possible errors and scales it.
#'
#' @param X The design matrix without intercept with the variables to be selected.
#' @param group A vector that specifies the group membership of each variable in X.
#'
#' @returns A list containing:
#' \describe{
#' \item{X}{The standardized design matrix X.}
#' \item{vars}{The variable names of the matrix.}
#' \item{center}{The center of the variables before the transformation.}
#' \item{scale}{The scale of the variables before the transformation.}
#' }
#'
process.X <- function(X, group) {

  # Validation and correction of X
  if(anyNA(X)) stop("Missing data detected in X. Please remove or impute cases with NA's.", call. = FALSE)
  if(!is.matrix(X)) {
    tmp <- try(X <- stats::model.matrix(~ 0 + ., X), silent = TRUE)
    if(inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call. = FALSE)
  }
  if ( dim(X)[2] <= 1) stop("X must be a matrix with more than one column")
  if(mode(X) == "integer") mode(X) <- "double"

  if(length(group) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)

  vars <- if(is.null(colnames(X))) paste0("Variable ", 1:ncol(X)) else colnames(X)

  Xsd <- scale(X)
  center <- attributes(Xsd)$'scaled:center'
  scale <- attributes(Xsd)$'scaled:scale'
  if (length(which(scale > 1e-6)) != ncol(X)) {
    stop ("Please remove constants (scale < 1e-6) from X.", call.=FALSE)
  }

  return(list(X = Xsd, vars = vars, center = center, scale = scale))
}
