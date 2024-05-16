#' Plots the coefficient path of an SGP object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted SGP object
#'
#' @param x A object that was generated with sgp.
#' @param alpha Tuning parameter for the alpha-blending.
#' @param legend.pos Coordinates or keyword for positioning the legend.
#' @param label A Boolean value that specifies whether the plot should be annotated.
#' @param log.l A Boolean value that specifies whether the horizontal axis should be on the log scale.
#' @param norm A Boolean value that specifies whether the norm of each group should be plotted.
#' @param \dots Other parameters of underlying basic functions.
#'
#' @return A plot object with the coefficient path of an SGP.
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
#' plot(lin_fit, legend.pos = "topright", label = TRUE)
#' plot(lin_fit,  label = TRUE, norm = TRUE)
#'
#' log_fit <- sgp(X, y_log, g, type = "logit")
#' plot(log_fit, legend.pos = "topright", label = TRUE)
#' plot(log_fit, label = TRUE, norm = TRUE)
#'
#' @export
#'
plot.sgp <- function(x, alpha = 1, legend.pos, label = FALSE, log.l = FALSE, norm = FALSE, ...) {

  if (norm) {
    y <- predict.sgp(x, extract = "norm")
    y <- y[which(apply(abs(y), 1, sum) != 0), ]
    g <- 1:nrow(y)
  } else {
    coef <- x$beta[-1, , drop = FALSE]
    y <- coef[which(apply(abs(coef), 1, sum) != 0), , drop = FALSE]
    g <- as.integer(as.factor(x$group[which(apply(abs(coef), 1, sum) != 0)]))
  }

  x_axis <- rev(x$lambda)
  xlab <- expression(lambda)
  if (log.l) {
    x_axis <- log(x_axis)
    xlab <- expression(log(lambda))
  }

  plot_agr <- list(x = x_axis, y = 1:length(x_axis), ylim = range(y), xlab = xlab, ylab = "",
                   type = "n", xlim = range(x_axis), las = 1, bty = "L")
  user_agr <- list(...)
  if (length(user_agr)) {
    user_plot_agr <- user_agr[names(user_agr) %in% c(names(graphics::par()), names(formals(graphics::plot.default)))]
    plot_agr[names(user_plot_agr)] <- user_plot_agr
  }
  do.call("plot", plot_agr)
  if (plot_agr$ylab == "") {
    ylab <- if (norm) expression("||"*hat(beta)*"||") else expression(hat(beta))
    graphics::mtext(ylab, 2, 3.5, las = 1, adj = 0)
  }

  cols <- grDevices::hcl(h = seq(50, 410, len = max(4, max(g)+1)), l = 50, c = 200, alpha = alpha)
  cols <- cols[1:max(g)]
  line_agr <- list(col = cols, lwd = 1 + 2 * exp(-nrow(y) / 20), lty = 1, pch = "", x = rev(x_axis), y = t(y))
  if (length(user_agr)) line_agr[names(user_agr)] <- user_agr
  do.call("matlines", line_agr)

  if (label) {
    ypos <- y[, ncol(y)]
    graphics::text(-0.001, ypos, names(ypos), xpd = NA, adj = c(0, NA))
  }

  if(!missing(legend.pos)) {
    legend_agr <- list(col = cols, lwd = line_agr$lwd, lty = line_agr$lty, legend = levels(x$group))
    if (length(user_agr)) {
      user_legend_agr <- user_agr[names(user_agr) %in% names(formals(graphics::legend))]
      legend_agr[names(user_legend_agr)] <- user_legend_agr
    }
    legend_agr$x <- legend.pos
    do.call("legend", legend_agr)
  }

  graphics::abline(h = 0, lwd = 0.5, col = "black")
}

#' Plots the cross-validation curve from a SGP object
#'
#' Plots the cross-validation curve as a function of the lambda values used.
#'
#' @param x A object that was generated with sgp.cv.
#' @param log.l A Boolean value that specifies whether the horizontal axis should be on the log scale.
#' @param highlight A Boolean value that specifies whether a vertical line should be added at the value where the cross-validation error is minimized.
#' @param col Controls the color of the dots.
#' @param \dots Other parameters of underlying basic functions.
#'
#' @return A plot object with the cross-validation curve of an SGP.
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
#' plot(lin_fit, col = "blue")
#'
#' log_fit <- sgp.cv(X, y_log, g, type = "logit")
#' plot(log_fit, col = "blue")
#'
#' @export
#'
plot.sgp.cv <- function(x, log.l = TRUE, highlight = TRUE, col = "firebrick3", ...) {

  y <- x$cve
  ylab <- "Cross-validation error"

  x_axis <- rev(x$lambdas)
  xlab <- expression(lambda)
  if (log.l) {
    x_axis <- log(x_axis)
    xlab <- expression(log(lambda))
  }

  Low <- x$cve - x$cvse
  Up <- x$cve + x$cvse
  ind_1 <- is.finite(x_axis[1:length(x$cve)])
  ylim <- range(c(Low[ind_1], Up[ind_1]))
  ind_2 <- ((Up-Low) / diff(ylim) > 1e-3) & ind_1

  plot_args = list(x = x_axis[ind_1], y = y[ind_1], ylim = ylim, xlab = xlab, ylab = ylab,
                   type = "n", xlim = range(x_axis[ind_1]), las = 1, bty = "L")
  new_args = list(...)
  if (length(new_args)) plot_args[names(new_args)] = new_args

  do.call("plot", plot_args)
  suppressWarnings(graphics::arrows(x0 = x_axis[ind_2], x1 = x_axis[ind_2], y0 = Low[ind_2], y1 = Up[ind_2],
                                    code = 3, angle = 90, col = "gray60", length = .05))
  graphics::points(x_axis[ind_1], y[ind_1], col = col, pch = 20, cex = 1)
  n_s <- sapply(predict.sgp(x$fit, lambda = x$lambdas, extract = "groups"), length)
  graphics::axis(3, at = x_axis, labels = n_s, tick = FALSE, line = -0.5)
  graphics::mtext("Groups selected", cex = 0.8, line = 1.5)

  if (highlight) graphics::abline(v = x_axis[x$min], lty = 2, lwd = 1, col = col)
}
