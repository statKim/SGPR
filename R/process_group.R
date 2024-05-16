#' Process groupings for a sparse group penalty
#'
#' A function that checks the group information for possible errors and processes it.
#'
#' @param group A vector that specifies the group membership of each variable in X.
#' @param group.weight A vector specifying weights that are multiplied by the group penalty to account for different group sizes.
#'
#' @returns A structure containing the prepared group structure and, as an attribute, its labels and group weights.
#'
process.group <- function(group, group.weight) {

  # Validation and correction of group and group.weight
  gf <- factor(group)
  label <- levels(gf)
  if (is.numeric(group) | is.integer(group)) {
    label <- paste0("Group ", label)
  }
  group <- as.integer(gf)

  if (any(levels(gf) == 'Group 0')) {
    stop('The group name "0" is not allowed: use Z to fit a model with unpenalized variables', call. = FALSE)
  }
  if (any(levels(gf) == 'Unpenalized')) {
    stop('The group name "Unpenalized" is not allowed: the name is reserved for the variables in Z', call. = FALSE)
  }

  if (missing(group.weight)) {
    group.weight <- rep(1, length(label))
    names(group.weight) <- label
  } else {
    if (any(group.weight < 0)) stop('A negative group.weight is not allowed', call. = FALSE)
    if (length(group.weight) != length(label)) stop("Groups and group.weight must be the same size", call. = FALSE)
    if (storage.mode(group.weight) != "double") storage.mode(group.weight) <- "double"
  }
  return(structure(group, levels = label, group.weight = group.weight))
}
