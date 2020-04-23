#' Assign age to group
#' @export
#' @importFrom stats reorder

group_age <- function(age, breaks) {
  res <- as.character(cut(age, breaks, right = FALSE,
                          labels = paste(breaks[-length(breaks)], c(breaks[-1]-1), sep="-")))
  reorder(res, age, min)
}
