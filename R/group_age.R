#' Assign age to group
#' 
#' @param age Vector of ages
#' @param breaks Ages that define strata
#' 
#' @return A factor that groups the ages into the age groups defined by `breaks`.
#' @export
#' @importFrom stats reorder

group_age <- function(age, breaks) {
  res <- as.character(cut(age, breaks, right = FALSE,
                          labels = paste(breaks[-length(breaks)], c(breaks[-1]-1), sep="-")))
  reorder(res, age, min)
}
