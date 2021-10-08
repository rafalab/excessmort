#' Callapse age groups into broader ones
#'
#' Collapse a count or demographic data frame into a broader age strata. 
#' 
#' @param demo A data frame with population sizes for different groups that will be collapsed.
#' @param counts A data frame with counts and population sizes for different groups that will be collapsed
#' @param breaks The new age breaks for the new, broader, age strata.
#' 
#' @return A age groups represented in `demo` or `counts` are grouped using the new age breaks
#' defined by breaks but containing the populations and counts, if applicable,
#' for age groups defined by `breaks`.
#' @examples
#' 
#' library(lubridate)
#' data(cook_records)
#' ## define smaller subset for example
#' cook_demographics_subset <- cook_demographics[year(cook_demographics$date)==2021, ]
#' demo <- collapse_age_dist(cook_demographics_subset, 
#'                           breaks = c(0, 20, 40, 60, 80, Inf))
#' @export
#' @importFrom stats reorder
#' @import dplyr
#' @importFrom tidyr separate drop_na

collapse_age_dist <- function(demo, breaks){
  ## assumes groups are start-end
  ## assumes agegroup contains group
  ## assumes column named population has population values

  demo <- separate(demo, agegroup, c("start", "end"), sep = "-", remove = FALSE, convert = TRUE)

  edges <- c(unique(demo$start), Inf)

  if(!all(breaks %in% edges)) stop("Age breaks must be a subset of", paste(edges, collapse = ","))

  g <- as.character(cut(demo$start, breaks, right = FALSE,
                        labels = paste(breaks[-length(breaks)], c(breaks[-1]-1), sep="-")))

  vars <- c(setdiff(names(demo), c("start", "end", "population")))

  return(demo %>% mutate(agegroup = g) %>%
           group_by_at(vars) %>%
           summarize(population = sum(population)) %>%
           ungroup() %>%
           separate(agegroup, c("start", "end"), sep = "-", remove = FALSE, convert = TRUE) %>%
           mutate(agegroup = reorder(agegroup, start, first)) %>%
           select(-start, -end) %>%
           arrange(agegroup) %>%
           drop_na(agegroup))
}




#' @rdname collapse_age_dist
#' @export


collapse_counts_by_age <- function(counts, breaks){
  ## assumes groups are start-end
  ## assumes agegroup contains group
  ## assumes column named population has population values
  ## assumes column names outcome has outcome values

  counts <- separate(counts, agegroup, c("start", "end"), sep = "-", remove = FALSE, convert = TRUE)

  edges <- c(unique(counts$start), Inf)

  if(!all(breaks %in% edges)) stop("Age breaks must be a subset of", paste(edges, collapse = ","))

  g <- as.character(cut(counts$start, breaks, right = FALSE,
                        labels = paste(breaks[-length(breaks)], c(breaks[-1]-1), sep="-")))

  vars <- c(setdiff(names(counts), c("start", "end", "outcome", "population")))

  return(counts %>% mutate(agegroup = g) %>%
           group_by_at(vars) %>%
           summarize(population = sum(population),
                     outcome = sum(outcome)) %>%
           ungroup() %>%
           separate(agegroup, c("start", "end"), sep = "-", remove = FALSE, convert = TRUE) %>%
           mutate(agegroup = reorder(agegroup, start, first)) %>%
           select(-start, -end) %>%
           drop_na(agegroup))
}
