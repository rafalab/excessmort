#' Callapse age groups into broader ones
#'
#' @export
#' @importFrom stats reorder
#' @import dplyr
collapse_age_dist <- function(demo, breaks){
  ## assumes groups are start-end
  ## assumes agegroup contains group
  ## assumes column named population has population values

  demo <- separate(demo, agegroup, c("start", "end"), sep = "-", remove = FALSE, convert = TRUE)

  edges <- c(unique(demo$start), Inf)

  if(!all(breaks %in% edges)) stop("Age breaks must be a subset", edges)

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
