#' Compute daily death counts for groups from individual records
#' @export
#' @import dplyr

compute_counts <- function(dat, group_by = NULL, demo = NULL,
                           date_column = "date",
                           age_column = "age",
                           agegroup_column = "agegroup",
                           age_breaks = NULL){

  if(!is.character(group_by) & !is.null(group_by)) stop("group_by needs to be a character verctor or NULL.")


  ## prepare demo if it is provided
  if(!is.null(demo) & (agegroup_column %in% group_by)){
    names(demo)[names(demo) == agegroup_column] <- "agegroup"
    if(is.null(age_breaks)){
      start <-grep("\\d+", demo$agegroup, value = TRUE) %>% unique() %>% as.numeric() %>% sort()
      age_breaks <- c(start, Inf)
    } else{
      demo <- collapse_age_dist(demo, age_breaks)
    }
  }

  if(agegroup_column %in% group_by){
    if(is.null(age_breaks)){
      stop("Need to provide age breaks or a demographics table")
    } else{
      if(age_column %in% names(dat)){
        dat$agegroup <-group_age(dat[[age_column]], age_breaks)
      } else stop(age_column, "not a column in dat")
    }
  }
  if(!is.null(group_by))
    if(!all(group_by %in% names(dat))) stop("group_by needs to be a subset of", setdiff(names(dat), date_column))

  ## Needs check that date_colum exists and has dates
  names(dat)[names(dat) == date_column] <- "date"

  group_by <- c("date", group_by)
  dat <- drop_na(dat, group_by)

  counts <- dat %>% filter(!(month(date) == 2 & day(date) == 29)) %>%
    group_by_at(group_by) %>%
    summarise(outcome = n()) %>%
    ungroup() %>%
    complete_(group_by, fill = list(outcome = 0)) %>%
    arrange(date)

  if(!is.null(demo)){
    by <- intersect(names(demo), names(counts))
    demo <- demo %>%
      group_by_at(by) %>%
      summarize(population = sum(population)) %>%
      ungroup()
    counts <- left_join(counts, demo, by = by)
  }
  return(counts)
}

