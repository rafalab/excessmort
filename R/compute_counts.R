#' Compute daily death counts for groups from individual records
#' @export
#' @import dplyr

compute_counts <- function(dat, group.by = NULL, demo = NULL,
                           date = "date",
                           age = "age",
                           agegroup = "agegroup",
                           breaks = NULL){

  if(!is.character(group.by) & !is.null(group.by)) stop("group.by needs to be a character verctor or NULL.")


  ## prepare demo if it is provided
  if(!is.null(demo) & (agegroup %in% group.by)){
    names(demo)[names(demo) == agegroup] <- "agegroup"
    if(is.null(breaks)){
      start <-grep("\\d+", demo$agegroup, value = TRUE) %>% unique() %>% as.numeric() %>% sort()
      breaks <- c(start, Inf)
    } else{
      demo <- collapse_age_dist(demo, breaks)
    }
  }

  if(agegroup %in% group.by){
    if(is.null(breaks)){
      stop("Need to provide age breaks or a demographics table")
    } else{
      if(age %in% names(dat)){
        dat$agegroup <-group_age(dat[[age]], breaks)
      } else stop(age, "not a column in dat")
    }
  }
  if(!is.null(group.by))
    if(!all(group.by %in% names(dat))) stop("group_by needs to be a subset of", setdiff(names(dat), date))

  ## Needs check that date_colum exists and has dates
  names(dat)[names(dat) == date] <- "date"

  group.by <- c("date", group.by)
  dat <- drop_na(dat, group.by)

  counts <- dat %>% filter(!(month(date) == 2 & day(date) == 29)) %>%
    group_by_at(group.by) %>%
    summarise(outcome = n()) %>%
    ungroup() %>%
    complete_(group.by, fill = list(outcome = 0)) %>%
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

