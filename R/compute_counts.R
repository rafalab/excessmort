#' Compute daily death counts for groups from individual records
#' @export
#' @import dplyr
#' @import rlang
#' @importFrom tidyr complete separate drop_na

compute_counts <- function(dat, by = NULL, demo = NULL,
                           date = "date",
                           age = "age",
                           agegroup = "agegroup",
                           breaks = NULL){

  if(!is.character(by) & !is.null(by)) stop("by needs to be a character verctor or NULL.")


  ## prepare demo if it is provided
  if(!is.null(demo) & (agegroup %in% by)){
    names(demo)[names(demo) == agegroup] <- "agegroup"
    if(is.null(breaks)){
      start <-grep("\\d+", demo$agegroup, value = TRUE) %>% unique() %>% as.numeric() %>% sort()
      breaks <- c(start, Inf)
    } else{
      demo <- collapse_age_dist(demo, breaks)
    }
  }

  if(agegroup %in% by){
    if(is.null(breaks)){
      stop("Need to provide age breaks or a demographics table")
    } else{
      if(age %in% names(dat)){
        dat$agegroup <-group_age(dat[[age]], breaks)
      } else stop(age, "not a column in dat")
    }
  }
  if(!is.null(by))
    if(!all(by %in% names(dat))) stop("by needs to be a subset of", setdiff(names(dat), date))

  ## Needs check that date_colum exists and has dates
  names(dat)[names(dat) == date] <- "date"

  by <- c("date", by)
  dat <- drop_na(dat, by)

  counts <- dat %>% filter(!(lubridate::month(date) == 2 & lubridate::day(date) == 29)) %>%
    group_by_at(by) %>%
    summarise(outcome = n()) %>%
    ungroup() %>%
    complete(!!!syms(by), fill = list(outcome = 0)) %>%
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

