#' Compute counts
#' 
#' Compute counts for groups from individual records
#' 
#' This is helper function that helps convert individual records data, in which each death is a row, 
#' to a count data frame where each row is a date. It is particulalry helpful for defining agegroups. If you
#' provide the `breaks` argument it will automaticall divided the data and provide the counts for each age strata. 
#' You can also select variables to group by using the `by` argument. 
#' One can provide a data frame with demogrpahic inform through the `demo` argument. This tabe must have the 
#' population size for each group for each data.
#'
#' 
#' @param dat The data frame with the individual records
#' @param by A character vector with the column names the define the groups for which we will compute counts 
#' @param demo A data frame with population sizes for each time point for each of the groups for which we will compute counts
#' @param date The column name of the column that contains dates
#' @param age The column name of the column that contains age
#' @param agegroup The column name of the column that contains agegroup
#' @param breaks The ages that define the agegroups
#' 
#' @return A data frame with counts for each group for each date with population sizes, if demo was provided.
#' 
#' @examples
#' library(lubridate)
#' data("cook_records")
#' the_breaks <- c(0, 20, 40, 60, 75, Inf)
#' 
#' ## take subset for example
#' cook_records_subset <- cook_records[year(cook_records$date)==2021, ]
#' cook_demographics_subset <- cook_demographics[year(cook_demographics$date)==2021, ]
#' 
#' cook_counts <- compute_counts(cook_records_subset, 
#'                        demo = cook_demographics_subset, 
#'                        by = c("agegroup", "race", "sex"), 
#'                        breaks = the_breaks)
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

