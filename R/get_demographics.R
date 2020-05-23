#' Get demographic data from Census
#' @param geography state or county
#' @param state name of the state
#' @param county name of the county
#' @param years vector of years for which we obtain data
#' @param vars names of variables that define strat of which we want population estimates
#' @export
#' @importFrom stats setNames
#' @import dplyr

get_demographics <- function(geography="state",
                             state,
                             county=NULL,
                             years=2018,
                             vars = c("SEX", "AGEGROUP", "RACE", "HISP"))
{
  ################## ----- PARAMETERS ----- ##################
  # 1. geography: Geographical granularity of interest (state, county)
  # 2. state    : State of interest
  # 3. county   : County of interest within the chosen state
  # 4. years    : A vector of year values of interest
  # 5. vars     : Vector of variable of interest

  requireNamespace("tidycensus")

  if(! "package:tidycensus" %in% search() ) stop("tidycensus must be loaded.\nLoad with library(tidycensus). Note that you will need to get a key.")
  if (Sys.getenv("CENSUS_API_KEY") != "") {
    key <- Sys.getenv("CENSUS_API_KEY")
  }
  else if (is.null(key)) {
    stop("A Census API key is required.  Obtain one at http://api.census.gov/data/key_signup.html, and then supply the key to the `census_api_key` function to use it throughout your tidycensus session.")
  }

  # -- Check state name
  if(nchar(state) > 2){

    # -- Correct spell for state name
    if(!state %in% datasets::state.name){ stop(paste0(state," is not a state"))}

    # -- Getting state abbreviation
    state <- datasets::state.abb[grep(state, datasets::state.name)]
  } else {

    # -- Correct spell for state abbreviation
    if(!state %in% datasets::state.name){ stop(paste0(state," is not a state"))}


  }

  # -- Getting data
  demographics <- lapply(years, function(x){

    tmp <- tidycensus::get_estimates(geography        = geography,
                         product          = "characteristics",
                         breakdown        = vars,
                         breakdown_labels = TRUE,
                         state            = state,
                         county           = county,
                         year             = x)

    mutate(tmp, year = x)
  })

  # -- Row binding and wrangling
  demographics <- do.call(rbind, demographics) %>%
    filter((grepl("to", AGEGROUP) | AGEGROUP == "Age 85 years and older") &
             SEX != "Both sexes" & HISP != "Both Hispanic Origins" &
             (grepl("alone$", RACE) | RACE == "Two or more races")) %>%
    mutate(AGEGROUP = droplevels(AGEGROUP)) %>%
    mutate(RACE = case_when(
      HISP == "Hispanic" ~ "Hispanic",
      grepl("White", RACE) ~ "White",
      grepl("Black", RACE) ~ "Black",
      grepl("Asian", RACE) ~ "Asian",
      TRUE ~ "Other")) %>%
    group_by(SEX, AGEGROUP, RACE, year) %>%
    summarize(value = sum(value)) %>%
    ungroup() %>%
    mutate(AGEGROUP = gsub(" to ", "-", AGEGROUP),
           AGEGROUP = gsub("[[:alpha:]]| ", "", AGEGROUP),
           AGEGROUP = ifelse(AGEGROUP=="85","85-Inf", AGEGROUP),
           RACE     = tolower(RACE),
           SEX      = tolower(SEX)) %>%
    setNames(c("sex", "agegroup", "race", "year", "population")) %>%
    select(year, sex, race, agegroup, population)

  return(demographics)
}
