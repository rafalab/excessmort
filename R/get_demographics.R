# -- Function to get demographic data from Census
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

  # -- Check state name
  if(nchar(state) > 2){

    # -- Correct spell for state name
    if(!state %in% state.name){ stop(paste0(state," is not a state"))}

    # -- Getting state abbreviation
    state <- state.abb[grep(state, state.name)]
  } else {

    # -- Correct spell for state abbreviation
    if(!state %in% state.abb){ stop(paste0(state," is not a state"))}


  }

  # -- Getting data
  demographics <- lapply(years, function(x){

    tmp <- get_estimates(geography        = geography,
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
