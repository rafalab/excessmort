#' Cook County Medical Examiner Records
#' 
#' A data frame with each row represening a death. The data includes 
#' death date and demographic information. When you load this dataset you also 
#' load `cook_demographics` which includes the population size for each demographice
#' group by date.
#' 
#' @usage data("cook_records")
#' 
#' @format A data frame with these columns
#' \describe{
#' \item{sex}{Sex of the deceased}
#' \item{age}{Age of the deceased}
#' \item{race}{Race of the deceased}
#' \item{residenceplace}{Residence placed of the deceased}
#' \item{date}{Date of the death}
#' \item{cause_1}{Reported cause of death}
#' \item{type_of_death}{Type of death}
#' }
#' 
#' @docType data
#'
#'@aliases cook_demographics
"cook_records"
