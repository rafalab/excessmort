#' Plot Expected Counts
#' 
#' Check if expected counts fit data
#' 
#' @param expected The output from `excess_counts`
#' @param title A title to add to plot
#' @param start First day to show
#' @param end Last day to show
#' @param ylim A vector with two numbers that determines the kimits for the y-axis
#' @param weekly Logical that determines if data should be summarized into weekly counts
#' @param color Color for the expected curve
#' @param alpha alpha blending for points
#' @examples
#' data(florida_counts)
#' exclude_dates <- as.Date("2017-09-10") + 0:180
#' e <- compute_expected(florida_counts, exclude = exclude_dates, weekday.effect = TRUE)
#' 
#' library(ggplot2)
#' expected_plot(e,
#' start = as.Date("2017-9-1"), 
#' end = as.Date("2018-9-1"), )
#' 

#' @export
#' @import dplyr
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line


expected_plot <- function(expected, 
                          title  = "",
                          start  = NULL,
                          end    = NULL,
                          ylim   = NULL,
                          weekly = TRUE,
                          color  = "#3366FF",
                          alpha=0.50){

  if(!"compute_expected" %in% class(expected)) stop("The expected argument needs to be the output of the computed_expected function.")
  if(!attr(expected, "frequency") %in% c(365, 52)) warning("This function assumes weekly or daily data. This dataset has ", attr(expected, "frequency"), "counts per year.")
  requireNamespace("ggplot2")

  if(is.null(start)) start <- min(expected$date)
  if(is.null(end)) end <- max(expected$date)

  dat <- with(expected,
              data.frame(date = date,
                         observed = outcome,
                         expected = expected)) %>%
    filter(date >= start & date <= end)
  if(weekly){
    dat <- dat %>%
      mutate(date = lubridate::floor_date(date, unit = "week")) %>%
      group_by(date) %>%
      summarize(expected = sum(expected),
                observed = sum(observed),
                n = n()) %>%
      ungroup()
    if(attr(expected, "frequency") == 365) dat <- filter(dat, n == 7) 
  }
  
  yl <- "Counts"
  if(weekly) yl <- paste("Weekly", yl)
  p <- dat %>%
    ggplot(aes(x = date)) +
    geom_point(aes(y = observed), alpha = alpha) +
    geom_line(aes(y = expected), size=0.70, color = color) +
    ylab(yl) +
    scale_y_continuous(labels = scales::comma) +
    scale_x_date(date_labels = "%b %Y") +
    ggtitle(title) +
    xlab("Date") +
    theme_bw()

  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  return(p)
}
