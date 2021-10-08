#' Plot results from fitted excess count model
#' 
#' @param fit The output from `excess_model`
#' @param title A title to add to plot
#' @param ylim A vector with two numbers that determines the kimits for the y-axis
#' @param show.data A logical that determines if the observed percent changes are shown
#' @param alpha  1 - `alpha` confidence intervals are shown
#' 
#' @return An ggplot object containing the plot.
#' 
#' @examples
#' data(new_jersey_counts)
#' exclude_dates <- as.Date("2012-10-29") + 0:180
#' control_dates <- seq(min(new_jersey_counts$date), min(exclude_dates) - 1, by = "day")
#' f <- excess_model(new_jersey_counts,
#' start = as.Date("2012-09-01"), 
#' end = as.Date("2013-09-01"), 
#' exclude = exclude_dates,
#' model = "correlated",
#' weekday.effect = TRUE,
#' control.dates = control_dates)
#' 
#' library(ggplot2)
#' excess_plot(f)
#' 
#' @export
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_hline ggtitle aes geom_ribbon ylab xlab geom_point coord_cartesian geom_line scale_y_continuous scale_y_continuous scale_x_date theme_bw 
#'

excess_plot <- function(fit, 
                        title     = "", 
                        ylim      = NULL,
                        show.data = TRUE,
                        alpha     = 0.05){

  requireNamespace("ggplot2")

  z <- qnorm(1 - alpha/2)


  p <- with(fit, data.frame(date = date,
                        y        = (observed - expected)/expected,
                        increase = fitted,
                        sd       = sd,
                        se       = se)) %>%
    ggplot(aes(date, y)) +
    geom_ribbon(aes(ymin = increase - z * se, ymax = increase + z * se), alpha = 0.5, fill = "#3366FF") +
    geom_hline(yintercept = 0) +
    ylab("% increase from expected death rate") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_date(date_labels = "%b %Y") +
    ggtitle(title) +
    xlab("Date") 

  if(show.data) p <- p + geom_point(alpha = 0.3)

  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  return(p + geom_line(aes(y = increase), size=0.70, col="#3366FF"))
}
