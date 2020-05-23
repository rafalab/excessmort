#' Plot results from fitted excess count model
#' 
#' @param fit The output from `excess_model`
#' @param title A title to add to plot
#' @param ylim A vector with two numbers that determines the kimits for the y-axis
#' @param show.data A logical that determines if the observed percent changes are shown
#' @param alpha  1 - `alpha` confidence intervals are shown
#' @examples
#' data(florida_counts)
#' exclude_dates <- as.Date("2017-09-10") + 0:180
#' f <- excess_model(florida_counts, 
#' start = as.Date("2017-9-1"), 
#' end = as.Date("2018-9-1"), 
#' exclude = exclude_dates)
#' 
#' library(ggplot2)
#' excess_plot(f)
#' 
#' @export
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_hline ggtitle aes geom_ribbon ylab xlab geom_point coord_cartesian geom_line
#'

excess_plot <- function(fit, title = "", ylim = NULL,
                        show.data = TRUE,
                        alpha = 0.05){

  requireNamespace("ggplot2")

  z <- qnorm(1 - alpha/2)


  p <- with(fit, data.frame(date = date,
                        y = 100 * (observed - expected)/expected,
                        increase = 100 * fitted,
                        sd = 100 * sd,
                        se = 100 * se)) %>%
    ggplot(aes(date, y)) +
    geom_ribbon(aes(ymin = increase - z * se, ymax = increase + z * se), alpha = 0.5) +
    geom_hline(yintercept = 0) +
    ylab("% increase from expected death rate") +
    ggtitle(title) +
    xlab("Date")

  if(show.data) p <- p + geom_point(alpha = 0.5)

  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  return(p + geom_line(aes(y = increase), col="#3366FF"))
}
