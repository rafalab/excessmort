#' Plot results from fitted excess count model
#' @export
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line
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
