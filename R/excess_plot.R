#' Plot results from fitted excess count model
#' @export
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line
#'

excess_plot <- function(fit, title = "", ylim = NULL,
                        center = c("zero", "fit"),
                        show.data = TRUE,
                        alpha = 0.01){

  requireNamespace("ggplot2")

  z <- qnorm(1 - alpha/2)

  if(match.arg(center) == "zero") center <- 0 else center <- fit$fitted

  p <- with(fit, data.frame(date = date,
                        y = 100 * (observed - expected)/expected,
                        increase = 100 * fitted,
                        lower = 100 * (center - z*se),
                        upper = 100 * (center + z*se))) %>%
    ggplot(aes(date, y, ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 0) +
    ylab("% increase from expected death rate") +
    ggtitle(title) +
    xlab("Date")

  if(show.data) p <- p + geom_point(alpha = 0.5)

  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  return(p + geom_line(aes(y = increase), col="#3366FF"))
}
