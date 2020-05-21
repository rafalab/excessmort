#' Check if expected counts fit data
#' @export
#' @import dplyr
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line


expected_plot <- function(expected, title = "",
                          start = NULL,
                          end = NULL,
                          ylim = NULL,
                          weekly = TRUE,
                          color = "#3366FF",
                          alpha = 1){

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
      summarize(expected = mean(expected)*7,
                observed = mean(observed)*7) %>%
      ungroup()
  }


  yl <- "Counts"
  if(weekly) yl <- paste("Weekly", yl)
  p <- dat %>%
    ggplot(aes(x = date)) +
    geom_point(aes(y = observed), alpha = alpha) +
    geom_line(aes(y = expected), color = color) +
    ylab(yl) +
    ggtitle(title) +
    xlab("Date")

  if(!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  return(p)
}
