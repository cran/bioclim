#' Function to plot bioclimatic balance
#'
#' @description Function to plot bioclimatic balance.
#' @param intens bioclimatic intensities in data.frame format from bioint() function.
#' @return Plot of bioclimatic balance
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_minimal theme guides xlab ylab scale_x_discrete ggtitle element_text element_blank guide_legend rel unit
#' @import reshape2
#' @examples
#' wb <- watbal(t = c(10, 11.5, 14, 16.5, 20, 24.5, 27.5, 28, 24.5, 19.5, 14.5, 11),
#'              p = c(55, 73, 84, 58, 33, 23, 2, 2, 28, 66, 94, 71), lat = 35, CC = 400)
#'bb <- biobal(wb, 400)
#'bi <- bioint(bb)
#'plotBiobal(bi)
#'
#' @export
#'

plotBiobal <- function(intens){
  #
  dbb <- intens[, rev(c(3, 5, 2, 7, 10, 12, 14, 9))]

  dbb$mon <- as.character(formatC(1:12, width = 2, flag = '0'))
  dbb <- melt(dbb, id.var="mon")

  ggplot(dbb[order(dbb$variable),], aes(x = mon, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("legend", values = c("IBPf_g" = 'azure',
                                           "IBSf" = "yellow",
                                           "IBLf_g" = "blue1",
                                           'IBCf' = 'cyan',
                                           'IBSc' = 'red',
                                           'IBPc_g' = 'cornsilk',
                                           'IBLc_g' = 'darkgreen',
                                           'IBCc' = 'chartreuse'),
                      labels = c('IBPf', 'IBSf', 'IBLf', 'IBCf',
                                 'IBSc', 'IBPc', 'IBLc', 'IBCc')) +
    theme_minimal() +
    theme(legend.position="bottom",
          legend.title = element_blank(),
          axis.title.y = element_text(size = rel(0.8)),
          plot.title = element_text(size = rel(0.9), face = 'bold'),
          legend.spacing.x = unit(0.15, 'cm')) +
    guides(fill = guide_legend(nrow = 1)) +
    xlab('') + ylab('') +
    scale_x_discrete(labels = month.abb) +
    ggtitle("Bioclimatic balance")

}
