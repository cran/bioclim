#' Computation of Bioclimatic Intensities
#'
#' @description Computes bioclimatic intensities from bioclimatic balance.
#' @param bb Bioclimatic balance.
#' @return data frame with 14 variables: "IBPc","IBPc_g","IBCc","IBLc","IBLc_g","IBRc","IBSc","IBPf","IBPf_g","IBCf","IBLf","IBLf_g","IBRf","IBSf".
#' @examples
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' bb <- biobal(wb, 400)
#' bi <- bioint(bb)
#' @export
#'

bioint <- function(bb){
  intens <- matrix(NA, ncol = 14, nrow = 12)
  colnames(intens) <-  c('IBPc', 'IBPc_g', 'IBCc', 'IBLc', 'IBLc_g', 'IBRc', 'IBSc',
                         'IBPf', 'IBPf_g', 'IBCf', 'IBLf', 'IBLf_g', 'IBRf', 'IBSf')
  rownames(intens) <- month.name
  intens <- as.data.frame(intens)

  for(i in 1:12){
    if(bb$B[i] > 0) intens$IBPc[i] <- bb$B[i] else intens$IBPc[i] <- 0
    if(bb$bc[i] > 0) intens$IBCc[i] <- bb$bc[i] else intens$IBCc[i] <-0
    if(bb$B[i] > 0 & bb$b[i] > 0){if(bb$bl[i] > 0) intens$IBLc[i] <- bb$bl[i] else intens$IBLc[i] <- bb$b[i]} else intens$IBLc[i] <- 0
    if(intens$IBLc[i] < intens$IBCc[i]) intens$IBLc_g[i] <- 0 else intens$IBLc_g[i] <- intens$IBLc[i] - intens$IBCc[i]
    intens$IBPc_g[i] <- intens$IBPc[i] - (intens$IBCc[i] + intens$IBLc_g[i])
    intens$IBRc[i] <- intens$IBLc[i] + intens$IBCc[i]
    if(bb$B[i] > 0 & bb$b[i] < 0) intens$IBSc[i] <- bb$b[i] else intens$IBSc[i] <- 0
    if(bb$B[i] < 0) intens$IBPf[i] <- bb$B[i] else intens$IBPf[i] <- 0
    if(bb$bc[3] < 0) intens$IBCf[i] <- bb$bc[3] else intens$IBCf[i] <- 0
    if(bb$B[i] < 0 & bb$b[i] < 0) {if(bb$bl[i] < 0) intens$IBLf[i] <- bb$bl[i] else intens$IBLf[i] <- bb$b[i]} else intens$IBLf[i] <- 0
    if(intens$IBCf[i] < intens$IBLf[i]) intens$IBLf_g[i] <- 0 else intens$IBLf_g[i] <- intens$IBLf[i] - intens$IBCf[i]
    intens$IBPf_g[i] <- intens$IBPf[i] - (intens$IBCf[i] + intens$IBLf_g[i])
    intens$IBRf[i] <- intens$IBLf[i] + intens$IBCf[i]
    if(bb$B[i] < 0) {if(bb$b[i] > 0) intens$IBSf[i] <- bb$b[i] else intens$IBSf[i] <- 0} else intens$IBSf[i] <- 0
  }
  intens
}
