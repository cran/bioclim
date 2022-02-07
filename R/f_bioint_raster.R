#' Computation of Bioclimatic Intensities (raster mode)
#'
#' @description Computes bioclimatic intensities from bioclimatic balance.
#' @param bb Bioclimatic balance in raster format.
#' @param wout name of the output raster. Must include ".tif".
#' @return SpatRaster with 120 layers corresponding to the 12 monthly values of "IBPc","IBCc","IBLc","IBRc","IBSc","IBPf","IBCf","IBLf","IBRf","IBSf".
#' @examples
#' \donttest{
#' bb <- terra::rast(bbRast)
#' bi <- biointRaster(bb, wout=NULL)
#' }
#' @export
#'

biointRaster <- function(bb, wout=NULL){

  B <- bb[[1:12]]
  b <- bb[[13:24]]
  bc <- bb[[25:36]]
  bl <- bb[[37:48]]
  # IBPc
  IBPc <- B * (B > 0)
  names(IBPc) <- paste0('IBPc', formatC(1:12,width = 2, flag = '0'))
  # IBCc
  IBCc <- bc * (bc > 0)
  names(IBCc) <- paste0('IBCc', formatC(1:12,width = 2, flag = '0'))
  # IBLc
  w1 <- (B > 0 & b > 0)
  IBLc <- (bl * w1) * (bl > 0)
  names(IBLc) <- paste0('IBLc', formatC(1:12,width = 2, flag = '0'))

  # IBLc_g
  IBLc_g <- (IBLc > IBCc) * (IBLc - IBCc)
  names(IBLc_g) <- paste0('IBLc_g', formatC(1:12,width = 2, flag = '0'))

  # IBPc_g
  IBPc_g <- IBPc - (IBCc + IBLc_g)
  names(IBPc_g) <- paste0('IBPc_g', formatC(1:12,width = 2, flag = '0'))

  # IBRc
  IBRc <- IBLc - IBCc
  names(IBRc) <- paste0('IBRc', formatC(1:12,width = 2, flag = '0'))

  # IBSc
  IBSc <- b * (B > 0 & b < 0)
  names(IBSc) <- paste0('IBSc', formatC(1:12,width = 2, flag = '0'))

  # IBPf
  IBPf <- B * (B < 0)
  names(IBPf) <- paste0('IBPf', formatC(1:12,width = 2, flag = '0'))

  # IBCf
  IBCf <- bc * (bc[[3]] < 0)
  names(IBCf) <- paste0('IBCf', formatC(1:12,width = 2, flag = '0'))

  # IBLf
  w1 <- (B < 0 & b < 0)
  IBLf <- ((bl * (bl < 0)) * w1) + ((b * (bl > 0)) * w1)
  names(IBLf) <- paste0('IBLf', formatC(1:12,width = 2, flag = '0'))

  # IBLf_g
  IBLf_g <- (IBCf > IBLf) * (IBLf - IBCf)
  names(IBLf_g) <- paste0('IBLf_g', formatC(1:12,width = 2, flag = '0'))

  # IBPf_g
  IBPf_g <- IBPf + (IBCf - IBLf_g)
  names(IBPf_g) <- paste0('IBPf_g', formatC(1:12,width = 2, flag = '0'))

  # IBRf
  IBRf <- IBLf + IBCf
  names(IBRf) <- paste0('IBRf', formatC(1:12,width = 2, flag = '0'))

  # IBSf
  w1 <- (B < 0)
  IBSf <- (b * (b > 0)) * w1
  names(IBSf) <- paste0('IBSf', formatC(1:12,width = 2, flag = '0'))

  # intens <- c(IBPc, IBPc_g, IBCc, IBLc, IBLc_g, IBRc, IBSc,
  #             IBPf, IBPf_g, IBCf, IBLf, IBLf_g, IBRf, IBSf)
  intens <- c(IBPc, IBCc, IBLc, IBRc, IBSc,
              IBPf, IBCf, IBLf, IBRf, IBSf)

  if(!is.null(wout)){
    terra::writeRaster(intens, wout, overwrite=TRUE)
  }
  return(intens)
}
