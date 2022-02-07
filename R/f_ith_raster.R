#' Function to calculate Thornthwaite’s index (raster format)
#'
#' @description This function calculates Thornthwaite’s index to refine the bioclimatic classification.
#' @param bh Water balance in SpatRaster format from watbalRaster() function.
#' @return Numeric, describing the humid characteristics of the climate. 1: 'HyperArid', 2: 'Arid', 3: 'Semiarid', 4: 'Dry humid', 5: 'Moist humid', 6 'Low humid', 7: 'Moderate humid', 8: 'Highly humid', 9: 'Very humid', 10: 'Perhumid'.
#' @examples
#' \donttest{
#' wb <- terra::rast(wbRast)
#' itr <- ithRaster(wb)
#' }
#' @export

ithRaster <- function(bh){
  num <- (sum(bh[[37:48]])/sum(bh[[25:36]]))*100

  m <- c(-99999999, -60, 1,
         -60, -40, 2,
         -40, -20, 3,
         -20, 0, 4,
         0, 20, 5,
         20, 40, 6,
         40, 60, 7,
         60, 80, 8,
         80, 100, 9,
         100, 99999999, 10)
  m <- matrix(m, ncol=3, byrow=TRUE)
  clf <- terra::classify(num, m)

  # 1: 'HyperArid'
  # 2: 'Arid',
  # 3: 'Semiarid',
  # 4: 'Dry humid',
  # 5: 'Moist humid',
  # 6 'Low humid',
  # 7: 'Moderate humid',
  # 8: 'Highly humid',
  # 9: 'Very humid',
  # 10: 'Perhumid'

  return(clf)
}
