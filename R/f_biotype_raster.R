#' Bioclimatic classification (raster mode)
#'
#' @description Calculates bioclimatic classification based on bioclimatic balance.
#' @param temp SpatRaster object with 12 layers representing temperature from January to December.
#' @param prec SpatRaster object with 12 layers representing precipitation from January to December.
#' @param CC Field capacity. It can be numeric (1 value) or a SpatRaster object.
#' @param PET Potential evapotranspiration. Optional. It must be a SpatRaster object.
#' @param bh Water balance. Optional. It must be a SpatRaster object.
#' @param wout Optional. Path and name of the output raster file.
#' @param ncpu number of cores to use in calculation. If not provided, sequential mode is used (1 core).
#' @return SpatRaster with 3 variables ("TBR": Types of Bioclimatic Regime; "zonal": zonal units; "sub": bioclimatic regime subtypes).
#' @examples
#' \donttest{
#' wb <- terra::rast(wbRast)
#' btr <- biotypeRaster(bh = wb)
#'}
#' @export
#'


biotypeRaster <- function(temp=NULL, prec=NULL, CC=NULL, wout=NULL, ncpu = 1, PET = NULL, bh = NULL){

  # If water balance is not provided, CC is computed
  if(is.null(bh)){
    # CC
    if(is.numeric(CC)){
      if(length(CC)>1) stop("CC must be a single number") else{
        # mask
        msk <- prec[[1]]/prec[[1]]
        CC <- msk*CC
      }
    } else if(class(CC) !='SpatRaster') stop("CC must be numeric or a SpatRaster")
  } else{
    if(class(bh) !='SpatRaster') stop("bh must be a SpatRaster")
    prec <- bh[[13:24]]
    temp <- bh[[1:12]]
    PET <- bh[[25:36]]
  }

  # check if PET is provided, then water balance, then computes it
  if(!is.null(PET)){
    if(class(PET) !='SpatRaster') stop("PET must be a SpatRaster")
  } else{
    if(is.null(bh)){
      message(paste0("Computing water balance ",'[',Sys.time(),']'))
      bh <- watbalRaster(temp, prec, CC, ncpu)
      PET <- bh[[25:36]]
    }
  }

    # message("Computing Thornthwaite’s index")
    subtype <- ithRaster(bh)

  # other required variables
  PVT <- sum(temp < 7.5)
  PVH <- sum((prec-0.2*PET) < 0)

  Tf <- min(temp)
  P <- sum(prec)


  # classification ----
  message(paste0("Computing Classification ",'[',Sys.time(),']'))
  zonal <- type <- P * 0
  # non tropical
  wnt <- (Tf < 18)
    wst <- wnt * (Tf > 7.5)
    zonal <- zonal + (wst * 1) # subtropical
      wPVH <- (PVH == 0) * wst
      type <- type + (((P > 2000) * wPVH) * 1) # Euritermo-Ombrophyllo
      type <- type + (((P <= 2000) * wPVH) * 2) # Euritermo-Mesophyllo
      wPVH <- (PVH != 0) * wst
        type <- type + (((PVH > 0 & PVH <= 4) * wPVH) * 3) # Euritermo-Tropophyllo
        type <- type + (((PVH > 4 & PVH < 8) * wPVH) * 4) # Euritermo-Xerophyllo
        type <- type + (((PVH >= 8) * wPVH) * 5) # Euritermo-Hyperxerophyllo
    wnt <- wnt * (Tf <= 7.5)
      wtm <- (PVT >= 1 & PVT < 6) * wnt
        zonal <- zonal + (wtm * 2) # Mid-latitudes
        w2 <- wtm * (PVH == 0)
          type <- type + (((P > 2000) * w2) * 6) # Cryo-Ombrophyllo
          type <- type + (((P <= 2000) * w2) * 7) # Cryo-Mesophyllo
        w2 <- wtm * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 8) # Cryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 9) # Cryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 10) # Cryo-Hyperxerophyllo
      wsp <- (PVT >= 6 & PVT < 11) * wnt
        zonal <- zonal + (wsp * 3) # Subpolar
        w2 <- wsp * (PVH == 0)
          type <- type + (((P > 1100) * w2) * 11) # Mesocryo-Ombrophyllo
          type <- type + (((P <= 1100) * w2) * 12) # Mesocryo-mesophyllo
        w2 <- wsp * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 13) # Mesocryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 14) # Mesocryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 15) # Mesocryo-Hyperxerophyllo
      wp <- (PVT >= 11) * wnt
        zonal <- zonal + (wp * 4) # Polar
        w2 <- wp * (PVH == 0)
          type <- type + (((P > 1100) * w2) * 16) # Hypercryo-Ombrophyllo
          type <- type + (((P <= 1100) * w2) * 17) # Hypercryo-mesophyllo
        w2 <- wp * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 18) # Hypercryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 19) # Hypercryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 20) # Hypercryo-Hyperxerophyllo
  wt <- (Tf >= 18)
    zonal <- zonal + (wt * 5) # Tropical
    w2 <- wt * (PVH == 0)
      type <- type + (((P > 2000) * w2) * 21) # Ombrophyllo
      type <- type + (((P <= 2000) * w2) * 22) # Mesophyllo
    w2 <- wt * (PVH > 0 & PVH <= 4)
      type <- type + (((P > 2000) * w2) * 23) # Ombro-Tropophyllo !
      type <- type + (((P <= 2000) * w2) * 24) # Tropophyllo
    w2 <- wt * (PVH > 4 & PVH < 8)
      type <- type + (((P > 2000) * w2) * 25) # Ombro-Xerophyllo !
      type <- type + (((P <= 2000) * w2) * 26) # Xerophyllo
    w2 <- wt * (PVH >= 8)
      type <- type + (((PVH >= 8) * w2) * 27) # Hyperxerophyllo

  # remove zeros
  zonal[zonal == 0] <- NA
  type[type == 0] <- NA


  zonalnames <- c('Subtropical', 'Mid latitudes', 'Subpolar', 'Polar', 'Tropical')
  subtypenames <- c('HyperArid', 'Arid', 'Semiarid', 'Dry humid', 'Moist humid',
                    'Low humid', 'Moderate humid', 'Highly humid', 'Very humid','Perhumid')
  typenames <- c('Euritermo-Ombrophyllo', 'Euritermo-Mesophyllo', 'Euritermo-Tropophyllo',
                 'Euritermo-Xerophyllo', 'Euritermo-Hyperxerophyllo', 'Cryo-Ombrophyllo',
                 'Cryo-Mesophyllo', 'Cryo-Tropophyllo', 'Cryo-Xerophyllo', 'Cryo-Hyperxerophyllo',
                 'Mesocryo-Ombrophyllo', 'Mesocryo-Mesophyllo', 'Mesocryo-Tropophyllo', 'Mesocryo-Xerophyllo',
                 'Mesocryo-Hyperxerophyllo', 'Hypercryo-Ombrophyllo', 'Hypercryo-Mesophyllo',
                 'Hypercryo-Tropophyllo', 'Hypercryo-Xerophyllo', 'Hypercryo-Hyperxerophyllo',
                 'Ombrophyllo', 'Mesophyllo', 'Ombro-Tropophyllo', 'Tropophyllo', 'Ombro-Xerophyllo',
                 'Xerophyllo', 'Hyperxerophyllo')
  df <- data.frame(ID = as.numeric(sapply(1:10, paste0, 1:27)),
                   category = as.character(sapply(subtypenames, paste, typenames)))

  levels(zonal) <- data.frame(ID = 1:5, category = zonalnames)
  # plot(zonal)

  levels(type) <- data.frame(ID = 1:27, category = typenames)
  # plot(type)

    levels(subtype) <- data.frame(ID = 1:10, category = subtypenames)

    funion <- function(uni){
      p <- as.numeric(paste0(uni[1], uni[2]))
      return(p)
    }
    ct <- suppressWarnings(app(c(subtype, type), funion))
    d <- df[match(unique(ct)[[1]], df$ID), 2]
    m <- cbind(unique(ct)[[1]], 1:length(unique(ct)[[1]]))
    ct <- terra::classify(ct, m)
    levels(ct) <- data.frame(ID = 1:length(unique(ct)[[1]]), category = d)
    # plot(ct, plg=list(x="bottom"))
    res <- c(zonal, type, ct)
    names(res) <- c('Zonal', 'TBR', 'Subtype')

  # zonal:
  # 1: Subtropical
  # 2: Mid latitudes
  # 3: Subpolar
  # 4: Polar
  # 5: Tropical

  # type:
  # (1: Subtropical)
    # 1: Euritermo-Ombrophyllo
    # 2: Euritermo-Mesophyllo
    # 3: Euritermo-Tropophyllo
    # 4: Euritermo-Xerophyllo
    # 5: Euritermo-Hyperxerophyllo
  # (2: Mid latitudes)
    # 6: Cryo-Ombrophyllo
    # 7: Cryo-Mesophyllo
    # 8: Cryo-Tropophyllo
    # 9: Cryo-Xerophyllo
    # 10: Cryo-Hyperxerophyllo
  # (3: Subpolar)
    # 11: Mesocryo-Ombrophyllo
    # 12: Mesocryo-Mesophyllo
    # 13: Mesocryo-Tropophyllo
    # 14: Mesocryo-Xerophyllo
    # 15: Mesocryo-Hyperxerophyllo
  # (4: Polar)
    # 16: Hypercryo-Ombrophyllo
    # 17: Hypercryo-Mesophyllo
    # 18: Hypercryo-Tropophyllo
    # 19: Hypercryo-Xerophyllo
    # 20: Hypercryo-Hyperxerophyllo
  # (5: Tropical)
    # 21: Ombrophyllo
    # 22: Mesophyllo
    # 23: Tropophyllo
    # 24: Xerophyllo
    # 25: Hyperxerophyllo



  if(!is.null(wout)){
    terra::writeRaster(res, wout, overwrite=TRUE)
  }

  message(paste0("End ",'[',Sys.time(),']'))
  return(res)
}
