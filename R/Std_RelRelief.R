#' Standard Relative Relief
#'
#' @param inRas Input DEM. (SpatRas)
#' @param sc Scale of window.
#'
#' @return Relative relief raster (SpatRast)
#' @export
#'
relRelief <- function(inRas, sc){

  # Focal min and max with moving window of size sc x sc
  mn <- terra::focal(inRas, w = sc, fun = "min", na.rm = TRUE,
                     na.policy = "All", fillvalue = NA)

  mx <- terra::focal(inRas, w = sc, fun = "max", na.rm = TRUE,
                     na.policy="All", fillvalue=NA)

  # Calc. RR
  relRelief <- (inRas - mn) / (mx - mn)

  # Return RR SpatRast
  return(relRelief)
}

#' Standard Topographic Position Index
#'
#' @param inRas Input DEM. (SpatRaster)
#' @param sc Scale of window.
#'
#' @return Relative relief raster (SpatRaster)
#' @export
#'
TPI <- function(inRas, sc){

  # Focal mean with moving window of size sc x sc
  mn <- terra::focal(inRas, w = sc, fun = "mean", na.rm = TRUE,
                     na.policy = "All", fillvalue = NA)

  # Calc. TPI
  TPI <- (inRas - mn)

  # Return TPI SpatRast
  return(TPI)
}
