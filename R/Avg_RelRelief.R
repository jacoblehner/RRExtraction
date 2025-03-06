#' Multi-scale average RR at a given cell
#'
#' @param inRas Input DEM. (SpatRast)
#' @param inMn Minimum scale. (default = 7)
#' @param inMx Maximum scale. (default = 21)
#' @param stp Scale step. (default = 2)
#' @param ii Row position of cell.
#' @param jj Column position of cell.
#'
#' @return Average relative relief at chosen cell. (numeric)
#' @export
#'
rr.AvgPoint <- function(inRas, inMn = 7, inMx = 21, stp = 2, ii, jj){
  x <- base::as.matrix(inRas, wide = TRUE)
  # Sequence of scales between min and max bounds
  seq.RR <- base::seq(inMn, inMx, stp)

  # offset from center based on scale
  m = (seq.RR - 1) / 2

  rr.here <- NA
  for(k in 1:base::length(m)){

    tmp <- c(x[(ii - m[k]):(ii + m[k]),(jj - m[k]):(jj + m[k])])
    rr.here[k] <- (x[ii, jj] - base::min(tmp, na.rm=TRUE)) /
      (base::max(tmp, na.rm=TRUE) - base::min(tmp, na.rm=TRUE))
  }
  return(base::mean(rr.here, na.rm=TRUE))
}

# -------------------------------------------------------------------------

#' Average RR Surface Raster
#'
#' @param inRas Input DEM. (SpatRast)
#' @param inMn Minimum scale. (default = 7)
#' @param inMx Maximum scale. (default = 21)
#' @param stp Scale step. (default = 2)
#'
#' @return Average relative relief surface raster. (SpatRast)
#' @export
#'
rr.Average <- function(inRas, inMn = 7, inMx = 21, stp = 2){

  x <- base::as.matrix(inRas, wide = TRUE)

  # Sequence of scales between min and max bounds
  seq.RR <- base::seq(inMn, inMx, stp)

  # Offset from window center based on scale
  m = (seq.RR - 1) / 2
  # Central cell offset
  mid = ((inMx - 1) / 2)

  nr = base::nrow(x) - mid
  nc = base::ncol(x) - mid

  out <- x * NA

  arr <- base::array(NA, dim = c(base::nrow(inRas), base::ncol(inRas), base::length(seq.RR)))
  for(i in 1:base::length(seq.RR)){
    base::cat(base::paste0('Scale: ', seq.RR[i]), '\r')
    tmp <- relRelief(inRas = inRas, sc = seq.RR[i])
    tmp <- base::as.matrix(tmp, wide = TRUE)
    arr[,,i] <- tmp
  }

  avg.RR <- base::apply(arr, c(1,2), mean, na.rm=F)

  arr.Out <- base::array(NA, dim = c(base::nrow(inRas), base::ncol(inRas), base::length(seq.RR) + 1))

  arr.Out[,,1:base::length(seq.RR)] <- arr
  arr.Out[,,base::length(seq.RR) + 1] <- avg.RR

  arr.Out <- terra::rast(arr.Out, extent = terra::ext(inRas), crs = terra::crs(inRas))

  return(arr.Out)
}

# -------------------------------------------------------------------------

#' Single scale RR at a cell/point
#'
#' @param inRas Input DEM. (matrix NOT raster)
#' @param in.S Scale.
#' @param ii i index.
#' @param jj j index.
#' @param thrsh Threshold elevation to calculate RR above.
#'
#' @return RR at cell/point for given scale.
#' @export
#'
rr.Point <- function(inRas, in.S = NA, ii = NA, jj = NA, thrsh = 0){
  # print(base::class(inRas))
  if(base::class(inRas)[1] == 'SpatRaster'){
    inRas <- base::as.matrix(inRas, wide = TRUE)
  }else{
    inRas <- inRas
  }

  # Using the scale input create the window parameters
  m = in.S - 1
  mid = m / 2

  # Check if the point is NA, if not then calculate RR otherwise return NA
  if(!base::is.na(inRas[ii, jj]) & inRas[ii, jj] >= thrsh){
    mat <- c(inRas[((ii-mid):(ii+mid)), ((jj-mid):(jj+mid))]) # Window creation

    return((inRas[ii, jj] - base::min(mat, na.rm=TRUE)) /
             (base::max(mat, na.rm=TRUE) - base::min(mat, na.rm=TRUE))) # RR Calc

  }else{
    return(NA)
  }

}

# -------------------------------------------------------------------------

#' RR of single scale along profile
#'
#' @param inRas Input DEM. (matrix NOT SpatRast)
#' @param in.Sc Scale.
#' @param in.ID ID value of profile.
#' @param thresh Threshold elevation to calculate RR above.
#' @param sea Direction/orientation of the sea/ocean. ("N", "E", "S", or "W")
#' @param chk Where to start calculating.
#'
#' @return RR vector at scale S for the given profile ID
#' @export
#'
rr.Vec <- function(inRas, in.Sc = 7, in.ID = NA, thresh = 0, sea = NA, chk = 0){
  if(base::class(inRas)[1] == 'SpatRaster'){
    inRas <- base::as.matrix(inRas, wide = TRUE)
  }else{
    inRas <- inRas
  }

  if(!base::is.na(sea)){
    if(sea == "N" || sea == "S"){
      out.V <- NA
      ct <- 1
      for(i in (chk + 1):(base::nrow(inRas) - chk)){

        out.V[ct] <- rr.Point(inRas = inRas,
                              in.S = in.Sc,
                              ii = i,
                              jj = in.ID + chk + 1,
                              thrsh = thresh)
        ct <- ct + 1
      }

    }else if(sea == "W" || sea == "E"){
      out.V <- NA
      ct <- 1
      for(i in (chk + 1):(base::ncol(inRas) - chk)){
        out.V[ct] <- rr.Point(inRas = inRas,
                              in.S = in.Sc,
                              ii = in.ID + chk + 1,
                              jj = i,
                              thrsh = thresh)
        ct <- ct + 1
      }

    }else{
      return(base::print("Error: 'sea' must be either N, S, E, or W"))
    }
  }else{
    return(base::print("Error: 'sea' must be either N, S, E, or W, cannot be NA"))
  }
  return(out.V)
}

# -------------------------------------------------------------------------

#' RR of random profiles at all scales
#'
#' @param inRas Input DEM. (matrix NOT SpatRast)
#' @param mn.S Min scale.
#' @param mx.S Max scale.
#' @param stp Scale sequence step. (default = 2)
#' @param samp.P Vector of the random profiles sampled.
#' @param thresh Threshold elevation to calculate RR above.
#' @param sea Direction/orientation of the sea/ocean. ("N", "E", "S", or "W")
#'
#' @return List containing 2 Lists:
#' Out1 = List of dataframes containing RR for all scales at random profiles
#' Out2 = List of dataframes containing "moving" Avg of RRs (moving avg. based
#' scale i.e., avg 3-5, 3-7, 3-9, etc.)
#' Out 3 = Sample profile.
#' @export
#'
rr.Prof.Vecs <- function(inRas, mn.S, mx.S, stp = 2, samp.P, thresh = 0, sea = NA){
  if(base::class(inRas) == 'SpatRaster'){
    inM <- base::as.matrix(inRas, wide = TRUE)
  }else{
    inM <- inRas
  }

  # Set warnings off because NA values return redundant warnings
  defaultW <- base::getOption("warn")
  base::options(warn=-1)
  # Create sequence of Scales
  seq.S <- base::seq(mn.S, mx.S, stp)
  # Create names for dataframe later
  names <- "z"
  names <- base::append(names, base::as.character(seq.S))

  # Create a new matrix with an offset of NA values on the border
  na.mat <- base::matrix(NA, nrow=base::nrow(inM)+(2 * mx.S), ncol=base::ncol(inM) + (2 * mx.S))

  na.mat[((mx.S+1):(base::nrow(na.mat)-mx.S)),((mx.S+1):(base::ncol(na.mat)-mx.S))] <- inM
  inEl2 <- na.mat

  # Initialize final list for RAW RR values ====================================
  out.List1 <- base::list()
  # Need to loop through all scales and then all profiles
  for(i in 1:base::length(samp.P)){ # Loop through profile IDs
    base::cat("\r", i, "/", base::length(samp.P), "   |   ", "Computing RRs for Profile: ",
              samp.P[i],  sep="")
    utils::flush.console()
    if(!base::is.na(sea)){
      # TODO: Check if S vs N matters
      if(sea == "N" || sea == "S"){
        # Initialize dataframe as a matrix that will fill the list later
        mat.S <- base::matrix(NA, nr = base::nrow(inM), nc = base::length(names))
        # Fill the first column with elevation values for this profile
        mat.S[,1] <- inM[,samp.P[i]]

        for(j in 1:base::length(seq.S)){ # Loop through scales
          # Fill the remainder of the columns with RR values corresponding to scale
          mat.S[, j + 1] <- rr.Vec(inRas=inEl2, in.Sc=seq.S[j], in.ID=samp.P[i], thresh=thresh, sea=sea, chk=mx.S)
        }
      }else if(sea == "W" || sea == "E"){
        # Initialize dataframe as a matrix that will fill the list later
        mat.S <- base::matrix(NA, nr=base::ncol(inM), nc=base::length(names))
        # Fill the first column with elevation values for this profile
        mat.S[,1] <- inM[samp.P[i],]

        for(j in 1:base::length(seq.S)){ # Loop through scales
          # Fill the remainder of the columns with RR values corresponding to scale
          mat.S[, j + 1] <- rr.Vec(inRas=inEl2, in.Sc=seq.S[j], in.ID=samp.P[i], thresh=thresh, sea=sea, chk=mx.S)
        }
      }else{
        return(base::print("Error: 'sea' must be either N, S, E, or W"))
      }
    }else{
      return(base::print("Error: 'sea' must be either N, S, E, or W, cannot be NA"))
    }

    # Convert to data frame
    mat.S <- base::as.data.frame(mat.S)
    # Set column names
    base::names(mat.S) <- names
    # Store data frame in the output list
    out.List1[[i]] <- mat.S
  }
  # out.List1[[1]]
  # ============================================================================
  # Initialize list for Averaged RR Values =====================================
  out.List2 <- base::list()

  base::cat("\n", "Calculating Average RR List...", sep="")
  utils::flush.console()
  # Need to loop through all scales and then all profiles
  for(i in 1:base::length(out.List1)){

    # Grab a data frame of RRs and convert to matrix
    mat.R <- base::as.matrix(out.List1[[i]])
    # Initialize temporary matrix
    tmp.Mat <- mat.R
    # Calculate "moving" average of RRs
    for(k in 3:base::ncol(tmp.Mat)){
      for(j in 1:base::nrow(tmp.Mat)){
        tmp.Mat[j,k] <- base::mean(mat.R[j,2:k],na.rm=T)
      }
    }
    # Convert to data frame
    tmp.Mat <- base::as.data.frame(tmp.Mat)
    # Set column names
    base::names(tmp.Mat) <- names

    # Store data frame in the output list
    out.List2[[i]] <- tmp.Mat
  }
  # out.List2[[1]]
  # ============================================================================

  base::options(warn=defaultW)
  return(base::list(out.List1, out.List2, samp.P))

}

# -------------------------------------------------------------------------

#' Plot multiple RR and moving avg. RR
#'
#' @param rr.L Input RR data. (List of lists)
#'
#' @return List of lists of local minima.
#' @export
#'
rr.Mins <- function(rr.L){

  # Initialize Out lists
  out.L1 <- base::list()
  out.L2 <- base::list()
  for(i in 1:base::length(rr.L[[1]])){
    # i = 1
    # Separate the two different lists
    rr.Raw <- rr.L[[1]][[i]] # Raw RR values
    rr.Avg <- rr.L[[2]][[i]] # Avg RR Values
    SS <- rr.L[[3]][[i]]

    # Crest inflection
    topsZ <- base::lapply(1:5, function(x) inflect(rr.Raw$z, threshold = x)$maxima)

    # Fill in the crest locations
    cc <- rr.Raw$z
    cc <- cc * NA
    cc[topsZ[[5]]] <- 1

    # Filter out crest locations below mean elevation
    cc[rr.Raw$z < (base::mean(rr.Raw$z, na.rm=T))] <- NA

    # Find distance locations for each crest location
    tmp.c.dist <- cc * c(1:base::length(cc))
    dd <- c(1:base::length(cc))

    # Set the most seaward location to -1 and filter the remaining locations out
    cc[dd == base::max(tmp.c.dist, na.rm=T)] <- -1
    cc[cc != -1 | base::is.na(cc)] <- 0
    cc <- base::abs(cc)
    cc[base::which(cc == 0)] <- NA


    # Crest Location
    C.x <- base::which(!base::is.na(cc))
    # print(C.x)
    # set values C.x to beginning of profile to NA
    if(base::length(C.x) > 0){
      nnn <- 0
      rr.Raw[1:(C.x-nnn),] <- NA
      rr.Avg[1:(C.x-nnn),] <- NA
    }

    # Initialize the output data frames
    raw.Mat <- rr.Raw
    avg.Mat <- rr.Avg
    # Loop through RRs and remove points that are not the min RR
    for(j in 2:base::ncol(raw.Mat)){
      # j = 2
      # Initialize vectors
      raw.V <- raw.Mat[,j]
      avg.V <- avg.Mat[,j]

      raw.V[base::which(raw.V != base::min(raw.V, na.rm=T))] <- NA
      avg.V[base::which(avg.V != base::min(avg.V, na.rm=T))] <- NA

      raw.V[!base::is.na(raw.V)] <- 1
      avg.V[!base::is.na(avg.V)] <- 1

      raw.Mat[,j] <- raw.V
      avg.Mat[,j] <- avg.V
    }
    # raw.Mat[151,] <- NA
    out.L1[[i]] <- raw.Mat
    out.L2[[i]] <- avg.Mat
  }

  return(base::list(out.L1, out.L2))
}

# -------------------------------------------------------------------------

#' Rework RR.mn DFs for Interpretation
#'
#' @param inMn Input RR.mn data. (List of lists)
#' @param inRR Input RR data. (List of lists)
#'
#' @return Simplified minima location. (List of lists)
#' @export
#'
rr.Mins.Analysis <- function(inMn, inRR){
  # Separate data
  rr.Raw <- inMn[[1]] # Raw RR minima values
  rr.Avg <- inMn[[2]] # Avg RR minima Values
  SS <- inRR[[3]]     # Profile ID


  # Names vector for later
  nms <- c("Scale", "Sx", "Sz", "Cx", "Cz")

  # Initialize Out lists
  out.L1 <- base::list()
  out.L2 <- base::list()

  # Initialize profile vector stuff
  p.out <- NA
  ct <- 1

  # Loop through the lists for each profile
  for(i in 1:base::length(rr.Raw)){
    # i = 2
    rr.Raw2 <- inRR[[1]][[i]] # Raw RR values
    # Crest inflection
    topsZ <- base::lapply(1:5, function(x) inflect(rr.Raw2$z, threshold = x)$maxima)

    # Fill in the crest locations
    cc <- rr.Raw2$z
    cc <- cc * NA
    cc[topsZ[[5]]] <- 1

    # Filter out crest locations below mean elevation
    cc[rr.Raw2$z < (base::mean(rr.Raw2$z, na.rm=T))] <- NA

    # Find distance locations for each crest location
    tmp.c.dist <- cc * c(1:base::length(cc))
    ddd <- c(1:base::length(cc))
    # Set the most seaward location to -1 and filter the remaining locations out
    cc[ddd == base::max(tmp.c.dist, na.rm=T)] <- -1
    cc[cc != -1 | base::is.na(cc)] <- 0
    cc <- base::abs(cc)
    cc[base::which(cc == 0)] <- NA
    ccz <- cc * rr.Raw2$z

    # Crest Location
    C.x <- base::which(!base::is.na(cc))
    if(base::length(C.x) > 0){
      C.x <- C.x
    }else{
      C.x <- 100
    }
    C.z <- ccz[base::which(!base::is.na(ccz))]
    if(base::length(C.z) > 0){
      C.z <- C.z
    }else{
      C.z <- 100
    }
    # # Grab the subset profile for the crest location and Elevation
    # tmp.P <- inP[inP$ID == SS[i],]
    # C.x <- tmp.P$c * tmp.P$dis
    # C.x <- C.x[which(!is.na(C.x))]
    # C.z <- tmp.P$c * tmp.P$z
    # C.z <- C.z[which(!is.na(C.z))]

    # Grab the RR minima for this profile
    tmp.Raw <- rr.Raw[[i]]
    tmp.Avg <- rr.Avg[[i]]

    scales <- base::colnames(tmp.Raw)
    scales <- scales[-1]
    scales <- base::as.numeric(scales)

    # Initialize new matrices for output
    out.Raw <- base::matrix(NA, nrow=(base::ncol(tmp.Raw) - 1), ncol = 5)
    out.Avg = out.Raw
    dd <- c(1:base::length(tmp.Raw[,1]))
    # Loop through all scales and find x and z values of minima
    for(j in 2:base::ncol(tmp.Raw)){
      # Find x and z of minima for scale j
      # j = 4
      S.x.Raw <- tmp.Raw[,j] * dd
      S.z.Raw <- tmp.Raw[,j] * tmp.Raw[,1]

      S.x.Raw <- base::min(S.x.Raw[base::which(!base::is.na(S.x.Raw))])
      S.z.Raw <- S.z.Raw[S.x.Raw]

      S.x.Avg <- tmp.Avg[,j] * dd
      S.z.Avg <- tmp.Avg[,j] * tmp.Avg[,1]

      S.x.Avg <- base::min(S.x.Avg[base::which(!base::is.na(S.x.Avg))])
      S.z.Avg <- S.z.Avg[S.x.Avg]

      dis.tmp <- tmp.Raw[,j] * dd
      dis.tmp <- dis.tmp[!base::is.na(dis.tmp)]

      # Create vector to fill the out matrix
      fill.Vec.Raw <- c(scales[j-1], S.x.Raw, S.z.Raw, C.x, C.z)
      fill.Vec.Avg <- c(scales[j-1], S.x.Avg, S.z.Avg, C.x, C.z)

      # Fill the out matrices by row
      out.Raw[j-1,] <- fill.Vec.Raw
      out.Avg[j-1,] <- fill.Vec.Avg

    }

    # Convert to data frame and set column names
    out.Raw <- base::as.data.frame(out.Raw)
    out.Avg <- base::as.data.frame(out.Avg)

    base::names(out.Raw) <- nms
    base::names(out.Avg) <- nms

    out.Raw$Dx <- base::abs(out.Raw$Sx - out.Raw$Cx)
    out.Raw$Dz <- base::abs(out.Raw$Sz - out.Raw$Cz)

    out.Avg$Dx <- base::abs(out.Avg$Sx - out.Avg$Cx)
    out.Avg$Dz <- base::abs(out.Avg$Sz - out.Avg$Cz)

    if(getmode(out.Raw$Sx)==base::length(dd)){
      # out.Raw$Sx[out.Raw$Sx==getmode(out.Raw$Sx)] <- NA
      p.out[ct] <- i
      ct <- ct + 1
    }

    # Fill output lists
    out.L1[[i]] <- out.Raw
    out.L2[[i]] <- out.Avg
  }

  return(base::list(out.L1, out.L2, p.out))

}

# -------------------------------------------------------------------------

#' Local extrema
#'
#' @param x Vector of values.
#' @param threshold Sensitivity threshold. (best results = 2-6)
#'
#' @return Extrema information.
#' @export
#'
inflect <- function(x, threshold = 1){
  up   <- base::sapply(1:threshold, function(n) c(x[-(base::seq(n))], base::rep(NA, n)))
  down <-  base::sapply(-1:-threshold, function(n) c(base::rep(NA,base::abs(n)),
                                                     x[-base::seq(base::length(x), base::length(x) - base::abs(n) + 1)]))
  a    <- base::cbind(x,up,down)
  base::list(minima = base::which(base::apply(a, 1, min) == a[,1]),
             maxima = base::which(base::apply(a, 1, max) == a[,1]))
}

# -------------------------------------------------------------------------

#' Get mode
#'
#' @param v Numeric vector.
#'
#' @return mode
#' @export
#'
getmode <- function(v) {
  uniqv <- base::unique(v)
  uniqv[base::which.max(base::tabulate(base::match(v, uniqv)))]
}

# -------------------------------------------------------------------------
#' RR Average
#'
#' @param inMn Minimum cale
#' @param inMx Maximum scale
#' @param inRas Input DEM
#'
#' @return Average RR surface
#' @export
#'
rr.AverageLayer <- function(inMn=7, inMx=Id.Sc, x=D, inRas=M, save=T){
  if(base::class(inRas) == 'SpatRaster'){
    x <- base::as.matrix(inRas, wide=TRUE)
  }else{
    x <- inRas
  }

  seq.RR <- base::seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = base::nrow(x) - mid
  nc = base::ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!base::is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:base::length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- (x[i,j] - base::min(tmp, na.rm=T))/(base::max(tmp, na.rm=T) - base::min(tmp, na.rm=T))
        }
        out[i,j] = base::mean(rr.here, na.rm=T)

      }
    }
  }
  out <- terra::rast(out, ext=terra::ext(inRas), crs=terra::crs(inRas))

  return(out)
}

# -------------------------------------------------------------------------

#' TPI Average
#'
#' @param inMn Minimum cale
#' @param inMx Maximum scale
#' @param inRas Input DEM
#'
#' @return Average TPI surface
#' @export
#'
tpi.Average <- function(inMn, inMx, inRas){
  if(base::class(inRas) == 'SpatRaster'){
    x <- base::as.matrix(inRas, wide=TRUE)
  }else{
    x <- inRas
  }

  seq.RR <- base::seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = base::nrow(x) - mid
  nc = base::ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!base::is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:base::length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- (x[i,j] - base::mean(tmp, na.rm=T))
        }
        out[i,j] = base::mean(rr.here, na.rm=T)

      }
    }
  }
  out <- terra::rast(out, ext=terra::ext(inRas), crs=terra::crs(inRas))

  return(out)
}

# -------------------------------------------------------------------------
#' TRI Average
#'
#' @param inMn Minimum cale
#' @param inMx Maximum scale
#' @param inRas Input DEM
#'
#' @return Average TRI surface
#' @export
#'
tri.Average <- function(inMn, inMx, inRas){
  if(base::class(inRas) == 'SpatRaster'){
    x <- base::as.matrix(inRas, wide=TRUE)
  }else{
    x <- inRas
  }
  seq.RR <- base::seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = base::nrow(x) - mid
  nc = base::ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!base::is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:base::length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- base::sum(base::abs(x[i,j] - tmp), na.rm=T)/(base::length(tmp)-1)
        }
        out[i,j] = base::mean(rr.here, na.rm=T)

      }
    }
  }

  out <- terra::rast(out, ext=terra::ext(inRas), crs=terra::crs(inRas))
  return(out)
}
