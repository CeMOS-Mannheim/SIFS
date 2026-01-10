#' Superimpose MSI pixels on an optical plot
#'
#' This functions plots MSI measurement locations (pixels) on top of the
#' corresponding (reference) optical image for visual inspection.
#'
#' @param msiTeachingImgPath  a character, path to the optical image used for MSI teaching. 
#' @param targetImgPath  a character, path to the target optical image (e.g., H&E). 
#' @param regions a `regions` object, see `?getRegions`.
#' @param mpm an `mpm` object. 
#' @param transMat a 3x3 transformation matrix for projecting points from `msiTeachingImg` (moving)
#' to `targetImg` (fixed).
#' @param InverseTransMat a logical, whether to use the inverse of `transMat`. 
#' @param outputPath a character, path for writing the resulting images. 
#' @param plotPoints logical, whether to plot the points and their intensities specified in `mpm`.
#' @param plotContours a logical (ignored for now), whether to plot the hotspot/coldspot contours of the
#' `mpm` object. 
#' @param whichContours a character, `c("hotspot", "coldspot", "both")`. 
#' @param smooth a logical, whether to fill the hotspot/coldspot contours with a smoothed
#' version of intensities computed via the KDE. 
#' @param col.hotspot character, the color to be used for hotspots contours
#' @param col.coldspot character, the color to be used for coldspot contours
#' @param lwd a numeric, line width of the contours
#' @param pch an integer, the points type when `plotPoints=TRUE`. 
#' @param point.cex a numeric, the size of the points when `plotPoints=TRUE`.
#' @param alpha a fraction, alpha chanel for the raste points when `plotPoints=TRUE`.
#' @param label.cex a numeric, the size of the region label text.
#' @param xlim a numeric vector specifying `xlim` in msi coordinates.
#' @param ylim a numeric vector specifying `ylim` in msi coordinates.
#' @param ... ignored.
#'
#'
#' @return
#' Generates plots and writes them to `outputPath`. 
#'
#' @export
#'
#' @author Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
#'


plotOptical3 <- function(msiTeachingImgPath,
                         targetImgPath,
                         regions,
                         mpm,
                         transMat,
                         InverseTransMat = FALSE,
                         outputPath,
                         plotPoints = TRUE,
                         plotContours = TRUE,
                         whichContours = "both",
                         smooth = FALSE,
                         col.hotspot = "red",
                         col.coldspot = "blue",
                         lwd = 2,
                         pch = 19,
                         point.cex = 0.1,
                         alpha = 0.5,
                         label.cex = 1,
                         xlim = NULL,
                         ylim = NULL,
                         colPal = "turbo",
                         ...){
  
  
  
  
  
  
  if(!file.exists(msiTeachingImgPath)){
    stop("msiTeachingImg file does not exist!\n")
  }
  
  if(!file.exists(targetImgPath)){
    stop("targetImg file does not exist!\n")
  }
  
  if(class(regions) != "regions"){
    stop("regions object is not of class 'regions'. Consider calling msiImporter::getRegions.\n")
  }
  
  if(class(mpm) != "molProbMap"){
    stop("mpm object is not of class 'molProbMap'\n")
  }
  
  if(prod(dim(transMat)) != 9){
    stop("transMat object must be a 3x3 matrix.\n")
  }
  
  
  
  # col <- col2rgb(col, TRUE) / 255
  # col <- apply(col, 2, FUN = function(x){
  #   rgb(red = x[1], 
  #       green = x[2], 
  #       blue = x[3], 
  #       alpha = alpha, 
  #       maxColorValue = 1)
  # })
  
  cat("reading msiTeachingImg .. \n")
  
  msiTeachingImg <- magick::image_read(msiTeachingImgPath)
  msiTeachingImg <- as.raster(msiTeachingImg)
  msiTeachingImg <- msiTeachingImg[nrow(msiTeachingImg):1, ]
  cat("msiTeachingImg dimensions = ", dim(msiTeachingImg), "\n")
  
  cat("reading targetImg .. \n")
  
  targetImg <- magick::image_read(targetImgPath)
  targetImg <- as.raster(targetImg)
  targetImg <- targetImg[nrow(targetImg):1, ]
  cat("targetImg dimensions = ", dim(targetImg), "\n")
  
  # specify xlim
  if(!is.null(xlim)){
    xlim <- sort(xlim)
    
    msiXmin <- data.frame(x = xlim[1], y = min(regions$coordsSummary$ymin))
    msiXmax <- data.frame(x = xlim[2], y = max(regions$coordsSummary$ymax))
    
    # transofrm
    opticalXmin <- regions$msiToOpticalTransFun(msiXmin)
    opticalXmax <- regions$msiToOpticalTransFun(msiXmax)
    
    opticalXlim <- c(opticalXmin$x, opticalXmax$x)
    
    #cat("optical pixel xlim = ", opticalXlim, "\n")
    
  } else{
    opticalXlim <- c(0, ncol(msiTeachingImg))
  }
  
  # specify ylim
  if(!is.null(ylim)){
    ylim <- sort(ylim)
    
    msiYmin <- data.frame(x = min(regions$coordsSummary$xmin), y = ylim[1])
    msiYmax <- data.frame(x = max(regions$coordsSummary$xmax), y = ylim[2])
    
    # transofrm
    opticalYmin <- regions$msiToOpticalTransFun(msiYmin)
    opticalYmax <- regions$msiToOpticalTransFun(msiYmax)
    
    opticalYlim <- c(opticalYmin$y, opticalYmax$y)
    
    #cat("optical pixel ylim = ", opticalYlim, "\n")
    
  } else{
    opticalYlim <- c(0, nrow(msiTeachingImg))
  }
  
  
  # extract the object of interest
  if(smooth){
    sppimg <- mpm$rhoMoi * spatstat.geom::as.im(mpm$hotspotMask)
  } else{
    #sppimg <- moleculaR::spp2im(mpm$hotspotpp, zero.rm = T)
    sppimg <- moleculaR::spp2im(mpm$sppMoi, zero.rm = F)
  }
  
  # transform to optical
  cat("applying the MSI -> msiTeachingImg transformation .. \n")
  
  cat("range original spp in x = ", range(mpm$sppMoi$x), "\n")
  cat("range original spp in y = ", range(mpm$sppMoi$y), "\n\n")
  cat("range original img in x = ", sppimg$xrange, "\n")
  cat("range original img in y = ", sppimg$yrange, "\n\n")
  
  spp <- .transformspp(mpm$sppMoi, transFun = misRegion$msiToOpticalTransFun)
  sppimg <- .transformImg(sppimg, misRegion$msiToOpticalTransFun)
  
  cat("range transformed to optical spp in x = ", range(spp$x), "\n")
  cat("range transformed to optical spp in y = ", range(spp$y), "\n\n")
  cat("range transformed to optical img in x = ", sppimg$xrange, "\n")
  cat("range transformed to optical img in y = ", sppimg$yrange, "\n\n")
  
  tiff(file.path(outputPath, "msiTeachingImg.tiff"), width = 15, height = 15, res = 300,
       compression = "lzw", units = "cm")

  plot(msiTeachingImg, xlim = opticalXlim, ylim = rev(opticalYlim))
  
  moleculaR::plotImg(sppimg, transpFactor = alpha, add = TRUE,
                     colourPal = colPal, smooth = !smooth)
  
  dev.off()


  cat("msiTeachingImg written to ", outputPath, "\n")


  # transform optical to targetImage
  cat("applying the msiTeachingImg -> targetImg transformation .. \n")

  if(InverseTransMat){
    transMat <- solve(transMat)
  }
  
  spp <- .transformspp(spp, transFun = .transmat2fun,  transMat = transMat)
  sppimg <- .transformImg(sppimg, transFun = .transmat2fun, transMat = transMat)
  #sppimg <- spatstat.geom::affine.im(sppimg, mat = t(tf[1:2,1:2]), vec = c(tf[2:1, 3]))
  
  cat("range transformed to optical spp in x = ", range(spp$x), "\n")
  cat("range transformed to optical spp in y = ", range(spp$y), "\n\n")
  cat("range transformed to optical img in x = ", sppimg$xrange, "\n")
  cat("range transformed to optical img in y = ", sppimg$yrange, "\n\n")

  tiff(file.path(outputPath, "targetImg.tiff"), width = 15, height = 15, res = 300,
       compression = "lzw", units = "cm")
  plot(targetImg, ylim = c(nrow(targetImg), 0))
  moleculaR::plotImg(sppimg, transpFactor = alpha, add = TRUE,
                     colourPal = colPal, smooth = !smooth)
  dev.off()

  cat("targetImg written to ", outputPath, "\n")


  
}

.transformWindow <- function(win, transFun, ...){
  
  xyRange <- data.frame(x = win$xrange, y = win$yrange)
  xyRange <- transFun(xyRange, ...)
  
  win$xrange <- xyRange$x
  win$yrange <- xyRange$y
  
  bdryCoords <-  lapply(win$bdry, FUN = function(x){
    tmp <- as.data.frame(x)
    tmp <- transFun(tmp, ...)
    
    as.list(tmp)
  })
  
  win$bdry <- bdryCoords
  
  return(win)
}

.transformspp <- function(spp, transFun, ...){
  
  xyCoords <- data.frame(x = spp$x, y = spp$y)
  xyCoords <- transFun(xyCoords, ...)
  
  spp$x <- xyCoords$x
  spp$y <- xyCoords$y
  
  spp$window <- .transformWindow(spp$window, transFun, ...)
  
  return(spp)
}



.transformImg <- function(img, transFun,  ...){



  
    xCoords <- data.frame(x = img$xcol , y = rep(100, length(img$xcol)))
    yCoords <- data.frame(x = rep(100 ,length(img$yrow)), y = img$yrow)
    
  
  
  #xyRange <- transFun(xyRange, ...)
  xCoords <- transFun(xCoords, ...)
  yCoords <- transFun(yCoords, ...)

  xrange <- range(xCoords$x)
  yrange <- range(yCoords$y)

  xstep <- diff(xrange)/length(xCoords$x)
  ystep <- diff(yrange)/length(yCoords$y)
  
  cat("calculated xstep = ", xstep, " | ystep = ", ystep, "\n")
  
  # 
  # img <- spatstat.geom::im(mat = img$v,
  #                          xcol = xCoords$x,
  #                          yrow = yCoords$y,
  #                          xrange = xrange,
  #                          yrange = yrange)#

  img$xrange <- xrange 
  img$yrange <- yrange 
  img$xstep  <- xstep
  img$ystep  <- ystep
  img$xcol   <- xCoords$x 
  img$yrow   <- yCoords$y 


  return(img)
}

.transmat2fun <- function(coords, transMat){
  
  # if(class(coords) != "data.frame" | dim(coords)[2] != 2 )
  # {
  #   stop("error in .transmat2fun. coords must be a data.frame with 2 columns 
  #        representing x and y coordinates, respectively.\n")
  # }
  # 
  # if(class(transMat) != "matrix" | dim(transMat)[2] != 3 | dim(transMat)[1] != 3)
  # {
  #   stop("error in .transmat2fun. transMat must be a 3x3 matrix.\n")
  # }
  
  #// convert to homogeneous coordinates
  coords[ , 3]      <- 1
  coords            <- t(coords)
  
  
  newCoords        <- round(transMat %*% coords)
  
  
  newCoords        <- data.frame(x = newCoords[1, ],
                                    y = newCoords[2, ])
  
  return(newCoords)
  
}


safeDevCapabilities <- function() {
  a <- try(dev.capabilities())
  if(!inherits(a, "try-error"))
    return(a)
  warning("dev.capabilities() caused an error!", call.=FALSE)
  return(NULL)
}