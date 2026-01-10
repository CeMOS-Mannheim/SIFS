#' Superimpose MSI pixels on an optical plot
#'
#' This functions plots MSI measurement locations (pixels) on top of the
#' corresponding (reference) optical image for visual inspection.
#'
#' @param msiTeachingImgPath  a character, path to the optical image used for MSI teaching.
#' @param targetImgPath  a character, path to the target optical image (e.g., H&E).
#' @param regions a `regions` object, see `?getRegions`.
#' @param coords a list of data.frames representing the different classes/regions
#' to be transfomed, each with 2 columns, x and y.
#' @param transMat a 3x3 transformation matrix for projecting points from `msiTeachingImg` (moving)
#' to `targetImg` (fixed).
#' @param InverseTransMat a logical, whether to use the inverse of `transMat`.
#' @param outputPath a character, path for writing the resulting images.
#' @param raster logical, whether to plot the points specified in `coords`.
#' @param polygons a logical (ignored for now), whether to compute a plot the boundaries of the point
#' specified in `coords`.
#' @param polygons.fill a logical (ignored for now), whether to fill the computed boundaries of the point
#' specified in `coords`.
#' @param col character, the color to be used for points/polygons.
#' @param lwd a numeric, line width of polygons.
#' @param pch an integer, the points type.
#' @param point.cex a numeric, the size of the points for the raster.
#' @param alpha a fraction, alpha chanel for the raste points.
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


plotOptical2 <- function(msiTeachingImgPath,
                         optical_imgName,
                        targetImgPath,
                        HE_imgName,
                        regions,
                        coords,
                        transMat,
                        InverseTransMat = FALSE,
                        outputPath,
                        raster = TRUE,
                        polygons = FALSE,
                        polygons.fill = FALSE,
                        col = "black",
                        lwd = 2,
                        pch = 19,
                        point.cex = 0.1,
                        alpha = 0.5,
                        label.cex = 1,
                        xlim = NULL,
                        ylim = NULL,
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

  if(!is.list(coords)){
    stop("coords object must be a list of data.frames, each with 2 columns, x and y.\n")
  }

  if(is.data.frame(coords)){
    stop("coords object must be a list of data.frames representing the different classes/regions to be transfomed, each with 2 columns, x and y.\n")
  }

  if(prod(dim(transMat)) != 9){
    stop("transMat object must be a 3x3 matrix.\n")
  }

  if((length(col) == 1) & (length(coords) > 1)){
    col <- rep(col, length(coords))
  }

  if((length(col) != 1) & (length(coords) != length(col))){
    stop("length of provided colors vector 'col' must have the same length as the 'coords' object.\n")
  }

  col <- col2rgb(col, TRUE) / 255
  col <- apply(col, 2, FUN = function(x){
    rgb(red = x[1],
        green = x[2],
        blue = x[3],
        alpha = alpha,
        maxColorValue = 1)
  })

  cat("reading msiTeachingImg .. \n")

  # msiTeachingImg <- magick::image_read(msiTeachingImgPath) # image size prevents this line to work, update it: 24.04.2024
  ## ___________________________________________________________________________
  # UPDATED: 24.04.2024
  # library(imager)
  #
  # # Function to read image without size limitations
  # readImage <- function(filepath) {
  #       img <- imager::load.image(filepath)
  #       return(img)
  # }
  # Example usage
  # msiTeachingImg <- readImage(msiTeachingImgPath)
  ## ___________________________________________________________________________

  # get the file extension:
  # fileExtension <- strsplit(basename(msiTeachingImgPath), split = ".", fixed = T)[[1]][2]
  # ifelse(fileExtension == "jpg", yes = msiTeachingImg <-  magick::image_read(msiTeachingImgPath), no <-  msiTeachingImg <- tiff::readTIFF(msiTeachingImgPath, native = F))

  # # choose between jpg or tif files:
  # if(fileExtension == "jpg"){
  #       msiTeachingImg <- magick::image_read(msiTeachingImgPath)
  # } else if(fileExtension == "tif"){
  #       # using tiff function directly: 25.04.2024
  #       msiTeachingImg <- tiff::readTIFF(msiTeachingImgPath, native = F)
  #
  #       }


  # only use original tiff files: 25.04.2024
  msiTeachingImg <- tiff::readTIFF(msiTeachingImgPath, native = F)

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


  # transform to optical
  cat("applying the MSI -> msiTeachingImg transformation .. \n")

  coordOpt <- lapply(coords, FUN = function(x){
    tmp <- data.frame(x = x$x, y = x$y)
    misRegion$msiToOpticalTransFun(tmp)
  })

  tiff(file.path(outputPath, paste0(optical_imgName,"msiTeachingImg.tiff")), width = 15, height = 15, res = 300,
       compression = "lzw", units = "cm")
  plot(msiTeachingImg, xlim = opticalXlim, ylim = rev(opticalYlim))
  for(i in 1:length(coordOpt)){
    points(coordOpt[[i]],
           pch = pch, col = col[i], cex = point.cex)

    # if(polygons){
    #   cat("creating polygonal regions: ", i, " of ", length(coordOpt) , "\n")
    #   pol <- spatstat.geom::owin(mask = data.frame(x = coordOpt[[i]]$x,
    #                                                y = coordOpt[[i]]$y))
    #   pol <- spatstat.geom::as.polygonal(W = pol)
    #
    #   if(polygons.fill){
    #     plot(pol, add = TRUE, col = col[i])
    #   } else {
    #     plot(pol, add = TRUE, border = col[i])
    #   }
    # }

  }
  dev.off()

  cat("msiTeachingImg written to ", outputPath, "\n")


  # transform optical to targetImage
  cat("applying the msiTeachingImg -> targetImg transformation .. \n")

  if(InverseTransMat){
    transMat <- solve(transMat)
  }

  coordTarget <- lapply(coordOpt, FUN = function(i){

    tmp <- transMat %*% rbind(i$x, i$y, rep(1, nrow(i)))
    data.frame(x = tmp[1, ], y = tmp[2, ])

  })



  tiff(file.path(outputPath, paste0(HE_imgName,"targetImg.tiff")), width = 15, height = 15, res = 300,
       compression = "lzw", units = "cm")
  plot(targetImg, ylim = c(nrow(targetImg), 0))
  for(i in 1:length(coordTarget)){
    points(coordTarget[[i]],
           pch = pch, col = col[i], cex = point.cex)

    # if(polygons){
    #   cat("creating polygonal regions: ", i, " of ", length(coordTarget) , "\n")
    #   pol <- spatstat.geom::owin(mask = data.frame(x = coordTarget[[i]]$x,
    #                                                y = coordTarget[[i]]$y))
    #   #pol <- spatstat.geom::as.polygonal(W = pol)
    #
    #   if(polygons.fill){
    #     plot(pol, add = TRUE, col = col[i])
    #   } else {
    #     plot(pol, add = TRUE, border = col[i])
    #   }
    # }


  }
  dev.off()

  cat("targetImg written to ", outputPath, "\n")



}


# ______________________________________________________________________________










