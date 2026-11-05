#' Compute binary HMCS values from spatial point patterns
#'
#' Run the SPP/MPM/DSC spatial workflow for a set of m/z features and combine
#' the resulting Dice similarity coefficients into a binary HMCS value.
#' With `nonvt` as the reference region, the calculation is:
#'
#' `HMCS = DSC(MPM hotspot, VT) - DSC(MPM hotspot, non-VT)`
#'
#' This is the Option A workflow: it uses only the VT labels and defines the
#' comparison region as the complement of VT, so no necrosis labels are read.
#'
#' @param spData A `moleculaR::sparseIntensityMatrix` object, as returned in
#'   `sifs_imzML_proc()$spData`.
#' @param mz Numeric vector of m/z values to test.
#' @param fwhm_obj An FWHM object returned by `readRDS()` of the FWHM file.
#' @param vt Logical vector marking VT pixels, ordered like `spData$spmat` rows.
#' @param nonvt Logical vector marking the reference (non-VT) pixels. If
#'   `NULL`, the complement of `vt` is used.
#' @param spwin Optional spatial window for the whole tissue. If `NULL`, it is
#'   created from `spData$coordinates`.
#' @param bw Bandwidth passed to `moleculaR::probMap()`. Default 2.5.
#' @param seed Integer seed passed to `moleculaR::probMap()`. Default 123.
#' @param min_points Minimum number of detected SPP points required before
#'   running `probMap()`. Features below this threshold return `NA`.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to `moleculaR::probMap()`.
#'
#' @return A data frame with one row per requested m/z and columns:
#'   `mzRequested`, `mzList`, `n_spp`, `hotspot_area`, `dsc_mpm_VT`,
#'   `dsc_mpm_nonVT`, and `HMCS_VT`.
#'
#' @export
sifs_compute_hmcs <- function(spData,
                              mz,
                              fwhm_obj,
                              vt,
                              nonvt = NULL,
                              spwin = NULL,
                              bw = 2.5,
                              seed = 123,
                              min_points = 100,
                              verbose = FALSE,
                              ...) {
  if (!inherits(spData, "sparseIntensityMatrix")) {
    stop(
      "`spData` must be a `sparseIntensityMatrix` from ",
      "`sifs_imzML_proc()$spData`.",
      call. = FALSE
    )
  }

  if (!inherits(fwhm_obj, "fwhm")) {
    stop("`fwhm_obj` must be an object of class `fwhm`.", call. = FALSE)
  }

  if (is.data.frame(mz)) {
    if (!"mzList" %in% names(mz)) {
      stop("`mz` data frame must contain an `mzList` column.", call. = FALSE)
    }
    mz <- mz$mzList
  }

  mz <- as.numeric(mz)
  mz <- unique(mz[is.finite(mz)])

  if (length(mz) == 0L) {
    stop("No finite m/z values were supplied.", call. = FALSE)
  }

  n_pixels <- nrow(spData$spmat)

  vt <- as.logical(vt)

  if (length(vt) != n_pixels || anyNA(vt)) {
    stop("`vt` must be a logical vector with one value per pixel.",
         call. = FALSE)
  }

  if (is.null(nonvt)) {
    nonvt <- !vt
  } else {
    nonvt <- as.logical(nonvt)
    if (length(nonvt) != n_pixels || anyNA(nonvt)) {
      stop("`nonvt` must be a logical vector with one value per pixel.",
           call. = FALSE)
    }
  }

  if (!any(vt)) {
    stop("No VT pixels are marked.", call. = FALSE)
  }

  if (!any(nonvt)) {
    stop("No non-VT reference pixels are marked.", call. = FALSE)
  }

  if (any(vt & nonvt)) {
    stop("`vt` and `nonvt` must not overlap.", call. = FALSE)
  }

  if (is.null(spwin)) {
    spwin <- moleculaR::createSpatialWindow(
      pixelCoords = spData$coordinates,
      plot = FALSE
    )
  }

  vt_win <- moleculaR::createSpatialWindow(
    pixelCoords = spData$coordinates[vt, , drop = FALSE],
    plot = FALSE
  )

  nonvt_win <- moleculaR::createSpatialWindow(
    pixelCoords = spData$coordinates[nonvt, , drop = FALSE],
    plot = FALSE
  )

  results <- vector("list", length(mz))

  for (i in seq_along(mz)) {
    mz_requested <- mz[[i]]
    mz_idx <- which.min(abs(spData$mzAxis - mz_requested))
    mz_matched <- spData$mzAxis[[mz_idx]]

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[%d/%d] requested m/z %.4f; matched %.4f",
          i, length(mz), mz_requested, mz_matched
        )
      )
    }

    result <- data.frame(
      mzRequested = mz_requested,
      mzList = mz_matched,
      n_spp = NA_integer_,
      hotspot_area = NA_real_,
      dsc_mpm_VT = NA_real_,
      dsc_mpm_nonVT = NA_real_,
      HMCS_VT = NA_real_,
      stringsAsFactors = FALSE
    )

    fwhm_mz <- moleculaR::getFwhm(fwhm_obj, mz_matched)

    if (!is.finite(fwhm_mz) || fwhm_mz <= 0) {
      results[[i]] <- result
      next
    }

    spp <- tryCatch(
      suppressWarnings(
        moleculaR::searchAnalyte(
          m = mz_matched,
          fwhm = fwhm_mz,
          spData = spData,
          spwin = spwin
        )
      ),
      error = function(e) NULL
    )

    if (is.null(spp)) {
      results[[i]] <- result
      next
    }

    result$n_spp <- spp$n

    if (spp$n < min_points) {
      results[[i]] <- result
      next
    }

    pm <- tryCatch(
      suppressWarnings(
        moleculaR::probMap(
          sppMoi = spp,
          bw = bw,
          seed = seed,
          verbose = FALSE,
          ...
        )
      ),
      error = function(e) NULL
    )

    if (is.null(pm)) {
      results[[i]] <- result
      next
    }

    dsc_vt <- tryCatch(
      suppressWarnings(moleculaR::dsc(pm$hotspotMask, vt_win)),
      error = function(e) NA_real_
    )

    dsc_nonvt <- tryCatch(
      suppressWarnings(moleculaR::dsc(pm$hotspotMask, nonvt_win)),
      error = function(e) NA_real_
    )

    result$hotspot_area <- tryCatch(
      summary(pm$hotspotMask)$area,
      error = function(e) NA_real_
    )
    result$dsc_mpm_VT <- dsc_vt
    result$dsc_mpm_nonVT <- dsc_nonvt
    result$HMCS_VT <- dsc_vt - dsc_nonvt

    results[[i]] <- result
  }

  do.call(rbind, results)
}
