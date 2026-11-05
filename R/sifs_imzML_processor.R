#' Load and process an imzML dataset for SIFS
#'
#' Read a processed imzML file, bin peaks using the supplied full width at
#' half maximum (FWHM) object, optionally remove low-frequency peaks, and
#' return a model-ready sparse intensity matrix together with aligned labels.
#'
#' The returned `spData` element is a `moleculaR::sparseIntensityMatrix` object
#' with `spmat`, `mzAxis`, and `coordinates` components. It is the object that
#' should be passed to `moleculaR::searchAnalyte()` and related functions.
#'
#' @param imzML_path Character string, path to the imzML file.
#' @param imzML_annotations_path Character string, path to the label table.
#' @param fwhm_path Character string, path to the FWHM object (`.rds`).
#' @param labelType Character string, one of `"rds"`, `"csv"`, or `"tsv"`.
#' @param sample_id Optional sample identifier in `pos_refNames`. If `NULL`,
#'   the sample whose coordinates match the imzML is selected automatically.
#' @param FWHM_setUp Numeric m/z value used as the reference for FWHM-based
#'   binning. Defaults to 400 (lipid-focused setting).
#' @param minFreq Numeric minimum peak frequency across spectra. Use `0` to
#'   skip frequency filtering. Defaults to 0.01.
#'
#' @return A list containing:
#'   \item{spmat}{The sparse intensity matrix (pixels x m/z).}
#'   \item{spData}{A `sparseIntensityMatrix` with `spmat`, `mzAxis`, and
#'   `coordinates`.}
#'   \item{coords}{Data frame with x, y coordinates for the imzML sample.}
#'   \item{lbls}{The selected sample labels, ordered to match `spmat` rows.}
#'   \item{spwin}{Spatial window object of the entire tissue sample.}
#'   \item{mzList}{The m/z axis.}
#'   \item{meta}{Processing metadata.}
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom stats setNames
#'
sifs_imzML_proc <- function(imzML_path,
                            imzML_annotations_path,
                            fwhm_path,
                            labelType = "rds",
                            sample_id = NULL,
                            FWHM_setUp = 400,
                            minFreq = 0.01) {
  imzML_path <- normalizePath(imzML_path, mustWork = TRUE)
  imzML_annotations_path <- normalizePath(
    imzML_annotations_path, mustWork = TRUE
  )
  fwhm_path <- normalizePath(fwhm_path, mustWork = TRUE)

  # Read imzML spectra
  imzML <- moleculaR::readCentrData(imzML_path)

  # FWHM-based binning. The tolerance is the FWHM at FWHM_setUp divided by
  # FWHM_setUp, i.e. a relative tolerance focused on the requested m/z region.
  fwhmObj <- readRDS(fwhm_path)
  bin_tol <- moleculaR::getFwhm(fwhmObj, FWHM_setUp) / FWHM_setUp

  if (!is.finite(bin_tol) || bin_tol <= 0) {
    stop(
      "Could not derive a valid FWHM-based binning tolerance from ",
      "`fwhm_path`. Check the FWHM object and `FWHM_setUp`.",
      call. = FALSE
    )
  }

  imzML <- MALDIquant::binPeaks(
    imzML,
    tolerance = bin_tol,
    method = "relaxed"
  )

  if (!is.null(minFreq) && isTRUE(minFreq > 0)) {
    imzML <- moleculaR::filterPeaks(imzML, minFreq = minFreq)
  }

  if (length(imzML) == 0L) {
    stop("No spectra remain after binning/filtering.", call. = FALSE)
  }

  # createSparseMat() is the single source of the sparseIntensityMatrix.
  spData <- moleculaR::createSparseMat(imzML)
  spmat  <- spData$spmat
  coords <- spData$coordinates
  mzList <- spData$mzAxis

  # Read labels
  lbls_raw <- switch(
    labelType,
    rds = readRDS(imzML_annotations_path),
    csv = data.table::fread(imzML_annotations_path),
    tsv = data.table::fread(imzML_annotations_path, sep = "\t"),
    stop("Unsupported labelType: ", labelType, call. = FALSE)
  )

  # Select the matching sample and align labels to matrix row order
  lbls <- .select_sifs_labels(lbls_raw, coords, sample_id = sample_id)

  if (nrow(lbls) != nrow(spmat)) {
    stop(
      "Label/imzML mismatch: ", nrow(lbls), " labels vs ",
      nrow(spmat), " spectra.",
      call. = FALSE
    )
  }

  spwin <- moleculaR::createSpatialWindow(
    pixelCoords = coords,
    plot = FALSE
  )

  list(
    spmat  = spmat,
    spData = spData,
    coords = coords,
    lbls   = lbls,
    spwin  = spwin,
    mzList = mzList,
    meta   = list(
      sample_id          = attr(lbls, "sifs_sample_id"),
      n_label_rows_total = attr(lbls, "sifs_label_rows_total"),
      n_spectra          = nrow(spmat),
      n_features         = length(mzList),
      fwhm_path          = fwhm_path,
      FWHM_setUp         = FWHM_setUp,
      fwhm_tolerance     = bin_tol,
      minFreq            = minFreq,
      imzML_path         = imzML_path,
      labels_path        = imzML_annotations_path,
      labelType          = labelType
    )
  )
}


# Select one sample from a label table and order it to match the imzML pixels.
#
# The label table may contain many samples (for example,
# `Lbls_df_Pos_all.rds`). When `sample_id` is supplied, that sample is used
# directly. Otherwise, each `pos_refNames` group is compared with the imzML
# coordinates and the unique group covering every coordinate is selected.
.select_sifs_labels <- function(lbls, coords, sample_id = NULL) {
  lbls <- as.data.frame(lbls, stringsAsFactors = FALSE)

  if (!all(c("xCoord", "yCoord") %in% names(lbls))) {
    stop("Labels must contain `xCoord` and `yCoord` columns.", call. = FALSE)
  }

  coord_key <- function(x, y) {
    paste(
      format(x, trim = TRUE, scientific = FALSE),
      format(y, trim = TRUE, scientific = FALSE),
      sep = "_"
    )
  }

  imz_key <- coord_key(coords$x, coords$y)

  if (anyDuplicated(imz_key)) {
    stop("The imzML coordinates contain duplicates.", call. = FALSE)
  }

  total_rows <- nrow(lbls)

  if (!is.null(sample_id)) {
    if (!"pos_refNames" %in% names(lbls)) {
      stop(
        "`sample_id` supplied, but labels have no `pos_refNames` column.",
        call. = FALSE
      )
    }

    lbls <- lbls[lbls$pos_refNames == sample_id, , drop = FALSE]

    if (nrow(lbls) == 0L) {
      stop(
        "No labels found for sample_id = '", sample_id, "'.",
        call. = FALSE
      )
    }
  } else if ("pos_refNames" %in% names(lbls)) {
    refs <- as.character(unique(lbls$pos_refNames))
    refs <- refs[!is.na(refs)]

    if (length(refs) == 0L) {
      stop("Labels have no usable `pos_refNames` values.", call. = FALSE)
    }

    lbl_key_all <- coord_key(lbls$xCoord, lbls$yCoord)

    counts <- vapply(refs, function(ref) {
      sum(imz_key %in% lbl_key_all[lbls$pos_refNames == ref])
    }, integer(1))

    full <- refs[counts == length(imz_key)]

    if (length(full) == 0L) {
      best <- refs[which.max(counts)]
      stop(
        "No label sample covers all imzML coordinates. Best match: ",
        best, " (", max(counts), "/", length(imz_key), "). ",
        "Supply `sample_id` explicitly.",
        call. = FALSE
      )
    }

    if (length(full) > 1L) {
      stop(
        "Multiple samples cover all imzML coordinates: ",
        paste(full, collapse = ", "),
        ". Supply `sample_id` explicitly.",
        call. = FALSE
      )
    }

    sample_id <- full[[1L]]
    lbls <- lbls[lbls$pos_refNames == sample_id, , drop = FALSE]
  }

  lbl_key <- coord_key(lbls$xCoord, lbls$yCoord)

  if (anyDuplicated(lbl_key)) {
    stop("Selected labels contain duplicate coordinates.", call. = FALSE)
  }

  idx <- match(imz_key, lbl_key)

  if (anyNA(idx)) {
    stop(
      "Selected labels miss ", sum(is.na(idx)),
      " imzML coordinates. Check `sample_id` and coordinate columns.",
      call. = FALSE
    )
  }

  lbls <- lbls[idx, , drop = FALSE]
  rownames(lbls) <- NULL

  attr(lbls, "sifs_sample_id") <- sample_id
  attr(lbls, "sifs_label_rows_total") <- total_rows

  lbls
}
