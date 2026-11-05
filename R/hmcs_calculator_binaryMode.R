#' HMCS calculator for binary datasets
#'
#' Calculates HMCS values for a binary comparison, for example viable tumor
#' (VT) versus a reference region.
#'
#' @param DSC_df A data frame containing `mzList` and the columns
#'   `dsc_mpm_<focusROI>` and `dsc_mpm_<referenceROI>`.
#' @param focusROI Character string specifying the tissue feature to focus on,
#'   for example `"VT"`.
#' @param referenceROI Character string specifying the comparison region.
#'   Defaults to `"Nec"` for backward compatibility. Use `"nonVT"` for a
#'   binary VT versus non-VT contrast that does not require necrosis labels.
#' @param ... Ignored.
#'
#' @return A numeric vector of HMCS values.
#'
#' @export
#'
#' @author Shad Arif Mohammed, \email{s.mohammed@doktoranden.hs-mannheim.de}
HMCS_calculatorBinary <- function(DSC_df,
                                  focusROI = "VT",
                                  referenceROI = "Nec",
                                  ...) {
  focus_col <- paste0("dsc_mpm_", focusROI)
  ref_col <- paste0("dsc_mpm_", referenceROI)

  if (!is.data.frame(DSC_df)) {
    stop("`DSC_df` must be a data frame.", call. = FALSE)
  }

  if (!all(c("mzList", focus_col, ref_col) %in% names(DSC_df))) {
    missing_cols <- setdiff(c("mzList", focus_col, ref_col), names(DSC_df))
    stop(
      "`DSC_df` is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  DSC_df[[focus_col]] - DSC_df[[ref_col]]
}
