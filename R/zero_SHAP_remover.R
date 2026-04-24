#' Zero_SHAP_remover
#'
#' a function to dynamically remove values from a dataframe based on a criterion 
#' when a given number of columns have zero values in that particular list of columns. 
#'It removes majority of irrelevant features in the m/z axis  
#'
#' @param df  a dataframe containing SHAP values for each tissue histology feature, this should contain the following vectors: 
#' 1) "mzValues" as ranked by SHAP, 2) Average SHAP value for VT, 3) SHAP value for Necrosis, ... etc     
#' @param col_idx  a range of columns, specifying which tissue histology features to consider, for example 2:3 for VT and Nec.
#' @param ... ignored.
#'
#'
#' @return
#' clean dataframe that filters SHAP values DataFrames.
#'
#' @export
#'
#' @author Shad Arif Mohammed, \email{s.mohammed@doktoranden.hs-mannheim.de}
#'


Zero_remover <- function(df, col_idx, ...) {
      UseDataTable <- is.data.table(df)
      if (!UseDataTable) {
        mat <- as.matrix(df[, col_idx, drop = FALSE])
        df[rowSums(mat) != 0, , drop = FALSE]
      } else {
        mat <- as.matrix(df[, ..col_idx])
        df[rowSums(mat) != 0]
      }
}
