#' HMCS calcualtor for binary datasets (Example: viable tumor vs all)
#'
#' This functions calculates HMCS values for binary comparisons.  
#' 
#'
#' @param DSC_df  a dataframe, containing the following vectors: mzList, dsc_mpm_VT. 
#' @param focusROI  a character, specifying which tissue histology feaure you want to focus on, for example: "VT".
#' @param regions a `regions` object, see `?getRegions`.
#' @param ... ignored.
#'
#'
#' @return
#' Calculates HMCS in binary mode, returns HMCS values from DSC.
#'
#' @export
#'
#' @author Shad Arif Mohammed, \email{s.mohammed@doktoranden.hs-mannheim.de}
#'


HMCS_calculatorBinary <- function(DSC_df, focusROI = "VT", ...){
      # equation
      if(focusROI == "VT"){
            hmcs_VT <- (DSC_df$dsc_mpm_VT)- (DSC_df$dsc_mpm_Nec )#+ DSC_df$dsc_PreNec
            return(hmcs_VT)
      } else if(focusROI == "Nec"){
            hmcs_Nec <- (DSC_df$dsc_mpm_Nec)- (DSC_df$dsc_mpm_VT)# + DSC_df$dsc_PreNec
            return(hmcs_Nec)
      }


}
