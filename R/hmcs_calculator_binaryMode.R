# [Date]   02.07.2025
# [Author] SHAD A. M.
# [Method] HMCS alculator binary mode.

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
