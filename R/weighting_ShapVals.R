#' Conditional SHAP weighing function
#'
#' a function to conditionally weigh SHAP values across the different tissue histology feaures. 
#' WARNING: only use this method if you have domain knowledge and you delibrately want to attend SIFS to a specific tissue histology feature of interest. 
#' when a given number of columns have zero values in that particular list of columns. 
#' It removes majority of irrelevant features in the m/z axis  
#' This method weighs SHAP values per class, and thereby decreases False Positive Rates (FPR) inherent to Spectral SHAP. EXPERIMENTAL as of 21.03.2024 and requires extensive domain knowlege. 
#' NOTE: outputs of this method are not included in Mohammed et al., 2026 manuscript.
#' Weighing criterion can be based on positive, negative, average, absolute, absolute average SHAP values, in combination with the different tissue histology features of interest.
#' outputs vary based on user defined parameter selections. 
#' @param df  a dataframe containing SHAP values for each tissue histology feature:    
#' @param focusColumn  a character, specifying which tissue histology feature you want SIFS to focus on and weigh based on it.
#' @param condition a preset character setting, like "VT" for viable tumor, to set the condition on viable tumor and penalize localization inside other tissue histolgy features.
#' @param weighted a boolean value, TRUE = run SHAP weighing, FALSE = return to classical HMCS calculation. 
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


# Function to weigh SHAP values based on different conditions
SHAP_weighing <- function(df, focusColumn = "VtShap", condition = which_shap_focus, weighted = FALSE, ...) {

      # Initialize a list to store weighted values
      weightedValues <- list()
      mz_perClass    <- list()
      mz_with_wt     <- list(mz= list(), WeightOfClass = list(), class = "")

      if (weighted == TRUE){
            # Determine the condition for weighting
            if (condition == "pos") {
                  # Weighting based on positive SHAP values
                  if (focusColumn == "VtShap") {
                        # subset the dataframe to account for the sign of the dedicated class, such that proper assignment is performed:

                        # index of the positive signed shap values focued on the column as selected in the focusColumn and the sign of the selection.
                        idx_signed = which(df[["Vt_sign"]]== "Positive")
                        #df <- df[idx_signed, ]


                        # Update : 21.03.2024 excluding this weighing criterion.

                        weightedValues <- df[["VtShap"]] - (df[["NecShap"]] + df[["PreNecShap"]]) # + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else if (focusColumn == "VtVascShap") {
                        idx_signed = which(df[["VtVasc_sign"]]== "Positive")
                        #df <- df[idx_signed, ]

                        weightedValues <- df[["VtVascShap"]]  - (df[["VtShap"]] + df[["NecShap"]] + df[["PreNecShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "NecShap") {
                        idx_signed = which(df[["Nec_sign"]]== "Positive")
                        #df <- df[idx_signed, ]
                        weightedValues <- df[["NecShap"]] - (df[["VtShap"]] + df[["PreNecShap"]]) # + df[["VtVascShap"]]

                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        idx_signed = which(df[["PreNec_sign"]]== "Positive")
                        #df <- df[idx_signed, ]

                        weightedValues <- df[["PreNecShap"]]  - (df[["VtShap"]] + df[["NecShap"]]) # + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")
                  }
            } else if (condition == "neg") {
                  # Weighting based on negative SHAP values
                  if (focusColumn == "VtShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["Vt_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["VtShap"]] + ( df[["NecShap"]] + df[["PreNecShap"]]) #df[["VtVascShap"]] +
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "VtVascShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["VtVasc_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["VtVascShap"]] + (df[["VtShap"]] + df[["NecShap"]] + df[["PreNecShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "NecShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["Nec_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["NecShap"]] + (df[["VtShap"]] + df[["PreNecShap"]]) # + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        # subsetting criterion:-
                        idx_signed = which(df[["PreNec_sign"]]== "Negative")
                        #df <- df[idx_signed, ]
                        weightedValues <- -df[["PreNecShap"]] + (df[["VtShap"]]  + df[["NecShap"]]) #+ df[["VtVascShap"]]

                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  }
            } else {
                  # Weighting based on absolute SHAP values
                  if (focusColumn == "VtShap") {
                        weightedValues <- abs(df[["VtShap"]]) - ( abs(df[["NecShap"]]) + abs(df[["PreNecShap"]])) #abs(df[["VtVascShap"]]) +
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "VtVascShap") {
                        weightedValues <- abs(df[["VtVascShap"]]) - (abs(df[["VtShap"]]) + abs(df[["NecShap"]]) + abs(df[["PreNecShap"]]))
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")
                  } else if (focusColumn == "NecShap") {
                        weightedValues <- abs(df[["NecShap"]]) - (abs(df[["VtShap"]]) + abs(df[["PreNecShap"]])) #+ abs(df[["VtVascShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        weightedValues <- abs(df[["PreNecShap"]]) - (abs(df[["VtShap"]])+ abs(df[["NecShap"]])) # + abs(df[["VtVascShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  }
            }

            return(mz_with_wt)







      } else{
            # Determine the condition for weighting
            if (condition == "pos") {
                  # Weighting based on positive SHAP values
                  if (focusColumn == "VtShap") {
                        # subset the dataframe to account for the sign of the dedicated class, such that proper assignment is performed:

                        # index of the positive signed shap values focued on the column as selected in the focusColumn and the sign of the selection.
                        idx_signed = which(df[["Vt_sign"]]== "Positive")
                        #df <- df[idx_signed, ]


                        # Update : 21.03.2024 excluding this weighing criterion.

                        weightedValues <- df[["VtShap"]] # - (df[["NecShap"]] + df[["PreNecShap"]])  + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else if (focusColumn == "VtVascShap") {
                        idx_signed = which(df[["VtVasc_sign"]]== "Positive")
                        #df <- df[idx_signed, ]

                        weightedValues <- df[["VtVascShap"]] # - (df[["VtShap"]] + df[["NecShap"]] + df[["PreNecShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "NecShap") {
                        idx_signed = which(df[["Nec_sign"]]== "Positive")
                        #df <- df[idx_signed, ]
                        weightedValues <- df[["NecShap"]] #- (df[["VtShap"]] + df[["PreNecShap"]]) # + df[["VtVascShap"]]

                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        idx_signed = which(df[["PreNec_sign"]]== "Positive")
                        #df <- df[idx_signed, ]

                        weightedValues <- df[["PreNecShap"]] # - (df[["VtShap"]] + df[["NecShap"]]) # + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")
                  }
            } else if (condition == "neg") {
                  # Weighting based on negative SHAP values
                  if (focusColumn == "VtShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["Vt_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["VtShap"]] # + ( df[["NecShap"]] + df[["PreNecShap"]]) #df[["VtVascShap"]] +
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "VtVascShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["VtVasc_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["VtVascShap"]] #+ (df[["VtShap"]] + df[["NecShap"]] + df[["PreNecShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "NecShap") {
                        # subsetting criterion:-
                        idx_signed = which(df[["Nec_sign"]]== "Negative")
                        #df <- df[idx_signed, ]

                        weightedValues <- -df[["NecShap"]] #+ (df[["VtShap"]] + df[["PreNecShap"]]) # + df[["VtVascShap"]]
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        # subsetting criterion:-
                        idx_signed = which(df[["PreNec_sign"]]== "Negative")
                        #df <- df[idx_signed, ]
                        weightedValues <- -df[["PreNecShap"]] #+ (df[["VtShap"]]  + df[["NecShap"]]) #+ df[["VtVascShap"]]

                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  }
            } else {
                  # Weighting based on absolute SHAP values
                  if (focusColumn == "VtShap") {
                        weightedValues <- abs(df[["VtShap"]]) - ( abs(df[["NecShap"]]) + abs(df[["PreNecShap"]])) #abs(df[["VtVascShap"]]) +
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")


                  } else if (focusColumn == "VtVascShap") {
                        weightedValues <- abs(df[["VtVascShap"]]) - (abs(df[["VtShap"]]) + abs(df[["NecShap"]]) + abs(df[["PreNecShap"]]))
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")
                  } else if (focusColumn == "NecShap") {
                        weightedValues <- abs(df[["NecShap"]]) - (abs(df[["VtShap"]]) + abs(df[["PreNecShap"]])) #+ abs(df[["VtVascShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  } else {
                        weightedValues <- abs(df[["PreNecShap"]]) - (abs(df[["VtShap"]])+ abs(df[["NecShap"]])) # + abs(df[["VtVascShap"]])
                        mz_perClass    <- df[["mzVals"]]
                        mz_with_wt     <- list(mz= mz_perClass, WeightOfClass= weightedValues, class = "Vt")

                  }
            }

            return(mz_with_wt)



      }



# final closing backet:
}

