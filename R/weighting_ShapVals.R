
# Function to weigh SHAP values based on different conditions
SHAP_weighing <- function(df, focusColumn = "VtShap", condition = which_shap_focus, weighted = FALSE) {

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

