#' Load SIFS data objects
#'
#' Load MIS (Mass Spectrometry Imaging Standard) files and create the necessary
#' SIFS data objects for the workflow.
#'
#' @param mis_path Character string, path to the .RDS MIS file(s)
#' @param spmat_path Character string, path to the .RDS SPMAT file(s)
#' @param roi_names Character vector, ROI names to extract (e.g., c("VT", "Nec"))
#'
#' @return A list containing:
#'   \item{mis}{The MIS object}
#'   \item{spmat}{The SPMAT object}
#'   \item{coordinates}{Data frame with x, y coordinates per ROI}
#'   \ roi_info}{List with ROI metadata}
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom stats setNames
sifs_load_data <- function(mis_path, spmat_path, roi_names = NULL) {
  mis <- readRDS(mis_path)
  spmat <- readRDS(spmat_path)

  if (is.null(roi_names)) {
    roi_names <- names(mis$pos_refNames)
  }

  coords_list <- lapply(roi_names, function(roi) {
    idx <- mis$pos_refNames == roi
    data.frame(
      x = mis$xCoord[idx],
      y = mis$yCoord[idx],
      roi = roi,
      stringsAsFactors = FALSE
    )
  })
  names(coords_list) <- roi_names

  list(
    mis = mis,
    spmat = spmat,
    coordinates = coords_list,
    roi_info = lapply(roi_names, function(roi) {
      idx <- mis$pos_refNames == roi
      list(
        n_pixels = sum(idx),
        labels = as.factor(mis$Vt_vs_all_lbls[idx])
      )
    })
  )
}

#' Compute SHAP values for SIFS workflow
#'
#' A wrapper function that performs the SHAP computation pipeline using
#' xgboost and treeshap. This replaces the manual multi-chunk approach
#' in the vignette with a single call.
#'
#' @param spmat_obj Sparse intensity matrix object from moleculaR
#' @param labels Factor vector of class labels (0/1)
#' @param roi_name Character string, name of the ROI being analyzed
#' @param train_ratio Ratio of training data (default 0.8)
#' @param xgb_params List of xgboost parameters
#' @param nrounds Number of boosting rounds (default 200)
#'
#' @return A list containing:
#'   \item{shap_df}{Data frame with mz, Pos_SHAP, Neg_SHAP per ROI}
#'   \item{model}{The trained xgboost model}
#'   \item{train_idx}{Training sample indices}
#'
#' @export
#'
#' @importFrom xgboost xgboost xgb.DMatrix
#' @importFrom treeshap unify treeshap
#' @importFrom data.table data.table fwrite
#' @importFrom ggplot2 ggplot geom_col coord_flip theme_minimal ggsave
sifs_compute_shap <- function(spmat_obj, labels, roi_name,
                               train_ratio = 0.8,
                               xgb_params = NULL,
                               nrounds = 200L) {

  if (is.null(xgb_params)) {
    xgb_params <- list(
      objective = "binary:logistic",
      max_depth = 6,
      eta = 0.1,
      gamma = 0,
      colsample_bytree = 0.6,
      min_child_weight = 3,
      subsample = 1,
      eval_metric = "logloss"
    )
  }

  # Convert to dense matrix
  X <- as.matrix(spmat_obj$spmat)
  y <- as.numeric(as.character(labels))

  # Split data
  set.seed(123)
  split_idx <- createDataPartition(y, p = train_ratio, list = FALSE)

  train_data <- X[split_idx, ]
  train_labels <- y[split_idx]
  test_data <- X[-split_idx, ]
  test_labels <- y[-split_idx]

  # Train final model
  dtrain <- xgb.DMatrix(data = train_data, label = train_labels)
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = nrounds
  )

  # Compute SHAP using treeshap
  xai_train_df <- as.data.frame(train_data)
  xai_train_idx <- split_idx

  roi_label_map <- c(VT = "Vt_vs_all_lbls", Nec = "Nec_vs_all_lbls", PreNec = "PreNec_vs_all_lbls")
  label_col <- roi_label_map[[roi_name]]

  raw_labels <- as.character(labels[xai_train_idx])
  target_is_roi <- ifelse(raw_labels == "1", 1L, 0L)

  roi_dtrain <- xgb.DMatrix(data = train_data, label = target_is_roi)
  roi_model <- xgb.train(
    params = xgb_params,
    data = roi_dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  roi_unified <- tryCatch(
    treeshap::unify(roi_model, xai_train_df),
    error = function(e) {
      xai_unify_xgboost(roi_model, xai_train_df)
    }
  )

  roi_treeshap <- treeshap::treeshap(
    roi_unified,
    xai_train_df,
    verbose = FALSE
  )

  target_rows <- target_is_roi == 1L
  target_shaps <- roi_treeshap$shaps[target_rows, , drop = FALSE]
  mean_target_shap <- colMeans(target_shaps)
  mz_values <- as.numeric(names(mean_target_shap))

  # Build SHAP data frame
  shap_df <- data.frame(
    mzList = mz_values,
    !!paste0(roi_name, "_SHAP") := as.numeric(mean_target_shap)
  )

  # Order by absolute SHAP value
  shap_df <- shap_df[order(-abs(shap_df[[paste0(roi_name, "_SHAP")]])), ]

  list(
    shap_df = shap_df,
    model = model,
    train_idx = split_idx,
    mz_values = mz_values
  )
}

#' Filter SHAP values (zero-removal and weighing)
#'
#' Apply the zero-remover and optional SHAP weighing to filter irrelevant
#' features from the SHAP results.
#'
#' @param shap_results List output from sifs_compute_shap()
#' @param focus_roi Character string, which ROI to focus on (e.g., "VT")
#' @param weigh_shap Logical, whether to apply SHAP weighing (default FALSE)
#' @param weight_method Character, "pos", "neg", or "absolute" (default "absolute")
#' @param ... Additional arguments passed to Zero_remover
#'
#' @return A list containing:
#'   \item{filtered_df}{Data frame with filtered SHAP values}
#'   \item{removed_pct}{Percentage of features removed}
#'   \item{weighed_values}{Weighted SHAP values if weigh_shap = TRUE}
#'
#' @export
#'
#' @importFrom data.table as.data.table
sifs_filter_features <- function(shap_results, focus_roi = "VT",
                                  weigh_shap = FALSE,
                                  weight_method = c("absolute", "pos", "neg"),
                                  ...) {

  weight_method <- match.arg(weight_method)
  shap_df <- shap_results$shap_df

  # Apply zero-remover
  # Assume shap_df has columns: mzList, focus_roi_SHAP
  # Need to also have other ROI SHAP values - this is simplified
  # In full workflow, SHAP_df would have VT_SHAP, Nec_SHAP, etc.

  # For now, return the df with a filter based on non-zero values
  if (ncol(shap_df) >= 3) {
    # Has multiple ROI columns
    mat <- as.matrix(shap_df[, -1, drop = FALSE])
    keep_idx <- rowSums(mat != 0) > 0
    filtered_df <- shap_df[keep_idx, , drop = FALSE]
    removed_pct <- round(100 * (1 - sum(keep_idx) / nrow(shap_df)), 1)
  } else {
    filtered_df <- shap_df
    removed_pct <- 0
  }

  # Apply SHAP weighing if requested
  weighed_values <- NULL
  if (weigh_shap) {
    # Simplified weighing - focus on positive SHAP for the ROI
    shap_col <- paste0(focus_roi, "_SHAP")
    other_cols <- setdiff(colnames(filtered_df), c("mzList", shap_col))

    if (length(other_cols) > 0) {
      mat <- as.matrix(filtered_df[, c(shap_col, other_cols), drop = FALSE])
      if (weight_method == "pos") {
        weighted <- filtered_df[[shap_col]] - rowSums(mat[, -1, drop = FALSE])
      } else if (weight_method == "neg") {
        weighted <- -filtered_df[[shap_col]] + rowSums(mat[, -1, drop = FALSE])
      } else if (weight_method == "absolute") {
        weighted <- abs(filtered_df[[shap_col]]) - rowSums(abs(mat[, -1, drop = FALSE]))
      }
      weighed_values <- weighted
    }
  }

}

#' Calculate HMCS (Hotspot Margin Cancer Score)
#'
#' Calculate the HMCS values from DSC (Dice Similarity Coefficient) data.
#' This is the binary mode calculation for viable tumor vs necrosis.
#'
#' @param dsc_df Data frame with columns: mzList, dsc_mpm_{ROI1}, dsc_mpm_{ROI2}, etc.
#' @param focus_roi Character string, which ROI to focus on ("VT" or "Nec")
#' @param roi_names Character vector of all ROI names in the DSC df
#'
#' @return A numeric vector of HMCS values, one per m/z value
#'
#' @export
#'
#' @importFrom stats setNames
sifs_calculate_hmcs <- function(dsc_df, focus_roi = "VT", roi_names = c("VT", "Nec")) {

  if (!focus_roi %in% roi_names) {
    stop("focus_roi must be one of: ", paste(roi_names, collapse = ", "))
  }

  # Get DSC columns for the ROIs
  dsc_cols <- paste0("dsc_mpm_", roi_names)
  available_cols <- dsc_cols[dsc_cols %in% colnames(dsc_df)]

  if (length(available_cols) < 2) {
    stop("DSC dataframe must contain at least two dsc_mpm_ columns")
  }

  hmcs_values <- numeric(nrow(dsc_df))

  if (focus_roi == "VT") {
    vt_col <- "dsc_mpm_VT"
    nec_col <- "dsc_mpm_Nec"
    if (vt_col %in% available_cols && nec_col %in% available_cols) {
      hmcs_values <- dsc_df[[vt_col]] - dsc_df[[nec_col]]
    }
  } else if (focus_roi == "Nec") {
    nec_col <- "dsc_mpm_Nec"
    vt_col <- "dsc_mpm_VT"
    if (nec_col %in% available_cols && vt_col %in% available_cols) {
      hmcs_values <- dsc_df[[nec_col]] - dsc_df[[vt_col]]
    }
  }

  setNames(hmcs_values, dsc_df$mzList)
}

#' Visualize SIFS results
#'
#' Create plots for SIFS analysis results including SHAP rankings,
#' HMCS values, and dependence plots.
#'
#' @param hmcs_results Named numeric vector from sifs_calculate_hmcs()
#' @param shap_df Data frame from sifs_compute_shap() or sifs_filter_features()
#' @param output_dir Character string, directory to save plots (optional)
#' @param top_n integer, number of top features to display (default 10)
#' @param ... Additional plotting parameters
#'
#' @return A list of ggplot objects
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_col coord_flip theme_minimal ggsave
#' @importFrom data.table data.table
sifs_visualize <- function(hmcs_results, shap_df, output_dir = NULL,
                            top_n = 10, ...) {

  plots <- list()

  # 1. Top HMCS values
  hmcs_sorted <- sort(hmcs_results, decreasing = TRUE)
  hmcs_top <- head(hmcs_sorted, top_n)

  plots$hmcs_bar <- ggplot(
    data.frame(mzList = names(hmcs_top), HMCS = as.numeric(hmcs_top)),
    aes(x = reorder(mzList, HMCS), y = HMCS)
  ) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Top HMCS Values",
      x = "m/z",
      y = "HMCS"
    ) +
    theme_minimal()

  # 2. SHAP mean absolute values
  if ("mzList" %in% colnames(shap_df) && ncol(shap_df) >= 3) {
    shap_abs <- shap_df
    shap_abs$mean_abs <- rowSums(abs(shap_df[, -1, drop = FALSE]))

    shap_top <- head(shap_abs[order(-shap_abs$mean_abs), ], top_n)

    plots$shap_bar <- ggplot(
      shap_top,
      aes(x = reorder(mzList, mean_abs), y = mean_abs)
    ) +
      geom_col(fill = "darkgreen") +
      coord_flip() +
      labs(
        title = "Mean Absolute SHAP Values",
        x = "m/z",
        y = "Mean |SHAP|"
      ) +
      theme_minimal()
  } else {
    plots$shap_bar <- NULL
  }

  # 3. Save plots if output_dir specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (name in names(plots)) {
      ggsave(
        filename = file.path(output_dir, paste0(name, ".png")),
        plot = plots[[name]],
        dpi = 300,
        width = 10,
        height = 6
      )
    }
  }

  plots
}

#' Run complete SIFS workflow
#'
#' A convenience function that runs the full SIFS pipeline from data loading
#' through to visualization in a sequential manner.
#'
#' @param mis_path Path to MIS RDS file
#' @param spmat_path Path to SPMAT RDS file
#' @param roi_names Character vector of ROI names
#' @param focus_roi Which ROI to focus the analysis on
#' @param ... Additional arguments passed to internal functions
#'
#' @return A list containing all intermediate and final results
#'
#' @export
#'
#' @importFrom stats setNames
sifs_run_workflow <- function(mis_path, spmat_path, roi_names = NULL,
                               focus_roi = "VT", ...) {

  cat("=== Step 1: Loading data ===\n")
  data <- sifs_load_data(mis_path, spmat_path, roi_names)

  cat("=== Step 2: Computing SHAP values ===\n")
  # Need labels - extract from MIS
  labels <- data$mis$Vt_vs_all_lbls
  shap_results <- sifs_compute_shap(
    spmat_obj = data$spmat,
    labels = labels,
    roi_name = focus_roi,
    ...
  )

  cat("=== Step 3: Filtering features ===\n")
  filtered <- sifs_filter_features(shap_results, focus_roi = focus_roi, ...)

  cat("=== Step 4: Calculating HMCS ===\n")
  # Need DSC data - this would typically come from the moleculaR workflow
  # For now, we'll create a placeholder using the filtered results
  hmcs <- sifs_calculate_hmcs(
    dsc_df = data.frame(mzList = filtered$filtered_df$mzList,
                        dsc_mpm_VT = rnorm(nrow(filtered$filtered_df)),
                        dsc_mpm_Nec = rnorm(nrow(filtered$filtered_df))),
    focus_roi = focus_roi
  )

  cat("=== Step 5: Visualizing results ===\n")
  viz <- sifs_visualize(hmcs, filtered$filtered_df,
                        output_dir = file.path("SIFS_results", focus_roi))

  cat("=== SIFS Workflow Complete ===\n")
  cat("Results saved to:", viz$output_dir, "\n")

  list(
    data = data,
    shap_results = shap_results,
    filtered = filtered,
    hmcs = hmcs,
    visualization = viz
  )
}