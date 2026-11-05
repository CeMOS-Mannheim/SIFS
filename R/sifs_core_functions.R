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
#'   \item{roi_info}{List with ROI metadata}
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
#' @param spmat_obj Sparse intensity matrix object from moleculaR, or a matrix.
#' @param labels Vector of class labels (0/1).
#' @param roi_name Character string, name of the ROI being analyzed.
#' @param train_ratio Ratio of training data (default 0.8).
#' @param xgb_params List of xgboost parameters.
#' @param nrounds Number of boosting rounds (default 200).
#' @param positive_class Value of `labels` that represents the ROI. The
#'   default is `1`. Internally the model is always trained with the ROI as 1.
#' @param seed Integer seed used for the stratified split (default 123).
#'
#' @return A list containing:
#'   \item{shap_df}{Data frame with m/z and mean SHAP for the ROI.}
#'   \item{model}{The trained xgboost model.}
#'   \item{train_idx}{Training sample indices.}
#'   \item{mz_values}{The m/z values corresponding to SHAP columns.}
#'   \item{shap_matrix}{SHAP values for the training pixels.}
#'   \item{positive_class}{The label value treated as the ROI.}
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
                               nrounds = 200L,
                               positive_class = 1L,
                               seed = 123L) {

  if (is.null(xgb_params)) {
    xgb_params <- list(
      objective = "binary:logistic",
      max_depth = 6,
      eta = 0.1,
      gamma = 0,
      colsample_bytree = 0.6,
      min_child_weight = 3,
      subsample = 1,
      eval_metric = "logloss",
      nthread = 2,
      verbosity = 0
    )
  }

  if (inherits(spmat_obj, "sparseIntensityMatrix")) {
    X <- spmat_obj$spmat
    mz_values <- spmat_obj$mzAxis
  } else {
    X <- spmat_obj
    mz_values <- if (!is.null(colnames(X))) {
      suppressWarnings(as.numeric(colnames(X)))
    } else {
      seq_len(ncol(X))
    }
    if (anyNA(mz_values)) {
      mz_values <- seq_len(ncol(X))
    }
  }

  if (!inherits(X, "dgCMatrix") && !is.matrix(X)) {
    stop("`spmat_obj` must be a sparseIntensityMatrix or a matrix.",
         call. = FALSE)
  }

  X <- as.matrix(X)
  mz_values <- as.numeric(mz_values)

  raw_labels <- as.character(labels)

  if (length(raw_labels) != nrow(X)) {
    stop(
      "`labels` length (", length(raw_labels), ") must equal the number of ",
      "pixels (", nrow(X), ").",
      call. = FALSE
    )
  }

  if (anyNA(raw_labels)) {
    stop("`labels` contains missing values.", call. = FALSE)
  }

  y_roi <- as.integer(raw_labels == as.character(positive_class))
  y <- y_roi

  if (length(unique(y)) < 2L) {
    stop("Both ROI and non-ROI classes are required.", call. = FALSE)
  }

  if (!is.numeric(train_ratio) || length(train_ratio) != 1L ||
      train_ratio <= 0 || train_ratio >= 1) {
    stop("`train_ratio` must be a single number between 0 and 1.",
         call. = FALSE)
  }

  # Stratified split without requiring the caret package.
  set.seed(seed)
  train_idx <- sort(unlist(lapply(split(seq_along(y), y), function(idx) {
    sample(idx, size = max(1L, floor(train_ratio * length(idx))))
  }), use.names = FALSE))

  if (length(unique(y[train_idx])) < 2L) {
    stop("The training split contains only one class. Increase `train_ratio`.",
         call. = FALSE)
  }

  train_data <- X[train_idx, , drop = FALSE]
  train_labels <- y[train_idx]

  colnames(train_data) <- as.character(mz_values)

  dtrain <- xgboost::xgb.DMatrix(
    data = train_data,
    label = train_labels,
    missing = NA
  )

  model <- xgboost::xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  # Compute SHAP on the training pixels. Feature names are preserved so that
  # downstream HMCS ranking can map SHAP columns back to m/z values.
  xai_train_df <- as.data.frame(train_data, check.names = FALSE)

  unified <- tryCatch(
    treeshap::unify(model, xai_train_df),
    error = function(e) {
      stop(
        "Could not unify the xgboost model for treeshap: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  shap_result <- treeshap::treeshap(
    unified,
    xai_train_df,
    verbose = FALSE
  )

  shap_matrix <- as.matrix(shap_result$shaps)

  if (ncol(shap_matrix) != length(mz_values)) {
    stop(
      "SHAP output has ", ncol(shap_matrix), " features but the m/z axis has ",
      length(mz_values), ".",
      call. = FALSE
    )
  }

  colnames(shap_matrix) <- as.character(mz_values)

  target_rows <- train_labels == 1L
  target_shaps <- shap_matrix[target_rows, , drop = FALSE]

  if (nrow(target_shaps) == 0L) {
    stop("No ROI pixels were present in the training split.", call. = FALSE)
  }

  mean_target_shap <- colMeans(target_shaps)

  shap_df <- data.frame(
    mzList = mz_values,
    shap = as.numeric(mean_target_shap),
    check.names = FALSE
  )
  names(shap_df)[2] <- paste0(roi_name, "_SHAP")

  shap_df <- shap_df[order(-abs(shap_df[[2]])), , drop = FALSE]
  rownames(shap_df) <- NULL

  list(
    shap_df = shap_df,
    model = model,
    train_idx = train_idx,
    mz_values = mz_values,
    shap_matrix = shap_matrix,
    positive_class = positive_class
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
#' @param ... Currently ignored.
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

  if (is.null(shap_df) || !"mzList" %in% names(shap_df)) {
    stop("`shap_results$shap_df` must contain an `mzList` column.",
         call. = FALSE)
  }

  # A feature is relevant if any of its SHAP contributions are non-zero.
  # Use absolute values so that large negative and positive contributions are
  # treated equally.
  if (ncol(shap_df) >= 2L) {
    mat <- as.matrix(shap_df[, -1, drop = FALSE])
    keep_idx <- rowSums(abs(mat), na.rm = TRUE) > 0

    if (anyNA(keep_idx)) {
      keep_idx[is.na(keep_idx)] <- FALSE
    }

    filtered_df <- shap_df[keep_idx, , drop = FALSE]
    removed_pct <- round(
      100 * (1 - nrow(filtered_df) / nrow(shap_df)),
      1
    )
  } else {
    filtered_df <- shap_df
    removed_pct <- 0
  }

  rownames(filtered_df) <- NULL

  # Optional SHAP weighing. The main vignette uses the unweighted HMCS path.
  weighed_values <- NULL
  if (isTRUE(weigh_shap)) {
    shap_col <- paste0(focus_roi, "_SHAP")

    if (!shap_col %in% names(filtered_df)) {
      stop("Column `", shap_col, "` was not found in the SHAP data frame.",
           call. = FALSE)
    }

    other_cols <- setdiff(names(filtered_df), c("mzList", shap_col))

    if (length(other_cols) > 0L) {
      mat <- as.matrix(filtered_df[, c(shap_col, other_cols), drop = FALSE])
      if (weight_method == "pos") {
        weighted <- filtered_df[[shap_col]] -
          rowSums(mat[, -1, drop = FALSE])
      } else if (weight_method == "neg") {
        weighted <- -filtered_df[[shap_col]] +
          rowSums(mat[, -1, drop = FALSE])
      } else {
        weighted <- abs(filtered_df[[shap_col]]) -
          rowSums(abs(mat[, -1, drop = FALSE]))
      }
      weighed_values <- weighted
    } else {
      weighed_values <- filtered_df[[shap_col]]
    }
  }

  list(
    filtered_df = filtered_df,
    removed_pct = removed_pct,
    weighed_values = weighed_values
  )
}

#' Calculate HMCS (Hotspot Margin Cancer Score)
#'
#' Calculate the HMCS values from DSC (Dice Similarity Coefficient) data.
#' This is the binary mode calculation for viable tumor vs necrosis.
#'
#' @param dsc_df Data frame with `mzList` and `dsc_mpm_*` columns.
#' @param focus_roi Character string, the ROI to focus on (for example
#'   `"VT"`).
#' @param reference_roi Character string, the comparison ROI. Defaults to
#'   `"Nec"` for backward compatibility. Use `"nonVT"` for a binary VT versus
#'   non-VT contrast without reading necrosis labels.
#' @param roi_names Optional character vector of ROI names. If supplied,
#'   `reference_roi` must be one of them.
#'
#' @return A named numeric vector of HMCS values, one per m/z value.
#'
#' @export
#'
#' @importFrom stats setNames
sifs_calculate_hmcs <- function(dsc_df,
                                focus_roi = "VT",
                                reference_roi = "Nec",
                                roi_names = NULL) {

  if (!"mzList" %in% names(dsc_df)) {
    stop("`dsc_df` must contain an `mzList` column.", call. = FALSE)
  }

  focus_col <- paste0("dsc_mpm_", focus_roi)
  ref_col <- paste0("dsc_mpm_", reference_roi)

  if (!is.null(roi_names)) {
    if (!focus_roi %in% roi_names) {
      stop("focus_roi must be one of: ",
           paste(roi_names, collapse = ", "), call. = FALSE)
    }
    if (!reference_roi %in% roi_names) {
      stop("reference_roi must be one of: ",
           paste(roi_names, collapse = ", "), call. = FALSE)
    }
  }

  missing_cols <- setdiff(c(focus_col, ref_col), names(dsc_df))

  if (length(missing_cols) > 0L) {
    stop(
      "DSC dataframe is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  hmcs_values <- dsc_df[[focus_col]] - dsc_df[[ref_col]]
  setNames(hmcs_values, dsc_df$mzList)
}

#' Visualize SIFS results
#'
#' Create plots for SIFS analysis results including SHAP rankings,
#' HMCS values, and dependence plots.
#'
#' @param hmcs_results Named numeric vector from `sifs_calculate_hmcs()` or
#'   a data frame with `mzList` and an HMCS column.
#' @param shap_df Data frame from `sifs_compute_shap()` or
#'   `sifs_filter_features()`.
#' @param output_dir Character string, directory to save plots (optional).
#' @param top_n Integer, number of top features to display (default 10).
#' @param ... Additional plotting parameters (currently unused).
#'
#' @return A list of ggplot objects. The output directory, if used, is stored
#'   in the `output_dir` attribute.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_col coord_flip theme_minimal ggsave aes labs
#' @importFrom data.table data.table
#' @importFrom utils globalVariables head
sifs_visualize <- function(hmcs_results, shap_df, output_dir = NULL,
                            top_n = 10, ...) {

  plots <- list()

  if (is.data.frame(hmcs_results)) {
    if (!"mzList" %in% names(hmcs_results)) {
      stop("HMCS data frame must contain `mzList`.", call. = FALSE)
    }
    hmcs_col <- setdiff(names(hmcs_results), "mzList")[[1L]]
    hmcs_named <- stats::setNames(
      as.numeric(hmcs_results[[hmcs_col]]),
      as.character(hmcs_results$mzList)
    )
  } else {
    hmcs_named <- hmcs_results
  }

  if (is.null(names(hmcs_named))) {
    stop("`hmcs_results` must be a named vector or a data frame.",
         call. = FALSE)
  }

  # 1. Top HMCS values
  hmcs_sorted <- sort(hmcs_named, decreasing = TRUE, na.last = TRUE)
  hmcs_top <- head(hmcs_sorted, top_n)

  plots$hmcs_bar <- ggplot2::ggplot(
    data.frame(
      mzList = names(hmcs_top),
      HMCS = as.numeric(hmcs_top)
    ),
    ggplot2::aes(x = stats::reorder(mzList, HMCS), y = HMCS)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top HMCS Values",
      x = "m/z",
      y = "HMCS"
    ) +
    ggplot2::theme_minimal()

  # 2. SHAP mean absolute values
  if (is.data.frame(shap_df) &&
      "mzList" %in% names(shap_df) &&
      ncol(shap_df) >= 2L) {
    shap_abs <- shap_df
    shap_mat <- as.matrix(shap_abs[, setdiff(names(shap_abs), "mzList"),
                                   drop = FALSE])
    shap_abs$mean_abs <- rowSums(abs(shap_mat), na.rm = TRUE)

    shap_top <- head(
      shap_abs[order(-shap_abs$mean_abs), , drop = FALSE],
      top_n
    )

    plots$shap_bar <- ggplot2::ggplot(
      shap_top,
      ggplot2::aes(x = stats::reorder(mzList, mean_abs), y = mean_abs)
    ) +
      ggplot2::geom_col(fill = "darkgreen") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Mean Absolute SHAP Values",
        x = "m/z",
        y = "Mean |SHAP|"
      ) +
      ggplot2::theme_minimal()
  } else {
    plots$shap_bar <- NULL
  }

  # 3. Save plots if output_dir specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    for (name in names(plots)) {
      if (!is.null(plots[[name]])) {
        ggplot2::ggsave(
          filename = file.path(output_dir, paste0(name, ".png")),
          plot = plots[[name]],
          dpi = 300,
          width = 10,
          height = 6
        )
      }
    }
  }

  attr(plots, "output_dir") <- output_dir
  plots
}

#' Run complete SIFS workflow
#'
#' A convenience function that runs the full SIFS pipeline from data loading
#' through to visualization in a sequential manner.
#'
#' @param mis_path Path to MIS RDS file.
#' @param spmat_path Path to SPMAT RDS file.
#' @param roi_names Character vector of ROI names.
#' @param focus_roi Which ROI to focus the analysis on.
#' @param dsc_df Data frame with `mzList` and `dsc_mpm_*` columns. This must
#'   be supplied by the caller because DSC depends on the spatial model
#'   workflow; the function no longer invents placeholder values.
#' @param reference_roi Character string, the comparison ROI passed to
#'   `sifs_calculate_hmcs()`. Defaults to `"Nec"`; use `"nonVT"` for the
#'   binary VT versus non-VT contrast.
#' @param positive_class Value of the MIS label column that represents the
#'   focus ROI. Defaults to `1`.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing all intermediate and final results.
#'
#' @export
#'
#' @importFrom stats setNames
sifs_run_workflow <- function(mis_path, spmat_path, roi_names = NULL,
                               focus_roi = "VT", dsc_df = NULL,
                               reference_roi = "Nec",
                               positive_class = 1L, ...) {

  if (is.null(dsc_df)) {
    stop(
      "`dsc_df` is required. Compute DSC values with the moleculaR spatial ",
      "workflow (SPP/MPM/DSC) and pass the resulting data frame.",
      call. = FALSE
    )
  }

  cat("=== Step 1: Loading data ===\n")
  data <- sifs_load_data(mis_path, spmat_path, roi_names)

  cat("=== Step 2: Computing SHAP values ===\n")
  labels <- data$mis$Vt_vs_all_lbls
  shap_results <- sifs_compute_shap(
    spmat_obj = data$spmat,
    labels = labels,
    roi_name = focus_roi,
    positive_class = positive_class,
    ...
  )

  cat("=== Step 3: Filtering features ===\n")
  filtered <- sifs_filter_features(shap_results, focus_roi = focus_roi, ...)

  cat("=== Step 4: Calculating HMCS ===\n")
  hmcs <- sifs_calculate_hmcs(
    dsc_df = dsc_df,
    focus_roi = focus_roi,
    reference_roi = reference_roi
  )

  cat("=== Step 5: Visualizing results ===\n")
  output_dir <- file.path("SIFS_results", focus_roi)
  viz <- sifs_visualize(
    hmcs,
    filtered$filtered_df,
    output_dir = output_dir
  )

  cat("=== SIFS Workflow Complete ===\n")
  cat("Results saved to:", attr(viz, "output_dir"), "\n")

  list(
    data = data,
    shap_results = shap_results,
    filtered = filtered,
    hmcs = hmcs,
    visualization = viz
  )
}

# Suppress R CMD check notes for ggplot2 non-standard evaluation symbols.
utils::globalVariables(c("mzList", "HMCS", "mean_abs"))