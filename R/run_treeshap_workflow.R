# Install treeshap if not available
if (!require("treeshap", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("ModelOriented/treeshap")
}

library(treeshap)
library(shapviz)
library(xgboost)
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)

abs_path <- "D:/Projects/Rpackages/SIFS/"

cat("Loading SPMAT objects...\n")
t3_spmat_pos <- readRDS(file = file.path(abs_path, "data/input/01_spmatObjects/trainset_tissue/tr_spmat_XIII_t3_83318.rds"))
t4_spmat_pos <- readRDS(file = file.path(abs_path, "data/input/01_spmatObjects/testset_tissue/tst_spmat_XV_t1_82558.rds"))

lbls_pos_All <- readRDS(file = file.path(abs_path, "data/input/04_y_lbls/spectraLblsAll/04_Lbl_df_giant_pos/Lbls_df_Pos_all.rds"))

ROI <- "VT"
if (ROI == "VT") {
  t3_lbls <- as.factor(lbls_pos_All$Vt_vs_all_lbls[lbls_pos_All$pos_refNames == "t3"])
  t4_lbls <- as.factor(lbls_pos_All$Vt_vs_all_lbls[lbls_pos_All$pos_refNames == "t4"])
}

y_train <- as.numeric(t3_lbls)
y_test <- as.numeric(t4_lbls)

X_train <- as.matrix(t3_spmat_pos[["spmat"]])
X_test <- as.matrix(t4_spmat_pos[["spmat"]])

MZ <- round(t3_spmat_pos[["mzAxis"]], 6)

model_path <- file.path(abs_path, "data/output/models/02_xgb_final_VT.rds")
if (file.exists(model_path)) {
  cat("Loading existing XGBoost model...\n")
  xgb_final <- readRDS(model_path)
} else {
  cat("Training new XGBoost model...\n")
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test, label = y_test)

  xgb_final <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 4,
      eta = 0.05,
      nthread = 4,
      seed = 42
    ),
    data = dtrain,
    nrounds = 60,
    verbose = FALSE
  )
  saveRDS(xgb_final, model_path)
}

cat("Preparing treeshap analysis...\n")
set.seed(123)
train_data_df <- as.data.frame(X_train)
MZ_cols <- MZ[1:ncol(train_data_df)]
colnames(train_data_df) <- MZ_cols
num_features <- ncol(train_data_df)

cat("Creating model with feature names...\n")
feature_names <- MZ_cols
attr(xgb_final, "feature_names") <- feature_names

cat("Unifying model...\n")
unified_model <- unify(xgb_final, train_data_df)

cat("Calculating SHAP values...\n")
shap_values_treshap <- treeshap(unified_model, train_data_df, verbose = 0)

cat("SHAP dimensions:", nrow(shap_values_treshap$shaps), "x", ncol(shap_values_treshap$shaps), "\n")

cat("Analyzing mean SHAP values...\n")
mean_shap_values <- colMeans(shap_values_treshap$shaps)
mean_shap_df <- data.frame(
  feature = names(mean_shap_values),
  mean_shap = mean_shap_values
)

top_mean_shap_df <- mean_shap_df %>%
  arrange(desc(mean_shap)) %>%
  head(10)

top_negative_shap_df <- mean_shap_df %>%
  filter(mean_shap < 0) %>%
  arrange(mean_shap) %>%
  head(10)

output_dir_shap <- file.path(abs_path, "data/output/SHAP_outputs", paste0(ROI, "_treeshap"))
dir.create(output_dir_shap, recursive = TRUE, showWarnings = FALSE)

p_pos <- ggplot(top_mean_shap_df, aes(x = reorder(feature, mean_shap), y = mean_shap)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Positive Mean SHAP Values", x = "Features", y = "Mean SHAP Value") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
ggsave(file.path(output_dir_shap, "top_mean_shap_values.png"), plot = p_pos, width = 10, height = 6, dpi = 300)

p_neg <- ggplot(top_negative_shap_df, aes(x = reorder(feature, mean_shap), y = mean_shap)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  labs(title = "Top 10 Negative Mean SHAP Values", x = "Features", y = "Mean SHAP Value") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
ggsave(file.path(output_dir_shap, "top_negative_mean_shap_values.png"), plot = p_neg, width = 10, height = 6, dpi = 300)

fwrite(top_mean_shap_df, file = file.path(output_dir_shap, "01_top10_positive_shapVals.csv"))
fwrite(top_negative_shap_df, file = file.path(output_dir_shap, "02_top10_negative_shapVals.csv"))

cat("\n=== Output Summary ===\n")
cat("Output directory:", output_dir_shap, "\n")
cat("\nTop 10 Positive Mean SHAP Values:\n")
print(top_mean_shap_df)
cat("\nTop 10 Negative Mean SHAP Values:\n")
print(top_negative_shap_df)

cat("\ntreeshap workflow completed successfully!\n")