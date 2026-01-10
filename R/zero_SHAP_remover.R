# a function to dynamically remove values from a dataframe based on a criterion when a given number of columns have zero values in that particular list of columns:
Zero_remover <- function(df, col_idx) {
      mat <- as.matrix(df[, ..col_idx])
      df[rowSums(mat) != 0]
}
