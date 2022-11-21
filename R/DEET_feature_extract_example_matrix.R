#' DEET_feature_extract_example_matrix
#'
#' An object of class data.frame where rows are genes and columns are comparisons.
#' The matrix is populated by the log2Fold-change of each gene within each study.
#' If the gene is not detected within that study, it is populated with 0 instead of the log2Fold-change.
#' This object is inputted into the `mat` input variable for the `DEET_feature_extract` function.
#' This example takes 1000 random genes and 200 random studies (seed = 1234s).
#'
#' @rdname DEET_feature_extract_example_matrix
#' @name DEET_feature_extract_example_matrix
#'
#' @usage data(DEET_feature_extract_example_matrix)
#'
#' @format  An object of class data.frame where rows are genes and columns are comparisons (1000 randomly selected genes and 200 randomly selected studies).
#'
#' @examples
#' data(DEET_feature_extract_example_matrix)
#' @keywords datasets
#'
NULL
