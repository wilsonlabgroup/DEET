testthat::test_that('Example, should run', {

  single1 <- DEET_feature_extract(DEET_feature_extract_example_matrix,
                                  DEET_feature_extract_example_response,"categorical")
  outRight <- all(names(single1) == c("elastic_net_coefficients", "elastic_net", "basic_features"))
  testthat::expect_true(outRight)

})

testthat::test_that('Response variable mismatch -- throw error', {

  single1 <- testthat::expect_error(DEET_feature_extract(DEET_feature_extract_example_matrix,
                                  DEET_feature_extract_example_response,"continuous"))

})



testthat::test_that('One gene, Error due to downstream analysis.', {
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  Signature <- odds_ratio_in

  genes <- rownames(Signature)[1]
  testthat::expect_error(tissue_scMappR_custom( genes, signature_matrix = Signature,
                                         output_directory =  "scMappR_test", toSave = FALSE))
})



