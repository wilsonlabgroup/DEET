testthat::test_that('Example, should run', {

  DEET_out <- DEET_enrich(example_DEET_enrich_input, DEET_dataset = DEET_example_data)

testthat::expect_true(all(names(DEET_out) == c("AP_INPUT_BP_output", "AP_INPUT_TF_output", "AP_DEET_DE_output", "AP_DEET_BP_output", "AP_DEET_TF_output", "DE_correlations")))
})

testthat::test_that('Unordered dataset, should run without correlations', {

  DEET_out <- DEET_enrich(example_DEET_enrich_input$gene_symbol, DEET_dataset = DEET_example_data)

  isTrue <- (DEET_out$DE_correlations == "No variance in coefs. Cannot proceed with correlation.")

  testthat::expect_true(isTrue)
})


testthat::test_that('Incorrect dataset file (just a gmt) -- throw error', {


  DEET_out <- testthat::expect_error(DEET_enrich(example_DEET_enrich_input$gene_symbol, DEET_dataset = DEET_example_data$DEET_gmt_DE))


})
