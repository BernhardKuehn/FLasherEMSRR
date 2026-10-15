test_that("custom SRR preserves the SSB derivative", {
  result <- FLasherEMSRR:::customSRR_AD_diagnostic(
    ssb = 3,
    model_name = "ssb^2 + params[1]",
    params = 5
  )
  
  expect_equal(unname(result[["value"]]), 14)
  expect_equal(unname(result[["derivative"]]), 6, tolerance = 1e-8)
})