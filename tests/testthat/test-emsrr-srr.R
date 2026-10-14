test_that("fwd evaluates a custom SRR using its environmental parameter", {
  library(FLCore)
  
  data("ple4", package = "FLCore")
  
  projection_year <- dims(ple4)$maxyear + 1L
  stock <- fwdWindow(ple4, end = projection_year)
  
  make_srr <- function(environment_value) {
    sr <- as.FLSR(stock)
    
    # A deliberately simple custom SRR:
    # recruitment = SSB * environmental value / 100
    model(sr) <- function() {
      list(model = rec ~ ssb * params / 100)
    }
    
    # Put the environmental value in the SRR params slot for the projection year.
    params(sr) <- FLPar(array(
      environment_value,
      dim = c(1L, 1L, 1L),
      dimnames = list(
        params = "environment",
        year = as.character(projection_year),
        iter = "1"
      )
    ))
    
    sr
  }
  
  control <- fwdControl(data.frame(
    year = projection_year,
    quant = "fbar",
    value = 0.2
  ))
  
  low_environment <- fwd(
    stock,
    sr = make_srr(10),
    control = control
  )
  
  high_environment <- fwd(
    stock,
    sr = make_srr(20),
    control = control
  )
  
  low_recruitment <- as.numeric(
    stock.n(low_environment)[1, ac(projection_year)]
  )
  high_recruitment <- as.numeric(
    stock.n(high_environment)[1, ac(projection_year)]
  )
  
  expect_true(is.finite(low_recruitment))
  expect_true(is.finite(high_recruitment))
  expect_gt(high_recruitment, low_recruitment)
})