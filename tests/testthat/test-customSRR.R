test_that("custom SRR parameters override the FLCore params generic", {
  
  gampredict_test <- function(ssb, params) {
    stopifnot(is.numeric(params))
    
    FLQuant(
      1000 + 0 * as.numeric(ssb) + as.numeric(params)
    )
  }
  
  old <- get0(
    "gampredict_test",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  
  assign(
    "gampredict_test",
    gampredict_test,
    envir = .GlobalEnv
  )
  
  on.exit({
    if (is.null(old)) {
      rm("gampredict_test", envir = .GlobalEnv)
    } else {
      assign("gampredict_test", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  data(ple4)
  
  srr <- as.FLSR(ple4)
  
  model(srr) <- rec ~ gampredict_test(ssb, params)
  
  params(srr) <- FLPar(
    sst = FLQuant(
      15,
      dimnames = list(
        year = dimnames(ple4)$year,
        iter = dimnames(ple4)$iter
      )
    )
  )
  
  # Deliberate collision with the FLCore accessor.
  params <- FLCore::params
  
  ctrl <- FLasherEMSRR::fwdControl(
    year = 2010,
    quant = "fbar",
    value = 0.2
  )
  
  result <- FLasherEMSRR::fwd(
    ple4,
    sr = srr,
    control = ctrl
  )
  
  expect_s4_class(result, "FLStock")
})