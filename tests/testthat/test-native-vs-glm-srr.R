library(FLCore)

test_that("fwd uses Beverton-Holt parameters for each iteration", {
  data("ple4", package = "FLCore")
  
  projection_year <- dims(ple4)$maxyear + 1L
  stock <- propagate(fwdWindow(ple4, end = projection_year), 3L)
  years <- dimnames(stock)$year
  
  # Parameter order is a, b. Keep every stock year in the array so that
  # conversion to FLBiolcpp cannot introduce a year-position mismatch.
  values <- array(
    NA_real_,
    dim = c(2L, length(years), 3L),
    dimnames = list(
      params = c("a", "b"),
      year = years,
      iter = as.character(1:3)
    )
  )
  values["a", , 1] <- 1e6
  values["b", , 1] <- 1e5
  values["a", , 2] <- 2e6  # Twice a, same b: twice recruitment.
  values["b", , 2] <- 1e5
  values["a", , 3] <- 1e6
  values["b", , 3] <- 2e5  # Same a, higher b: lower recruitment.
  
  sr <- as.FLSR(stock)
  model(sr) <- bevholt()
  params(sr) <- FLPar(values)
  
  control <- FLasherEMSRR::fwdControl(data.frame(
    year = projection_year, quant = "fbar", value = 0.2
  ))
  projected <- FLasherEMSRR::fwd(
    stock, sr = sr, control = control
  )
  
  recruitment <- as.numeric(
    stock.n(projected)[1, ac(projection_year)]
  )
  
  expect_length(recruitment, 3L)
  expect_true(all(is.finite(recruitment)))
  expect_true(all(recruitment > 0))
  expect_equal(recruitment[2] / recruitment[1], 2, tolerance = 1e-4)
  expect_lt(recruitment[3], recruitment[1])
})

test_that("fwd applies both GLM covariates by year and iteration", {
  data("ple4", package = "FLCore")
  
  # Deterministic synthetic observations: one positive and one negative
  # environmental effect. Fit an actual model rather than hard-coding
  # the coefficients in the prediction callback.
  training <- expand.grid(
    ssb = c(1e5, 5e5, 1e6, 2e6),
    warm = c(0, 1),
    stress = c(0, 1)
  )
  training$rec <- exp(
    log(1e6) + 1e-7 * training$ssb +
      0.3 * training$warm - 0.5 * training$stress
  )
  fit <- glm(
    rec ~ ssb + warm + stress,
    data = training,
    family = gaussian(link = "log")
  )
  expect_gt(unname(coef(fit)["warm"]), 0)
  expect_lt(unname(coef(fit)["stress"]), 0)
  
  # The C++ callback evaluates its expression in an environment that
  # inherits from .GlobalEnv. Temporarily install this closure there,
  # preserving and restoring any pre-existing binding.
  function_name <- "predict_test_glm_srr"
  previously_present <- exists(
    function_name, envir = .GlobalEnv, inherits = FALSE
  )
  if (previously_present) {
    previous_value <- get(
      function_name, envir = .GlobalEnv, inherits = FALSE
    )
  }
  on.exit({
    if (previously_present) {
      assign(function_name, previous_value, envir = .GlobalEnv)
    } else {
      rm(list = function_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  predict_test_glm_srr <- function(ssb, params) {
    stopifnot(length(params) == 2L)
    as.numeric(predict(
      fit,
      newdata = data.frame(
        ssb = as.numeric(ssb),
        warm = params[1],
        stress = params[2]
      ),
      type = "response"
    ))
  }
  assign(function_name, predict_test_glm_srr, envir = .GlobalEnv)
  
  projection_year <- dims(ple4)$maxyear + 1L
  stock <- propagate(fwdWindow(ple4, end = projection_year), 4L)
  years <- dimnames(stock)$year
  
  # Iterations: baseline; positive covariate; negative covariate; both.
  # All earlier years are zero; only the projection year differs.
  covariates <- array(
    0,
    dim = c(2L, length(years), 4L),
    dimnames = list(
      params = c("warm", "stress"),
      year = years,
      iter = as.character(1:4)
    )
  )
  covariates["warm", ac(projection_year), 2] <- 1
  covariates["stress", ac(projection_year), 3] <- 1
  covariates["warm", ac(projection_year), 4] <- 1
  covariates["stress", ac(projection_year), 4] <- 1
  
  sr <- as.FLSR(stock)
  model(sr) <- rec ~ predict_test_glm_srr(ssb, params)
  params(sr) <- FLPar(covariates)
  
  control <- FLasherEMSRR::fwdControl(data.frame(
    year = projection_year, quant = "fbar", value = 0.2
  ))
  projected <- FLasherEMSRR::fwd(
    stock, sr = sr, control = control
  )
  
  recruitment <- as.numeric(
    stock.n(projected)[1, ac(projection_year)]
  )
  
  expect_length(recruitment, 4L)
  expect_true(all(is.finite(recruitment)))
  expect_true(all(recruitment > 0))
  expect_gt(recruitment[2], recruitment[1])
  expect_lt(recruitment[3], recruitment[1])
  expect_equal(
    recruitment[2] / recruitment[1],
    exp(unname(coef(fit)["warm"])),
    tolerance = 1e-3
  )
  expect_equal(
    recruitment[3] / recruitment[1],
    exp(unname(coef(fit)["stress"])),
    tolerance = 1e-3
  )
  expect_equal(
    recruitment[4] / recruitment[1],
    exp(sum(coef(fit)[c("warm", "stress")])),
    tolerance = 1e-3
  )
})