# Simple custom SRR for testing the FLasherEMSRR callback.
# Recruitment is constant across SSB and set by params[1].
test_custom_srr <- function(stock) {
  rate <- mean(rec(stock), na.rm = TRUE)
  
  predictModel(
    model = rec ~ params[1],
    params = FLPar(recruitment = rate)
  )
}