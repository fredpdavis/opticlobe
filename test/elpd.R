source("test/environment.R")

# Set up dummy data for InferState()
set.seed(2022)
dummy_data <- data.frame(
  value = c(rnorm(501, mean = 7, sd = 5), rnorm(501, mean = 10, sd = 5)),
  class = c(rep("A", 501), rep("B", 501)),
  driver = "driver"
)

bim.elpd_i <- numeric(length = 1002)

kpart1 <- inferState.partitionK(
  n = nrow(dummy_data),
  k = 10,
  labels = dummy_data$driver
)

# Only run one iteration as a proof of principle
for (k in seq(1)) {
  hList <- which(kpart1 == k)
  tList <- which(kpart1 != k)
  
  # This is modified from line 2401 - 2416 of src/R/analyzeOpticLobeExpr.R
  stan_bimodal <- list(
    file        = "test/level1_ordered_CV.stan",
    # From line 2439 - 2442
    data        = list(
      nSamples_t = length(which(kpart1 != k)),
      nSamples_h = length(which(kpart1 == k)),
      logE_t     = dummy_data$value[which(kpart1 != k)],
      logE_h     = dummy_data$value[which(kpart1 == k)]
    ),
    iter        = 100, # Just a quick test
    chains      = 4,
    cores       = 4,
    verbose     = FALSE,
    control     = list(adapt_delta = 0.98)
  )
  bimod_fit <- do.call(stan, stan_bimodal)
  
  ## Extract elpd: From line 2525 - 2536
  bim.kll <- extract_log_lik(bimod_fit,
                             parameter_name="log_lik_h")
  
  for (h in 1:length(hList)) {
    bim.elpd_i[hList[h]] <- logSumExp(bim.kll[h,]) -
      log(nrow(bim.kll)) ;
  }
}

## Structure of log likelihood table
dim(bim.kll)
# [1] 200 100
# Row: 2 * number of iterations
# Column: Size of a fold that is left-out in each iteration

## Original calculation of elpd
sum(bim.elpd_i)
# [1] -361.8318

## Expected value calculated by implementation in the LOO package
loo::elpd(bim.kll)
# Estimate   SE
# elpd   -311.0  9.3
# ic      622.0 18.6

bim.elpd_i_pr <- numeric(length = 1002)
# Updated calculation method in PR (take col sum)
for (h in 1:length(hList)) {
  bim.elpd_i_pr[hList[h]] <- logSumExp(bim.kll[, h]) -
    log(nrow(bim.kll)) ;
}
sum(bim.elpd_i_pr)
# [1] -310.9898