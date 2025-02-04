#' Simulate data set for Szenario 2:
#' - Many influential Variables
#' - Small main effects
#' @param data        Number of observations
#' @param n_test      Number of obseravtions in test data set
#' @param n_cal       Number of observations in calibration data set
#' @param p           Number of covariables in total
#' @param p_eff       Number of binary influential variables
#' @param beta_eff    Effect coefficients for influential variables
#' @param eps         standard deviation of noise term
#' @param num.trees   Number of trees for ranger
#' @param mtry        Mtry for ranger

simulate_rf_setting_2 <- function(data, n_test, n_cal, p, p_eff, beta_eff, eps, num.trees, mtry, min_node_size, ...){
  
  n <- data
  
  params <- data.frame(min_node_size = min_node_size,
                       setting     = "Setting 2"
  )
  
  ## Simulate effect variables
  dat_bin <- lapply(1:(p_eff), function(x){
    as.data.frame(t(sample(c(0,1), size = n + n_test + n_cal, replace = TRUE)))
  })
  dat_bin <- as.data.frame(t(data.table::rbindlist(dat_bin)))
  
  ## Simulate continous dependent variable
  ## Calculate main effects
  y <- rowSums(dat_bin * beta_eff) + rnorm((n + n_test + n_cal), mean = 0, sd = eps)
  
  # transform y to Values 0-1
  y <- (y-min(y))/(max(y)-min(y))
  
  ## Fill up with random noise variables
  dat_noise <- lapply(1:(p - p_eff), function(x){
    as.data.frame(t(sample(c(0,1), size = n + n_test + n_cal, replace = TRUE)))
  })
  dat_noise <- as.data.frame(t(data.table::rbindlist(dat_noise)))
  
  ## Merge data sets
  dat_all <- cbind(y, dat_bin, dat_noise)
  names(dat_all) <- c("y", paste("V", 1:p, sep = ""))
  
  ## Split in training, testing and calibration
  train_dat <- dat_all[1:n,]
  test_dat  <- dat_all[(n+1):(n + n_test),]
  cal_dat   <- dat_all[(n + n_test + 1):nrow(dat_all),]
  
  ## Train ranger object
  rf <- ranger(y ~ ., data = train_dat, mtry = mtry, num.trees = num.trees, min.node.size = min_node_size, importance = "permutation")
  
  ## Set effect_var_ids
  effect_var_ids <- names(dat_all)[2:(p_eff+1)]
  noise_var_ids <- names(dat_all)[(p_eff+2):p]
  
  ## Return 
  return(list(test_dat, rf, params, cal_dat, effect_var_ids, noise_var_ids, train_dat, "y"))
}