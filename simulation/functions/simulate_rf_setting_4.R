#' Simulation and code from Papers: Laabs, B.-H., Westenberger, A., and König, I. R. (2024). Identification of representative trees in random forests based on a new tree-based distance measure. Advances in Data Analysis and Classification, 18, 363–380.
#' Simulate data set for Szenario 4:
#' - Few influential Variables
#' - Large main effects
#' - Additional interaction effects
#' @param data        Number of observations
#' @param n_test      Number of obseravtions in test data set
#' @param n_val       Number of observations in validation data set
#' @param p           Number of covariables in total
#' @param p_eff       Number of binary influential variables
#' @parma p_int       Number of pair-wise interactions
#' @param beta_eff    Effect coefficients for influential variables
#' @param beta_int    Effect coefficient for pair-wise interactions
#' @param eps         standard deviation of noise term
#' @param num.trees   Number of trees for ranger
#' @param mtry        Mtry for ranger

simulate_rf_setting_4 <- function(data, n_test, n_val, p, p_eff, p_int, beta_eff, beta_int, eps, num.trees, mtry, min_node_size, 
                                  ...){
  
  n <- data
  
  params <- data.frame(min_node_size = min_node_size,
                       setting     = "Setting 4"
  )
  
  ## Check if p large enough
  if (p < p_eff + 2*p_int){
    stop("p is too small")
  }
  
  ## Simulate effect variables
  dat_bin <- lapply(1:(p_eff), function(x){
    as.data.frame(t(sample(c(0,1), size = n + n_test + n_val, replace = TRUE)))
  })
  dat_bin <- as.data.frame(t(data.table::rbindlist(dat_bin)))
  
  ## Simulate interaction variables
  dat_int <- lapply(1:(2*p_int), function(x){
    as.data.frame(t(sample(c(0,1), size = n + n_test + n_val, replace = TRUE)))
  })
  dat_int <- as.data.frame(t(data.table::rbindlist(dat_int)))
  
  ## Simulate continous dependent variable
  ## Calculate main effects
  y <- rowSums(dat_bin * beta_eff) 
  
  ## Calculate interaction effects
  for (i in 1:p_int){
    y <- y + dat_int[,2*(i-1) + 1] * dat_int[,2*(i-1) + 2] * beta_int
  }
  
  ## Add random noise
  y <- y + rnorm((n + n_test + n_val), mean = 0, sd = eps)
  
  
  # transform y to Values 0-1
  y <- (y-min(y))/(max(y)-min(y))
  
  ## Fill up with random noise variables
  dat_noise <- lapply(1:(p - (p_eff + 2*p_int)), function(x){
    as.data.frame(t(sample(c(0,1), size = n + n_test + n_val, replace = TRUE)))
  })
  dat_noise <- as.data.frame(t(data.table::rbindlist(dat_noise)))
  
  ## Merge data sets
  dat_all <- cbind(y, dat_bin, dat_int, dat_noise)
  names(dat_all) <- c("y", paste("V", 1:p, sep = ""))
  
  ## Split in training, testing and validation
  train_dat <- dat_all[1:n,]
  test_dat  <- dat_all[(n+1):(n + n_test),]
  val_dat   <- dat_all[(n + n_test + 1):nrow(dat_all),]
  
  ## Train ranger object
  rf <- ranger(y ~ ., data = train_dat, mtry = mtry, num.trees = num.trees, min.node.size = min_node_size, importance = "permutation")
  
  ## Set effect_var_ids
  effect_var_ids <- names(dat_all)[2:(p_eff+2*p_int+1)]
  noise_var_ids <- names(dat_all)[(p_eff+2*p_int+2):p]
  
  # prob RF
  train_data_0.5 <- train_dat %>% 
    mutate(y = ifelse(y>=0.5,1,0))
  rf_prob0.5 <- ranger(y~., data = train_data_0.5, min.node.size = min_node_size, probability = TRUE)
  
  
  
  ## Return 
  return(list(train_dat, test_dat, val_dat, rf, params, rf_prob0.5, effect_var_ids, noise_var_ids, "y"))
}