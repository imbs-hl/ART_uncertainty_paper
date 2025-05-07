calculate_art_uncertainty <- function(data, instance, metric, imp.num.var, probs_quantiles, epsilon, min.bucket = 0, 
                                      num.splits = NULL, significance_level = c(0.1), ...){

  
  # Exctract data from instance 
  test_dat       <- instance[[1]]
  rf             <- instance[[2]]
  params         <- instance[[3]]
  cal_dat        <- instance[[4]]
  effect_var_ids <- instance[[5]]
  noise_var_ids  <- instance[[6]]
  train_dat      <- instance[[7]]
  dependent_varname <- instance[[8]]
  
  metric_test_data = NULL
  if(metric == "prediction"){
    # Divide calibration data set into 80% calibration, 20% extra test data for distance metric prediction 
    ids_cal <- sample(1:nrow(cal_dat), nrow(cal_dat)*0.8, replace = FALSE)
    
    metric_test_data <- cal_dat[-ids_cal,]
    cal_dat <- cal_dat[ids_cal,]
  }
  
  significance_level <- unlist(significance_level)
  probs_quantiles <- unlist(probs_quantiles)
  start <- proc.time()
  # Generate artificial representative tree (ART)
  rf_rep <- generate_tree(rf = rf, metric = metric, test_data = metric_test_data, train_data = train_dat, 
                          importance.mode = TRUE, imp.num.var = imp.num.var, dependent_varname = dependent_varname,
                          probs_quantiles = probs_quantiles, epsilon = epsilon,
                          min.bucket = min.bucket, num.splits = num.splits)
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # Calculate how many effect and noise variables are used as split variables in ART
  covered_effect_vars <- sum(effect_var_ids %in% treeInfo(rf_rep)$splitvarName)/length(effect_var_ids)
  covered_noise_vars  <- sum(noise_var_ids %in% treeInfo(rf_rep)$splitvarName)/length(noise_var_ids)
  effect_var_fdr      <- 1 - sum(treeInfo(rf_rep)$splitvarName %in% effect_var_ids)/(nrow(treeInfo(rf_rep)) - sum(is.na(treeInfo(rf_rep)$splitvarName)))
  
  # Prediction accuracy on test data set
  mse_test_dat_rf <- 1/nrow(test_dat) * sum((test_dat$y - predict(rf, data = test_dat[,-1])$predictions)^2)
  mse_test_dat_rep_tree <- 1/nrow(test_dat) * sum((test_dat$y - predict(rf_rep, data = test_dat[,-1])$predictions)^2)
  
  # Prediction accuracy on forest prediction on test data
  mse_rf_pred <- 1/nrow(test_dat) * sum((predict(rf, data = test_dat[,-1])$predictions - predict(rf_rep, data = test_dat[,-1])$predictions)^2)
  
  # Number of splits
  tree_info <- treeInfo(rf_rep)
  number_splits <- nrow(tree_info %>% filter(!terminal))
  
  # Calculate predictions for validation and test data
  y_cal <- cal_dat[,dependent_varname]
  y_test <- test_dat[,dependent_varname]
  y_cal_pred <- predict(rf_rep, cal_dat)$predictions
  y_test_pred <- predict(rf_rep, test_dat)$predictions
  
  # Calculate uncertainty for each significance level
  mean_error_icp_vec = c()
  mean_errors_mondrian_icp_vec = c()
  mean_error_cps_two_tailed_vec = c()
  mean_error_cps_left_tailed_vec = c()
  mean_error_cps_right_tailed_vec = c()
  mean_error_cps_two_tailed_mondrian_vec = c()
  mean_error_cps_left_tailed_mondrian_vec = c()
  mean_error_cps_right_tailed_mondrian_vec = c()
  mean_interval_size_icp_vec = c()
  mean_interval_size_mondrian_icp_vec = c()
  mean_interval_width_cps_two_tailed_vec = c()
  mean_interval_width_cps_two_tailed_mondrian_vec = c()
  for(alpha in significance_level){
    # Calibration using inductive conformal prediction (ICP)
    calibrated_predictions <- get_calibrated_prediction_regression(y_cal_pred = y_cal_pred, 
                                                                   y_cal = y_cal, 
                                                                   y_test_pred = y_test_pred, 
                                                                   significance_level = alpha)
    
    # Calibration per terminal node using Mondrian ICP
    # For Mondrian ICP the number of observations in the terminal nodes might be too small, in this case NA is returned
    calibrated_mondrian_predictions <- get_calibrated_prediction_regression_mondrian(y_cal_pred = y_cal_pred, 
                                                                                     y_cal = y_cal, 
                                                                                     y_test_pred = y_test_pred, 
                                                                                     significance_level = alpha,
                                                                                     tree = rf_rep,
                                                                                     cal_data = cal_dat, 
                                                                                     test_data = test_dat,
                                                                                     dependent_varname = dependent_varname,
                                                                                     show_node_id = TRUE
    )
    
    # Calibration using conformal predictive systems (CPS)
    # one- and two-tailed
    calibrates_cps_two_tailed <- get_calibrated_predictive_system(y_cal_pred, y_cal, y_test_pred, alpha,interval_type = "two-tailed", direction = NULL)
    calibrates_cps_left_tailed <- get_calibrated_predictive_system(y_cal_pred, y_cal, y_test_pred, alpha,interval_type = "one-tailed", direction = "left-tailed")
    calibrates_cps_right_tailed <- get_calibrated_predictive_system(y_cal_pred, y_cal, y_test_pred, alpha,interval_type = "one-tailed", direction = "right-tailed")
    
    # Calibration per terminal node using Mondrian CPS
    # For Mondrian CPS the number of observations in the terminal nodes might be too small, in this case NA is returned
    # one- and two-tailed
    calibrates_cps_two_tailed_mondrian <- NA
    calibrates_cps_left_tailed_mondrian <- NA
    calibrates_cps_right_tailed_mondrian <- NA
    
    try(calibrates_cps_two_tailed_mondrian <- get_calibrated_predictive_system_mondrian(y_cal_pred, y_cal, y_test_pred, tree = rf_rep, cal_data = cal_dat, test_data = test_dat, dependent_varname, significance_level = alpha,interval_type = "two-tailed", direction = NULL))
    try(calibrates_cps_left_tailed_mondrian <- get_calibrated_predictive_system_mondrian(y_cal_pred, y_cal, y_test_pred, tree = rf_rep, cal_data = cal_dat, test_data = test_dat, dependent_varname, significance_level = alpha,interval_type = "one-tailed", direction = "left-tailed"))
    try(calibrates_cps_right_tailed_mondrian <- get_calibrated_predictive_system_mondrian(y_cal_pred, y_cal, y_test_pred, tree = rf_rep, cal_data = cal_dat, test_data = test_dat, dependent_varname, significance_level = alpha,interval_type = "one-tailed", direction = "right-tailed"))
    
    # Error rates for 4 methods
    # ICP
    error_icp <- (calibrated_predictions$upper_bound < y_test | calibrated_predictions$lower_bound > y_test)
    mean_error_icp <- sum(error_icp)/length(y_test)*100
    
    # Mondrian ICP
    errors_mondrian_icp <- (calibrated_mondrian_predictions$upper_bound < y_test | calibrated_mondrian_predictions$lower_bound > y_test)
    mean_errors_mondrian_icp <- round(sum(errors_mondrian_icp)/length(y_test)*100, 4)
    
    # CPS
    mean_error_cps_two_tailed <- round(sum(calibrates_cps_two_tailed$lower_bound>y_test|y_test>calibrates_cps_two_tailed$upper_bound)/length(y_test)*100, 4)
    mean_error_cps_left_tailed <- round(sum(calibrates_cps_left_tailed$lower_bound>y_test)/length(y_test)*100, 4)
    mean_error_cps_right_tailed <- round(sum(y_test>calibrates_cps_right_tailed$upper_bound)/length(y_test)*100, 4)
    
    # Mondrian CPS
    mean_error_cps_two_tailed_mondrian <- NA
    mean_error_cps_left_tailed_mondrian <- NA
    mean_error_cps_right_tailed_mondrian <- NA
    try(mean_error_cps_two_tailed_mondrian <- round(sum(calibrates_cps_two_tailed_mondrian$lower_bound>y_test|y_test>calibrates_cps_two_tailed_mondrian$upper_bound)/length(y_test)*100, 4))
    try(mean_error_cps_left_tailed_mondrian <- round(sum(calibrates_cps_left_tailed_mondrian$lower_bound>y_test)/length(y_test)*100, 4))
    try(mean_error_cps_right_tailed_mondrian <- round(sum(y_test>calibrates_cps_right_tailed_mondrian$upper_bound)/length(y_test)*100, 4))
    
    
    # Interval width for 4 methods, for CPS and Mondrian CPS only two-tailed intervals are considered
    # ICP
    mean_interval_size_icp <- mean(calibrated_predictions$upper_bound-calibrated_predictions$lower_bound)
    
    # Mondrian ICP
    mean_interval_size_mondrian_icp <- mean(calibrated_mondrian_predictions$upper_bound-calibrated_mondrian_predictions$lower_bound)
    
    # CPS
    mean_interval_width_cps_two_tailed <- mean(calibrates_cps_two_tailed$upper_bound-calibrates_cps_two_tailed$lower_bound)
    
    # Mondrian CPS
    mean_interval_width_cps_two_tailed_mondrian <- NA
    try(mean_interval_width_cps_two_tailed_mondrian <- mean(calibrates_cps_two_tailed_mondrian$upper_bound-calibrates_cps_two_tailed_mondrian$lower_bound))
    
    # Add all values to list with results of all loops
    mean_error_icp_vec = c(mean_error_icp_vec, mean_error_icp)
    mean_errors_mondrian_icp_vec = c(mean_errors_mondrian_icp_vec,  mean_errors_mondrian_icp)
    mean_error_cps_two_tailed_vec = c(mean_error_cps_two_tailed_vec, mean_error_cps_two_tailed)
    mean_error_cps_left_tailed_vec = c(mean_error_cps_left_tailed_vec, mean_error_cps_left_tailed)
    mean_error_cps_right_tailed_vec = c(mean_error_cps_right_tailed_vec, mean_error_cps_right_tailed)
    mean_error_cps_two_tailed_mondrian_vec = c(mean_error_cps_two_tailed_mondrian_vec, mean_error_cps_two_tailed_mondrian)
    mean_error_cps_left_tailed_mondrian_vec = c(mean_error_cps_left_tailed_mondrian_vec, mean_error_cps_left_tailed_mondrian)
    mean_error_cps_right_tailed_mondrian_vec = c(mean_error_cps_right_tailed_mondrian_vec, mean_error_cps_right_tailed_mondrian)
    
    mean_interval_size_icp_vec = c(mean_interval_size_icp_vec, mean_interval_size_icp)
    mean_interval_size_mondrian_icp_vec = c(mean_interval_size_mondrian_icp_vec, mean_interval_size_mondrian_icp)
    mean_interval_width_cps_two_tailed_vec = c(mean_interval_width_cps_two_tailed_vec, mean_interval_width_cps_two_tailed)
    mean_interval_width_cps_two_tailed_mondrian_vec = c(mean_interval_width_cps_two_tailed_mondrian_vec, mean_interval_width_cps_two_tailed_mondrian)
    
  }
  
  return(results = data.frame(metric= metric, 
                    method               = "artificial tree with uncertainty",
                    setting              = params$setting,
                    min_node_size        = params$min_node_size,
                    covered_effect_vars  = covered_effect_vars,
                    covered_noise_vars   = covered_noise_vars,
                    effect_var_fdr       = effect_var_fdr,
                    mse_test_dat_rf      = mse_test_dat_rf,
                    mse_test_dat_tree    = mse_test_dat_rep_tree,
                    mse_rf_pred          = mse_rf_pred,
                    imp.num.var          = imp.num.var,
                    time                 = time,
                    probs_quantiles      = paste0(probs_quantiles, collapse = ","),
                    number_splits        = number_splits,
                    epsilon              = epsilon,
                    min.bucket           = min.bucket,
                    num.splits           = paste0(num.splits, collapse = ","),
                    significance_level   = significance_level,
                    mean_error_icp                  = mean_error_icp_vec,
                    mean_errors_mondrian_icp         = mean_errors_mondrian_icp_vec,
                    mean_error_cps_two_tailed = mean_error_cps_two_tailed_vec,
                    mean_error_cps_left_tailed = mean_error_cps_left_tailed_vec,
                    mean_error_cps_right_tailed = mean_error_cps_right_tailed_vec,
                    mean_error_cps_two_tailed_mondrian = mean_error_cps_two_tailed_mondrian_vec,
                    mean_error_cps_left_tailed_mondrian = mean_error_cps_left_tailed_mondrian_vec,
                    mean_error_cps_right_tailed_mondrian = mean_error_cps_right_tailed_mondrian_vec,
                    mean_interval_size_icp          = mean_interval_size_icp_vec,
                    mean_interval_size_mondrian_icp = mean_interval_size_mondrian_icp_vec,
                    mean_interval_width_cps_two_tailed = mean_interval_width_cps_two_tailed_vec,
                    mean_interval_width_cps_two_tailed_mondrian = mean_interval_width_cps_two_tailed_mondrian_vec,
                    n_test                                  = nrow(test_dat),
                    n_train                                 = nrow(train_dat),
                    n_cal                                   = nrow(cal_dat)
                    )
  )
}