calculate_art_rep_tree <- function(data, instance, metric, imp.num.var, probs_quantiles, epsilon, min.bucket = 0, num.splits = NULL, significance_level = 0.1, ...){

  # Exctract data from instance 
  test_dat        <- instance[[1]]
  rf              <- instance[[2]]
  params          <- instance[[3]]
  train_dat       <- instance[[4]]
  dependent_varname <- instance[[5]]
  information_df  <- instance[[6]]
  cal_dat  <- instance[[7]]
  
  probs_quantiles <- unlist(probs_quantiles)
  start <- proc.time()
  # Generate artificial representative tree (ART)
  rf_rep <- generate_tree_reimplementation(rf = rf, metric = metric, test_data = test_dat, train_data = train_dat, 
                                           importance.mode = TRUE, imp.num.var = imp.num.var, dependent_varname = dependent_varname,
                                           probs_quantiles = probs_quantiles, epsilon = epsilon,
                                           min.bucket = min.bucket, num.splits = num.splits)
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])

  # Prediction accuracy on test data set
  mse_test_dat_rf <- 1/nrow(test_dat) * sum((test_dat$y - predict(rf, data = test_dat)$predictions)^2)
  mse_test_dat_rep_tree <- 1/nrow(test_dat) * sum((test_dat$y - predict(rf_rep, data = test_dat)$predictions)^2)
  
  # Prediction accuracy on forest prediction on test data
  mse_rf_pred <- 1/nrow(test_dat) * sum((predict(rf, data = test_dat)$predictions - predict(rf_rep, data = test_dat)$predictions)^2)
  
  # number of splits
  tree_info <- treeInfo(rf_rep)
  number_splits <- nrow(tree_info %>% filter(!terminal))
  # Calculate uncertainty
  # Calculate predictions for calibration and test data
  y_cal <- cal_dat[,dependent_varname]
  y_test <- test_dat[,dependent_varname]
  y_cal_pred <- predict(rf_rep, cal_dat)$predictions
  y_test_pred <- predict(rf_rep, test_dat)$predictions
  
  # Calibration using conformal prediction
  # Get calibrated prediction for test data
  calibrated_predictions <- get_calibrated_prediction_regression(y_cal_pred = y_cal_pred, 
                                                                 y_cal = y_cal, 
                                                                 y_test_pred = y_test_pred, 
                                                                 significance_level = significance_level)
  
  # Calculate error rate and mean intercal size
  error_calibration <- (calibrated_predictions$upper_bound < y_test | calibrated_predictions$lower_bound > y_test)
  mean_error_calibration <- sum(error_calibration)/length(y_test)*100
  
  # Mean interval size
  mean_interval_size_calibration <- mean(calibrated_predictions$upper_bound-calibrated_predictions$lower_bound)
  
  # Calibration per class using mondrian conformal prediction
  # Get calibrated prediction for test data
  calibrated_mondrian_predictions <- get_calibrated_prediction_regression_mondrian(y_cal_pred = y_cal_pred, 
                                                                                   y_cal = y_cal, 
                                                                                   y_test_pred = y_test_pred, 
                                                                                   significance_level = significance_level,
                                                                                   tree = rf_rep,
                                                                                   cal_data = cal_dat, 
                                                                                   test_data = test_dat,
                                                                                   dependent_varname = dependent_varname,
                                                                                   show_node_id = TRUE
  )
  
  # Calculate error rate and mean intercal size
  errors_mondrian <- (calibrated_mondrian_predictions$upper_bound < y_test | calibrated_mondrian_predictions$lower_bound > y_test)
  
  mean_error_mondrian_calibration <- round(sum(errors_mondrian)/length(y_test)*100, 4)
  
  # Mean interval size
  mean_interval_size_mondrian_calibration <- mean(calibrated_mondrian_predictions$upper_bound-calibrated_mondrian_predictions$lower_bound)
  
  return(data.frame(metric               = metric, 
                    method               = "artificial tree with uncertainty",
                    min_node_size        = params$min_node_size,
                    mse_test_dat_rf      = mse_test_dat_rf,
                    mse_test_dat_tree    = mse_test_dat_rep_tree,
                    mse_rf_pred          = mse_rf_pred,
                    imp.num.var          = imp.num.var,
                    runtime              = time,
                    probs_quantiles      = paste0(probs_quantiles, collapse = ","),
                    number_splits        = number_splits,
                    epsilon              = epsilon,
                    r.squared_ranger     = rf$r.squared,
                    min.bucket_art       = min.bucket,
                    num.splits           = paste0(num.splits, collapse = ","),
                    significance_level   = significance_level,
                    mean_error_calibration                  = mean_error_calibration,
                    mean_error_mondrian_calibration         = mean_error_mondrian_calibration,
                    mean_interval_size_mondrian_calibration = mean_interval_size_mondrian_calibration,
                    mean_interval_size_calibration          = mean_interval_size_calibration
                    ) %>% cbind(., information_df)
  )
}