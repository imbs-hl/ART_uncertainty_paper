# Function to build three decision trees (DTs) as input for parameter "algorithm" using batchtools
# Prediction of continuous outcome: glycohemoglobin
# Prediction of probability outcome: prediabetes and diabetes
# One DT for each outcome
get_seperate_dt <- function(data, instance, min.bucket = 0, ...){
  
  # Extract components from instance
  #----------------
  train_data     <- instance[[1]]  # training data
  test_data      <- instance[[2]]  # test data
  cal_data       <- instance[[3]]  # calibration data
  rf             <- instance[[4]]  # regression random forest
  information_df <- instance[[5]]  # meta information
  
  # RF mean squared error on test data
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$glycohemoglobin - pred_regr_rf)^2)
  
  # Build DT for regression
  #----------------
  start <- proc.time()  # start runtime measurement
  
  # Train regression decision tree
  art_regr <- ranger(
    glycohemoglobin ~ .,
    data = train_data,
    min.bucket = min.bucket,
    mtry = (ncol(train_data) - 1),
    sample.fraction = 1,
    replace = FALSE,
    num.trees = 1
  )
  
  # Predictions from regression DT
  pred_regr_art <- predict(art_regr, test_data)$predictions
  
  # Prepare calibration and test values
  #----------------
  y_cal <- cal_data$glycohemoglobin
  y_test <- test_data$glycohemoglobin
  y_cal_pred <- predict(art_regr, cal_data)$predictions
  y_test_pred <- predict(art_regr, test_data)$predictions
  
  # Train probability DT for prediabetes (>= 5.7)
  train_data_5.7 <- train_data %>%
    mutate(glycohemoglobin = ifelse(glycohemoglobin >= 5.7, 1, 0))
  
  art_prob5.7 <- ranger(
    glycohemoglobin ~ .,
    data = train_data_5.7,
    min.bucket = min.bucket,
    mtry = (ncol(train_data) - 1),
    sample.fraction = 1,
    replace = FALSE,
    num.trees = 1,
    probability = TRUE
  )
  
  # Train probability DT for diabetes (>= 6.5)
  train_data_6.5 <- train_data %>%
    mutate(glycohemoglobin = ifelse(glycohemoglobin >= 6.5, 1, 0))
  
  art_prob6.5 <- ranger(
    glycohemoglobin ~ .,
    data = train_data_6.5,
    min.bucket = min.bucket,
    mtry = (ncol(train_data) - 1),
    sample.fraction = 1,
    replace = FALSE,
    num.trees = 1,
    probability = TRUE
  )
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])  # total runtime
  
  # Performance measures
  #----------------
  # Regression DT MSE
  mse_regr_art <- mean((test_data$glycohemoglobin - pred_regr_art)^2)
  
  # Average number of splits across trees
  num_splits_regr_art <- mean(c(
    nrow(treeInfo(art_regr) %>% filter(!terminal)),
    nrow(treeInfo(art_prob5.7) %>% filter(!terminal)),
    nrow(treeInfo(art_prob6.5) %>% filter(!terminal))
  ))
  
  # Brier scores per test observation
  #----------------
  # Assign observations to terminal nodes
  test_data_leafs  <- predict(art_regr, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(art_regr, train_data, type = "terminalNodes")$predictions
  cal_data_leafs   <- predict(art_regr, cal_data, type = "terminalNodes")$predictions
  
  brier_score_df <- data.frame(
    y = y_test,
    nodeID = test_data_leafs,
    prob5.7 = predict(art_prob5.7, test_data)$predictions[, 2],
    prob6.5 = predict(art_prob6.5, test_data)$predictions[, 2]
  ) %>%
    mutate(
      over5.7 = ifelse(y >= 5.7, 1, 0),
      over6.5 = ifelse(y >= 6.5, 1, 0)
    )
  
  brier_score5.7 <- BrierScore(brier_score_df$over5.7, brier_score_df$prob5.7)
  brier_score6.5 <- BrierScore(brier_score_df$over6.5, brier_score_df$prob6.5)
  
  # Maximum tree depth
  #--------
  