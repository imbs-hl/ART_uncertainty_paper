# Function to construct three decision trees (DTs)
# Builds one DT for continuous outcome prediction (glycohemoglobin)
# Builds two separate DTs for probability prediction of prediabetes and diabetes
# Uses pre-trained random forests as reference for regression performance
# Returns all DTs together with performance metrics and meta information

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
  regression_tree <- ranger(
    glycohemoglobin ~ .,
    data = train_data,
    min.bucket = min.bucket,
    mtry = (ncol(train_data) - 1),
    sample.fraction = 1,
    replace = FALSE,
    num.trees = 1
  )
  
  # Predictions from regression DT
  pred_regr_dt <- predict(regression_tree, test_data)$predictions
  
  # Prepare calibration and test values
  #----------------
  y_cal <- cal_data$glycohemoglobin
  y_test <- test_data$glycohemoglobin
  y_cal_pred <- predict(regression_tree, cal_data)$predictions
  y_test_pred <- predict(regression_tree, test_data)$predictions
  
  # Train probability DT for prediabetes (>= 5.7)
  train_data_5.7 <- train_data %>%
    mutate(glycohemoglobin = as.factor(ifelse(glycohemoglobin >= 5.7, 1, 0)))
  
  probability_tree5.7 <- ranger(
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
    mutate(glycohemoglobin = as.factor(ifelse(glycohemoglobin >= 6.5, 1, 0)))
  
  probability_tree6.5 <- ranger(
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
  mse_regr_dt <- mean((test_data$glycohemoglobin - pred_regr_dt)^2)
  
  # Average number of splits across trees
  num_splits_regr_dt <- mean(c(
    nrow(treeInfo(regression_tree) %>% filter(!terminal)),
    nrow(treeInfo(probability_tree5.7) %>% filter(!terminal)),
    nrow(treeInfo(probability_tree6.5) %>% filter(!terminal))
  ))
  
  # Brier scores per test observation
  #----------------
  # Assign observations to terminal nodes
  test_data_leafs  <- predict(regression_tree, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(regression_tree, train_data, type = "terminalNodes")$predictions
  cal_data_leafs   <- predict(regression_tree, cal_data, type = "terminalNodes")$predictions
  
  brier_score_df <- data.frame(
    y = y_test,
    nodeID = test_data_leafs,
    prob5.7 = predict(probability_tree5.7, test_data)$predictions[, 2],
    prob6.5 = predict(probability_tree6.5, test_data)$predictions[, 2]
  ) %>%
    mutate(
      over5.7 = ifelse(y >= 5.7, 1, 0),
      over6.5 = ifelse(y >= 6.5, 1, 0)
    )
  
  brier_score5.7 <- BrierScore(brier_score_df$over5.7, brier_score_df$prob5.7)
  brier_score6.5 <- BrierScore(brier_score_df$over6.5, brier_score_df$prob6.5)
  
  # Maximum tree depth
  #----------------
  get_max_depth <- function(rf_object, tree = 1) {
    
    # Extract tree structure
    ti <- ranger::treeInfo(rf_object, tree = tree)
    
    # Compute node depths iteratively
    depths <- setNames(integer(nrow(ti)), ti$nodeID)
    depths["0"] <- 0  # root node
    
    for (i in seq_len(nrow(ti))) {
      node <- ti$nodeID[i]
      
      # Skip terminal nodes
      if (ti$terminal[i]) next
      
      # Read left and right child nodes
      left  <- ti$leftChild[i]
      right <- ti$rightChild[i]
      
      depths[as.character(left)]  <- depths[as.character(node)] + 1
      depths[as.character(right)] <- depths[as.character(node)] + 1
    }
    
    # Return maximum depth
    max(depths)
  }
  
  max_depth <- mean(
    max(get_max_depth(regression_tree)),
    max(get_max_depth(probability_tree5.7)),
    max(get_max_depth(probability_tree6.5))
  )
  
  # Number of terminal nodes
  num_leaves_regr <- nrow(treeInfo(regression_tree) %>% filter(terminal))
  num_leaves_prediabetes <- nrow(treeInfo(probability_tree5.7) %>% filter(terminal))
  num_leaves_diabetes <- nrow(treeInfo(probability_tree6.5) %>% filter(terminal))
  num_leaves <- mean(num_leaves_regr, num_leaves_prediabetes, num_leaves_diabetes)
  
  # Sparsity: number of splitting variables
  num_vars_split_regr <- length(na.omit(unique(treeInfo(regression_tree)$splitvarID)))
  num_vars_split_prediabetes <- length(na.omit(unique(treeInfo(probability_tree5.7)$splitvarID)))
  num_vars_split_diabetes <- length(na.omit(unique(treeInfo(probability_tree6.5)$splitvarID)))
  num_vars_split <- mean(num_vars_split_regr, num_vars_split_prediabetes, num_vars_split_diabetes)
  
  # Return results
  #----------------
  return(list(
    info_df = data.frame(
      method = "Mult. DTs",
      metric = NA,
      mse_test_dat_tree = mse_regr_dt,
      mse_test_dat_rf = mse_regr_rf,
      probs_quantiles = NA,
      epsilon = NA,
      min.bucket = min.bucket,
      significance_level = NA,
      r.squared_ranger = rf$r.squared,
      num.splits = paste0(num_splits_regr_dt, collapse = ","),
      runtime = time,
      n_test = nrow(test_data),
      n_train = nrow(train_data),
      n_cal = nrow(cal_data),
      brier_score5.7 = brier_score5.7,
      brier_score6.5 = brier_score6.5,
      interval_width = NA,
      coverage = NA,
      max_depth = max_depth,
      num_leaves = num_leaves,
      num_vars_split = num_vars_split
    ) %>% cbind(information_df),
    regression_trees = regression_tree,
    probability_trees5.7 = probability_tree5.7,
    probability_trees6.5 = probability_tree6.5
  ))
}
