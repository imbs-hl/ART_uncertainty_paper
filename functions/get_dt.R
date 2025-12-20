# Function to construct an decision tree (DT)
# Predicts a continuous outcome (glycohemoglobin) using a pre-trained random forest
# Estimates leaf-wise probabilities for prediabetes and diabetes via conformal predictive systems (CPS)
# Provides uncertainty quantification through calibrated prediction intervals
# Returns the DT, probability estimates per leaf, and multiple performance metrics

get_dt <- function(data, instance, min.bucket = 0, significance_level = 0.1, ...){
  
  # Extract components from instance
  #----------------
  train_data     <- instance[[1]]  # training data
  test_data      <- instance[[2]]  # test data
  cal_data       <- instance[[3]]  # calibration data
  rf             <- instance[[4]]  # trained random forest
  information_df <- instance[[5]]  # meta information
  
  # RF mean squared error on test data
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$glycohemoglobin - pred_regr_rf)^2)
  
  # Build decision tree (DT)
  #----------------
  start <- proc.time()  # start runtime measurement
  
  # Train regression decision tree (single tree ranger)
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
  
  # Build CPS
  #----------------
  # Prepare calibration and test values
  y_cal <- cal_data$glycohemoglobin
  y_test <- test_data$glycohemoglobin
  y_cal_pred <- predict(art_regr, cal_data)$predictions
  y_test_pred <- predict(art_regr, test_data)$predictions
  
  # Compute calibrated predictive system (Mondrian CPS)
  cpd_art <- get_calibrated_predictive_system_mondrian(
    y_cal_pred = y_cal_pred,
    y_cal = y_cal,
    y_test_pred = y_test_pred,
    significance_level = significance_level,
    interval_type = "two-tailed",
    tree = art_regr,
    cal_data = cal_data,
    test_data = test_data,
    show_node_id = TRUE,
    dependent_varname = "glycohemoglobin"
  ) %>%
    mutate(nodeID = leaf - 1)
  
  # Merge tree structure with CPS results
  tree_info_df <- treeInfo(art_regr) %>%
    left_join(unique(cpd_art)) %>%
    mutate(
      lower_bound = round(lower_bound, 2),
      upper_bound = round(upper_bound, 2),
      prediction  = round(prediction, 2)
    )
  
  # Predicted vs true probabilities per node
  #----------------
  # Vector of terminal nodes
  leafs <- tree_info_df %>% filter(terminal) %>% select(nodeID) %>% unlist() %>% unname()
  
  # Assign observations to terminal nodes
  test_data_leafs  <- predict(art_regr, test_data,  type = "terminalNodes")$predictions
  train_data_leafs <- predict(art_regr, train_data, type = "terminalNodes")$predictions
  cal_data_leafs   <- predict(art_regr, cal_data,   type = "terminalNodes")$predictions
  
  probs_df <- list()
  
  for(i in 1:length(leafs)){
    
    # Calibration data within leaf
    cal_dat_leaf <- cal_data[cal_data_leafs == leafs[i], ]
    
    # Single test observation from leaf
    test_dat_leaf <- test_data[test_data_leafs == leafs[i], ][1, ]
    
    # Predictive distribution for leaf
    cpd_dist <- get_predictive_distribution(
      y_cal = cal_dat_leaf$glycohemoglobin,
      y_cal_pred = predict(art_regr, cal_dat_leaf)$predictions,
      y_test_pred = predict(art_regr, test_dat_leaf)$predictions
    )[[1]]
    
    # Construct cumulative probability grid
    p <- seq(0, 1 - 1/length(cpd_dist), length.out = length(cpd_dist))
    data_cps <- data.frame(
      cpd = c(cpd_dist, max(cpd_dist + 0.000001), max(cal_data$glycohemoglobin)),
      p   = c(p, 1, 1)
    )
    
    # Probability to exceed clinical thresholds
    p_5.7 <- data_cps$p[which(data_cps$cpd > 5.7)[1]]
    p_6.5 <- data_cps$p[which(data_cps$cpd > 6.5)[1]]
    
    probs_df[[i]] <- c(
      nodeID = leafs[i],
      prob5.7 = 1 - p_5.7,
      prob6.5 = 1 - p_6.5
    )
    
    # Two-sided prediction interval bounds
    lower_index <- tail(which(data_cps$p <= 0.025), 1)
    upper_index <- which(data_cps$p >= 0.975)[1]
    
    low_percentile  <- data_cps$cpd[lower_index]
    high_percentile <- data_cps$cpd[upper_index]
  }
  
  probs_df <- left_join(tree_info_df, bind_rows(probs_df))
  
  # Empirical frequencies of test data exceeding thresholds per leaf
  prob_df_i <- data.frame(
    nodeID = predict(art_regr, test_data, type = "terminalNodes")$prediction,
    glycohemoglobin = test_data$glycohemoglobin
  ) %>%
    group_by(nodeID) %>%
    dplyr::summarise(
      n_total = n(),                                   # total observations
      n_over_5_7 = sum(glycohemoglobin > 5.7, na.rm = TRUE),
      proportion_over_5_7_test_data = mean(glycohemoglobin > 5.7, na.rm = TRUE),
      n_over_6_5 = sum(glycohemoglobin > 6.5, na.rm = TRUE),
      proportion_over_6_5_test_data = mean(glycohemoglobin > 6.5, na.rm = TRUE)
    ) %>%
    left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound, prob5.7, prob6.5))
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])  # total runtime
  
  # Performance measures
  #----------------
  # DT mean squared error
  mse_regr_art <- mean((test_data$glycohemoglobin - pred_regr_art)^2)
  
  # Number of splits
  num_splits_regr_art <- nrow(treeInfo(art_regr) %>% filter(!terminal))
  
  # Brier scores
  brier_score_df <- data.frame(
    y = y_test,
    nodeID = test_data_leafs
  ) %>%
    left_join(prob_df_i %>% select(nodeID, prob5.7, prob6.5)) %>%
    mutate(
      over5.7 = ifelse(y >= 5.7, 1, 0),
      over6.5 = ifelse(y >= 6.5, 1, 0)
    )
  
  brier_score5.7 <- BrierScore(brier_score_df$over5.7, brier_score_df$prob5.7)
  brier_score6.5 <- BrierScore(brier_score_df$over6.5, brier_score_df$prob6.5)
  
  # Coverage of prediction intervals
  coverage_df <- data.frame(
    y = y_test,
    nodeID = test_data_leafs
  ) %>%
    left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound)) %>%
    mutate(inside_interval = y >= lower_bound & y <= upper_bound)
  
  coverage <- sum(coverage_df$inside_interval) / nrow(coverage_df)
  
  # Average interval width
  interval_width_df <- cpd_art %>%
    select(nodeID, lower_bound, upper_bound) %>%
    unique()
  
  interval_width <- mean(interval_width_df$upper_bound - interval_width_df$lower_bound)
  
  # Tree depth calculation
  ti <- treeInfo(art_regr)
  
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
  
  max_depth <- max(depths)
  
  # Number of terminal nodes
  num_leaves <- nrow(treeInfo(art_regr) %>% filter(terminal))
  
  # Sparsity: number of splitting variables
  num_vars_split <- length(na.omit(unique(treeInfo(art_regr)$splitvarID)))
  
  # Return results
  #----------------
  return(list(
    info_df = data.frame(
      method = "Regression DT + CPS",
      metric = NA,
      mse_test_dat_tree = mse_regr_art,
      mse_test_dat_rf = mse_regr_rf,
      probs_quantiles = NA,
      epsilon = NA,
      min.bucket = min.bucket,
      significance_level = significance_level,
      r.squared_ranger = rf$r.squared,
      num.splits = paste0(num_splits_regr_art, collapse = ","),
      runtime = time,
      n_test = nrow(test_data),
      n_train = nrow(train_data),
      n_cal = nrow(cal_data),
      brier_score5.7 = brier_score5.7,
      brier_score6.5 = brier_score6.5,
      interval_width = interval_width,
      coverage = coverage,
      max_depth = max_depth,
      num_leaves = num_leaves,
      num_vars_split = num_vars_split
    ) %>% cbind(information_df),
    tree = art_regr,
    test_data = test_data,
    cal_data = cal_data,
    prob_df = prob_df_i %>% select(nodeID, prob5.7, prob6.5)
  ))
}