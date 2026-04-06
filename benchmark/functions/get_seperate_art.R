#' Multiple Artificial Representative Tree (ART) for continuous and threshold outcomes
#'
#' Trains two ARTs:
#' 1. A regression ART for the continuous outcome.
#' 2. A probability ART for a binary threshold (here y >= 0.5).
#' Calculates performance metrics including MSE, Brier score, tree depth, number of leaves, and sparsity.
#'
#' @param data                Original dataset (for reference)
#' @param instance            List containing: 
#'                              [[1]] training data, 
#'                              [[2]] test data, 
#'                              [[3]] calibration data, 
#'                              [[4]] trained random forest (ranger object), 
#'                              [[5]] additional information data frame,
#'                              [[6]] trained random forest for probability tree (rf_prob0.5)
#' @param metric              Performance metric used for tree building
#' @param probs_quantiles     Quantiles for ART splitting (you can save time here)
#' @param epsilon             Continue adding more nodes to the ART if the similarity remains the same but the prediction improves by 1 - epsilon
#' @param min.bucket          Minimum number of training observations in a leaf node (default 0)
#' @param significance_level  Significance level for conformal prediction intervals (default 0.1)
#' @param ...                 Additional parameters passed to underlying functions

get_seperate_art <- function(data, instance, metric, probs_quantiles, epsilon, min.bucket = 0, significance_level = 0.1, ...){
  
  # Extract data from instance
  #----------------
  train_data     <- instance[[1]]
  test_data      <- instance[[2]]
  cal_data       <- instance[[3]]
  rf             <- instance[[4]]
  information_df <- instance[[5]]
  rf_prob0.5     <- instance[[6]]
  
  # Compute MSE of random forest
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$y - pred_regr_rf)^2)
  
  # Build regression ART (continuous outcome)
  #----------------
  probs_quantiles <- unlist(probs_quantiles)
  start <- proc.time()
  art_regr <- generate_tree(
    rf = rf, metric = metric, train_data = train_data, test_data = train_data, 
    dependent_varname = "y", probs_quantiles = probs_quantiles, min.bucket=min.bucket
  )
  
  # Predictions from regression ART
  pred_regr_art <- predict(art_regr, test_data)$predictions
  
  # Calibration predictions for CPS (not used here but can be used for intervals)
  y_cal <- cal_data$y
  y_test <- test_data$y
  y_cal_pred <- predict(art_regr, cal_data)$predictions
  y_test_pred <- predict(art_regr, test_data)$predictions
  
  # Build probability ART (for threshold outcome, e.g., y >= 0.5)
  train_data_0.5 <- train_data %>% mutate(y = ifelse(y >= 0.5, 1, 0))
  probability_tree <- generate_tree(
    rf = rf_prob0.5, metric = metric, train_data = train_data_0.5, test_data = train_data_0.5, 
    dependent_varname = "y", probs_quantiles = probs_quantiles, min.bucket=min.bucket
  )
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # Performance metrics
  #----------------
  mse_regr_art <- mean((test_data$y - pred_regr_art)^2)
  
  # Average number of splits across both trees
  num_splits_regr_art <- mean(
    c(
      nrow(treeInfo(art_regr) %>% filter(!terminal)),
      nrow(treeInfo(probability_tree) %>% filter(!terminal))
    )
  )
  
  # Brier score for threshold 0.5
  test_data_leafs <- predict(art_regr, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(art_regr, train_data, type = "terminalNodes")$predictions
  cal_data_leafs <- predict(art_regr, cal_data, type = "terminalNodes")$predictions
  
  brier_score_df <- data.frame(
    y = y_test,
    nodeID = test_data_leafs,
    prob0.5 = predict(probability_tree, test_data)$predictions[,2]
  ) %>% mutate(over0.5 = ifelse(y >= 0.5, 1, 0))
  
  brier_score0.5 <- BrierScore(brier_score_df$over0.5, brier_score_df$prob0.5)
  
  # Function to calculate maximum tree depth
  get_max_depth <- function(rf_object, tree = 1) {
    ti <- ranger::treeInfo(rf_object, tree = tree)
    depths <- setNames(integer(nrow(ti)), ti$nodeID)
    depths["0"] <- 0
    for (i in seq_len(nrow(ti))) {
      node <- ti$nodeID[i]
      if (ti$terminal[i]) next
      left  <- ti$leftChild[i]
      right <- ti$rightChild[i]
      depths[as.character(left)]  <- depths[as.character(node)] + 1
      depths[as.character(right)] <- depths[as.character(node)] + 1
    }
    max(depths)
  }
  
  # Average max depth across both trees
  max_depth <- mean(max(get_max_depth(art_regr)), max(get_max_depth(probability_tree)))
  
  # Average number of terminal nodes and variables used for splits
  num_leaves <- mean(
    nrow(treeInfo(art_regr) %>% filter(terminal)),
    nrow(treeInfo(probability_tree) %>% filter(terminal))
  )
  
  num_vars_split <- mean(
    length(na.omit(unique(treeInfo(art_regr)$splitvarID))),
    length(na.omit(unique(treeInfo(probability_tree)$splitvarID)))
  )
  
  # Return results: summary info, regression and probability trees
  #----------------
  return(list(
    info_df = data.frame(
      method              = "Mult. ARTs",
      metric              = metric,
      mse_test_dat_tree   = mse_regr_art,
      mse_test_dat_rf     = mse_regr_rf,
      probs_quantiles     = paste0(probs_quantiles, collapse = ","),
      epsilon             = epsilon,
      min.bucket          = min.bucket,
      significance_level  = NA,
      r.squared_ranger    = rf$r.squared,
      num.splits          = paste0(num_splits_regr_art, collapse = ","),
      runtime             = time,
      n_test              = nrow(test_data),
      n_train             = nrow(train_data),
      n_cal               = nrow(cal_data),
      brier_score0.5      = brier_score0.5,
      interval_width      = NA,
      coverage            = NA,
      max_depth           = max_depth,
      num_leaves          = num_leaves,
      num_vars_split      = num_vars_split
    ) %>% cbind(information_df),
    regression_trees = art_regr,
    probability_trees = probability_tree
  ))
}