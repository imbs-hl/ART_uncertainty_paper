#' Multiple Decision Trees (DT) for continuous and threshold outcomes
#'
#' Trains two separate DTs:
#' 1. A regression DT for the continuous outcome.
#' 2. A probability DT for a binary threshold outcome (e.g., y >= 0.5).
#' Calculates performance metrics including MSE, Brier score, tree depth, number of leaves, and sparsity.
#'
#' @param data        Original dataset (for reference)
#' @param instance    List containing: 
#'                      [[1]] training data, 
#'                      [[2]] test data, 
#'                      [[3]] calibration data, 
#'                      [[4]] trained random forest (ranger object), 
#'                      [[5]] additional information data frame
#' @param min.bucket  Minimum number of training observations in a leaf (default 0)
#' @param ...         Additional parameters passed to underlying functions

get_seperate_dt <- function(data, instance, min.bucket = 0, ...){
  
  # Extract data from instance
  #----------------
  train_data     <- instance[[1]]
  test_data      <- instance[[2]]
  cal_data       <- instance[[3]]
  rf             <- instance[[4]]
  information_df <- instance[[5]]
  
  # Compute MSE of random forest
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$y - pred_regr_rf)^2)
  
  # Build regression DT for continuous outcome
  #----------------
  start <- proc.time()
  art_regr <- ranger(y ~ ., data = train_data, 
                     min.bucket = min.bucket, 
                     mtry = ncol(train_data) - 1,
                     sample.fraction = 1, replace = FALSE, num.trees = 1)
  
  # Predictions from regression DT
  pred_regr_art <- predict(art_regr, test_data)$predictions
  
  # Build probability DT for threshold outcome (y >= 0.5)
  train_data_0.5 <- train_data %>% mutate(y = as.factor(ifelse(y >= 0.5, 1, 0)))
  
  probability_tree <- ranger(y ~ ., data = train_data_0.5, 
                             min.bucket = min.bucket, 
                             mtry = ncol(train_data) - 1,
                             sample.fraction = 1, replace = FALSE, num.trees = 1,
                             probability = TRUE)
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # Calculate performance measures
  #----------------
  mse_regr_art <- mean((test_data$y - pred_regr_art)^2)
  
  # Average number of splits across both trees
  num_splits_regr_art <- mean(
    c(nrow(treeInfo(art_regr) %>% filter(!terminal)),
      nrow(treeInfo(probability_tree) %>% filter(!terminal)))
  )
  
  # Brier score for threshold 0.5
  test_data_leafs <- predict(art_regr, test_data, type = "terminalNodes")$predictions
  
  brier_score_df <- data.frame(
    y = test_data$y, 
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
  
  max_depth <- mean(max(get_max_depth(art_regr)), max(get_max_depth(probability_tree)))
  
  # Average number of terminal nodes and sparsity
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
      method               = "Mult. DTs",
      metric               = NA,
      mse_test_dat_tree    = mse_regr_art,
      mse_test_dat_rf      = mse_regr_rf,
      probs_quantiles      = NA,
      epsilon              = NA,
      min.bucket           = min.bucket,
      significance_level   = NA,
      r.squared_ranger     = rf$r.squared,
      num.splits           = paste0(num_splits_regr_art, collapse = ","),
      runtime              = time,
      n_test               = nrow(test_data),
      n_train              = nrow(train_data),
      n_cal                = nrow(cal_data),
      brier_score0.5       = brier_score0.5,
      interval_width       = NA,
      coverage             = NA,
      max_depth            = max_depth,
      num_leaves           = num_leaves,
      num_vars_split       = num_vars_split
    ) %>% cbind(information_df),
    regression_trees = art_regr,
    probability_trees = probability_tree
  ))
}