#' Artificial Representative Tree (ART) with Conformal Predictive System (CPS)
#' ---------------
#' This function builds an ART based on a trained random forest and calculates 
#' calibrated predictive intervals using a Mondrian CPS. 
#' Computes performance metrics including MSE, Brier score, coverage, interval width, tree depth, 
#' and sparsity.
#' 
#' @param data                Input parameter of package batchtools that is needed from the package (but not used here)
#' @param instance            List containing: 
#'                              [[1]] training data, 
#'                              [[2]] test data, 
#'                              [[3]] calibration data, 
#'                              [[4]] trained random forest model (ranger object), 
#'                              [[5]] additional information of simulated data as data.frame
#' @param metric              Performance metric used for tree building
#' @param probs_quantiles     Quantiles for ART splitting (you can save time here)
#' @param epsilon             Continue adding more nodes to the ART if the similarity remains the same but the prediction improves by 1 - epsilon
#' @param min.bucket          Minimum number of training observations in a leaf node (default 0)
#' @param significance_level  Significance level for conformal prediction intervals (default 0.1)
#' @param ...                 Additional parameters passed to underlying functions
get_art <- function(data, instance, metric, probs_quantiles, epsilon, min.bucket = 0, significance_level = 0.1, ...){
  
  # Extract data from instance
  #----------------
  train_data     <- instance[[1]]  # training dataset
  test_data      <- instance[[2]]  # test dataset
  cal_data       <- instance[[3]]  # calibration dataset
  rf             <- instance[[4]]  # random forest model
  information_df <- instance[[5]]  # additional information/data frame
  
  # Calculate MSE of random forest predictions
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$y - pred_regr_rf)^2)
  
  # Build ART (Adaptive Regression Tree)
  #----------------
  probs_quantiles <- unlist(probs_quantiles)
  start <- proc.time()
  
  # Train regression ART on training data
  regression_tree <- generate_tree(rf = rf, metric = metric, train_data = train_data, test_data = train_data, 
                                   dependent_varname = "y", 
                                   probs_quantiles = probs_quantiles, min.bucket=min.bucket)
  
  # Predictions from regression ART
  pred_regr_art <- predict(regression_tree, test_data)$predictions
  
  
  # Build CPS (Conformal Predictive System)
  #----------------
  # Calculate probabilities for CPS
  y_cal <- cal_data$y
  y_test <- test_data$y
  y_cal_pred <- predict(regression_tree, cal_data)$predictions
  y_test_pred <- predict(regression_tree, test_data)$predictions
  
  # Generate calibrated predictive distributions using Mondrian conformal prediction
  cpd_art <- get_calibrated_predictive_system_mondrian(y_cal_pred = y_cal_pred, y_cal = y_cal, y_test_pred = y_test_pred, 
                                                       significance_level = significance_level, interval_type = "two-tailed", 
                                                       tree = regression_tree, 
                                                       cal_data = cal_data, test_data = test_data, show_node_id = TRUE, 
                                                       dependent_varname = "y") %>% 
    mutate(nodeID = leaf-1)
  
  # Join tree info with calibrated predictions and round values
  tree_info_df <- treeInfo(regression_tree) %>% 
    left_join(unique(cpd_art)) %>% 
    mutate(lower_bound = round(lower_bound, 2),
           upper_bound = round(upper_bound, 2),
           prediction = round(prediction, 2))
  
  # Predictions vs true probability per leaf node
  # Generate calibrated predictive distributions per terminal node
  leafs <- tree_info_df %>% filter(terminal) %>% select(nodeID) %>% unlist() %>% unname()
  
  # Determine which leaf each observation belongs to
  test_data_leafs <- predict(regression_tree, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(regression_tree, train_data, type = "terminalNodes")$predictions
  cal_data_leafs <- predict(regression_tree, cal_data, type = "terminalNodes")$predictions
  
  probs_df <- list()
  
  for(i in 1:length(leafs)){
    # All calibration data in this leaf
    cal_dat_leaf <- cal_data[cal_data_leafs==leafs[i],]
    # A single test observation in this leaf
    test_dat_leaf <- test_data[test_data_leafs==leafs[i],][1,]
    # Generate predictive distribution for the test observation
    cpd_dist <- get_predictive_distribution(y_cal = cal_dat_leaf$y, 
                                            y_cal_pred = predict(regression_tree, cal_dat_leaf)$predictions, 
                                            y_test_pred = predict(regression_tree, test_dat_leaf)$predictions)[[1]]
    # Calculate "p-value"
    p <- seq(0, 1 - 1/length(cpd_dist), length.out = length(cpd_dist))  
    data_cps <- data.frame(cpd = c(cpd_dist,max(cpd_dist+0.000001), max(cal_data$y)), p = c(p,1,1))
    
    # Probability to exceed the threshold 0.5
    p_0.5 <- data_cps$p[which(data_cps$cpd > 0.5)[1]]
    probs_df[[i]] <- c(nodeID=leafs[i], prob0.5=1-p_0.5)
    
    # Upper and lower bounds of prediction intervals for two-tailed
    lower_index <- tail(which(data_cps$p <= 0.025), 1)
    upper_index <- which(data_cps$p >= 0.975)[1]
    
    low_percentile <-  data_cps$cpd[lower_index]
    high_percentile <-  data_cps$cpd[upper_index]
  }
  
  probs_df <- left_join(tree_info_df, bind_rows(probs_df))
  
  # Relative and absolute frequency of test data exceeding thresholds per leaf
  prob_df_i <- data.frame(nodeID = predict(regression_tree, test_data, type = "terminalNodes")$prediction,
                          y = test_data$y) %>% 
    group_by(nodeID)  %>%
    dplyr::summarise(
      n_total = n(),  # Total number of observations per leaf
      n_over_5_7 = sum(y > 0.5, na.rm = TRUE),  # Number of observations over 0.5
      proportion_over_5_7_test_data = mean(y > 0.5, na.rm = TRUE),  # Proportion over 0.5
    ) %>% left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound, prob0.5))
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # Calculate performance measures
  #----------------
  # MSE of regression tree
  mse_regr_art <- mean((test_data$y - pred_regr_art)^2)
  
  # Number of splits in the regression ART
  num_splits_regr_art <- nrow(treeInfo(regression_tree) %>% filter(!terminal))
  
  # Brier Score per test observation
  brier_score_df <- data.frame(y = y_test, 
                               nodeID = test_data_leafs) %>% 
    left_join(prob_df_i %>% select(nodeID, prob0.5)) %>% 
    mutate(over0.5 = ifelse(y>= 0.5,1,0))
  
  brier_score0.5 = BrierScore(brier_score_df$over0.5, brier_score_df$prob0.5)
  
  # Coverage of test data (proportion of true y within prediction interval)
  coverage_df <- data.frame(y = y_test, 
                            nodeID = test_data_leafs) %>% 
    left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound)) %>% 
    mutate(inside_interval = y >= lower_bound & y <= upper_bound)
  
  coverage <- sum(coverage_df$inside_interval)/nrow(coverage_df)
  
  # Average interval width
  interval_width_df <- cpd_art %>% 
    select(nodeID, lower_bound, upper_bound) %>% 
    unique()
  
  interval_width <- mean(interval_width_df$upper_bound - interval_width_df$lower_bound)
  
  # Maximum tree depth
  ti <- treeInfo(regression_tree)
  
  # Ranger returns tree as node list with child nodes
  depths <- setNames(integer(nrow(ti)), ti$nodeID)
  depths["0"] <- 0  # Root node depth = 0
  
  # Calculate depth for each node
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
  num_leaves <- nrow(treeInfo(regression_tree) %>% filter(terminal))
  
  # Sparsity: number of variables used for splitting
  num_vars_split <- length(na.omit(unique(treeInfo(regression_tree)$splitvarID)))
  
  # Return model, test & calibration data, and performance summary
  #----------------
  return(list(info_df = data.frame(method               = "ART + CPS",
                                   metric               = metric, 
                                   mse_test_dat_tree    = mse_regr_art,
                                   mse_test_dat_rf      = mse_regr_rf,
                                   probs_quantiles      = paste0(probs_quantiles, collapse = ","),
                                   epsilon              = epsilon,
                                   min.bucket           = min.bucket,
                                   significance_level   = significance_level,
                                   r.squared_ranger     = rf$r.squared,
                                   num.splits           = paste0(num_splits_regr_art, collapse = ","),
                                   runtime              = time,
                                   n_test               = nrow(test_data),
                                   n_train              = nrow(train_data),
                                   n_cal                = nrow(cal_data),
                                   brier_score0.5       = brier_score0.5,
                                   interval_width       = interval_width,
                                   coverage             = coverage,
                                   max_depth            = max_depth,
                                   num_leaves           = num_leaves,
                                   num_vars_split       = num_vars_split
  ) %>% cbind(information_df),
  regression_trees = regression_tree,
  prob_df = prob_df_i %>% select(nodeID, prob0.5)
  ))
}