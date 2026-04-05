# Function to compute variable usage in a random forest (RF)
# Calculates the relative frequency with which each variable is used for splitting
# Variable usage is averaged across all trees in the forest
# Results are returned together with fold-specific meta information

get_rf_variable_usage <- function(data, instance, ...){
  
  # Extract data from instance
  #----------------
  train_data     <- instance[[1]]
  test_data      <- instance[[2]]
  cal_data       <- instance[[3]]
  rf             <- instance[[4]]
  information_df <- instance[[5]]
  
  # Calculate percentage usage of each variable for splitting
  # Variables not used in a tree appear as NA and are later set to 0
  df <- list()
  
  for(i in 1:rf$num.trees){
    tree_info_df <- treeInfo(rf, i)
    
    # Extract split variable names (remove terminal nodes)
    split_var_names <- tree_info_df$splitvarName %>% na.omit()
    num_splits <- length(split_var_names)
    
    # Relative frequency of splits per variable
    table <- data.frame(table(split_var_names) / num_splits)
    df_i <- data.frame(t(table$Freq))
    colnames(df_i) <- table$split_var_names
    
    df[[i]] <- df_i
  }
  
  # Combine all trees
  df <- bind_rows(df)
  
  # Replace NA with 0 (variable not used in a tree)
  df <- df %>% dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))
  
  # Average variable usage across trees
  variable_usage <- colMeans(df)
  
  # Return results
  #----------------
  return(list(
    info_df = data.frame(
      method             = "RF",
      metric             = NA, 
      mse_test_dat_tree  = NA,
      mse_test_dat_rf    = NA,
      probs_quantiles    = NA,
      epsilon            = NA,
      min.bucket         = NA,
      significance_level = NA,
      r.squared_ranger   = rf$r.squared,
      num.splits         = NA,
      runtime            = NA,
      brier_score5.7     = NA,
      brier_score6.5     = NA,
      interval_width     = NA,
      coverage           = NA,
      max_depth          = NA,
      num_leaves         = NA,
      num_vars_split     = NA,
      t(variable_usage)
    ) %>% cbind(information_df),
    
    tree      = rf,
    test_data = test_data,
    cal_data  = cal_data,
    prob_df   = NA
  ))
}
