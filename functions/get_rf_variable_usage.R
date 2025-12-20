
get_rf_variable_usage <- function(data, instance, ...){
  # Exctract data from instance 
  #----------------
  train_data     <- instance[[1]]
  test_data      <- instance[[2]]
  cal_data       <- instance[[3]]
  rf             <- instance[[4]]
  information_df <- instance[[5]]
  
  # Wie oft wurde welche Variable prozentual zum splitten benutzt
  # wurde eine Variable nicht benutzt, wird sie automatisch zu NA, also tauchte sie 0 mal auf
  df <- list()
  for(i in 1:rf$num.trees){
    tree_info_df <- treeInfo(rf, i)
    split_var_names <- tree_info_df$splitvarName %>% na.omit()
    num_splits <- length(split_var_names)
    table <- data.frame(table(split_var_names)/num_splits)
    df_i <- data.frame(t(table$Freq))
    colnames(df_i) <- table$split_var_names
    df[[i]] <- df_i
  }
  df <- bind_rows(df)
  df <- df%>% dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))
  
  variable_usage <- colMeans(df)
  
  # return model, test_& cal data, performance df
  #----------------
    return(list(info_df = data.frame(method     = "RF",
                                   metric               = NA, 
                                   mse_test_dat_tree    = NA,
                                   mse_test_dat_rf      = NA,
                                   probs_quantiles      = NA,
                                   epsilon              = NA,
                                   min.bucket           = NA,
                                   significance_level   = NA,
                                   r.squared_ranger     = rf$r.squared,
                                   num.splits           = NA,
                                   runtime              = NA,
                                   #n_test = nrow(test_data),
                                   #n_train = nrow(train_data),
                                   #n_cal = nrow(cal_data),
                                   brier_score5.7 = NA,
                                   brier_score6.5 = NA, # ab hier
                                   interval_width = NA,
                                   coverage = NA,
                                   max_depth = NA,
                                   num_leaves = NA,
                                   num_vars_split = NA,
                                   t(variable_usage)
  ) %>% cbind(information_df),
  tree = rf,
  test_data = test_data,
  cal_data = cal_data,
  prob_df = NA
  ))


}
