

simulate_benchmark_data <- function(data, sim_id, current_fold, num_folds, num.trees, ...){
  # Get benchmark data set
  task = getOMLTask(sim_id)
  
  # Save data and needed information
  dat = task$input$data.set$data
  target_variable = task$input$data.set$target.features
  task_id = task$task.id
  task_name = task$input$data.set$desc$name
  task_version = task$input$data.set$desc$version
  
  # Remove NAs from data and all variables with only one value for all observations save dimension
  dat_na_free <- na.omit(dat)

  task_rows <- nrow(dat)
  task_cols <- ncol(dat)
  
  task_rows_na <- nrow(dat_na_free)
  task_cols_na <- ncol(dat_na_free)
  
  # Function for removing variables with constant values in all rows
  remove_constant_vars <- function(data) {
    constant_vars <- sapply(data, function(col) length(unique(col)) == 1)
    data <- data[, !constant_vars, drop = FALSE]
    return(data)
  }
  dat_na_free <- remove_constant_vars(dat_na_free)
  
  task_rows_na_unique <- nrow(dat_na_free)
  task_cols_na_unique <- ncol(dat_na_free)
  
  # Remove outliers: Median +/- 5* IQR
  median_y <- median(dat_na_free[,target_variable])
  iqr_y <- IQR(dat_na_free[,target_variable])
  
  lower_border <- median_y - 5 * iqr_y
  upper_border <- median_y + 5 * iqr_y
  
  dat_na_free <- dat_na_free[dat_na_free[,target_variable] > lower_border & 
                                dat_na_free[,target_variable] < upper_border,]
  
  task_rows_na_unique_outlier <- nrow(dat_na_free)
  task_cols_na_unique_outlier <- ncol(dat_na_free)
  
  # target in Bereich 0,1 transformieren
  dat_na_free[,target_variable] <- (dat_na_free[,target_variable]-min(dat_na_free[,target_variable]))/(max(dat_na_free[,target_variable])-min(dat_na_free[,target_variable]))
  
  # split data in train, calibration and test data with 7x cross validation
  # create folds for cv
  set.seed(sim_id)
  folds <- createFolds(dat_na_free[,target_variable], k = num_folds, list = TRUE, returnTrain = TRUE)
  
  # get current fold ids and separate data
  current_fold_ids <- folds[[current_fold]]
  
  # Test data
  test_dat <- dat_na_free[-current_fold_ids,]
  test_dat_y <- dat_na_free[-current_fold_ids,target_variable]
  
  # 2/3 Train, 1/3 Calibration data set
  train_ids <- sample(current_fold_ids, ((2/3)*length(current_fold_ids)))
  train_dat <- dat_na_free[train_ids,]
  
  calibration_ids <- current_fold_ids[which(!(current_fold_ids %in% train_ids))]
  cal_dat <- dat_na_free[calibration_ids,]
  
  x = train_dat[,!(colnames(train_dat) %in% target_variable)]
  
  # Train ranger object
  rf <- ranger(y = train_dat[,target_variable],
               x = x, 
               num.trees = num.trees, min.node.size = (nrow(train_dat)/10), importance = "permutation")
  
  # Parameters of rf
  params <- data.frame(min_node_size = nrow(train_dat)/10,
                       mtry = floor(sqrt(ncol(train_dat))), 
                       num.trees)
  
  # Collect information for results
  information_df <- data.frame(task_id, task_name, task_version, 
                               "n" = task_rows, "p" = task_cols, 
                               "n_na_free" = task_rows_na, 
                               "p_na_free" = task_cols_na,
                               "n_na_free_unique" = task_rows_na_unique, 
                               "p_na_free_unique" = task_cols_na_unique,
                               "n_na_free_unique_outlier" = task_rows_na_unique_outlier, 
                               "p_na_free_unique_outlier" = task_cols_na_unique_outlier,
                               fold = current_fold, target_variable,
                               var_test = var(test_dat_y),
                               mean_test = mean(test_dat_y),
                               task_id = task_id,
                               n_test = length(test_dat_y),
                               var_pop = sum((test_dat_y-mean(test_dat_y))^2),
                               sim_id)
  
  
  # Return 
  return(list(test_dat, rf, params, train_dat, target_variable, information_df, cal_dat))
  
  
}
