# Function to generate cross-validation (CV) data for one fold
# Splits the data into training, calibration, and test sets
# The selected fold is used as test data, remaining data is split into train and calibration sets
# Trains a regression random forest and two probability random forests (prediabetes and diabetes)
# Returns all data splits and trained models for artificial representative tree (ART) construction

get_cv_data <- function(data, current_fold, num_folds, repitition, ...){
  
  # Load data
  #----------------
  data_dir <- "/imbs/home/kronziel/ART_uncertainty/data/nhanes_prepandemic_complete.csv"
  data_imp <- read.csv(data_dir)
  
  # Remove ID and prediabetes column
  data <- data_imp %>% 
    select(-SEQN, -prediabetes)
  
  # Cross-validation setup
  #----------------
  set.seed(123 + repitition)
  folds <- createFolds(
    data[, "glycohemoglobin"],
    k = num_folds,
    list = TRUE,
    returnTrain = FALSE
  )
  
  set.seed(current_fold)
  
  # Split into training, test, and calibration data
  cv_ids_i <- folds[[current_fold]]
  test_data_i <- data[cv_ids_i, ]
  data_i <- data[-cv_ids_i, ]
  
  n_i <- nrow(data_i)
  
  train_ids_i <- sample(
    c(1:(n_i - length(cv_ids_i))),
    n_i * (2 / 3),
    replace = FALSE
  )
  
  train_data_i <- data_i[train_ids_i, ]
  cal_data_i <- data_i[-train_ids_i, ]
  n_train_i <- nrow(train_data_i)
  
  # Train RF and probability forests
  #----------------
  # Regression random forest
  rf <- ranger(
    glycohemoglobin ~ .,
    data = train_data_i,
    min.node.size = (n_train_i * 0.1)
  )
  
  # Probability RF for prediabetes (>= 5.7)
  train_data_5.7 <- train_data_i %>% 
    mutate(glycohemoglobin = ifelse(glycohemoglobin >= 5.7, 1, 0))
  
  rf_prob5.7 <- ranger(
    glycohemoglobin ~ .,
    data = train_data_5.7,
    min.node.size = (n_train_i * 0.1),
    probability = TRUE
  )
  
  # Probability RF for diabetes (>= 6.5)
  train_data_6.5 <- train_data_i %>% 
    mutate(glycohemoglobin = ifelse(glycohemoglobin >= 6.5, 1, 0))
  
  rf_prob6.5 <- ranger(
    glycohemoglobin ~ .,
    data = train_data_6.5,
    min.node.size = (n_train_i * 0.1),
    probability = TRUE
  )
  
  # Collect information for results
  #----------------
  y_test <- test_data_i$glycohemoglobin
  
  information_df <- data.frame(
    fold = current_fold,
    var_test_y = var(y_test),
    mean_test_y = mean(y_test),
    var_pop_y = sum((y_test - mean(y_test))^2),
    repitition = repitition
  )
  
  # Return objects
  #----------------
  return(list(
    train_data_i,
    test_data_i,
    cal_data_i,
    rf,
    information_df,
    rf_prob5.7,
    rf_prob6.5
  ))
}
