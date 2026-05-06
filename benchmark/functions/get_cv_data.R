#' Load external dataset and create cross-validation splits
#' 
#' This function loads a specified external benchmark dataset, preprocesses it,
#' and creates train, test, and calibration splits based on cross-validation.
#' Additionally, it trains a regression random forest and a probability random forest.
#' 
#' @param data            (unused) placeholder for compatibility
#' @param current_fold    Index of the current cross-validation fold
#' @param num_folds       Total number of folds for cross-validation
#' @param repitition      Repetition index (used for reproducibility)
#' @param dataset_name    Name of the dataset to load
#' @param ...             Additional arguments (not used)
#' 
#' @return A list containing:
#' \itemize{
#'   \item train_data_i     Training dataset
#'   \item test_data_i      Test dataset (current fold)
#'   \item cal_data_i       Calibration dataset
#'   \item rf               Trained regression random forest
#'   \item information_df   Data frame with meta information about the split
#'   \item rf_prob0.5       Trained probability random forest (binary outcome)
#' }


get_cv_data <- function(data, current_fold, num_folds, repitition, dataset_name, ...){
  
  #------------------------------------------------------------------
  # Load dataset depending on selected dataset_name
  # Outcome variable is renamed to "y"
  #------------------------------------------------------------------
  
  data_ext_dir <- file.path("external_data")
  if(!file.exists(data_ext_dir)){
    stop("Please set path of data_ext_dir in function get_cv_data.R to folder with downloaded external data.")
  }
  
  if(dataset_name == "abalone"){
    data <- read.table(file.path(data_ext_dir, "abalone/abalone.data"),
                       sep = ",", header = FALSE) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "airfoil"){
    data <- read.table(file.path(data_ext_dir, "airfoil/airfoil_self_noise.dat"),
                       sep = "\t", header = FALSE) %>% 
      rename(y = V6)
    
  }else if(dataset_name == "deltaA"){
    data <- read.table(file.path(data_ext_dir, "deltaA/delta_ail.dat"),
                       sep = ",", header = FALSE, skip = 10) %>% 
      rename(y = V6)
    
  }else if(dataset_name == "deltaE"){
    data <- read.table(file.path(data_ext_dir, "deltaE/delta_elv.dat"),
                       sep = ",", header = FALSE, skip = 11) %>% 
      rename(y = V7)
    
  }else if(dataset_name == "friedm"){
    data <- read.table(file.path(data_ext_dir, "friedm/friedman.dat"),
                       sep = ",", header = FALSE, skip = 10) %>% 
      rename(y = V6)
    
  }else if(dataset_name == "mortgage"){
    data <- read.table(file.path(data_ext_dir, "mortgage/mortgage.dat"),
                       sep = ",", header = FALSE, skip = 20) %>% 
      rename(y = V16)
    
  }else if(dataset_name == "wineW"){
    data <- read.csv(file.path(data_ext_dir, "wineW/winequality-white.csv"),
                     header = TRUE, sep = ";") %>% 
      rename(y = quality)
    
  }else if(dataset_name == "wizmir"){
    data <- read.table(file.path(data_ext_dir, "wizmir/wizmir.dat"),
                       sep = ",", skip = 14) %>% 
      rename(y = V10)
    
  }else if(dataset_name == "puma8fh"){
    data <- read.table(file.path(data_ext_dir, "puma8fh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "puma8fm"){
    data <- read.table(file.path(data_ext_dir, "puma8fm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "puma8nh"){
    data <- read.table(file.path(data_ext_dir, "puma8nh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "puma8nm"){
    data <- read.table(file.path(data_ext_dir, "puma8nm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "kin8fh"){
    data <- read.table(file.path(data_ext_dir, "kin8fh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "kin8fm"){
    data <- read.table(file.path(data_ext_dir, "kin8fm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "kin8nh"){
    data <- read.table(file.path(data_ext_dir, "kin8nh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "kin8nm"){
    data <- read.table(file.path(data_ext_dir, "kin8nm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "bank8fh"){
    data <- read.table(file.path(data_ext_dir, "bank8fh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "bank8fm"){
    data <- read.table(file.path(data_ext_dir, "bank8fm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "bank8nh"){
    data <- read.table(file.path(data_ext_dir, "bank8nh/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "bank8nm"){
    data <- read.table(file.path(data_ext_dir, "bank8nm/Dataset.data")) %>% 
      rename(y = V9)
    
  }else if(dataset_name == "comp"){
    data <- read.table(file.path(data_ext_dir, "comp/Prototask.data")) %>% 
      rename(y = V13)
  }
  
  
  #------------------------------------------------------------------
  # Normalize outcome variable y to range [0, 1]
  #------------------------------------------------------------------
  data$y <- (data$y - min(data$y)) / (max(data$y) - min(data$y))
  
  
  #------------------------------------------------------------------
  # Create cross-validation splits
  #------------------------------------------------------------------
  set.seed(123 + repitition)
  folds <- createFolds(data[,"y"], k = num_folds, list = TRUE, returnTrain = FALSE)
  
  set.seed(current_fold)
  
  # Define test set based on current fold
  cv_ids_i <- folds[[current_fold]]
  test_data_i <- data[cv_ids_i,]
  
  # Remaining data used for training and calibration
  data_i <- data[-cv_ids_i,]
  n_i <- nrow(data_i)
  
  # Randomly split into training (2/3) and calibration (1/3)
  train_ids_i <- sample(c(1:(n_i - length(cv_ids_i))), n_i * (2/3), replace = FALSE)
  train_data_i <- data_i[train_ids_i,]
  cal_data_i <- data_i[-train_ids_i,]
  n_train_i <- nrow(train_data_i)
  
  
  #------------------------------------------------------------------
  # Train models
  #------------------------------------------------------------------
  
  # Train regression random forest
  rf <- ranger(y ~ ., data = train_data_i, min.node.size = (n_train_i * 0.1))
  
  # Train probability random forest (binary outcome y >= 0.5)
  train_data_0.5 <- train_data_i %>% 
    mutate(y = as.factor(ifelse(y >= 0.5, 1, 0)))
  
  rf_prob0.5 <- ranger(y ~ ., data = train_data_0.5, 
                       min.node.size = (n_train_i * 0.1), 
                       probability = TRUE)
  
  
  #------------------------------------------------------------------
  # Collect meta information about current split
  #------------------------------------------------------------------
  y_test <- test_data_i$y
  
  information_df <- data.frame(
    fold         = current_fold,
    var_test_y   = var(y_test),
    mean_test_y  = mean(y_test),
    var_pop_y    = sum((y_test - mean(y_test))^2),
    repitition   = repitition,
    dataset_name = dataset_name
  )
  
  
  #------------------------------------------------------------------
  # Return datasets, models, and metadata
  #------------------------------------------------------------------
  return(list(train_data_i, test_data_i, cal_data_i, rf, information_df, rf_prob0.5))
}