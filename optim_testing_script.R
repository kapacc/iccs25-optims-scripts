source('datasets_prep.R')
source('scoring_methods.R')
source('_optims_wrapper.R')

usf <- c("Nelder-Mead"
         ,"BFGS"
         ,"CG"
         ,"SANN"
         ,"NLOPT_LD_SLSQP"
         ,"NLOPT_LD_TNEWTON_PRECOND"
         ,"NLOPT_LN_PRAXIS"
         , "naive")

datasets <- load_datasets(more_then_n_unique_values = 2, dataset_type = 'n')
list2env(datasets, envir = .GlobalEnv)


param_grid <- expand.grid(optims = usf,
                          dsets = c(
                                    'bank_marketing',
                                    'wine_quality',
                                    'breastc',
                                    'diabetes',
                                    'artif',
                                    'banknote',
                                    'credit_a',
                                    'adult' ,
                                    'credit_g',
                                    'heart_c',
                                    'spambase',
                                    'dhfr'
                                    ),
                          c = c(0.3, 0.5, 0.7, 0.9),
                          n_resamp = 1:100, stringsAsFactors = FALSE)

param_grid_addon <- expand.grid(optims = usf,
                                dsets = c('kdd_data'),
                                c = c(0.3, 0.5, 0.7, 0.9),
                                n_resamp = 1:10, stringsAsFactors = FALSE)

param_grid <- rbind(param_grid, param_grid_addon)


library(caret)
library(pROC)


evaluate_function <- function(optims, dsets, c, n_resamp) {
  
  set.seed(n_resamp)
  
  current_frame <- get(dsets)
  
  current_frame$Y <- as.numeric(current_frame$Y)
  
  current_frame$S <- current_frame$Y * rbinom(nrow(current_frame), 1, as.numeric(c))
  current_frame <- as.data.frame(sapply(current_frame, as.numeric))
  
  rowtrain <- createDataPartition(current_frame$Y, p = 0.8)[[1]]
  
  train <- current_frame[rowtrain, ]
  test <- current_frame[-rowtrain, ]
  
  record_start_dt <- Sys.time()
  
   if (optims == "naive") {
    {
      model <- glm(S ~ ., data = subset(train, select = -Y), family = gaussian)
    
      predictions <- predict(model, newdata = subset(test, select = -Y))
      predictions <- factor(ifelse(predictions > 0.5, 1, 0), levels = c(0, 1)) 
      
      } %>% system.time(.) -> time
    
    time <- time["elapsed"] %>% unname
  } else {
    x_train <- as.matrix(subset(train, select = -c(Y, S)))
    s_train <- train$S
    x_test <- as.matrix(subset(test, select = -c(Y, S)))
    y_test <- test$Y
    
    {
      predictions <- logistic_joint(x_train, s_train, x_test, optim_method = optims)
      predictions <- factor(predictions, levels = c(0, 1))
    } %>% system.time(.) -> time
    
    time <- time["elapsed"] %>% unname
  }
  
  y_test <- factor(test$Y, levels = c(0, 1))
  
  acc <- mean(predictions == y_test)
  
  confusion <- confusionMatrix(predictions, y_test)
  precision <- confusion$byClass["Pos Pred Value"] %>% as.numeric()
  recall <- confusion$byClass["Sensitivity"] %>% as.numeric()
  f1 <- 2 * ((precision * recall) / (precision + recall)) %>% as.numeric()
  
  roc_curve <- roc(y_test, as.numeric(predictions), quiet = TRUE)
  auc <- auc(roc_curve) %>% as.numeric()
  
  record_end_dt <- Sys.time()
  
  return(list(
    accuracy = acc,
    precision = precision,
    recall = recall,
    f1_score = f1,
    auc = auc,
    time = time,
    ncol_train = ncol(train)
  ))
}



results <- list()  
short_commit_hash <- system("git rev-parse --short HEAD", intern = TRUE)
output_filename <- paste0('res_optims_',short_commit_hash,'.RDS')

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i, ] %>% unlist   
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ":", "Parametry:", row, '\n')
  
  
  result <- evaluate_function(row["optims"], row["dsets"], row["c"], row["n_resamp"])
  results[[i]] <- result  # Dodanie wyniku do listy
  
  results_df <- bind_rows(results)
  
  res <- cbind.data.frame(param_grid[1:length(results),], results_df)
  
  saveRDS(object = res, file = output_filename)
}
  



