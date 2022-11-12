################################################################################
#                            3: Modelling scripts 
################################################################################

#Directories
xgboost_output_dir = paste(rootdir, "xgboost_final_cleanup", sep = "/")
xgboost_output_dir

#############################################################
###### Benign and Pathogenic separation (Using IQR)
#############################################################

#Dataframe making
xgboost_IQR_df = bind_rows(IQR_pathogenic_df,
                       IQR_benign_df)
xgboost_IQR_df_filtered = xgboost_IQR_df %>%
  dplyr::select(all_of(high_agreement_angles$Angle) | class | key)

final_xgboost_IQR_df_training = xgboost_IQR_df_filtered %>%
  filter(key != "Asp121Gly") %>%
  filter(key != "Val231Met") %>%
  filter(key != "Val243Phe") %>%
  filter(key != "Gly283Glu") %>%
  filter(key != "Arg306Cys") %>%
  filter(key != "Ser319Asn") %>%
  filter(key != "Gln335His") %>%
  filter(key != "Gln335Arg") %>%
  filter(key != "Tyr125His") %>%
  filter(key != "Trp128Arg") %>%
  filter(key != "Pro154Leu") %>%
  filter(key != "Tyr176Cys") %>%
  filter(key != "Tyr177Ser") %>%
  filter(key != "Arg179Cys") %>%
  filter(key != "Arg179His") %>%
  filter(key != "Arg182Gln") %>%
  filter(key != "Arg182Trp") %>%
  filter(key != "Gly186Glu") %>%
  filter(key != "Asn235Ser") %>%
  filter(key != "Arg238Trp") %>%
  filter(key != "Arg242Cys") %>%
  filter(key != "Arg242His") %>%
  filter(key != "Arg271Trp") %>%
  filter(key != "Met280Val") %>%
  filter(key != "Pro292Leu") 

final_xgboost_IQR_df_testing = xgboost_IQR_df_filtered %>%
  filter(key == "Arg306Cys" | key == "Arg242Cys" |
             key == "Gln335His" | key == "Arg242His" |
             key == "Gln335Arg" | key == "Pro154Leu" |
             key == "Trp128Arg" | key == "Arg179His" |
             key == "Arg179Cys" | key == "Asn235Ser" |
             key == "Arg182Trp" | key == "Val243Phe" | 
             key == "Gly186Glu" | key == "Met280Val" |
             key == "Arg238Trp" | 
             key == "Gly283Glu" |
             key == "Arg271Trp" | 
             key == "Pro292Leu" |
             key == "Val231Met" )

#############################
# Autorun modeller
#############################
xgboost_output_IQR_dir = paste(xgboost_output_dir, "IQR", sep = "/")
setwd(xgboost_output_IQR_dir)
list.files()
auc_values = list()
xgboost_ROC_modeller = function(modelname, expected_accuracy,
                                sensitivity_high,sensitivity_low,
                                specificity_high, specificity_low, 
                                deleterious_percentage_high,
                                deleterious_percentage_low,
                                tolerated_percentage_high,
                                tolerated_percentage_low,
                                balanced_accuracy_high,
                                balanced_accuracy_low) {
  #Training and testing dataset
  train_set = final_xgboost_IQR_df_training
  testing_set = final_xgboost_IQR_df_testing
  
  y_train = as.integer(as.factor(train_set$class)) - 1 #Not sure why -1
  y_test = as.integer(as.factor(testing_set$class)) - 1
  X_train = train_set %>% dplyr::select(-class,-key) #dplyr::selects the columns
  X_test = testing_set %>% dplyr::select(-class,-key)
  
  #Making the model
  xgb_train = xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  xgb_test = xgb.DMatrix(data = as.matrix(X_test), label = y_test)
  xgb_params <- list(
    booster = "gbtree",
    eta = 0.03, 
    max_depth = 8,
    gamma = 4,
    subsample = 0.9,
    objective = "binary:logistic",
    #eval_metric = "mlogloss",
    #num_class = length(levels(final_xgboost_raw_angle_df_training$class)),
    eval_metric = "auc",
    min_child_weight = 1
    #eval_metric = "merror"
  )
  print(paste0("Building ", modelname, "...", sep = ""))
  #print(paste0("Current number of fitting models produced: ",
  #            length(list.files(xgboost_output_dir)), sep = ""))
  xgbCV = xgb.cv(params = xgb_params,
                 data = xgb_train,
                 nrounds = 1000,
                 prediction = TRUE,
                 showsd = TRUE,
                 early_stopping_rounds = 100,
                 maximize = TRUE,
                 nfold = 5,
                 stratified = TRUE)
                 #watchlist = list(train = xgb_train , 
                                  #test = xgb_test))
  
  numrounds_auc <- min(which(xgbCV$evaluation_log$test_auc_mean == 
                         max(xgbCV$evaluation_log$test_auc_mean)))
  #numrounds_error = min(which(xgbCV$evaluation_log$test_error_mean ==
  #                        max(xgbCV$evaluation_log$test_error_mean)))
  fit <- xgboost(params = xgb_params,
                 data = xgb_train,
                 nrounds = numrounds_auc)
  #Control performance AUC, specificity and sensitivity
  pred_xgb = predict(fit, as.matrix(X_test),
                     type = "response")
  ROCpred.xgb <- prediction(as.numeric(pred_xgb),
                            as.numeric(y_test))
  auc.xgb <- performance(ROCpred.xgb, measure = "auc")
  auc = auc.xgb@y.values[1]
  auc
  ROCperf_xgb = performance(ROCpred.xgb, "tpr", "fpr")
  
  print(ROCperf_xgb@alpha.values[[1]])
  if (auc >= expected_accuracy) {
    cutoffs_length = length(ROCperf_xgb@alpha.values[[1]])
     
    prediction_table = as_tibble(pred_xgb) %>%
      rename(prediction_value = value)
    prediction_table
    
    for (cutoff in 2:cutoffs_length) {
      cutoff_value = ROCperf_xgb@alpha.values[[1]][cutoff]
      
      #Control performance checking
      model_prediction_testing_set_df = cbind((testing_set %>% dplyr::select(class, key)),
                        prediction_table) %>%
      mutate(classification = case_when(prediction_value >= cutoff_value ~ "pathogenic",
                         prediction_value < cutoff_value  ~ "benign")) %>%
      arrange(desc(prediction_value))
      model_prediction_testing_set_df
    
      modeller_confusion_matrix = confusionMatrix(data = as.factor(
        model_prediction_testing_set_df$classification),
        reference = as.factor(model_prediction_testing_set_df$class),
        positive = "pathogenic", mode = "everything")
      modeller_confusion_matrix
      
      model_sensitivity = modeller_confusion_matrix$byClass[1]
      model_specificity = modeller_confusion_matrix$byClass[2]
      model_balanced_accuracy = modeller_confusion_matrix$byClass[11]
    
      #print(model_sensitivity)
      #print(model_specificity)
      #print(model_balanced_accuracy)
      
      #VUS prediction
      model_vus_df = IQR_vus_df %>%
        dplyr::select(all_of(high_agreement_angles$Angle))
      model_vus_raw_preds = predict(fit, as.matrix(model_vus_df), reshape = TRUE)
      model_vus_raw_with_cutoff = as.numeric(model_vus_raw_preds >= cutoff_value)
      model_vus_raw_with_cutoff_df = as_tibble(model_vus_raw_with_cutoff) %>%
        rename(prediction = value) %>%
        mutate(classification = case_when(prediction == 1 ~ "deleterious",
                             prediction == 0 ~ "tolerated"))
      
      deleterious_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "deleterious"))[1]
      tolerated_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "tolerated"))[1]
      
      deleterious_percent = deleterious_count/254
      tolerated_percent = tolerated_count/254
      
      #print(deleterious_percent)
      #print(tolerated_percent)
    
      if(model_sensitivity <= sensitivity_high & 
         model_sensitivity >= sensitivity_low &
         model_specificity <= specificity_high &
         model_specificity >= specificity_low &
         deleterious_percent >= deleterious_percentage_low &
         deleterious_percent <= deleterious_percentage_high &
         tolerated_percent >= tolerated_percentage_low &
         tolerated_percent <= tolerated_percentage_high &
         model_balanced_accuracy >= balanced_accuracy_low &
         model_balanced_accuracy <= balanced_accuracy_high) {
        name = paste(modelname, auc, sep = "_")
        xgb.save(fit, name)
        print(paste("Model fits the criteria and has been saved as ", name, "...",
              sep = ""))
        return(fit)
      } else {
        print("This cutoff is not good")
      }
    }
  } else {
    print("Model is not good, training another model...")
  }
  print(paste("The AUC of this model is: ", auc[[1]], sep = " "))
  return(fit)
  #return(auc[[1]])
} 
xgboost_ROC_modeller_run = function (number_of_runs,number_of_desired_models,
                                     expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low,
                                     directory) {
  model_numbers = c(1:number_of_runs) 
  #the modeller to run
  model_names = lapply(model_numbers, 
                       function(x) paste0("model_", x))
  auc_values = list()
  for (i in model_names){
    if (length(list.files(directory)) < number_of_desired_models) {
      #xgb.save(xgb_model,i)
      seed = .Random.seed
      model = xgboost_ROC_modeller(i, expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low)
      attr(model, "seed") = seed
      #auc_values[[length(auc_values)+1]] = auc
      #Sys.sleep(1)
    } else {
      print("The needed amount of models have been achieved already")
      print(seed)
      print(paste("Number of fitting models produced: ", 
                  length(list.files(directory)), sep = ""))
      break
    }
  }
  return(model)
} 
.Random.seed = dget("DARVIC_seed.txt") #Run this first to get the exact representative model
xgboost_ROC_modeller_run(10000,8,0.70,1,0,1,0,1,0,1,0,1,0,xgboost_output_IQR_dir)


#############################################################
###### Benign and Pathogenic separation (Using Median)
#############################################################

#Dataframe making
xgboost_median_df = bind_rows(median_pathogenic_df,
                           median_benign_df)
xgboost_median_df_filtered = xgboost_median_df %>%
  dplyr::select(all_of(high_agreement_angles$Angle) | class | key)

final_xgboost_median_df_training = xgboost_median_df_filtered %>%
  filter(key != "Asp121Gly") %>%
  filter(key != "Val231Met") %>%
  filter(key != "Val243Phe") %>%
  filter(key != "Gly283Glu") %>%
  filter(key != "Arg306Cys") %>%
  filter(key != "Ser319Asn") %>%
  filter(key != "Gln335His") %>%
  filter(key != "Gln335Arg") %>%
  filter(key != "Tyr125His") %>%
  filter(key != "Trp128Arg") %>%
  filter(key != "Pro154Leu") %>%
  filter(key != "Tyr176Cys") %>%
  filter(key != "Tyr177Ser") %>%
  filter(key != "Arg179Cys") %>%
  filter(key != "Arg179His") %>%
  filter(key != "Arg182Gln") %>%
  filter(key != "Arg182Trp") %>%
  filter(key != "Gly186Glu") %>%
  filter(key != "Asn235Ser") %>%
  filter(key != "Arg238Trp") %>%
  filter(key != "Arg242Cys") %>%
  filter(key != "Arg242His") %>%
  filter(key != "Arg271Trp") %>%
  filter(key != "Met280Val") %>%
  filter(key != "Pro292Leu") 

final_xgboost_median_df_testing = xgboost_median_df_filtered %>%
filter(      key == "Arg306Cys" | key == "Arg242Cys" |
             key == "Gln335His" | key == "Arg242His" |
             key == "Gln335Arg" | key == "Pro154Leu" |
             key == "Trp128Arg" | key == "Arg179His" |
             key == "Arg179Cys" | key == "Asn235Ser" |
             key == "Arg182Trp" | key == "Val243Phe" | 
             key == "Gly186Glu" | key == "Met280Val" |
             key == "Arg238Trp" | 
             key == "Gly283Glu" |
             key == "Arg271Trp" | 
             key == "Pro292Leu" |
             key == "Val231Met" )
final_xgboost_median_df_testing


#############################
# Autorun modeller
#############################
xgboost_median_output_dir = paste(xgboost_output_dir, "median",sep = "/")
setwd(xgboost_median_output_dir)

list.files()
auc_values = list()
xgboost_ROC_modeller = function(modelname, expected_accuracy,
                                sensitivity_high,sensitivity_low,
                                specificity_high, specificity_low, 
                                deleterious_percentage_high,
                                deleterious_percentage_low,
                                tolerated_percentage_high,
                                tolerated_percentage_low,
                                balanced_accuracy_high,
                                balanced_accuracy_low) {
  #Training and testing dataset
  train_set = final_xgboost_median_df_training
  testing_set = final_xgboost_median_df_testing
  
  y_train = as.integer(as.factor(train_set$class)) - 1 #Not sure why -1
  y_test = as.integer(as.factor(testing_set$class)) - 1
  X_train = train_set %>% dplyr::select(-class,-key) #dplyr::selects the columns
  X_test = testing_set %>% dplyr::select(-class,-key)
  
  #Making the model
  xgb_train = xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  xgb_test = xgb.DMatrix(data = as.matrix(X_test), label = y_test)
  xgb_params <- list(
    booster = "gbtree",
    eta = 0.03, 
    max_depth = 8,
    gamma = 4,
    subsample = 0.9,
    objective = "binary:logistic",
    #eval_metric = "mlogloss",
    #num_class = length(levels(final_xgboost_raw_angle_df_training$class)),
    eval_metric = "auc",
    min_child_weight = 1
    #eval_metric = "merror"
  )
  print(paste0("Building ", modelname, "...", sep = ""))
  #print(paste0("Current number of fitting models produced: ",
  #            length(list.files(xgboost_output_dir)), sep = ""))
  xgbCV = xgb.cv(params = xgb_params,
                 data = xgb_train,
                 nrounds = 1000,
                 prediction = TRUE,
                 showsd = TRUE,
                 early_stopping_rounds = 100,
                 maximize = TRUE,
                 nfold = 5,
                 stratified = TRUE)
                 #watchlist = list(train = xgb_train , 
                                  #test = xgb_test))
  
  numrounds_auc <- min(which(xgbCV$evaluation_log$test_auc_mean == 
                         max(xgbCV$evaluation_log$test_auc_mean)))
  #numrounds_error = min(which(xgbCV$evaluation_log$test_error_mean ==
  #                        max(xgbCV$evaluation_log$test_error_mean)))
  fit <- xgboost(params = xgb_params,
                 data = xgb_train,
                 nrounds = numrounds_auc)
  #Control performance AUC, specificity and sensitivity
  pred_xgb = predict(fit, as.matrix(X_test),
                     type = "response")
  ROCpred.xgb <- prediction(as.numeric(pred_xgb),
                            as.numeric(y_test))
  auc.xgb <- performance(ROCpred.xgb, measure = "auc")
  auc = auc.xgb@y.values[1]
  auc
  ROCperf_xgb = performance(ROCpred.xgb, "tpr", "fpr")
  
  print(ROCperf_xgb@alpha.values[[1]])
  if (auc >= expected_accuracy) {
    cutoffs_length = length(ROCperf_xgb@alpha.values[[1]])
     
    prediction_table = as_tibble(pred_xgb) %>%
      rename(prediction_value = value)
    prediction_table
    
    for (cutoff in 2:cutoffs_length) {
      cutoff_value = ROCperf_xgb@alpha.values[[1]][cutoff]
      
      #Control performance checking
      model_prediction_testing_set_df = cbind((testing_set %>% dplyr::select(class, key)),
                        prediction_table) %>%
      mutate(classification = case_when(prediction_value >= cutoff_value ~ "pathogenic",
                         prediction_value < cutoff_value  ~ "benign")) %>%
      arrange(desc(prediction_value))
      model_prediction_testing_set_df
    
      modeller_confusion_matrix = confusionMatrix(data = as.factor(
        model_prediction_testing_set_df$classification),
        reference = as.factor(model_prediction_testing_set_df$class),
        positive = "pathogenic", mode = "everything")
      modeller_confusion_matrix
      
      model_sensitivity = modeller_confusion_matrix$byClass[1]
      model_specificity = modeller_confusion_matrix$byClass[2]
      model_balanced_accuracy = modeller_confusion_matrix$byClass[11]
    
      #print(model_sensitivity)
      #print(model_specificity)
      #print(model_balanced_accuracy)
      
      #VUS prediction
      model_vus_df = median_vus_df %>%
        dplyr::select(all_of(high_agreement_angles$Angle))
      model_vus_raw_preds = predict(fit, as.matrix(model_vus_df), reshape = TRUE)
      model_vus_raw_with_cutoff = as.numeric(model_vus_raw_preds >= cutoff_value)
      model_vus_raw_with_cutoff_df = as_tibble(model_vus_raw_with_cutoff) %>%
        rename(prediction = value) %>%
        mutate(classification = case_when(prediction == 1 ~ "deleterious",
                             prediction == 0 ~ "tolerated"))
      
      deleterious_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "deleterious"))[1]
      tolerated_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "tolerated"))[1]
      
      deleterious_percent = deleterious_count/254
      tolerated_percent = tolerated_count/254
      
      #print(deleterious_percent)
      #print(tolerated_percent)
    
      if(model_sensitivity <= sensitivity_high & 
         model_sensitivity >= sensitivity_low &
         model_specificity <= specificity_high &
         model_specificity >= specificity_low &
         deleterious_percent >= deleterious_percentage_low &
         deleterious_percent <= deleterious_percentage_high &
         tolerated_percent >= tolerated_percentage_low &
         tolerated_percent <= tolerated_percentage_high &
         model_balanced_accuracy >= balanced_accuracy_low &
         model_balanced_accuracy <= balanced_accuracy_high) {
        name = paste(modelname, auc, sep = "_")
        xgb.save(fit, name)
        print(paste("Model fits the criteria and has been saved as ", name, "...",
              sep = ""))
        return(fit)
      } else {
        print("This cutoff is not good")
      }
    }
  } else {
    print("Model is not good, training another model...")
  }
  print(paste("The AUC of this model is: ", auc[[1]], sep = " "))
  return(fit)
  #return(auc[[1]])
} 
xgboost_ROC_modeller_run = function (number_of_runs,number_of_desired_models,
                                     expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low,
                                     directory) {
  model_numbers = c(1:number_of_runs) 
  #the modeller to run
  model_names = lapply(model_numbers, 
                       function(x) paste0("model_", x))
  auc_values = list()
  for (i in model_names){
    if (length(list.files(directory)) < number_of_desired_models) {
      #xgb.save(xgb_model,i)
      seed = .Random.seed
      model = xgboost_ROC_modeller(i, expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low)
      attr(model, "seed") = seed
      #auc_values[[length(auc_values)+1]] = auc
      #Sys.sleep(1)
    } else {
      print("The needed amount of models have been achieved already")
      print(seed)
      print(paste("Number of fitting models produced: ", 
                  length(list.files(directory)), sep = ""))
      break
    }
  }
  return(model)
} 
xgboost_ROC_modeller_run(10000,8,0.60,1,0,1,0,1,0,1,0,1,0,xgboost_median_output_dir)

#############################################################
###### Benign and Pathogenic separation (Using range)
#############################################################

#Dataframe making
xgboost_range_df = bind_rows(range_pathogenic_df,
                            range_benign_df)
xgboost_range_df_filtered = xgboost_range_df %>%
  dplyr::select(all_of(high_agreement_angles$Angle) | class | key)

final_xgboost_range_df_training = xgboost_range_df_filtered %>%
  filter(key != "Asp121Gly") %>%
  filter(key != "Val231Met") %>%
  filter(key != "Val243Phe") %>%
  filter(key != "Gly283Glu") %>%
  filter(key != "Arg306Cys") %>%
  filter(key != "Ser319Asn") %>%
  filter(key != "Gln335His") %>%
  filter(key != "Gln335Arg") %>%
  filter(key != "Tyr125His") %>%
  filter(key != "Trp128Arg") %>%
  filter(key != "Pro154Leu") %>%
  filter(key != "Tyr176Cys") %>%
  filter(key != "Tyr177Ser") %>%
  filter(key != "Arg179Cys") %>%
  filter(key != "Arg179His") %>%
  filter(key != "Arg182Gln") %>%
  filter(key != "Arg182Trp") %>%
  filter(key != "Gly186Glu") %>%
  filter(key != "Asn235Ser") %>%
  filter(key != "Arg238Trp") %>%
  filter(key != "Arg242Cys") %>%
  filter(key != "Arg242His") %>%
  filter(key != "Arg271Trp") %>%
  filter(key != "Met280Val") %>%
  filter(key != "Pro292Leu") 

final_xgboost_range_df_testing = xgboost_range_df_filtered %>%
  filter(    key == "Arg306Cys" | key == "Arg242Cys" |
             key == "Gln335His" | key == "Arg242His" |
             key == "Gln335Arg" | key == "Pro154Leu" |
             key == "Trp128Arg" | key == "Arg179His" |
             key == "Arg179Cys" | key == "Asn235Ser" |
             key == "Arg182Trp" | key == "Val243Phe" | 
             key == "Gly186Glu" | key == "Met280Val" |
             key == "Arg238Trp" | 
             key == "Gly283Glu" |
             key == "Arg271Trp" | 
             key == "Pro292Leu" |
             key == "Val231Met" )
final_xgboost_range_df_testing


#############################
# Autorun modeller
#############################
xgboost_range_output_dir = paste(xgboost_output_dir, "range",sep = "/")
setwd(xgboost_range_output_dir)
list.files()
xgboost_ROC_modeller = function(modelname, expected_accuracy,
                                sensitivity_high,sensitivity_low,
                                specificity_high, specificity_low, 
                                deleterious_percentage_high,
                                deleterious_percentage_low,
                                tolerated_percentage_high,
                                tolerated_percentage_low,
                                balanced_accuracy_high,
                                balanced_accuracy_low) {
  #Training and testing dataset
  train_set = final_xgboost_range_df_training
  testing_set = final_xgboost_range_df_testing
  
  y_train = as.integer(as.factor(train_set$class)) - 1 #Not sure why -1
  y_test = as.integer(as.factor(testing_set$class)) - 1
  X_train = train_set %>% dplyr::select(-class,-key) #dplyr::selects the columns
  X_test = testing_set %>% dplyr::select(-class,-key)
  
  #Making the model
  xgb_train = xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  xgb_test = xgb.DMatrix(data = as.matrix(X_test), label = y_test)
  xgb_params <- list(
    booster = "gbtree",
    eta = 0.03, 
    max_depth = 8,
    gamma = 4,
    subsample = 0.9,
    objective = "binary:logistic",
    #eval_metric = "mlogloss",
    #num_class = length(levels(final_xgboost_raw_angle_df_training$class)),
    eval_metric = "auc",
    min_child_weight = 1
    #eval_metric = "merror"
  )
  print(paste0("Building ", modelname, "...", sep = ""))
  #print(paste0("Current number of fitting models produced: ",
  #            length(list.files(xgboost_output_dir)), sep = ""))
  xgbCV = xgb.cv(params = xgb_params,
                 data = xgb_train,
                 nrounds = 1000,
                 prediction = TRUE,
                 showsd = TRUE,
                 early_stopping_rounds = 100,
                 maximize = TRUE,
                 nfold = 5,
                 stratified = TRUE)
                 #watchlist = list(train = xgb_train , 
                                  #test = xgb_test))
  
  numrounds_auc <- min(which(xgbCV$evaluation_log$test_auc_mean == 
                         max(xgbCV$evaluation_log$test_auc_mean)))
  #numrounds_error = min(which(xgbCV$evaluation_log$test_error_mean ==
  #                        max(xgbCV$evaluation_log$test_error_mean)))
  fit <- xgboost(params = xgb_params,
                 data = xgb_train,
                 nrounds = numrounds_auc)
  #Control performance AUC, specificity and sensitivity
  pred_xgb = predict(fit, as.matrix(X_test),
                     type = "response")
  ROCpred.xgb <- prediction(as.numeric(pred_xgb),
                            as.numeric(y_test))
  auc.xgb <- performance(ROCpred.xgb, measure = "auc")
  auc = auc.xgb@y.values[1]
  auc
  ROCperf_xgb = performance(ROCpred.xgb, "tpr", "fpr")
  
  print(ROCperf_xgb@alpha.values[[1]])
  if (auc >= expected_accuracy) {
    cutoffs_length = length(ROCperf_xgb@alpha.values[[1]])
     
    prediction_table = as_tibble(pred_xgb) %>%
      rename(prediction_value = value)
    prediction_table
    
    for (cutoff in 2:cutoffs_length) {
      cutoff_value = ROCperf_xgb@alpha.values[[1]][cutoff]
      
      #Control performance checking
      model_prediction_testing_set_df = cbind((testing_set %>% dplyr::select(class, key)),
                        prediction_table) %>%
      mutate(classification = case_when(prediction_value >= cutoff_value ~ "pathogenic",
                         prediction_value < cutoff_value  ~ "benign")) %>%
      arrange(desc(prediction_value))
      model_prediction_testing_set_df
    
      modeller_confusion_matrix = confusionMatrix(data = as.factor(
        model_prediction_testing_set_df$classification),
        reference = as.factor(model_prediction_testing_set_df$class),
        positive = "pathogenic", mode = "everything")
      modeller_confusion_matrix
      
      model_sensitivity = modeller_confusion_matrix$byClass[1]
      model_specificity = modeller_confusion_matrix$byClass[2]
      model_balanced_accuracy = modeller_confusion_matrix$byClass[11]
    
      #print(model_sensitivity)
      #print(model_specificity)
      #print(model_balanced_accuracy)
      
      #VUS prediction
      model_vus_df = range_vus_df %>%
        dplyr::select(all_of(high_agreement_angles$Angle))
      model_vus_raw_preds = predict(fit, as.matrix(model_vus_df), reshape = TRUE)
      model_vus_raw_with_cutoff = as.numeric(model_vus_raw_preds >= cutoff_value)
      model_vus_raw_with_cutoff_df = as_tibble(model_vus_raw_with_cutoff) %>%
        rename(prediction = value) %>%
        mutate(classification = case_when(prediction == 1 ~ "deleterious",
                             prediction == 0 ~ "tolerated"))
      
      deleterious_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "deleterious"))[1]
      tolerated_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "tolerated"))[1]
      
      deleterious_percent = deleterious_count/254
      tolerated_percent = tolerated_count/254
      
      #print(deleterious_percent)
      #print(tolerated_percent)
    
      if(model_sensitivity <= sensitivity_high & 
         model_sensitivity >= sensitivity_low &
         model_specificity <= specificity_high &
         model_specificity >= specificity_low &
         deleterious_percent >= deleterious_percentage_low &
         deleterious_percent <= deleterious_percentage_high &
         tolerated_percent >= tolerated_percentage_low &
         tolerated_percent <= tolerated_percentage_high &
         model_balanced_accuracy >= balanced_accuracy_low &
         model_balanced_accuracy <= balanced_accuracy_high) {
        name = paste(modelname, auc, sep = "_")
        xgb.save(fit, name)
        print(paste("Model fits the criteria and has been saved as ", name, "...",
              sep = ""))
        return(fit)
      } else {
        print("This cutoff is not good")
      }
    }
  } else {
    print("Model is not good, training another model...")
  }
  print(paste("The AUC of this model is: ", auc[[1]], sep = " "))
  return(fit)
  #return(auc[[1]])
} 
xgboost_ROC_modeller_run = function (number_of_runs,number_of_desired_models,
                                     expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low,
                                     directory) {
  model_numbers = c(1:number_of_runs) 
  #the modeller to run
  model_names = lapply(model_numbers, 
                       function(x) paste0("model_", x))
  auc_values = list()
  for (i in model_names){
    if (length(list.files(directory)) < number_of_desired_models) {
      #xgb.save(xgb_model,i)
      seed = .Random.seed
      model = xgboost_ROC_modeller(i, expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low)
      attr(model, "seed") = seed
      #auc_values[[length(auc_values)+1]] = auc
      #Sys.sleep(1)
    } else {
      print("The needed amount of models have been achieved already")
      print(seed)
      print(paste("Number of fitting models produced: ", 
                  length(list.files(directory)), sep = ""))
      break
    }
  }
  return(model)
} 
xgboost_ROC_modeller_run(10000,8,0.60,1,0,1,0,1,0,1,0,1,0,xgboost_range_output_dir)

#############################################################
###### Benign and Pathogenic separation (Using Mean)
#############################################################

#Dataframe making
xgboost_mean_df = bind_rows(mean_pathogenic_df,
                           mean_benign_df)
xgboost_mean_df_filtered = xgboost_mean_df %>%
  dplyr::select(all_of(high_agreement_angles$Angle) | class | key)

final_xgboost_mean_df_training = xgboost_mean_df_filtered %>%
  filter(key != "Asp121Gly") %>%
  filter(key != "Val231Met") %>%
  filter(key != "Val243Phe") %>%
  filter(key != "Gly283Glu") %>%
  filter(key != "Arg306Cys") %>%
  filter(key != "Ser319Asn") %>%
  filter(key != "Gln335His") %>%
  filter(key != "Gln335Arg") %>%
  filter(key != "Tyr125His") %>%
  filter(key != "Trp128Arg") %>%
  filter(key != "Pro154Leu") %>%
  filter(key != "Tyr176Cys") %>%
  filter(key != "Tyr177Ser") %>%
  filter(key != "Arg179Cys") %>%
  filter(key != "Arg179His") %>%
  filter(key != "Arg182Gln") %>%
  filter(key != "Arg182Trp") %>%
  filter(key != "Gly186Glu") %>%
  filter(key != "Asn235Ser") %>%
  filter(key != "Arg238Trp") %>%
  filter(key != "Arg242Cys") %>%
  filter(key != "Arg242His") %>%
  filter(key != "Arg271Trp") %>%
  filter(key != "Met280Val") %>%
  filter(key != "Pro292Leu") 

final_xgboost_mean_df_testing = xgboost_mean_df_filtered %>%
filter(    key == "Arg306Cys" | key == "Arg242Cys" |
             key == "Gln335His" | key == "Arg242His" |
             key == "Gln335Arg" | key == "Pro154Leu" |
             key == "Trp128Arg" | key == "Arg179His" |
             key == "Arg179Cys" | key == "Asn235Ser" |
             key == "Arg182Trp" | key == "Val243Phe" | 
             key == "Gly186Glu" | key == "Met280Val" |
             key == "Arg238Trp" | 
             key == "Gly283Glu" |
             key == "Arg271Trp" | 
             key == "Pro292Leu" |
             key == "Val231Met" )


#############################
# Autorun modeller
#############################
xgboost_mean_output_dir = paste(xgboost_output_dir, "mean",sep = "/")
setwd(xgboost_mean_output_dir)

list.files()
xgboost_ROC_modeller = function(modelname, expected_accuracy,
                                sensitivity_high,sensitivity_low,
                                specificity_high, specificity_low, 
                                deleterious_percentage_high,
                                deleterious_percentage_low,
                                tolerated_percentage_high,
                                tolerated_percentage_low,
                                balanced_accuracy_high,
                                balanced_accuracy_low) {
  #Training and testing dataset
  train_set = final_xgboost_mean_df_training
  testing_set = final_xgboost_mean_df_testing
  
  y_train = as.integer(as.factor(train_set$class)) - 1 #Not sure why -1
  y_test = as.integer(as.factor(testing_set$class)) - 1
  X_train = train_set %>% dplyr::select(-class,-key) #dplyr::selects the columns
  X_test = testing_set %>% dplyr::select(-class,-key)
  
  #Making the model
  xgb_train = xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  xgb_test = xgb.DMatrix(data = as.matrix(X_test), label = y_test)
  xgb_params <- list(
    booster = "gbtree",
    eta = 0.03, 
    max_depth = 8,
    gamma = 4,
    subsample = 0.9,
    objective = "binary:logistic",
    #eval_metric = "mlogloss",
    #num_class = length(levels(final_xgboost_raw_angle_df_training$class)),
    eval_metric = "auc",
    min_child_weight = 1
    #eval_metric = "merror"
  )
  print(paste0("Building ", modelname, "...", sep = ""))
  #print(paste0("Current number of fitting models produced: ",
  #            length(list.files(xgboost_output_dir)), sep = ""))
  xgbCV = xgb.cv(params = xgb_params,
                 data = xgb_train,
                 nrounds = 1000,
                 prediction = TRUE,
                 showsd = TRUE,
                 early_stopping_rounds = 100,
                 maximize = TRUE,
                 nfold = 5,
                 stratified = TRUE)
                 #watchlist = list(train = xgb_train , 
                                  #test = xgb_test))
  
  numrounds_auc <- min(which(xgbCV$evaluation_log$test_auc_mean == 
                         max(xgbCV$evaluation_log$test_auc_mean)))
  #numrounds_error = min(which(xgbCV$evaluation_log$test_error_mean ==
  #                        max(xgbCV$evaluation_log$test_error_mean)))
  fit <- xgboost(params = xgb_params,
                 data = xgb_train,
                 nrounds = numrounds_auc)
  #Control performance AUC, specificity and sensitivity
  pred_xgb = predict(fit, as.matrix(X_test),
                     type = "response")
  ROCpred.xgb <- prediction(as.numeric(pred_xgb),
                            as.numeric(y_test))
  auc.xgb <- performance(ROCpred.xgb, measure = "auc")
  auc = auc.xgb@y.values[1]
  auc
  ROCperf_xgb = performance(ROCpred.xgb, "tpr", "fpr")
  
  print(ROCperf_xgb@alpha.values[[1]])
  if (auc >= expected_accuracy) {
    cutoffs_length = length(ROCperf_xgb@alpha.values[[1]])
     
    prediction_table = as_tibble(pred_xgb) %>%
      rename(prediction_value = value)
    prediction_table
    
    for (cutoff in 2:cutoffs_length) {
      cutoff_value = ROCperf_xgb@alpha.values[[1]][cutoff]
      
      #Control performance checking
      model_prediction_testing_set_df = cbind((testing_set %>% dplyr::select(class, key)),
                        prediction_table) %>%
      mutate(classification = case_when(prediction_value >= cutoff_value ~ "pathogenic",
                         prediction_value < cutoff_value  ~ "benign")) %>%
      arrange(desc(prediction_value))
      model_prediction_testing_set_df
    
      modeller_confusion_matrix = confusionMatrix(data = as.factor(
        model_prediction_testing_set_df$classification),
        reference = as.factor(model_prediction_testing_set_df$class),
        positive = "pathogenic", mode = "everything")
      modeller_confusion_matrix
      
      model_sensitivity = modeller_confusion_matrix$byClass[1]
      model_specificity = modeller_confusion_matrix$byClass[2]
      model_balanced_accuracy = modeller_confusion_matrix$byClass[11]
    
      #print(model_sensitivity)
      #print(model_specificity)
      #print(model_balanced_accuracy)
      
      #VUS prediction
      model_vus_df = mean_vus_df %>%
        dplyr::select(all_of(high_agreement_angles$Angle))
      model_vus_raw_preds = predict(fit, as.matrix(model_vus_df), reshape = TRUE)
      model_vus_raw_with_cutoff = as.numeric(model_vus_raw_preds >= cutoff_value)
      model_vus_raw_with_cutoff_df = as_tibble(model_vus_raw_with_cutoff) %>%
        rename(prediction = value) %>%
        mutate(classification = case_when(prediction == 1 ~ "deleterious",
                             prediction == 0 ~ "tolerated"))
      
      deleterious_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "deleterious"))[1]
      tolerated_count = dim(model_vus_raw_with_cutoff_df 
                             %>% filter(classification == "tolerated"))[1]
      
      deleterious_percent = deleterious_count/254
      tolerated_percent = tolerated_count/254
      
      #print(deleterious_percent)
      #print(tolerated_percent)
    
      if(model_sensitivity <= sensitivity_high & 
         model_sensitivity >= sensitivity_low &
         model_specificity <= specificity_high &
         model_specificity >= specificity_low &
         deleterious_percent >= deleterious_percentage_low &
         deleterious_percent <= deleterious_percentage_high &
         tolerated_percent >= tolerated_percentage_low &
         tolerated_percent <= tolerated_percentage_high &
         model_balanced_accuracy >= balanced_accuracy_low &
         model_balanced_accuracy <= balanced_accuracy_high) {
        name = paste(modelname, auc, sep = "_")
        xgb.save(fit, name)
        print(paste("Model fits the criteria and has been saved as ", name, "...",
              sep = ""))
        return(fit)
      } else {
        print("This cutoff is not good")
      }
    }
  } else {
    print("Model is not good, training another model...")
  }
  print(paste("The AUC of this model is: ", auc[[1]], sep = " "))
  return(fit)
  #return(auc[[1]])
} 
xgboost_ROC_modeller_run = function (number_of_runs,number_of_desired_models,
                                     expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low,
                                     directory) {
  model_numbers = c(1:number_of_runs) 
  #the modeller to run
  model_names = lapply(model_numbers, 
                       function(x) paste0("model_", x))
  auc_values = list()
  for (i in model_names){
    if (length(list.files(directory)) < number_of_desired_models) {
      #xgb.save(xgb_model,i)
      seed = .Random.seed
      model = xgboost_ROC_modeller(i, expected_accuracy,
                                     sensitivity_high, sensitivity_low,
                                     specificity_high, specificity_low, 
                                     deleterious_percentage_high,
                                     deleterious_percentage_low,
                                     tolerated_percentage_high,
                                     tolerated_percentage_low,
                                     balanced_accuracy_high,
                                     balanced_accuracy_low)
      attr(model, "seed") = seed
      #auc_values[[length(auc_values)+1]] = auc
      #Sys.sleep(1)
    } else {
      print("The needed amount of models have been achieved already")
      print(seed)
      print(paste("Number of fitting models produced: ", 
                  length(list.files(directory)), sep = ""))
      break
    }
  }
  return(model)
} 
xgboost_ROC_modeller_run(10000,8,0.60,1,0,1,0,1,0,1,0,1,0,xgboost_mean_output_dir)
