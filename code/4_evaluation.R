
#############################
# Evaluation Scripts
#############################

#Directory
xgboost_output_IQR_dir
setwd(xgboost_output_IQR_dir) #Change to whichever directory has the "DARVIC_model" or your model of interest
list.files()

#Dataframe setup

xgboost_IQR_df = bind_rows(IQR_pathogenic_df,
                       IQR_benign_df)
xgboost_IQR_df_filtered = xgboost_IQR_df %>%
  dplyr::select(all_of(high_agreement_angles$Angle) | class | key)

final_xgboost_IQR_df_testing = xgboost_IQR_df_filtered %>%
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

testing_set = final_xgboost_IQR_df_testing
y_test = as.integer(as.factor(testing_set$class)) - 1
X_test = testing_set %>% dplyr::select(-class,-key)

#Model 1
model = xgb.load("DARVIC_model") #Load the name of the model desired to be evaluated

model_pred_xgb = predict(model, as.matrix(X_test),
                   type = "response")
model_ROCpred.xgb <- prediction(as.numeric(model_pred_xgb),
                          as.numeric(y_test))
model_auc.xgb <- performance(model_ROCpred.xgb, measure = "auc")
model_auc = model_auc.xgb@y.values[[1]]
names(model_auc) = c("XGBoost AUC")
model_auc
model_ROCperf_xgb = performance(model_ROCpred.xgb, "tpr", "fpr")
model_roc = data.frame(FalsePositive = c(model_ROCperf_xgb@x.values[[1]]),
                        TruePositive = c(model_ROCperf_xgb@y.values[[1]]))
model_roc
length(model_ROCperf_xgb@alpha.values[[1]])
alphavalue = model_ROCperf_xgb@alpha.values[[1]][13] #
alphavalue

xgb.plot.tree(model = model)

#Making a dataframe with the key and predictions
#With experimental testing set
model_predictions_table = as_tibble(model_pred_xgb) %>%
  rename(prediction_value = value)
model_predictions_table

prediction_testing_set_df = cbind((testing_set %>% dplyr::select(class, key)),
                      model_predictions_table) %>%
 mutate(classification = case_when(prediction_value >= alphavalue ~ "Deleterious",
                       prediction_value < alphavalue  ~ "Tolerated")) %>%
  arrange(desc(prediction_value))
prediction_testing_set_df = prediction_testing_set_df %>%
  arrange(class)

model_performance = confusionMatrix(data = as.factor(
  prediction_testing_set_df$classification),
  reference = as.factor(prediction_testing_set_df$class),
  positive = "Deleterious", mode = "everything")
model_performance # This would tell you the performance metrics of the model given the dataset

prediction_testing_set_df


#Importance matrix

importance_matrix = xgb.importance(colnames(testing_set %>%
                                              dplyr::select(-class,-key)), model = model)
gg = xgb.ggplot.importance(importance_matrix, rel_to_first = FALSE, xlab = "importance")
gg + ggplot2::ylab("Frequency")

xgb.plot.tree(model = model, feature_names = colnames(X_test))


#############################
# VUS Prediction
#############################
vus_df = IQR_vus_df %>% 
  arrange(key) %>%
  dplyr::select(all_of(high_agreement_angles$Angle)) 

model_vus_raw_preds = predict(model, as.matrix(vus_df), reshape = TRUE)
model_vus_raw_preds_df = as_tibble(model_vus_raw_preds)

##Vus prediction
model_vus_preds_withcutoff = as.numeric(model_vus_raw_preds > alphavalue )

key_names = IQR_vus_df %>%
  arrange(key) %>%
  dplyr::select(key)

model_vus_preds_withcutoff_df = as_tibble(model_vus_preds_withcutoff) %>%
  dplyr::rename(prediction = value) %>%
  mutate(classification = case_when(prediction == 1 ~ "Deleterious",
                       prediction == 0 ~ "Tolerated")) %>%
  mutate(raw_probability = model_vus_raw_preds_df$value) %>%
  mutate(key = key_names$key)

model_vus_preds_withcutoff_df

model_vus_plot = ggplot() +
  geom_bar(data = model_vus_preds_withcutoff_df, aes(classification, fill = as.factor(classification))) + theme_cowplot() +
  ylim(0,250)
model_vus_plot 

##exporting
write.csv(model_vus_preds_withcutoff_df, paste(raw_data_dir, "vus_predictions.csv", sep = "/"),
          row.names = FALSE)
write.csv(prediction_testing_set_df, paste(raw_data_dir, "testing_predictions.csv", sep = "/"),
          row.names = FALSE)
