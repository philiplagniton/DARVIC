################################################################################
#               2: Data manipulation and cleaning
################################################################################

#################################### 
# Calculating mean for each residue
####################################


mean_byresidue = function(dataframe,variant,
                         position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "mean...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  mean = mean(wt_df3)
  return(mean)
}


#Function for loop for processing pvalues
mean_calculator = function(dataframe, refgroup,
                          angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_mean = list()
  for (testvar in list) {
    output_mean = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_mean[[position]] = mean_byresidue(
        dataframe,testvar,position,angle)
      #pvalue = pvalue_byresidue(refdf, refvar, 
      #                testdf, testvar, position, angle)
      #print(pvalue)
    }
    #print(output_mean)
    pervariant_mean[[testvar]] = output_mean
  }
  return(pervariant_mean)
}

mean_benign_phi = mean_calculator(benign_by_residue, "benign", "phi")
mean_benign_psi = mean_calculator(benign_by_residue, "benign", "psi")
mean_pathogenic_phi = mean_calculator(pathogenic_by_residue, "pathogenic", "phi")
mean_pathogenic_psi = mean_calculator(pathogenic_by_residue, "pathogenic", "psi")
mean_wtdiffsim_phi = mean_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
mean_wtdiffsim_psi = mean_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")
mean_vus_phi = mean_calculator(vus_by_residue, "vus","phi")
mean_vus_psi = mean_calculator(vus_by_residue, "vus","psi")

#Matrix maker function

aminoacid_residues
aminoacid_residues_phi = lapply(aminoacid_residues, 
                                function(x) paste0(x, "_phi"))
aminoacid_residues_psi = lapply(aminoacid_residues, 
                                function(x) paste0(x, "_psi"))

matrixmaker = function(dataframe,variant,angle){
  data1 = dataframe[[variant]]
  data2 = as.numeric(data1[80:353])
  aalength = length(aminoacid_residues)
  aanames = 
    if (angle == "phi") {
      aminoacid_residues_phi  
    } else if (angle == "psi") {
      aminoacid_residues_psi
    }
  dim(data2) = c(1,aalength)
  dimnames(data2) = list(variant, aanames)
  return(data2)
}

#Dataframe maker
mean_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }
  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = group) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

mean_benign_df = bind_rows(mean_dataframe_maker(
  mean_benign_phi,mean_benign_psi,"benign"))
mean_pathogenic_df = bind_rows(mean_dataframe_maker(
  mean_pathogenic_phi,mean_pathogenic_psi,"pathogenic"))
mean_wtdiffsim_df = bind_rows(mean_dataframe_maker(
  mean_wtdiffsim_phi,mean_wtdiffsim_psi,"wtdiffsim"))
mean_vus_df = bind_rows(mean_dataframe_maker(
  mean_vus_phi,mean_vus_psi,"vus"))


####################################
# Calculating IQR for each residue
####################################

IQR_byresidue = function(dataframe,variant,
                           position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "IQR...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  IQR = IQR(wt_df3)
  return(IQR)
}

#Function for loop for processing pvalues
IQR_calculator = function(dataframe, refgroup,
                            angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_IQR = list()
  for (testvar in list) {
    output_IQR = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_IQR[[position]] = IQR_byresidue(
        dataframe,testvar,position,angle)
    }
    pervariant_IQR[[testvar]] = output_IQR
  }
  return(pervariant_IQR)
}

IQR_benign_phi = IQR_calculator(benign_by_residue, "benign", "phi")
IQR_benign_psi = IQR_calculator(benign_by_residue, "benign", "psi")
IQR_pathogenic_phi = IQR_calculator(pathogenic_by_residue, "pathogenic", "phi")
IQR_pathogenic_psi = IQR_calculator(pathogenic_by_residue, "pathogenic", "psi")
IQR_wtdiffsim_phi = IQR_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
IQR_wtdiffsim_psi = IQR_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")
IQR_vus_phi = IQR_calculator(vus_by_residue, "vus","phi")
IQR_vus_psi = IQR_calculator(vus_by_residue, "vus","psi")

#Dataframe maker
IQR_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }

  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = group) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

IQR_benign_df = bind_rows(IQR_dataframe_maker(
  IQR_benign_phi,IQR_benign_psi,"benign"))
IQR_pathogenic_df = bind_rows(IQR_dataframe_maker(
  IQR_pathogenic_phi,IQR_pathogenic_psi,"pathogenic"))
IQR_wtdiffsim_df = bind_rows(IQR_dataframe_maker(
  IQR_wtdiffsim_phi,IQR_wtdiffsim_psi,"wtdiffsim"))
IQR_vus_df = bind_rows(IQR_dataframe_maker(
  IQR_vus_phi,IQR_vus_psi,"vus"))

####################################
# Calculating median for each residue
####################################


median_byresidue = function(dataframe,variant,
                            position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "median...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  median = median(wt_df3)
  return(median)
}

#Function for loop for processing pvalues
median_calculator = function(dataframe, refgroup,
                             angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_median = list()
  for (testvar in list) {
    output_median = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_median[[position]] = median_byresidue(
        dataframe,testvar,position,angle)

    }
    pervariant_median[[testvar]] = output_median
  }
  return(pervariant_median)
}

median_benign_phi = median_calculator(benign_by_residue, "benign", "phi")
median_benign_psi = median_calculator(benign_by_residue, "benign", "psi")
median_pathogenic_phi = median_calculator(pathogenic_by_residue, "pathogenic", "phi")
median_pathogenic_psi = median_calculator(pathogenic_by_residue, "pathogenic", "psi")
median_wtdiffsim_phi = median_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
median_wtdiffsim_psi = median_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")
median_vus_phi = median_calculator(vus_by_residue, "vus","phi")
median_vus_psi = median_calculator(vus_by_residue, "vus","psi")

#Dataframe maker
median_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }
  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = as.factor(group)) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

median_benign_df = bind_rows(median_dataframe_maker(
  median_benign_phi,median_benign_psi,"benign"))
median_pathogenic_df = bind_rows(median_dataframe_maker(
  median_pathogenic_phi,median_pathogenic_psi,"pathogenic"))
median_vus_df = bind_rows(median_dataframe_maker(
  median_vus_phi,median_vus_psi,"vus"))
median_wtdiffsim_df = bind_rows(median_dataframe_maker(
  median_wtdiffsim_phi,median_wtdiffsim_psi,"wtdiffsim"))

####################################
# Calculating range for each residue
####################################

range_byresidue = function(dataframe,variant,
                        position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "range...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  range = diff(range(wt_df3))
  return(range)
}

#Function for loop for processing pvalues
range_calculator = function(dataframe, refgroup,
                         angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_range = list()
  for (testvar in list) {
    output_range = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_range[[position]] = range_byresidue(
        dataframe,testvar,position,angle)
    }
    pervariant_range[[testvar]] = output_range
  }
  return(pervariant_range)
}

range_benign_phi = range_calculator(benign_by_residue, "benign", "phi")
range_benign_psi = range_calculator(benign_by_residue, "benign", "psi")
range_pathogenic_phi = range_calculator(pathogenic_by_residue, "pathogenic", "phi")
range_pathogenic_psi = range_calculator(pathogenic_by_residue, "pathogenic", "psi")
range_wtdiffsim_phi = range_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
range_wtdiffsim_psi = range_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")
range_vus_phi = range_calculator(vus_by_residue, "vus","phi")
range_vus_psi = range_calculator(vus_by_residue, "vus","psi")

#Dataframe maker
range_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }
  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = group) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

range_benign_df = bind_rows(range_dataframe_maker(
  range_benign_phi,range_benign_psi,"benign"))
range_pathogenic_df = bind_rows(range_dataframe_maker(
  range_pathogenic_phi,range_pathogenic_psi,"pathogenic"))
range_wtdiffsim_df = bind_rows(range_dataframe_maker(
  range_wtdiffsim_phi,range_wtdiffsim_psi,"wtdiffsim"))
range_vus_df = bind_rows(range_dataframe_maker(
  range_vus_phi,range_vus_psi,"vus"))

####################################
# Calculating sd for each residue
####################################

sd_byresidue = function(dataframe,variant,
                            position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "sd...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  sd = sd(wt_df3)
  return(sd)
}

#Function for loop for processing pvalues
sd_calculator = function(dataframe, refgroup,
                             angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_sd = list()
  for (testvar in list) {
    output_sd = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_sd[[position]] = sd_byresidue(
        dataframe,testvar,position,angle)
    }
    pervariant_sd[[testvar]] = output_sd
  }
  return(pervariant_sd)
}

sd_benign_phi = sd_calculator(benign_by_residue, "benign", "phi")
sd_benign_psi = sd_calculator(benign_by_residue, "benign", "psi")
sd_pathogenic_phi = sd_calculator(pathogenic_by_residue, "pathogenic", "phi")
sd_pathogenic_psi = sd_calculator(pathogenic_by_residue, "pathogenic", "psi")
sd_wtdiffsim_phi = sd_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
sd_wtdiffsim_psi = sd_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")

#Dataframe maker
sd_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }
  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = group) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

sd_benign_df = bind_rows(sd_dataframe_maker(
  sd_benign_phi,sd_benign_psi,"benign"))
sd_pathogenic_df = bind_rows(sd_dataframe_maker(
  sd_pathogenic_phi,sd_pathogenic_psi,"pathogenic"))
sd_wtdiffsim_df = bind_rows(sd_dataframe_maker(
  sd_wtdiffsim_phi,sd_wtdiffsim_psi,"wtdiffsim"))

str(sd_pathogenic_df)

####################################
# Calculating mad for each residue
####################################

mad_byresidue = function(dataframe,variant,
                        position,angle){
  current_values = paste("Calculating", variant, 
                         position, angle, "mad...", sep = " ")
  print(current_values)
  wt_df1 = dataframe[[position]]
  wt_df2 = wt_df1[[variant]]
  wt_df3 = wt_df2[[angle]]
  mad = mad(wt_df3)
  return(mad)
}

#Function for loop for processing pvalues
mad_calculator = function(dataframe, refgroup,
                         angle){
  if (refgroup == "pathogenic") {
    list = pathogenic_var_list
  } else if (refgroup == "benign") {
    list = benign_var_list 
  } else if (refgroup == "vus") {
    list = vus_var_list
  } else if (refgroup == "wtdiffsim") {
    list = wtdiffsim_var_list
  }
  
  pervariant_mad = list()
  for (testvar in list) {
    output_mad = list() 
    for (position in aminoacid_residues) {
      print(paste("Processing ", position, " in "
                  ,testvar , "...", sep = ""))
      output_mad[[position]] = mad_byresidue(
        dataframe,testvar,position,angle)
    }
    pervariant_mad[[testvar]] = output_mad
  }
  return(pervariant_mad)
}

mad_benign_phi = mad_calculator(benign_by_residue, "benign", "phi")
mad_benign_psi = mad_calculator(benign_by_residue, "benign", "psi")
mad_pathogenic_phi = mad_calculator(pathogenic_by_residue, "pathogenic", "phi")
mad_pathogenic_psi = mad_calculator(pathogenic_by_residue, "pathogenic", "psi")
mad_wtdiffsim_phi = mad_calculator(wtdiffsim_by_residue, "wtdiffsim", "phi")
mad_wtdiffsim_psi = mad_calculator(wtdiffsim_by_residue, "wtdiffsim", "psi")


#Dataframe maker
mad_dataframe_maker = function(dataframe_phi, dataframe_psi, group){
  if (group == "pathogenic") {
    list = pathogenic_var_list
  } else if (group == "benign") {
    list = benign_var_list 
  } else if (group == "wtdiffsim") {
    list = wtdiffsim_var_list
  } else if (group == "vus") {
    list = vus_var_list
  }
  
  dataframe_per_variant = list()
  for (variant in list) {
    print(paste("Processing ",variant, " from ", 
                group, " group...", sep = ""))
    df_phi = as_tibble(matrixmaker(dataframe_phi,variant,"phi"))
    df_psi = as_tibble(matrixmaker(dataframe_psi,variant,"psi"))
    df_phi_psi = bind_cols(df_phi,df_psi)
    final_df = df_phi_psi %>%
      mutate(class = group) %>%
      mutate(key = variant)
    dataframe_per_variant[[variant]] = final_df
  }
  return(dataframe_per_variant)
}

mad_benign_df = bind_rows(mad_dataframe_maker(
  mad_benign_phi,mad_benign_psi,"benign"))
mad_pathogenic_df = bind_rows(mad_dataframe_maker(
  mad_pathogenic_phi,mad_pathogenic_psi,"pathogenic"))
mad_wtdiffsim_df = bind_rows(mad_dataframe_maker(
  mad_wtdiffsim_phi,mad_wtdiffsim_psi,"wtdiffsim"))

####################################
# Finding wildtype comaparable angles
####################################


df_col_names = colnames(median_wtdiffsim_df[1:548])
df_col_names

sd_by_angle = function(dataframe,angle){
   current_values = paste("Calculating",angle,
                           "standard deviation...", sep = " ")
    print(current_values)
    df1 = dataframe[[angle]]
    sd = sd(df1)
    return(sd) 
  
}
mad_by_angle = function(dataframe,angle){
  current_values = paste("Calculating",angle,
                         "standard deviation...", sep = " ")
  print(current_values)
  df1 = dataframe[[angle]]
  mad = mad(df1)
  return(mad) 
  
} #Mean absolute Deviation
summary_by_angle = function(dataframe,angle){
  current_values = paste("Calculating",angle,
                         "standard deviation...", sep = " ")
  print(current_values)
  df1 = dataframe[[angle]]
  summary = summary(df1)
  return(summary) 
  
} 

#Statistics summary
sd_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
   df = sd_by_angle(dataframe, angle)
   df_sd_by_angle[[angle]] = df
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Std_dev")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
mad_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = mad_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Median_Std_Dev")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
min_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[1]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Min.")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
firstQ_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[2]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "1st Qu.")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
median_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[3]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Median")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
mean_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[4]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Mean")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
thirdQ_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[5]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "3rd Qu.")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}
max_by_angle_df = function(dataframe) {
  df_sd_by_angle = c()
  for (angle in df_col_names) {
    df = summary_by_angle(dataframe, angle)
    df_sd_by_angle[[angle]] = df[6]
  }
  print("Making the dataframe...")
  dataframe = as_tibble(df_sd_by_angle)
  dataframe = pivot_longer(dataframe, cols = 1:548, names_to = "Angle",
                           values_to = "Max")
  print("Dataframe amalgamation completed!")
  return(dataframe)
}

#Final dataframe for analysis
summary_by_angle_df = function(dataframe) {
  sd_df = sd_by_angle_df(dataframe)
  mad_df = mad_by_angle_df(dataframe)
  min_df = min_by_angle_df(dataframe)
  firstQ_df = firstQ_by_angle_df(dataframe)
  med_df = median_by_angle_df(dataframe)
  mean_df = mean_by_angle_df(dataframe)
  thirdQ_df = thirdQ_by_angle_df(dataframe)
  max_df = max_by_angle_df(dataframe)
  #Joining the dataframes pair 1
  summary_df1 = left_join(sd_df, mad_df, by = "Angle")
  summary_df2 = left_join(min_df, firstQ_df, by = "Angle")
  summary_df3 = left_join(summary_df1, summary_df2, by = "Angle")
  #Pair 2
  summary_df4 = left_join(med_df,mean_df, by = "Angle")
  summary_df5 = left_join(thirdQ_df,max_df, by = "Angle")
  summary_df6 = left_join(summary_df4,summary_df5, by = "Angle")
  #Final joining
  final_df = left_join(summary_df3,summary_df6, by = "Angle")
  
  return(final_df)
}
summary_wtdiffsim_df = summary_by_angle_df(median_wtdiffsim_df)
summary_wtdiffsim_df

#High agreement angles in different simulations
high_agreement_angles = summary_wtdiffsim_df %>%
  filter(Std_dev <= 1) %>%
  arrange(desc(Std_dev))
hist(high_agreement_angles$Std_dev)
tail(high_agreement_angles)
high_agreement_angles
