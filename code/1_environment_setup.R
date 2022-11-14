################################################################################
#                        1: Environment setup script 
################################################################################

library(tidyr)
library(ggplot2)
library(cowplot)
library(plyr) #Always load this before dplyr
library(dplyr)
library(readr)
library(purrr) #Needed to filter list of dataframes
library(jpeg)
library(RColorBrewer)
library(datasets)
library(questionr)
library(stringr) #Good for text processing
library(pheatmap)
library(moments) #Checking skewness and kurtosis
library(reshape2)
library(car)
library(cluster)
library(ClusterR)
library(factoextra)
library(tidyverse)
library(extrafont) #For more font sizes
library(nonpar)
library(NSM3)

#xgboost related
# https://www.r-bloggers.com/2021/02/machine-learning-with-r-a-complete-guide-to-gradient-boosting-and-xgboost/

library(xgboost)
library(caTools)
library(caret)
library(Ckmeans.1d.dp) #plotting importance
library(FactoMineR) #For PCA analysis
library(DiagrammeR) #For decision tree making
library(ROCR)
library(ParBayesianOptimization)
library(Boruta)
library(rehh) #For ROC curve evaluations
library(cowplot)
library(doParallel)
#############################
#Listing directories
#############################

#Change the root directory and create the file paths below as folders and put the data from the repository 
rootdir = "Insert your root directory here and make the subsequent folder as shown by the path below"
setwd(rootdir)
getwd()
list.files()
raw_data_dir = paste(rootdir,"raw_data",sep = "/")
benign_dir = paste(raw_data_dir,"benign", sep = "/")
pathogenic_dir = paste(raw_data_dir,"pathogenic", sep = "/")
vus_dir = paste(raw_data_dir,"vus", sep = "/")
wt_dir = paste(raw_data_dir,"wildtype", sep = "/")
wt_diff_sim_dir = paste(wt_dir, "wildtype_different_simulations", sep = "/")
plot_dir = paste(rootdir, "Plots", sep = "/")

xgboost_output_dir = paste(rootdir, "xgboost_output", sep = "/")
list.files()

#############################
#Making variables from list
#############################
varmaker.txt = function(inputfilename,ext) {
  file.ext=paste(inputfilename,ext,sep = ".")
  list=as_tibble(read.delim2(file.ext))
  var_list=list$variant
}

#Normal PDBs
setwd(pathogenic_dir)
pathogenic_var=varmaker.txt("pathogenic_list","txt")
setwd(benign_dir)
benign_var=varmaker.txt("benign_list","txt") 
setwd(vus_dir)
vus_var=varmaker.txt("vus_list","txt")
setwd(wt_diff_sim_dir)
wtdiffsim_var=varmaker.txt("wildtype_different_simulation_list","txt")


#############################
#Inputting data frames
#############################


#https://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
#https://stackoverflow.com/questions/15561331/release-memory-in-r

setwd(pathogenic_dir)
list.files(".")
pathogenic_var_data = list()
for (variant in pathogenic_var) {
  print(paste0("Processing ",variant, "..." ))
  #Try catch is for showing the error but will not stop the loop
  tryCatch({
    txtfile = paste(variant,"txt", sep = ".")
    pathogenic_var_data[[variant]] = 
      as_tibble(read.delim(txtfile, sep = "")) %>%
      separate(aminoacid, into = c("aminoacid", "residue"),
               sep = "([-])", convert = TRUE) %>%
      mutate(varianttype = "pathogenic") %>%
      mutate(variantname = variant)},
    error=function(e){cat("ERROR :", conditionMessage(e), "/n")})
}

setwd(benign_dir)
list.files(".")
benign_var_data = list()
for (variant in benign_var) {
  print(paste0("Processing ",variant, "..." ))
  tryCatch({
    txtfile = paste(variant,"txt", sep = ".")
    benign_var_data[[variant]] = 
      as_tibble(read.delim(txtfile, sep = "")) %>%
      separate(aminoacid, into = c("aminoacid", "residue"),
               sep = "([-])", convert = TRUE) %>%
      mutate(varianttype = "benign") %>%
      mutate(variantname = variant)},
    error=function(e){cat("ERROR :", conditionMessage(e), "/n")})
} 

setwd(vus_dir)
list.files(".")
vus_var_data = list()
for (variant in vus_var) {
  print(paste0("Processing ",variant, "..." ))
  tryCatch({
    txtfile = paste(variant,"txt", sep = ".")
    vus_var_data[[variant]] = 
      as_tibble(read.delim(txtfile, sep = "")) %>%
      separate(aminoacid, into = c("aminoacid", "residue"),
               sep = "([-])", convert = TRUE) %>%
      mutate(varianttype = "VUS") %>%
      mutate(variantname = variant)},
    error=function(e){cat("ERROR :", conditionMessage(e), "/n")})
}

setwd(wt_diff_sim_dir)
list.files(".")
wtdiffsim_var_data = list()
for (variant in wtdiffsim_var) {
  print(paste0("Processing ",variant, "..." ))
  #Try catch is for showing the error but will not stop the loop
  tryCatch({
    txtfile = paste(variant,"txt", sep = ".")
    wtdiffsim_var_data[[variant]] = 
      as_tibble(read.delim(txtfile, sep = "")) %>%
      separate(aminoacid, into = c("aminoacid", "residue"),
               sep = "([-])", convert = TRUE) %>%
      mutate(varianttype = "wildtype") %>%
      mutate(variantname = variant) %>%
      mutate(residue = residue + 11)},
    error=function(e){cat("ERROR :", conditionMessage(e), "/n")})
}


#############################
#Sorting by residue/position
#############################

######Listing available variants per group

#normal pdb
pathogenic_var_list = names(pathogenic_var_data)
benign_var_list = names(benign_var_data)
vus_var_list = names(vus_var_data)
wtdiffsim_var_list = names(wtdiffsim_var_data)


#####Sorting the list of dataframes containing variants by residue
by_residue_sort = function(list.df, residuenum) {
  residue = map(list.df, ~filter(.x, residue == residuenum))
}

#####listing all the residue positions

#Normal pdb
aminoacid_residues = unique(benign_var_data$Arg291His$residue)
aminoacid_residues

benign_by_residue = list()
for (residue in aminoacid_residues){
  print(paste0("Processing ", residue, " position..."))
  benign_by_residue[[residue]] = 
    by_residue_sort(benign_var_data,residue)
}
rm(benign_var_data)

pathogenic_by_residue = list()
for (residue in aminoacid_residues){
  print(paste0("Processing ", residue, " position..."))
  pathogenic_by_residue[[residue]] = 
    by_residue_sort(pathogenic_var_data,residue)
}
rm(pathogenic_var_data)

vus_by_residue = list()
for (residue in aminoacid_residues){
  print(paste0("Processing ", residue, " position..."))
  vus_by_residue[[residue]] = 
    by_residue_sort(vus_var_data,residue)
}
rm(vus_var_data)
gc()

wtdiffsim_by_residue = list()
for (residue in aminoacid_residues){
  print(paste0("Processing ", residue, " position..."))
  wtdiffsim_by_residue[[residue]] = 
    by_residue_sort(wtdiffsim_var_data,residue)
}
rm(wtdiffsim_var_data)
gc()
