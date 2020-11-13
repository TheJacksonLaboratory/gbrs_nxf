#!/usr/bin/env Rscript

################## Gathering counts from GBRS output file ######################

# This script goes to a directory containing the "gbrs_quantified_diploid_genes_
# expected_read_counts" for each sample and convert them into
# a matrix (samples x total gene counts per gene).

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 01_29_2020
# Modified for automation: Mike Lloyd
# Date: 20.02.10

################################################################################
############ 

library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Input, output, and output name must be supplied", call.=FALSE)
} else if (length(args)==3) {
  path_to_input <-  args[1]
  
  path_to_output <- args[2]
  
  output_name <- args[3]
  
  path_to_input <- file_path_as_absolute(path_to_input)
  path_to_output <- file_path_as_absolute(path_to_output)

  if(!dir.exists(path_to_input)) {
    stop("Input directory path does not exist as typed") 
  }

  if(!dir.exists(path_to_output)) {
    stop("Output directory path does not exist as typed") 
  }

}

################################################################################
############ loading libraries

library(tidyverse)
library(stringr)

################################################################################


################################################################################

# create list of all .gene.counts.txt files in folder

file_list <- list.files(path_to_input, pattern="*diploid.genes.expected_read_counts") 

# function read in each .txt file in file_list and create a data frame with the same name as the file

list_counts <- function (x) {
  x <- read.table(paste0(path_to_input,"/",x),header=TRUE,sep = "\t")
  colnames(x) <- c("Target_ID",LETTERS[1:8],"Total","Notes") 
  return(x) }

# Applying to the list

df_counts_list <- sapply(file_list, list_counts,simplify = FALSE, USE.NAMES = TRUE)

# Adjusting the names of the samples
# that all the txt files are named: DO#_SEX_Generation

names(df_counts_list) <- str_split(file_list,pattern = "_",simplify = TRUE)[,1]

names(df_counts_list) <- gsub("DO","DO.",names(df_counts_list))

rm(list_counts,file_list)

#Checking if all the genes ids (Target_ID column) for all the samples are the same

temp <- df_counts_list
for (i in 1:length(df_counts_list)){
  temp[[i]] <- select(temp[[i]],Target_ID)
}; rm(i)

b <- bind_cols(temp)

b <- apply(b, 2, function (x) factor(x,labels = c(1:nrow(b))))

b <- apply(b, 2, as.numeric)

compare_genes <- function(x) {
  
  c <- apply(x, 1, sd)
  
  if (sum(c) > 0) {
    print(paste0("Some ids don't match for the row index: ",which(sum(c) > 0)))
    stopifnot(sum(c) > 0)
  } else {
    print("All ids match!")
  }
  
}

compare_genes(b)

rm(compare_genes,temp,b)

# Creating input data with the total gene counts for all the samples

df_counts <- df_counts_list
for (i in 1:length(df_counts_list)){
  df_counts[[i]] <- select(df_counts_list[[i]],Target_ID,Total)
}; rm(i)

df_counts <- bind_rows(df_counts,.id = "Sample")

df_counts <- df_counts %>% 
  spread(Target_ID,Total) %>% 
  column_to_rownames("Sample")
  

# Saving the final expression matrix

saveRDS(df_counts,file = paste0(path_to_output,"/",output_name,".RDS"))