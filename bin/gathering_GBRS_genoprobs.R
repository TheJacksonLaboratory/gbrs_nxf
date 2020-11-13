#!/usr/bin/env Rscript
################## Gathering genoprobs from GBRS output file ###################

# This script goes to a directory containing the "gbrs.interpolated.genoprobs.tsv
# for each sample and convert them into the genoprobs format for QTL2.
# This script also save the output map to use it on QTL2.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 01_29_2020
# Modified for automation: Mike Lloyd
# Date: 20.02.10

################################################################################

# Set these variables acording to your data

library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Input, markergrid, output path, geoprob_name, and map_name must be supplied", call.=FALSE)
} else if (length(args)==5) {
  
  path_to_input <-  args[1]
  
  path_to_markergrid <- args[2]
  
  path_to_output <- args[3]
  
  genoprobs_name <- args[4]
  
  map_name <- args[5]
  
  
  path_to_input <- file_path_as_absolute(path_to_input)
  path_to_output <- file_path_as_absolute(path_to_output)
  path_to_markergrid <- file_path_as_absolute(path_to_markergrid)

  if(!dir.exists(path_to_input)) {
    stop("Input directory path does not exist as typed") 
  }
  
  if(!dir.exists(path_to_output)) {
    stop("Output directory path does not exist as typed") 
  }
  
  if(!file.exists(path_to_markergrid)) {
    stop("Markergrid file does not exist as typed") 
  }
  
}

################################################################################
############ loading libraries

library(stringr)
library(qtl2convert)
library(tidyverse)

################################################################################

# create list of all .gene.counts.txt files in folder

file_list <- list.files(path_to_input, pattern="*.tsv") # create list of all .tsv files in folder

# read in each .txt file in file_list and create a data frame with the same name as the file

list_geno <- function (x) {
  x <- read.table(paste0(path_to_input,"/",x),header=FALSE,sep = "\t")
  colnames(x) <- LETTERS[1:8]
  return(x) }

df_geno_list <- sapply(file_list, list_geno,simplify = FALSE, USE.NAMES = TRUE)

# Adjusting the names of the samples
# that all the txt files are named: DO#_SEX_Generation

names(df_geno_list) <- str_split(file_list,pattern = "_",simplify = TRUE)[,1]

names(df_geno_list) <- gsub("DO","DO.",names(df_geno_list))

rm(list_geno,file_list)

# Taking markers information

markers <- read.table(path_to_markergrid,header = TRUE,sep = "\t")

# Adjusting the markers object

markers <- markers %>% 
  mutate(chr = factor(chr,levels = c(1:19,"X")),
         marker = as.character(marker))

# Transforming the marker dataframe into QTL2 format (list)

map <- map_df_to_list(markers,pos_column = "pos")

# Making an empty array (mice, haplotypes, markers) and filling it by mouse

arr <- array(NA,dim=c(length(df_geno_list),8,69005),dimnames = list(names(df_geno_list),
                            LETTERS[1:8],markers$marker))

for (i in 1:length(df_geno_list)){
  arr[i,,] <- t(df_geno_list[[i]])
}; rm(i)

dim(arr)

if (dim(arr)[3] > dim(markers)[1]) {
  print(paste0("Some ids don't match for the row index: ",which(sum(c) > 0)))
  stopifnot(sum(c) > 0)
} else {
  print("Genotype array and grid file match")
}

# [1]   samples     8 69005

# Transforming the array to match QTL2 function

genoprobs <- probs_doqtl_to_qtl2(arr,map = markers,pos_column = "pos")

rm(arr)

# Saving the gathered/modified objects

saveRDS(genoprobs,file = paste0(path_to_output,"/",genoprobs_name,".RDS"))

saveRDS(map,file = paste0(path_to_output,"/",map_name,".RDS"))

