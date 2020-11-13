#!/usr/bin/env Rscript
################## QTLviewer formatting ###################

# This script transform the genoprobs and expression matrices created previously
# on the "gathering" scripts into a list called "dataset<name>" to be used 
# on QTL viewer.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 02_02_2020
# Modified for 

################ Set these variables acording to your data #####################
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=8) {
  stop("SOMETHING IS WRONG", call.=FALSE)
} else if (length(args)==8) {

  path_to_genoprobs <-  args[1]
  
  path_to_map <- args[2]
  
  path_to_samples_annotation <- args[3]
  
  path_to_expression <- args[4]
  
  path_to_output <- args[5]
  
  tissue <- args[6]
  
  display_name <- args[7]
  
  dataset_name <- args[8]  

  
  path_to_genoprobs <- file_path_as_absolute(path_to_genoprobs)
  path_to_map <- file_path_as_absolute(path_to_map)
  path_to_samples_annotation <- file_path_as_absolute(path_to_samples_annotation)
  path_to_expression <- file_path_as_absolute(path_to_expression)
  path_to_output <- file_path_as_absolute(path_to_output)
  
  if(!file.exists(path_to_genoprobs)) {
    stop("Genoprobs file/path does not exist as typed") 
  }
  
  if(!file.exists(path_to_map)) {
    stop("Map file/path does not exist as typed") 
  }
  
  if(!file.exists(path_to_samples_annotation)) {
    stop("Sample annotation file/path does not exist as typed") 
  }
  
  if(!file.exists(path_to_expression)) {
    stop("Expressoin file/path does not exist as typed") 
  }
  
  if(!dir.exists(path_to_output)) {
    stop("Output directory path does not exist as typed") 
  }
  
}


################################################################################
############ loading libraries

library(qtl2)
library(qtl2convert)
library(DESeq2)
library(ensimplR)
library(lme4)
library(tidyverse)
library(miqtl)



# 1) Reading variables

genoprobs <- readRDS(path_to_genoprobs)

map <- readRDS(path_to_map)

annot.samples <- read_csv(path_to_samples_annotation)

df_counts <- readRDS(path_to_expression)



# 2) Setting annotations

annot.samples <- annot.samples %>% 
  filter(Tissue == tissue) %>% 
  mutate(mouse.id = str_split(Customer_Sample_Name,pattern = "_",simplify = TRUE)[,1]) %>% 
  mutate(Sex = str_split(Customer_Sample_Name,pattern = "_",simplify = TRUE)[,2]) %>% 
  mutate(Generation = str_split(Customer_Sample_Name,pattern = "_",simplify = TRUE)[,3]) %>% 
  filter(grepl("DO",mouse.id)) %>% #There are non DO mice in the annotation dataframe
  mutate(mouse.id = gsub("DO","DO.",mouse.id)) %>% 
  rename(Batch = GT_Batch) %>% 
  select(mouse.id, Sex, Generation,Batch) %>% 
  filter(mouse.id %in% rownames(df_counts)) #Keeping just the mice that are in the
# expression matrix

### SETTING THE ANNOTATION MATRIX IS AN ISSUE WHICH NEEDS TO BE RESOLVED FOR AUTOMATED RUNNING.
  

# 3) Setting expression matrix

# Note: Filtering dataset, keeping only genes that have at least 1 read for at least half 
# of the samples

genes_names <- list()
for (i in 1:nrow(df_counts)) {
  genes_names[[i]] <- colnames(df_counts[,df_counts[i,]>=1])
}; rm(i)

genes_names <- table(unlist(genes_names))

filtered_genes_names <- genes_names[genes_names >= nrow(df_counts)/2] #~ half of the samples

df_counts <- df_counts[,colnames(df_counts) %in% names(filtered_genes_names)]

rm(filtered_genes_names,genes_names)

# 3.1) Saving raw matrix
raw <- df_counts

# 3.2) Normalizing data by vst and applying batch correction by lmer

norm <- floor(df_counts)

norm <- t(norm)

norm <- DESeq2::vst(norm)

norm <- t(norm)

norm <- norm %>% 
  as.data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  left_join(annot.samples,by = "mouse.id")

rm(df_counts)

debatch_lmer <- function(data) {
  
  genes <- grep(colnames(data), pattern = "ENSMUS", value = TRUE)
  
  resid_dat <- matrix(NA, nrow = nrow(data), ncol = length(genes))
  colnames(resid_dat) <- genes
  rownames(resid_dat) <- data$mouse.id
  
  
  for (i in 1:length(genes)) {
    
    lmer_fit <- suppressMessages(lmer(formula(paste(genes[i], "~ Sex + Generation + (1 | Batch)")), data = data))
    
    resid_dat[,i] <- data[,genes[i]] - ranef(lmer_fit)$Batch[data$Batch, 1]
    
    }
  
  return(resid_dat)
  
}


debatch_norm <- debatch_lmer(norm)



rm(debatch_lmer)

# 3.3) Creating normalized (norm) matrix

norm <- norm %>% 
  column_to_rownames("mouse.id") %>% 
  select(contains("ENSMUSG")) %>% 
  as.matrix()

# 3.4) Creating rankzed normalized (rz) matrix

rz <- apply(debatch_norm, 2, function(x) miqtl::rankZ(x))


######################### Generating checkup plots #############################
################################################################################

pdf(paste0(path_to_output,"/",tissue,"_checkup_boxplot_expression.pdf"))

# Norm before batch correction

y <- rowMeans(norm)

rownames(norm) %>% identical(names(y))
# [1] TRUE

expr_mean_rna <- cbind.data.frame(Expr = y, Sex=annot.samples$Sex)

#quartz()
expr_mean_rna %>% 
  ggplot(aes(x=as.factor(Sex),y=Expr)) + 
  geom_boxplot() +
  ggtitle("Expression - before batch correction") +
  labs(x = "Sex",
        y = "Vst transform expression")

# Norm after batch correction

y <- rowMeans(debatch_norm)

rownames(debatch_norm) %>% identical(names(y))
# [1] TRUE

expr_mean_rna <- cbind.data.frame(Expr = y, Sex=annot.samples$Sex)

#quartz()
expr_mean_rna %>% 
  ggplot(aes(x=as.factor(Sex),y=as.numeric(Expr))) + 
  geom_boxplot()+
  ggtitle("Expression - after batch correction") +
  labs(x = "Sex",
       y = "Vst transform expression")

# Rankz after batch correction

y <- rowMeans(rz)

rownames(rz) %>% identical(names(y))
# [1] TRUE

expr_mean_rna <- cbind.data.frame(Expr = y, Sex=annot.samples$Sex)

#quartz()
expr_mean_rna %>% 
  ggplot(aes(x=as.factor(Sex),y=as.numeric(Expr))) + 
  geom_boxplot()+
  ggtitle("Expression - after batch correction") +
  labs(x = "Sex",
       y = "Rankz transform expression")

dev.off()

rm(expr_mean_rna)


pdf(paste0(path_to_output,"/",tissue,"_checkup_density_expression.pdf"))

# Norm before batch correction

expr_density <- as.data.frame(norm)

expr_density$Mouse_ID <- rownames(expr_density)

expr_density <- gather(expr_density,"genes_ID","Expression",c(1:ncol(raw)))

expr_density$Mouse_ID <- as.factor(expr_density$Mouse_ID)

#quartz()
g <- expr_density %>% ggplot(aes(x=Expression,colour=Mouse_ID))+
  geom_line(stat = 'density', alpha = 0.05, aes(color = Mouse_ID))+
  ggtitle("Expression - before batch correction") +
  labs(x = "Vst transform expression",
       y = "Density")

g+theme(legend.position = "none")

# Norm after batch correction

expr_density <- as.data.frame(debatch_norm)

expr_density$Mouse_ID <- rownames(expr_density)

expr_density <- gather(expr_density,"genes_ID","Expression",c(1:ncol(debatch_norm)))

expr_density$Mouse_ID <- as.factor(expr_density$Mouse_ID)

#quartz()
g <- expr_density %>% ggplot(aes(x=Expression,colour=Mouse_ID))+
  geom_line(stat = 'density', alpha = 0.05, aes(color = Mouse_ID))+
  ggtitle("Expression - after batch correction") +
  labs(x = "Vst transform expression",
       y = "Density")

g+theme(legend.position = "none")

# Rankz after batch correction

expr_density <- as.data.frame(rz)

expr_density$Mouse_ID <- rownames(expr_density)

expr_density <- gather(expr_density,"genes_ID","Expression",c(1:ncol(rz)))

expr_density$Mouse_ID <- as.factor(expr_density$Mouse_ID)

#quartz()
g <- expr_density %>% ggplot(aes(x=Expression,colour=Mouse_ID))+
  geom_line(stat = 'density', alpha = 0.05, aes(color = Mouse_ID))+
  ggtitle("Expression - after batch correction") +
  labs(x = "Rankz transform expression",
       y = "Density")

g+theme(legend.position = "none")

dev.off()

rm(expr_density,g,norm)

################################################################################
################################################################################

# 4) Sorting genoprobs mouse.id according to the annotations dataframe

for (i in 1:length(genoprobs)){
  genoprobs[[i]] <- genoprobs[[i]][annot.samples$mouse.id,,]
};rm(i)




# 5) Creating the covariates variables

covar.matrix <- model.matrix(~ Sex + Generation, data = annot.samples)[,-1]

rownames(covar.matrix) <- annot.samples$mouse.id

class(covar.matrix) #matrix, done. 

covar.info <- tibble(sample.column = c("Sex","Generation"),
                     covar.column = c("Sex","Generation"),
                     display.name = c("Sex","Generation"),
                     interactive = as.factor(c("TRUE","TRUE")),
                     primary = as.factor(c("TRUE","FALSE")),
                     lod.peaks = c("sex_int","gen_int"))

covar.info <- covar.info %>% 
  mutate(interactive = as.logical(interactive),
         primary = as.logical(primary))



# 6) Creating the kinship, map and ensemble version variables

K <- calc_kinship(genoprobs, type = "loco", cores = 2)

markers <- map_list_to_df(map)               
str(markers) #It should be: data.frame(chr,num,chr) - done.

ensembl.version <- 94



# 7) Creating the annotations for RNA (using Matt function - BatchGenes)

annot.mrna <- batchGenes(colnames(debatch_norm))

# 7.1) Just keeping genes that have annotation

raw <- raw[,colnames(raw) %in% annot.mrna$id]

#dim(raw)

debatch_norm <- debatch_norm[,colnames(debatch_norm) %in% annot.mrna$id]

#dim(debatch_norm)

rz <- rz[,colnames(rz) %in% annot.mrna$id]

#dim(rz)

entrez_id <- str_split(annot.mrna$external_ids,",",simplify = TRUE)[,1]

entrez_id <- regmatches(entrez_id, gregexpr("[[:digit:]]+", entrez_id))

names(entrez_id) <- annot.mrna$id

test_list <- function(x) {
  if (length(x) > 1 | length(x) == 0) {
    x <- "NA"
  } else {
    x
  }
}

entrez_id <- lapply(entrez_id, test_list)

entrez_id <- bind_rows(entrez_id)

annot.mrna <- annot.mrna %>% 
  select(id,symbol,chromosome,start,end,strand) %>% 
  mutate(entrez_id = as.character(entrez_id),
         middle = (((end-start)/2)+start)/1000000,
         start = start/1000000,
         end = end/1000000) %>% 
  rename(gene.id = id,
         chr = chromosome) %>% 
  select(gene.id, symbol, entrez_id,chr,start,end,middle,strand)


nearest_marker <- character(length = nrow(annot.mrna))
for(i in 1:nrow(annot.mrna)){
  cm <- annot.mrna$chr[i]
  if (cm %in% c(1:19,"X")){
    sub <- subset(markers, chr == cm)
    nearest_marker[i] <- sub$marker[which.min(abs(sub$pos - annot.mrna$start[i]))]
    
  }};rm(i,sub,cm)

annot.mrna$nearest.marker.id <- nearest_marker 



# 8) Confirming that mice order is identical

identical(annot.samples$mouse.id,rownames(raw))
identical(annot.samples$mouse.id,rownames(debatch_norm))
identical(annot.samples$mouse.id,rownames(rz))
identical(rownames(K$`1`),rownames(genoprobs$`1`))



# 9) Creating the variables for mRNA QTLviewer

if (startsWith(dataset_name, 'dataset.')) {
  assign(dataset_name,
         list(annot.mrna = as_tibble(annot.mrna),
              covar.matrix = covar.matrix,
              covar.info = covar.info,
              datatype = "mrna",
              display.name = display_name,
              data = list(raw = as.matrix(raw),
                          norm = as.matrix(debatch_norm),
                          rz = as.matrix(rz)),
              annot.samples = as_tibble(annot.samples)))
  
  save(list = c("genoprobs","K","map","markers","ensembl.version",dataset_name),
       file = paste0(path_to_output,"/",dataset_name,".RData"))
  
} else {
  dataset_mod_name <- paste0('dataset.', dataset_name)
  
  assign(dataset_mod_name,
         list(annot.mrna = as_tibble(annot.mrna),
              covar.matrix = covar.matrix,
              covar.info = covar.info,
              datatype = "mrna",
              display.name = display_name,
              data = list(raw = as.matrix(raw),
                          norm = as.matrix(debatch_norm),
                          rz = as.matrix(rz)),
              annot.samples = as_tibble(annot.samples)))
  
  save(list = c("genoprobs","K","map","markers","ensembl.version",dataset_mod_name),
       file = paste0(path_to_output,"/",dataset_name,".RData"))
  
}

