#   设置工作目录
setwd("D:/Study the R 22.5.14/R TCGA test/Bladder cancer TCGA/Bladder TCGA")

#批量读取
# stage 1: clean the data

rm(list = ls())     
library(dplyr)
library(readxl)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)
library(survMisc)
library(tidyr)
library(stringi)


data_dir <- "D:/Study the R 22.5.14/R TCGA test/Lung cancer TCGA/Lung TCGA"
sig_dir <- "D:/Study the R 22.5.14/R TCGA test/Lung cancer TCGA/Lung TCGA"
star_counts_data_dir <- "D:/Study the R 22.5.14/R TCGA test/Lung cancer TCGA/Lung TCGA/Lung cancer data"
#Read the sample sheet(read.table)
fgss1 <- read.table("gdc_sample_sheet.2023-01-15.tsv",header = T, sep = "\t",
                    check.names = F, stringsAsFactors = F) 
fgss <- fgss1[which(fgss1$'Sample Type' %in% "Primary Tumor"),]
colnames(fgss)[2] <- "file_name"
metadata <- jsonlite::fromJSON("metadata.cart.2023-01-15.json")
View(fgss) #fgss包含File ID	，File Name等8项对应关系

files <- list.files(star_counts_data_dir,recursive = T)
files <- files[! str_detect(files,"parcel")]
files <- files[str_detect(files,"rna_seq")]
matrix <- data.frame()
i=0
for (file in files) {
  i=i+1
  print(i)
  filename <- file %>% str_remove(".*/")
  NewName <- metadata %>% filter(file_name == filename)
  NewName <- (NewName$associated_entities[[1]])$entity_submitter_id # %>% substr(1,16)
  print(NewName)
  middata <- read.table(file.path(star_counts_data_dir,file),header = T,row.names = 1,sep = "\t",skip = 5)
  midmatrix <- middata %>% select(X.4)
  colnames(midmatrix) = NewName
  matrix <-  midmatrix %>% t.data.frame() %>% rbind(matrix)
}
#save(matrix,file = "Loaded lung matrix.RData")
load(file="Loaded lung matrix.RData")

matrix <- matrix %>% t()
colnames(matrix) <- substr(colnames(matrix),1,16)##只取1-16位
a <- as.data.frame(gsub('[-]', '.', fgss$'Sample ID'))#把-替换成点
colnames(a) <- c("Sample ID")
fzt <- intersect(colnames(matrix),a$'Sample ID')#取交集
rownames(matrix) <- substr(rownames(matrix),1,17)#ENSG取小数点后一位便于匹配
matrix_pt <- matrix[,fzt]
matrix=as.data.frame(matrix)#$ operator is invalid for atomic vectors
#matrix_pt <- cbind(matrix$X,matrix_pt)
#colnames(matrix_pt)[1] <- "Ensembl ID"
dim(matrix) #检索matrix的维度
matrix[1:4,1:4]
#View(matrix)