# -*- coding: utf-8 -*-
##### TCGA data processing
options(stringsAsFactors = F)
if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr")
library(readr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(biomaRt)
library(tibble)
library(igraph)
library(readr)
#####
library(SummarizedExperiment)
library(dplyr)
library(lubridate)
library(tidyverse)
library(EDASeq)
############
library(Rtsne)
library(ggplot2)
# install.packages("plotly")
library(plotly)
library(colorspace)
#####################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("EDASeq", quietly = TRUE))
  BiocManager::install("EDASeq")



# ####################################################workround###############################################################
# https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/627
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

#############################################################################################################################
#### Collect gene expression data#################################################################
setwd("/projappl/project_2010541/code/")
# dataset_TCGA <- c("BLCA", "BRCA", "COAD", "ESCA", "KICH", "KIRC", 
#                   "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", 
#                   "READ", "SKCM", "STAD", "THCA", "THYM", "UCEC")
# 
# tried_groups<-c("BLCA",  "COAD",  "KICH", "KIRC", "KIRP","SKCM","LIHC","LUAD","LUSC","UCEC","READ","PAAD","BRCA")
# still_groups<-c( )
# not_installed<-c("ESCA","THYM","THCA","STAD")
# 
# #print(cancer_type)
# in_pool<-"THYM"
# doing_groups<-c("THYM")
##############################################################################
cancer_type<-"BRCA"
query_exp <- GDCquery(project = paste0("TCGA-",cancer_type), 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")#
GDCdownload(query_exp) 
data_exp <- GDCprepare(query_exp) 
print(cancer_type)
##########################################################################################
#https://www.biostars.org/p/9544856/
fpkm_uq_unstrand<-assays(data_exp)[["fpkm_uq_unstrand"]]
exp_assay<-as.data.frame(fpkm_uq_unstrand)
exp_rowRanges <- as.data.frame(rowRanges(data_exp)) # Gene annotation
saveRDS(exp_assay,paste0("../data/tcga_data/",cancer_type,"_exp_assay.RData"))
exp_assay<-readRDS(paste0("../data/tcga_data/",cancer_type,"_exp_assay.RData"))
saveRDS(exp_rowRanges, paste0("../data/tcga_data/",cancer_type,"_exp_rowRangess.RData"))
#exp_rowRanges <- as.data.frame(rowRanges(exp_assay))
############################################################################################################
#### Collect DNA methylation data
query_mty <- GDCquery(project = paste0("TCGA-",cancer_type), 
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.type = "Methylation Beta Value")
GDCdownload(query_mty)
data_mty <- GDCprepare(query_mty)
mty_assay <- as.data.frame(assay(data_mty)) # DNA methylation matrix
mty_colData <- as.data.frame(colData(data_mty)) # Patient annotation(475 patients)
#View(mty_colData)
mty_rowRanges <- as.data.frame(rowRanges(data_mty)) # cg probe annotation
#######################Save files##########################################
saveRDS(mty_colData, paste0("../data/tcga_data/",cancer_type,"_mty_colData.RData"))
saveRDS(mty_assay, paste0("../data/tcga_data/",cancer_type,"_mty_assay.RData"))
saveRDS(mty_rowRanges, paste0("../data/tcga_data/",cancer_type,"_mty_rowRanges.RData"))

mty_colData<-readRDS(paste0("../data/tcga_data/",cancer_type,"_mty_colData.RData"))
mty_assay<-readRDS(paste0("../data/tcga_data/",cancer_type,"_mty_assay.RData"))
mty_rowRanges<-readRDS(paste0("../data/tcga_data/",cancer_type,"_mty_rowRanges.RData"))
#########################################Masked Somatic Mutation ################################################
query <- GDCquery(
  project = paste0("TCGA-", cancer_type),
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  #workflow.type = "MuSE Variant Aggregation and Masking"
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
print(cancer_type)
######################################################################################
# file_paths <- getResults(query)$file_id
# # #file_path<-"0b0de7da-a900-4595-bf65-805513874ecf"
# # library(readr)
# # 
#  data <- read_tsv("./GDCdata/TCGA-ESCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/1cdbcfc2-2cfc-435b-b1f3-ea1824983b1b/0695bd65-8748-4396-b47b-f3fa238895fb.wxs.aliquot_ensemble_masked.maf.gz", comment = "#")
# # 
#  View(data)
# # data$Tumor_Seq_Allele2 <- as.character(data$Tumor_Seq_Allele2)
# # #View(data)
# # #View(data)
# # 
# gz_files <- c()
# for (file_path in file_paths) {
# #   # 使用 list.files 函数找到所有 .gz 文件
#   gz_files <- c(gz_files, list.files(path = paste0("./GDCdata/TCGA-ESCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/",file_path), pattern = "\\.gz$", full.names = TRUE))
#  }
# # 
# print(gz_files)
# # 
# # 
# # for (gz_file in gz_files) {
# #   df <- read_tsv(gz_file,comment = "#",show_col_types = FALSE)
# #   #df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE,comment = "#"))
# #   df$Tumor_Seq_Allele2 <- as.character(df$Tumor_Seq_Allele2)
# #   #write.table(df, file_path, row.names = FALSE)
# #   temp_file <- tempfile(fileext = ".gz")
# #   write.table(df, gzfile(temp_file), row.names = FALSE, sep = "\t", quote = FALSE)
# #   #gzip(temp_file, destname = "dbbd21e0-3f6a-43cd-868e-5d3ed3e7a024.wxs.aliquot_ensemble_masked.maf.gz", overwrite = TRUE)
# #   file.copy(temp_file, gz_file, overwrite = TRUE)
# # }
###################################################################################
#View(as.data.frame(file_paths))
maf <- GDCprepare(query)
#View(maf)
#View(clinicalInfo)
filename1 <- paste0("../data/tcga_data/", cancer_type, "_maf_test06.RData")
#print(filename1)
saveRDS(maf, filename1)
maf<-readRDS(filename1)
#############################################################"Copy Number Variation"###################################################################################
print(cancer_type)
query_cnv <- GDCquery(project = paste0("TCGA-", cancer_type),
                      data.category = "Copy Number Variation",
                      data.type = "Gene Level Copy Number")

GDCdownload(query_cnv)
potential_df <- query_cnv$results[[1]]
potential_df$created_datetime <- ymd_hms(potential_df$created_datetime)
potential_df$timestamp <- as.numeric(potential_df$created_datetime)

#######################################################################
unique_cases_df <- potential_df %>%
  distinct(cases, .keep_all = TRUE) %>%
  group_by(cases) %>%
  mutate(datetime = ymd_hms(created_datetime)) %>%  # 将created_datetime转换为日期时间对象
  filter(datetime == max(datetime)) %>%            # 过滤出最新的日期时间
  ungroup() %>%
  dplyr::select(-datetime)%>%
  dplyr::select(-timestamp)
#####################################################
query_cnv$results[[1]] <- unique_cases_df
##########################################################################
saveRDS(query_cnv,paste0("../data/tcga_data/",cancer_type,"_query_cnv_test06.RData"))
query_cnv<-readRDS(paste0("../data/tcga_data/",cancer_type,"_query_cnv_test06.RData"))
write.csv(query_cnv$results[[1]], file = paste0("../data/tcga_data/",cancer_type,"_query_results_test06.csv", row.names = TRUE))
############################################################
query_cnv_copy<-query_cnv
query_cnv_copy$results[[1]] <- query_cnv_copy$results[[1]] %>%
  distinct(sample.submitter_id, .keep_all = TRUE)

GDCdownload(query_cnv_copy)
data_cnv <- GDCprepare(query_cnv_copy)
data_cnv_savage<-data_cnv
saveRDS(data_cnv_savage,"../data/tcga_data/data_cnv_savage_GDCprepare_test.RData")
data_temp<-data_cnv_savage#
View(data_temp)
data_temp <- as.data.frame(assay(data_temp))
#View(data_temp)
dim(data_temp)
saveRDS(data_temp,paste0("../data/tcga_data/",cancer_type,"_data_temp_savage_test.csv"))
data_temp<-readRDS(paste0("../data/tcga_data/",cancer_type,"_data_temp_savage_test.csv"))
#################################################### Collect clinical data################################################################################################
clinical <- GDCquery_clinic(project = paste0("TCGA-",cancer_type), type = "clinical")
#View(clinical)
####################################################Collect clinical radiation and drug therapy data###############################################
#print(cancer_type)
query <- GDCquery(project = paste0("TCGA-",cancer_type), 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab"
)
GDCdownload(query)
sckm.tab.all <- GDCprepare(query)
column_name0 <- paste0("clinical_drug_", tolower(cancer_type))
therapy <- sckm.tab.all[[column_name0]] # clinical drug therapy ->NULL
#View(therapy)
saveRDS(therapy, paste0("../data/tcga_data/",cancer_type,"_therapy_test.RData"))
therapy<-readRDS(paste0("../data/tcga_data/",cancer_type,"_therapy_test.RData"))
therapy$pharmaceutical_therapy_type # Therapy types
column_name1<-paste0("clinical_radiation_", tolower(cancer_type))
radiation <- sckm.tab.all[[column_name1]] # Clinical radiation therapy
#View(radiation)
saveRDS(radiation, paste0("../data/tcga_data/",cancer_type,"_radiation_test.RData"))
radiation<-readRDS(paste0("../data/tcga_data/",cancer_type,"_radiation_test.RData"))
##########################################################################################################################

# ########################################################## Process gene expression data#################################

exp_assay_frame<-exp_assay
dim(exp_assay_frame)
colnames(exp_assay_frame) <- substr(colnames(exp_assay_frame), 1, 16)
true_index<-!(is.na(exp_rowRanges[row.names(exp_assay_frame), "gene_name"]))
exp_assay_frame$SYMBOL <- exp_rowRanges[row.names(exp_assay_frame), "gene_name"]
exp_assay_frame<-exp_assay_frame[true_index,]
#which(is.na(exp_assay_frame$SYMBOL))
saveRDS(exp_assay_frame,paste0("../data/tcga_data/",cancer_type,"_exp_assay_frame01.RData"))
####################################################################################################################################################
num_rows <- nrow(exp_assay_frame)
num_cols <- ncol(exp_assay_frame)

exp_assay_test1 <- exp_assay_frame[,c(num_cols,1:num_cols-1)]

exp_assay_test2 <- as.matrix(exp_assay_test1)

row.names(exp_assay_test2) <- exp_assay_test2[,1]

exp_assay_test2 <- exp_assay_test2[,2:num_cols] # Convert Ensemble ID to corresponding gene symbols

row_names_tcga<-rownames(as.data.frame(exp_assay_test2)) 

exp_assay_test2 <- matrix(as.numeric(exp_assay_test2), nrow = nrow(exp_assay_test2), dimnames = list(row.names(exp_assay_test2), colnames(exp_assay_test2)))
any(is.na(row.names(exp_assay_test2)))
exp_assay_test2 <- avereps(exp_assay_test2) # If the gene corresponds to multiple gene expression values, take the average value,remove duplicated gene name

exp_assy_trial<-exp_assay_test2
#############################################################################################################################

exp_intgr <- t(exp_assy_trial) # Gene expression data matrix
exp_intgr <- log10(exp_intgr + 1) # 1og10 conversion
saveRDS(exp_intgr, paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all01.RData"))
write.csv(exp_intgr,paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all01.csv"))
##############################################################################################################################
#### Process DNA methylation data

data_mty <- subset(data_mty, subset = (rowSums(is.na(assay(data_mty)))==0)) 

data_mty<- subset(data_mty, subset =!is.na(as.data.frame(rowRanges(data_mty))$gene_HGNC))
dim(data_mty)

# ############################################################

mty_assay <- as.data.frame(assay(data_mty))
mty_rowRanges <- as.data.frame(rowRanges(data_mty))
mty_symbol <- strsplit(mty_rowRanges$gene_HGNC, split = ";")
names(mty_symbol) <- row.names(mty_rowRanges)
mty_symbol <- lapply(mty_symbol, FUN = function(x){x<-unique(x)})
mty_symbol <- as.matrix(unlist(mty_symbol))
row.names(mty_symbol) <- substr(row.names(mty_symbol), 1, 10)
mty_symbol <- data.frame("probe" = row.names(mty_symbol),"SYMBOL" = mty_symbol[,1])#cg14528386.1 ->H19_1
mty_assay$probe <- row.names(mty_assay)
mty_assay <- merge(mty_assay, mty_symbol, by.x = "probe", by.y = "probe") 
num_rows <- nrow(mty_assay)

# 获取列数
num_cols <- ncol(mty_assay)
mty_mat <- as.matrix(mty_assay[,2:num_cols])
colnames(mty_mat) <- substr(colnames(mty_mat), 1, 16)
row.names(mty_mat) <- mty_assay$SYMBOL # Convert probe ID to corresponding gene symbols
mty_mat <- matrix(as.numeric(mty_mat), nrow = nrow(mty_mat), dimnames = list(row.names(mty_mat), colnames(mty_mat)))
mty_mat <- avereps(mty_mat) # If the gene corresponds to multiple methylation values, take the average value
dim(mty_mat)
write.csv(mty_mat,paste0("../data/",cancer_type,"_mty_mat_all.csv"))
saveRDS(mty_mat,paste0("../data/",cancer_type,"_mty_mat_all.RData"))
mty_mat_trial<-readRDS(paste0("../data/",cancer_type,"_mty_mat_all.RData"))
mty_mat_trial <- t(mty_mat_trial) # Gene expression data matrix
saveRDS(mty_mat_trial, paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all01.RData"))
write.csv(mty_mat_trial,paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all01.csv"))
samples <- intersect(colnames(exp_intgr), colnames(mty_mat_trial)) 
print(samples)
#dim(samples)
saveRDS(samples, paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))

#### Integrate the genomic data into the network
exp_intgr <- exp_intgr[,samples]
mty_mat_trial <- mty_mat_trial[,samples]
rownamesSamples<-intersect(rownames(exp_intgr),rownames(mty_mat_trial))
length(rownamesSamples)
#class(exp_intgr)# "matrix" "array
exp_intgr_trial<-exp_intgr[rownamesSamples,]
mty_mat_trial<-mty_mat_trial[rownamesSamples,]
dim(exp_intgr_trial)
#View(exp_intgr_trial)
################################################################################################################################################################

clinicalInfo <- mty_colData
clinicalInfo<-clinicalInfo[!duplicated(clinicalInfo$sample),]
indices <- which(clinicalInfo$shortLetterCode == "NT",clinicalInfo$sample)
#print(indices)
#View(clinicalInfo)

NT_remove_name<-(rownames(clinicalInfo[indices,]))
print(NT_remove_name)
# View(NT_remove_name)
# new_clinicalInfo<-clinicalInfo[samples,]
# length(samples)
# #dim(new_clinicalInfo)
# new_clinicalInfo<-clinicalInfo[-indices,]
# View(new_clinicalInfo)
# saveRDS(new_clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
############################################################
new_exp_intgr<-exp_intgr_trial[!(rownames(exp_intgr_trial) %in% NT_remove_name),]

new_mty_intgr<-mty_mat_trial[!(rownames(mty_mat_trial) %in% NT_remove_name),]

saveRDS(new_exp_intgr, paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all.RData"))
write.csv(new_exp_intgr,paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all.csv"))

write.csv(new_mty_intgr,paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.csv"))
saveRDS(new_mty_intgr,paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.RData"))

# dim(new_exp_intgr)
# dim(new_mty_intgr)

# #######################################################MAF start for data_cnv######################################################################

maf <- maf[,c("Tumor_Sample_Barcode","Hugo_Symbol","Gene","Variant_Classification")]

rnames <- unique(maf$Hugo_Symbol)
cnames <- unique(maf$Tumor_Sample_Barcode)
snv_count <- matrix(data = 0, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames,cnames)) 

# Calculate the frequency of genes' variants
for(i in 1:nrow(maf)){
  rname <- maf[i,]$Hugo_Symbol
  cname <- maf[i,]$Tumor_Sample_Barcode
  snv_count[rname, cname] <- snv_count[rname,cname] + 1
}
dim(snv_count)
#https://www.notion.so/57782b3efe144e91a93ae602f85a1cbd

#View(snv_count)
colnames(snv_count) <- substr(colnames(snv_count), 1, 16)

# ##############################################################################################!!!!!!!!!!!!!!!

data_temp$ensembl<-rownames(data_temp)
#View(as.data.frame(rownames(data_temp)))

rownames(data_temp)<-seq_along(rownames(data_temp))

data_temp$ensembl<-substr((data_temp$ensembl),1, 15)

na_count <- function(x) sum(is.na(x))

data_temp$na_count <- apply(data_temp, 1, na_count)

data_temp <- data_temp[order(data_temp$ensembl, data_temp$na_count), ]

data_temp <- data_temp[!duplicated(data_temp$ensembl), ]
#remove na_count
data_temp$na_count <- NULL

row.names(data_temp)<-data_temp[,ncol(data_temp)]

#
#dim(data_temp)
data_cnv_tmp<-data_temp
#dim(data_cnv)
colnames(data_cnv_tmp) <- substr(colnames(data_cnv_tmp), 1, 16)

# View(data_cnv_tmp)

###############################################################################################!!!!!!!!!!!!!!!!!!!!
col_data_cnv<-ncol(data_cnv_tmp)
#View(as.data.frame(col_data_cnv))#146
##########################################################################################################
# Convert Ensemble ID to corresponding gene symbols
# Using biomaRt for gene ID conversion will cause some loss
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")


# class(ensembl)

# View(ensembl)
# View(listFilters(ensembl))

# 	hgnc_symbol  HGNC symbol(s) [e.g. A1BG]
# ensembl_gene_id  Gene stable ID(s) [e.g. ENSG00000000003]

# print(row.names(data_cnv))
# View(data_cnv)

# View(as.data.frame(row.names(data_cnv_tmp)))


cnv_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = c("ensembl_gene_id"),
                values = row.names(data_cnv_tmp),
                mart = ensembl)

#View(cnv_df)
##############################################################
cnv_df <- cnv_df[which(cnv_df$hgnc_symbol !=''),]
#View(cnv_df)

cnv_df<-unique(cnv_df)

colnames(data_cnv_tmp)[colnames(data_cnv_tmp)=="ensembl"]<-"emsembl"

unique_cols<-!duplicated(t(data_cnv_tmp))
data_cnv_tmp<-data_cnv_tmp[,unique_cols]
#############################################################################################################

data_cnv_tmp <- merge(x = data_cnv_tmp, y = cnv_df, by.x = "emsembl", by.y = "ensembl_gene_id")
data_cnv_tmp <- as.matrix(data_cnv_tmp)
dim(data_cnv_tmp)
#############################################################################################################
##########################################################################################################

row.names(data_cnv_tmp) <- data_cnv_tmp[,ncol(data_cnv_tmp)]
#View(data_cnv_tmp)
data_cnv_tmp <- data_cnv_tmp[,2:ncol(data_cnv_tmp)-2]  # Convert Ensemble ID to corresponding gene symbols
data_cnv_tmp <- matrix(as.numeric(data_cnv_tmp), nrow = nrow(data_cnv_tmp), dimnames = list(row.names(data_cnv_tmp), colnames(data_cnv_tmp)))
data_cnv_tmp <- data_cnv_tmp[!duplicated(row.names(data_cnv_tmp)),] # Only one gene PRAMEF7 has duplicate copy number variation value, and the duplicate value is the same
data_cnv_tmp<-data_cnv_tmp[,-1]

print(rownamesSamples)
samples_f1<-intersect(rownamesSamples, colnames(snv_count))
print(samples_f1)

print(colnames(snv_count))
samples_f2 <- intersect(samples_f1, colnames(data_cnv_tmp))
print(samples_f2)
print(cancer_type)
saveRDS(samples_f2,paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))

sample_f2_check<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))
print(sample_f2_check)
#length(samples_f2)

# print(samples_f2)
#
# class(new_exp_intgr)
# class(samples_f2)
# new_exp_intgr_tmp1 <- new_exp_intgr[,samples_f2]
# new_mty_intgr_tmp1 <- new_mty_intgr[,samples_f2]
# snv_count_tmp1 <- snv_count[,samples_f2]
# data_cnv_tmp1 <- data_cnv_tmp[,samples_f2]
#
# print(samples_f2[2])

# dim(snv_count)
# length(samples)
# print(samples[ZZZ3])
# View(samples)
# dim(data_cnv_tmp)
# View(as.data.frame(samples))
# View(data_cnv_tmp)
# max(data_cnv_tmp)
# class(samples)
# rows_to_keep <- row.names(data_cnv_tmp) %in% samples
# class(rows_to_keep)
# data_cnv_tmp <- data_cnv_tmp[rows_to_keep, ]
# cols_to_keep<-colnames(data_cnv_tmp) %in% rownamesSamples


# class(cols_to_keep)
# data_cnv_samfilter<-data_cnv_tmp[,cols_to_keep]
# data_cnv_samfilter<-data_cnv_tmp[samples,]
# dim(data_cnv_tmp)
# max(samples)
# length(samples)
# dim(data_cnv_samfilter)



# View(as.data.frame(samples))
# View(data_cnv_tmp)
# snv_count_samfilter<-snv_count[samples,]
# snv_count_samfilter<-snv_count[,cols_to_keep] 
# dim(snv_count)
# data_cnv_samfilter<-t(data_cnv_samfilter)
# snv_count_samfilter<-t(snv_count_samfilter)

# ### Keep patient samples with four genomic profiles

# View(data_cnv)

# data_cnv_samfilter<- data_cnv[,samples]

# dim(data_cnv_samfilter)


print(samples_f2)

print(colnames(new_exp_intgr))

expSample_to_keep<-rownames(new_exp_intgr) %in% samples_f2
mytSample_to_keep<-rownames(new_mty_intgr) %in% samples_f2
snvSample_to_keep<-colnames(snv_count) %in% samples_f2 
dncSample_to_keep<-colnames(data_cnv_tmp) %in% samples_f2
print(samples_f2)
#View(snv_count)
#View(data_cnv_tmp)

#View(new_exp_intgr)
# length(expSample_to_keep)
# print(expSample_to_keep)
# print(dncSample_to_keep)
#print(expSample_to_keep)
###########################################
new_exp_intgr_tmp01<-new_exp_intgr[expSample_to_keep,]

new_mty_intgr_tmp01<-new_mty_intgr[mytSample_to_keep,]
# dim(new_mty_intgr_tmp01)
###########################################

snv_count_tmp01<-snv_count[,snvSample_to_keep]
data_cnv_tmp01<-data_cnv_tmp[,dncSample_to_keep]

# dim(snv_count_tmp01)
# dim(data_cnv_tmp01)

snv_count_tmp01<-t(snv_count_tmp01)
data_cnv_tmp01<-t(data_cnv_tmp01)


#View(NT_remove_name)
#class(clinicalInfo) data.frame
new_clinicalInfo<-clinicalInfo[samples_f2,]
View(as.data.frame(samples_f2))
#View(new_clinicalInfo)
#length(samples_f2)
#dim(new_clinicalInfo)
#print(indices)
#length(indices)

if (length(indices) > 0) {
  new_clinicalInfo <- clinicalInfo[-indices,]
} else {
  print("The length of NT indices is 0")  # 或者进行其他操作
}

#View(indices)
#View(new_clinicalInfo)
saveRDS(new_clinicalInfo, paste0("../data/tcga_data_processed/",cancer_type,"_clinical_info_test01.RData"))

View(new_clinicalInfo)
