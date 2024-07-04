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

library(dplyr)
library(lubridate)
library(tidyverse)

############
library(Rtsne)
library(ggplot2)
# install.packages("plotly")
library(plotly)
library(colorspace)
#####################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# if (!requireNamespace("EDASeq", quietly = TRUE))
#   BiocManager::install("EDASeq")


library(EDASeq)
# ####################################################workround###############################################################
# https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/627
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

#############################################################################################################################
#### Collect gene expression data#################################################################
#"SKCM"
dataset_TCGA <- c("BLCA", "BRCA", "COAD", "ESCA", "KICH", "KIRC",
                  "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
                   "READ", "STAD", "SKCM","THCA", "THYM", "UCEC")

# tried_groups<-c("BLCA",  "COAD",  "KICH", "KIRC", "KIRP","SKCM","LIHC","LUAD","LUSC","UCEC","READ","PAAD","BRCA")
# still_groups<-c( )
# not_installed<-c("ESCA","THYM","THCA","STAD")
#
# #print(cancer_type)
# in_pool<-"THYM"
# doing_groups<-c("THYM")
##############################################################################
for (cancer_type in dataset_TCGA){
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data21/code/")
#cancer_type<-"SKCM"


# #####################################section2 singlass for communities##########################################################

filename1 <- paste0("../data/tcga_data/", cancer_type, "_maf_test06.RData")


maf<-readRDS(filename1)



new_exp_intgr<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all.RData"))
new_mty_intgr<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.RData"))

data_temp<-readRDS(paste0("../data/tcga_data/",cancer_type,"_data_temp_savage_test.csv"))


therapy<-readRDS(paste0("../data/tcga_data/",cancer_type,"_therapy_test.RData"))
radiation<-readRDS(paste0("../data/tcga_data/",cancer_type,"_radiation_test.RData"))


exp_intgr<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all01.RData"))



mty_mat_trial<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all01.RData"))

samples<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))

sample_f2_check<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))


samples_f2<-sample_f2_check

new_exp_intgr<-readRDS( paste0("../data/tcga_data_processed/",cancer_type,"_exp_intgr_all.RData"))

new_mty_intgr<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.RData"))



#exp_intgr <- exp_intgr[,samples]#?
#mty_mat_trial <- mty_mat_trial[,samples]#?
#rownamesSamples<-intersect(rownames(exp_intgr),rownames(mty_mat_trial))

#nrow(rownamesSamples)

#exp_intgr_trial<-exp_intgr[rownamesSamples,]
#mty_mat_trial<-mty_mat_trial[rownamesSamples,]

#############################################################################################################################################################

maf <- maf[,c("Tumor_Sample_Barcode","Hugo_Symbol","Gene","Variant_Classification")]

rnames <- unique(maf$Hugo_Symbol)
cnames <- unique(maf$Tumor_Sample_Barcode)
#View(cnames)
snv_count <- matrix(data = 0, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames,cnames))

# Calculate the frequency of genes' variants
for(i in 1:nrow(maf)){
  rname <- maf[i,]$Hugo_Symbol
  cname <- maf[i,]$Tumor_Sample_Barcode
  snv_count[rname, cname] <- snv_count[rname,cname] + 1
}
#dim(snv_count)
#https://www.notion.so/57782b3efe144e91a93ae602f85a1cbd

#View(snv_count)
colnames(snv_count) <- substr(colnames(snv_count), 1, 16)

#View(snv_count)
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

###############################################################################################!!!!!!!!!!!!!!!!!!!!
col_data_cnv<-ncol(data_cnv_tmp)
#View(as.data.frame(col_data_cnv))#146
##########################################################################################################
# Convert Ensemble ID to corresponding gene symbols
# Using biomaRt for gene ID conversion will cause some loss
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

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


expSample_to_keep<-rownames(new_exp_intgr) %in% samples_f2
mytSample_to_keep<-rownames(new_mty_intgr) %in% samples_f2
snvSample_to_keep<-colnames(snv_count) %in% samples_f2
dncSample_to_keep<-colnames(data_cnv_tmp) %in% samples_f2


new_clinicalInfo<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_clinical_info_test01.RData"))

#View(new_clinicalInfo)

###############################################################################################################################################################################

##### The melanoma network and network community detection
options(stringsAsFactors = F)
library(igraph)
library(visNetwork)

library(OmnipathR)#Get interactions data
library(igraph)#network analysis
library(ggraph)#network visualization
library(RColorBrewer)#Network visualization, heat map drawing

library(factoextra)#principal component analysis
library(FactoMineR)#principal component analysis

library(DESeq2)#differential expression analysis
#library(SANTA)#Node Score Calculation

library(ComplexHeatmap)#heatmap drawing
library(stringr)#regular expression

library(randomForest)#Random Forest
library(caret)#Cross-validation
library(pROC)#AUC

library(clusterProfiler)#enrichment analysis
library(org.Hs.eg.db)
library(enrichplot)#enrichment analysis
library(msigdbr)#enrichment analysis
library(dplyr)

library(survival)#survival analysis
library(survminer)#survival analysis
library(glmnet)# cox regression
library(timeROC)#AUC

library(reshape2)#ggplot2 drawing

library(coin)#permutation test

library(networkD3)#sankey diagram

library(janitor)# clean data
library(here)# home directory setting
library(ggplot2) # draw figures
library(tidyr) # data processing
if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
library(visNetwork)
library(igraph)
### The melanoma network
#cancer_type<-"BLCA"
print(cancer_type)
setwd("C:/Users/gklizh/Documents/Workspace/fan_project/code_data/")
edgetest<-read.csv(paste0("Raw data/network_merge/edge/", cancer_type, ".csv"))
vertextest<-read.csv(paste0("Raw data/network_merge/vertex/", cancer_type, ".csv"))
interactions_E <- edgetest
interactions_V <- vertextest
graph_comp <- graph_from_data_frame(interactions_E, directed = TRUE, vertices = interactions_V)
graph_comp <- decompose(graph_comp, min.vertices = 15)[[1]]
graph_comp <- set_vertex_attr(graph_comp, "degree", index = V(graph_comp)$name, degree(graph_comp))

print(class(degree(graph_comp)))
print(length(degree(graph_comp)))
print(class(V(graph_comp)$log2FoldChange))
print(length(V(graph_comp)$log2FoldChange))
#print(graph_comp)
#################################################################

#setwd("/projappl/project_2010541/code/")
# spg01 <- list() # List of spinglass simulations
# spg_mod01 <- numeric() # List of modularity simulations
# count01<-1
#for (k in 1:100){
#  print(Sys.time())
#  print(k)
#  spins_value <- sample(20:30, 1)
#  spg01[[k]] <- cluster_spinglass(graph_comp,
#                                  spins_value = 25,
#                                  weights = E(graph_comp)$weight,
#                                  implementation = "neg")
#
#  spg_mod01[k] <- spg01[[k]]$modularity
#  #count = count + 1
#  print(Sys.time())
#}
# ##################################################################################
# spins_value <- sample(20:30, 1)
# print(spins_value)
# for (k in 1:21) {
#   print(Sys.time())
#   print(k)
# 
#   # 在20到30之间随机选择一个spins值
#   #spins_value <- sample(20:30, 1)
# 
#   # 运行社区检测
#   spg_result <- cluster_spinglass(graph_comp,
#                                   spins = spins_value,
#                                   weights = E(graph_comp)$weight,
#                                   implementation = "neg")
# 
#   # 提取membership并过滤社区
#   membership <- membership(spg_result)
#   filtered_membership <- membership
#   communities_sizes <- sizes(spg_result)
# 
#   for (i in seq_along(communities_sizes)) {
#     if (communities_sizes[i] < 5) {
#       filtered_membership[membership == i] <- NA
#     }
#   }
# 
#   # 删除NA值
#   filtered_membership <- filtered_membership[!is.na(filtered_membership)]
# 
#   # 如果过滤后的社区不为空，则保存结果
#   if (length(filtered_membership) > 0) {
#     spg01[[k]] <- make_clusters(graph_comp, membership = filtered_membership)
#     spg_mod01[k] <- modularity(graph_comp, membership = filtered_membership)
#   } else {
#     spg_mod01[k] <- -Inf # 设置为负无穷大，以便在选择最大模块度时被排除
#   }
# 
#   #print(Sys.time())
# }

setwd("C:/Users/gklizh/Documents/Workspace/code_and_data21/code/")
graph_comp_undirected <- as.undirected(graph_comp, mode = "collapse")
absolute_weights <- abs(E(graph_comp_undirected)$weight)
all_data <- list()
optimize_modularity <- function(graph_comp_undirected, res_low = 0.1, res_high = 5.0, step = 0.01) {
  best_modularity <- -Inf
  best_community <- NULL
  best_resolution <- NULL
  
  for (resolution in seq(res_low, res_high, by = step)) {
    community <- cluster_louvain(graph_comp_undirected,weights = absolute_weights, resolution = resolution)
    modularity_score <- modularity(community)
    
    all_data[[as.character(resolution)]] <- list(modularity = modularity_score, community = community)
    
    if (modularity_score > best_modularity) {
      best_modularity <- modularity_score
      best_community <- community
      best_resolution <- resolution
    }
  }
  
  saveRDS(all_data, paste0("../data/spinglass/",cancer_type,"_all_Modularity_Data.RData"))
  list(
    community = best_community,
    modularity = best_modularity,
    resolution = best_resolution
  )
}

# 对最大连通图应用优化函数
results <- optimize_modularity(graph_comp_undirected)

# 打印最优结果
print(results$modularity)#$modularity [1] 0.9012516

print(results$resolution)#resolution [1] 0.1

print(results$community)
View(results)
saveRDS(results, paste0("../data/spinglass/",cancer_type,"_best_Modularity_Results.RData"))
# #################################################################################
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data21/code/")

# max_index01 <- which.max(spg_mod01)
# length(spg_mod01)
# best_spg01 <- spg01[[max_index01]]
# best_spg01_copy<-best_spg01
# #print(best_spg01_copy)
# saveRDS(spg01, paste0("../data/spinglass",cancer_type,"_spg01.RData"))
# saveRDS(spg_mod01, paste0("../data/spinglass",cancer_type,"_spg_mod01.RData"))
# best_spg01<-as.list(communities(best_spg01))
# saveRDS(best_spg01, paste0("../data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))
# 
# print(best_spg01)
#class(best_spg01) array
genes <- as.character(unlist(vertex.attributes(graph_comp_undirected)))
#######################################################################################################
new_exp_intgr_tmp01<-new_exp_intgr[expSample_to_keep,]
#View(as.data.frame(expSample_to_keep))
new_mty_intgr_tmp01<-new_mty_intgr[mytSample_to_keep,]
snv_count_tmp01<-snv_count[,snvSample_to_keep]
data_cnv_tmp01<-data_cnv_tmp[,dncSample_to_keep]

# View(new_exp_intgr_tmp01)

y1<-which(colnames(new_exp_intgr_tmp01) %in% genes)
y2<-which(colnames(new_mty_intgr_tmp01) %in% genes)

# print(y2)


y3<-which(colnames(snv_count_tmp01) %in% genes)
y4<-which(colnames(data_cnv_tmp01) %in% genes)


new_exp_intgr_tmp02<-new_exp_intgr_tmp01[,y1]
new_mty_intgr_tmp02<-new_mty_intgr_tmp01[,y2]
snv_count_tmp02<-snv_count_tmp01[,y3]
data_cnv_tmp02<-data_cnv_tmp01[,y4]
# dim(new_exp_intgr_tmp02)
# dim(new_mty_intgr_tmp02)
# dim(snv_count_tmp02)
# dim(data_cnv_tmp02)

# print(y4)

# write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
# ###################################################################################################


#View(data_cnv_tmp02)
#print(data_cnv_tmp02)

# ###################################################################################################
# #################################section3 community mapping to TCGA genes###############################################

options(stringsAsFactors = F)
#setwd("/projappl/project_2010541/code/")
#exp_intgr <- readRDS("./data/tcga_data_processed/exp_intgr.RData")
#mty_intgr <- readRDS("./data/tcga_data_processed/mty_intgr.RData")
#error solved:mty_intgr <- readRDS("./data/tcga_data_processed/exp_intgr.RData")

hfeat <- colnames(new_exp_intgr_tmp02) # The features of gene expression data
#print(length(hfeat))
mfeat <- colnames(new_mty_intgr_tmp02) # The features of DNA methylation data
#print(length(mfeat))

#print(nrow(new_exp_intgr_tmp02))
#melanet_cmt<-readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))
melanet_cmt<-(results$community)

selected_features <- matrix(ncol = 3, nrow=2*ncol(new_exp_intgr_tmp02)+2) # Matrix of 3 columns; column1: community, column2: genomic type, column3: mapped component
#View(selected_features)
len <- 0
j <- 1
indx <- NULL
for(i in 1:length(melanet_cmt)){
  cmt = melanet_cmt[[i]]
  # Start mapping
  # Gene expression
  j = j + length(indx)
  indx = NULL
  indx = which(hfeat %in% cmt)
  if (length(indx) != 0){
    len = len +length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 1
    selected_features[j:len,3] = indx
  }

  # DNA methylation
  j = j+length(indx)
  indx = NULL
  indx = which(mfeat %in% cmt)
  if (length(indx) != 0){
    len = len + length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 2
    selected_features[j:len,3] = indx
  }
}
#print(cancer_type)
selected_features <- na.omit(selected_features)
write.csv(selected_features, paste0("../data/python_related/data/",cancer_type,"_selected_features01.csv"), row.names = F)
write.csv(colnames(new_exp_intgr_tmp02), paste0("../data/python_related/data/",cancer_type,"_exp_feature_names01.csv"), row.names = F)
write.csv(colnames(new_mty_intgr_tmp02), paste0("../data/python_related/data/",cancer_type,"_mty_feature_names01.csv"), row.names = F)
write.csv(new_exp_intgr_tmp02, paste0("../data/python_related/data/",cancer_type,"_exp_intgr01.csv"))
write.csv(new_mty_intgr_tmp02, paste0("../data/python_related/data/",cancer_type,"_mty_intgr01.csv"))
#View(new_exp_intgr_tmp02)
#check_new_mty_intgr_tmp02<-read.csv(paste0("../data/python_related/data/",cancer_type,"_mty_intgr01.csv"))
#print(cancer_type)
#nrow(check_new_mty_intgr_tmp02)
##print(best_spg01)
##print(cancer_type)
#print(spins_value)
}
