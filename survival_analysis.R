###############data analysis#################
##### Patient stratification
options(stringsAsFactors = F)

library(NbClust)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(maftools)
library(htmlwidgets)

library(pdftools)
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
library(NbClust)

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



dataset_TCGA <- c("BLCA","BRCA", "COAD", "ESCA", "KICH", "KIRC","KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD","READ", "SKCM", "STAD", "THCA", "THYM", "UCEC")


for (cancer_type in dataset_TCGA){
#cancer_type<-"SKCM"
print("----------------------------------------------------------------------------------------------------")
print(cancer_type)


setwd("/users/zhanglia/ondemand/testSpace/code_and_data21/code")

print(cancer_type)
#samples <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))


#length(samples_f2)
cmtScores <- read.csv(paste0("../data/python_related/result/community/combine_",cancer_type,"_score_profile_test01.csv"), check.names = F, header = F)

setwd("/projappl/project_2010541/code/")
new_clinicalInfo<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_clinical_info_test01.RData"))
#print(melanet_cmt)
#dim(samples_f2)
#print(samples_f2)
#View(clinicalInfo)
#View(new_clinicalInfo)
clinicalInfo_tmp<-new_clinicalInfo
#View(clinicalInfo_tmp)

samples_f2<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))


#View(new_clinicalInfo)
row.names(cmtScores) <- samples_f2
colnames(cmtScores) <- paste0("cmt", 1:ncol(cmtScores))
saveRDS(cmtScores, paste0("../data/",cancer_type,"_community_scores.RData"))


# #############################k-means#####################################

n_rows<-nrow(cmtScores)
#print(n_rows)
normalized_data<-scale(cmtScores)

samples_name<-readRDS(paste0("../data/tcga_data_processed/",cancer_type,"_samples.RData"))
rownames(cmtScores)<-samples_name
#View(cmtScores)
rownames(normalized_data)<-samples_name
#View(normalized_data)
# set.seed(42)  # 为了可重复性设置随机种子
# 
# tsne_results <- Rtsne(normalized_data, dims = 2, perplexity = ceiling(sqrt(n_rows)), check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
# 
# 
# tsne_data <- as.data.frame(tsne_results$Y)
# 
# # 使用 ggplot2 进行可视化
# 
# ggplot(tsne_data, aes(x = V1, y = V2)) +
#   geom_point(alpha = 0.6) +  # 点的透明度
#   theme_minimal() +
#   labs(title = "t-SNE Visualization of cmtScores")
set.seed(123)


#tsne_results01 <- Rtsne(normalized_data, dims = 3, perplexity = 30,max_iter=10000, check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
#View(data.frame(tail(tsne_results01$itercosts,1)))
#############################################################################################
run_tsne <- function(data, perplexity) {
  tsne_result <- Rtsne(data, dims = 3, perplexity = perplexity, max_iter = 1000, verbose = FALSE)
  # 返回最后一个itercosts的值
  mean_iter_cost <- mean(tsne_result$itercosts)
  return(mean_iter_cost)
}

###############################################################################################
find_optimal_perplexity <- function(data, start = 20, end = n_rows) {
  best_perplexity <- start
  min_error <- Inf

  for (p in seq(start, end, by = 1)) {
    #cat("Testing perplexity:", p, "\n")

    # 使用 tryCatch 来处理可能出现的错误
    mean_iter_cost <- tryCatch({
      run_tsne(data, p)
    }, error = function(e) {
      #cat("Error at perplexity", p, ": ", e$message, "\n")
      return(Inf)  # 当错误发生时返回Inf
    })

    if (is.infinite(mean_iter_cost)) {
      # 如果出现错误，终止循环
      break
    }

    #cat("Last iter cost:", mean_iter_cost, "\n")

    if (mean_iter_cost < min_error) {
      min_error <- mean_iter_cost
      best_perplexity <- p
    }
  }

  return(best_perplexity)
}


#########################################################################################
# find_optimal_perplexity <- function(data, start = 20, end = n_rows) {
#   best_perplexity <- start
#   min_error <- Inf
#   error_occurred <- FALSE
#   
#   for (p in seq(start, end, by = 2)) {
#     # 使用 tryCatch 来处理可能出现的错误
#     mean_iter_cost <- tryCatch({
#       run_tsne(data, p)
#     }, error = function(e) {
#       error_occurred <<- TRUE  # 设置错误标志
#       return(Inf)  # 当错误发生时返回 Inf
#     })
#     
#     if (error_occurred) {
#       # 如果出现错误，退出循环
#       break
#     }
#     
#     if (mean_iter_cost < min_error) {
#       min_error <- mean_iter_cost
#       best_perplexity <- p
#     }
#   }
#   
#   # 检查是否发生错误
#   if (error_occurred) {
#     cat("Error occurred during t-SNE computation. Exiting loop.\n")
#   } else {
#     cat("Optimal perplexity found:", best_perplexity, "\n")
#   }
#   
#   return(best_perplexity)
# }


##############################################################################################
b_perplexity<-find_optimal_perplexity(normalized_data,start = 20, end = n_rows)

#b_perplexity<-154
#print(b_perplexity)


########################################################################################################

tsne_results01 <- Rtsne(normalized_data, dims = 3, perplexity = b_perplexity,max_iter=5000, check_duplicates = FALSE, pca = FALSE, verbose = FALSE)


##################################################################################################
tsne_data01 <- as.data.frame(tsne_results01$Y)
tsne_data_knn<-as.data.frame(tsne_results01$Y)

rownames(tsne_data_knn)<-samples_name

tsne_data_knn_test<-tsne_data_knn

tsne_data_numeric <- tsne_data_knn[, c("V1", "V2", "V3")]

# 确保所有列都是数值型
tsne_data_numeric[] <- lapply(tsne_data_numeric, as.numeric)

tsne_matrix<-as.matrix(tsne_data_numeric)




#########################################

# gap_stat <- clusGap(tsne_matrix, FUN = kmeans, nstart = 500, K.max = 15, B = n_rows,verbose=FALSE)
# 
# # 打印Gap Statistic结果
# print(gap_stat, method = "firstmax")
# 
# # 绘制Gap Statistic图
# fviz_gap_stat(gap_stat)
# View(gap_stat)
# 
# 
# gap_values <- gap_stat$Tab[, "gap"]
# max_gap_index <- which.max(gap_values)  # 找到最大Gap值的索引
# 
# # 提取最优聚类数
# optimal_clusters <- max_gap_index
# 
# # 打印最优聚类数
# print(paste("Optimal number of clusters based on Gap statistic:", optimal_clusters))



###############################################
nc_result <- NbClust(tsne_matrix, distance = "euclidean", min.nc = 2, max.nc = 12, method = "complete", index = "all")
best_clusters_kmeans <- nc_result$Best.nc[1, -1]  # 获取所有指数的推荐集群数量
cluster_counts_kmeans <- table(best_clusters_kmeans)  # 计算每个集群数量的出现频率
print(cluster_counts_kmeans)
View(cluster_counts_kmeans)
freq_cluster_means <- as.integer(names(which.max(cluster_counts_kmeans)))  # 找到频率最高的集群数量
most_freq_cluster<-freq_cluster_means

#View(cluster_counts_kmeans)

print(most_freq_cluster)

######################################################################################################################################################
k_num<-most_freq_cluster

kmeans_result <- kmeans(tsne_data_knn, centers = k_num)  # 假定我们有 4 个聚类中心most_freq_cluster
#View(as.data.frame(kmeans_result))
# 将聚类结果添加到数据框中

tsne_data_knn$cluster <- as.factor(kmeans_result$cluster)

centroids <- kmeans_result$centers
centroids_df <- as.data.frame(centroids)
centroids_df$cluster <- as.factor(1:nrow(centroids))  # 为每个中心点创建一个聚类标签



#View(tsne_data_knn)
#write.csv(tsne_data_knn,paste0("../figure/",cancer_type,"_tableClusters.csv"))
#View(tsne_data_knn)
# 使用 plot_ly 创建 3D 散点图来展示聚类结果
plot_kmeans_3d <- plot_ly(data = tsne_data_knn, x = ~V1, y = ~V2, z = ~V3,
                          type = 'scatter3d', mode = 'markers',
                          marker = list(size = 5, opacity = 0.6),
                          color = ~cluster,  # 根据聚类结果上色
                          hoverinfo = 'text+name')  # 显示行名称和聚类

# 添加聚类中心到图形
plot_kmeans_3d <- plot_kmeans_3d %>%
  add_trace(data = centroids_df, x = ~V1, y = ~V2, z = ~V3,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 6, color = 'black', symbol = 'star'),  # 更大更显眼的中心点
            name = 'Centroids') %>%
  layout(title = "3D Visualization of k-Means Clustering with Centroids",
         scene = list(xaxis = list(title = 'Component 1'),
                      yaxis = list(title = 'Component 2'),
                      zaxis = list(title = 'Component 3')))

# 显示图形
#plot_kmeans_3d
saveWidget(plot_kmeans_3d, file = paste0("/users/zhanglia/ondemand/testSpace/code_and_data21/figure",cancer_type,"_plot_kmeans_3d.html"))

# View(tsne_data_knn)

# -----------------------------------------------
# install.packages("dbscan")
# library(dbscan)
#
# library(rgl)
#
# # 进行 DBSCAN 聚类
# dbscan_result <- dbscan(tsne_dbscan, eps = 5, MinPts = 5)
#
# # 获取聚类结果
# clusters <- dbscan_result$cluster
#
# print(clusters)
#
# # 创建3D散点图
# plot3d(tsne_dbscan, col = clusters, size = 5)
#
# # 添加标签
# text3d(tsne_dbscan, texts = rownames(tsne_dbscan), adj = c(-0.5,0), cex = 0.7)
#
# # 设置3D场景
# rglwidget()  # 打开3D视图窗口，可以旋转和缩放
# ##############################################################
# print(cancer)
#
# best_clusters <- nc$Best.nc[1, ]  # 获取所有指数的推荐集群数量
#
# # 计算每个集群数量的出现频率
# cluster_counts <- table(best_clusters)

# View(cluster_counts)
# # 找到频率最高的集群数量
# most_frequent_cluster <- as.numeric(names(which.max(cluster_counts)))
#
# # 打印最频繁推荐的集群数量
# print(most_frequent_cluster)
#
# pdf(paste0("./figure/",cancer_type,"_best_number_of_clusters01.pdf", width = 7, height = 7))
# ggplot(data.frame(cluster = factor(nc$Best.nc[1,])), aes(x = cluster)) + 
#   geom_bar(stat = "count", fill = "#C1BFBF") + 
#   labs(x = "Number of clusters", y = "Number of criteria", title = "Number of clusters chosen by 26 criteria") + 
#   theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) + 
#   scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14))
# dev.off()

#View(clinicalInfo)
### Community scores of clustered patients
tumorType <- clinicalInfo_tmp[,"shortLetterCode"]

# View(clinicalInfo_tmp)
#View(clinicalInfo_tmp)
# tumorStage <- clinicalInfo_tmp[,"tumor_stage"]

mty_colData<-readRDS(paste0("../data/tcga_data/",cancer_type,"_mty_colData.RData"))
clinicalInfo <- mty_colData
clinicalInfo<-clinicalInfo[!duplicated(clinicalInfo$sample),]


clinicalInfo_col<-as.data.frame(colnames(clinicalInfo))



matching_columns <- grep("ajcc_pathologic_stage|ajcc_pathologic_t|figo_stage|masaoka_stage", colnames(clinicalInfo_tmp), value = TRUE)
survivalInfo <- clinicalInfo_tmp[,c("shortLetterCode", "vital_status","days_to_death","days_to_last_follow_up",matching_columns)]
# View(matching_columns)
# print(matching_columns)
#View(clinicalInfo_tmp)

#View(survivalInfo)
tsne_data_knn$SampleID<-rownames(tsne_data_knn)
#View(tsne_data_knn)
#View(survivalInfo)
survivalInfo_tmp_test<-survivalInfo

new_column_names <- setNames(rep("ajcc_pathologic_stage", length(matching_columns)), matching_columns)

colnames(survivalInfo_tmp_test)[colnames(survivalInfo_tmp_test) %in% names(new_column_names)] <- new_column_names
#View(survivalInfo_tmp_test)
survivalInfo_tmp_test$PatientID <- rownames(survivalInfo)
#View(as.data.frame(substr((survivalInfo_tmp_test$PatientID),1,15)))
survivalInfo_tmp_test$PatientID<-substr((survivalInfo_tmp_test$PatientID),1,16)
#View(survivalInfo_tmp_test)
#View(tsne_data01)
combined_data <- merge(survivalInfo_tmp_test, tsne_data_knn, by.x = "PatientID", by.y = "SampleID")
#View(survivalInfo_tmp_test)
#View(tsne_data_knn)
combined_data$V1<-NULL
combined_data$V2<-NULL
combined_data$V3<-NULL

#colnames(combined_data)[colnames(data) == "cluster"] <- "hc"
#View(combined_data)
#View(survivalInfo_tmp_test)
###################################Section 5################################################################
##### Survival analysis
options(stringsAsFactors = F)
library(TCGAbiolinks)
#cancer_type<-"BLCA"
setwd("/projappl/project_2010541/code/")
print(cancer_type)
#setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19")
#clinicalInfo <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_new_clinicalInfo_test06.RData"))



survivalInfo_tmp_test_record<-survivalInfo_tmp_test

# View(combined_data)

#View(survivalInfo_tmp_test_record)
#View(samplePartition)
#View(survivalInfo_tmp)
#View(survivalInfo_tmp_test)
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage 0")] <- 0
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage i")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ia")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ib")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ic")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ii")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iia")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iib")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iic")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iii")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiia")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiib")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiic")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iv")] <- 4
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iva")] <- 4
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ivb")] <- 4
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ivc")] <- 4
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "not reported")] <- "NA"
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "i/ii nos")] <- "NA"
survivalInfo_tmp_test$ajcc_pathologic_stage[which(is.na(survivalInfo_tmp_test$ajcc_pathologic_stage))] <- "NA"


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TCGAanalyze_survival(combined_data, clusterCol = "cluster", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565","#FB9A99", 
#                                                                       "#CAB2D6", "#FF7F00", "#C2C2F0", "#B0E2FF"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update01.pdf"), conf.int = F, width = 7, height = 7)

color_palette <- rainbow_hcl(k_num)

# 查看生成的颜色
print(color_palette)

# 调用 TCGAanalyze_survival 函数

setwd("/users/zhanglia/ondemand/testSpace/code_and_data21/code/")

View(combined_data)

try({
  TCGAanalyze_survival(
  combined_data, 
  clusterCol = "cluster", 
  color = color_palette, 
  filename = paste0("../figure/surv_analysis/", cancer_type, "_survival_analysis01_update06.pdf"), 
  conf.int = FALSE, 
  width = 7, 
  height = 7
)},silent = TRUE)




#####################################



try({
  TCGAanalyze_survival(survivalInfo_tmp_test, clusterCol = "ajcc_pathologic_stage", filename = paste0("../figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage_update06.pdf"), conf.int = F, width = 7, height = 7)
},silent = TRUE)


ajcc_pathologic_stage_count <- length(unique(survivalInfo_tmp_test$ajcc_pathologic_stage))
#print(paste("Number of unique ajcc_pathologic_stage categories:", ajcc_pathologic_stage_count))

#####################################


try({
TCGAanalyze_survival(combined_data, clusterCol = "shortLetterCode", filename = paste0("../figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update01.pdf"), conf.int = F, width = 7, height = 7)
},silent = TRUE)


#####################################


#################################################################################









extract_p_value <- function(pdf_path) {
  # 读取 PDF 文件的文本内容
  pdf_text_content <- pdf_text(pdf_path)
  full_text <- paste(pdf_text_content, collapse = " ")
  
  # 定义正则表达式模式
  p_value_pattern_equal <- "p = ([0-9.]+)"
  p_value_pattern_less <- "p < ([0-9.]+)"
  
  # 尝试匹配 "p = 0.xxxx" 格式
  p_value_match_equal <- regexpr(p_value_pattern_equal, full_text)
  
  if (p_value_match_equal[1] != -1) {
    # 匹配 "p = 0.xxxx" 格式
    p_value <- regmatches(full_text, p_value_match_equal)
    p_value <- sub("p = ", "", p_value)
    return(as.numeric(p_value))
  } else {
    # 尝试匹配 "p < 0.xxxx" 格式
    p_value_match_less <- regexpr(p_value_pattern_less, full_text)
    if (p_value_match_less[1] != -1) {
      p_value <- regmatches(full_text, p_value_match_less)
      p_value <- sub("p < ", "", p_value)
      return(as.numeric(p_value))
    } else {
      return(NA)
    }
  }
}

# 示例：提取 PDF 文件中的 p 值
pdf_path_cluster <- paste0("../figure/surv_analysis/", cancer_type, "_survival_analysis01_update01.pdf")
p_value_cluster <- extract_p_value(pdf_path_cluster)
print(paste("p_value_cluster:", p_value_cluster))

pdf_path_stage <- paste0("../figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage_update03.pdf")
p_value_stage <- extract_p_value(pdf_path_stage)
print(paste("p_value_stage:", p_value_stage))



results <- data.frame(
  cancer_type = cancer_type,
  num_clusters = k_num,
  p_value_cluster = p_value_cluster,
  ajcc_pathologic_stage_count = ajcc_pathologic_stage_count,
  p_value_stage = p_value_stage
 
)

write.csv(results, paste0("../figure/surv_analysis/", cancer_type, "_survival_analysis_results.csv"), row.names = FALSE)
#print(results)


}




