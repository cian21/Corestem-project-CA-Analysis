# Corestem project - CA Analysis.proj
# 20201102 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for CA data analysis

#### 1. source Loading ####
load("./Data/CA analysis files.image")
source("./Script/Functions/CA analysis - functions.R")

#### 2. cpm cleaning ####
row.names(expr_file) = expr_file$gene_id
expr_file$gene_id = NULL

second_samples_num = which(substr(colnames(expr_file),4,4) =="2")
expr_file = expr_file[,-second_samples_num] #KYS 2nd removed

group_name = substr(colnames(expr_file),1,4)
group_count = make_group_count(expr_file)

save.image("./Data/CA analysis expr.image")

expr_file_full_cpm = cpm_calculation(expr_file)
expr_file_full_cpm = anno_DEG(as.data.frame(expr_file_full_cpm))

expr_file_030 = expr_file[,which(group_count %in% c("d0", "d30") == T)]
expr_file_030_cpm = cpm_calculation(expr_file_030)
expr_file_030_cpm = anno_DEG(as.data.frame(expr_file_030_cpm))

expr_file_090 = expr_file[,which(group_count %in% c("d0", "d90") == T)]
expr_file_090_cpm = cpm_calculation(expr_file_090)
expr_file_090_cpm = anno_DEG(as.data.frame(expr_file_090_cpm))

group_name_030_cpm = substr(colnames(expr_file_030_cpm),1,4)

expr_file_030_SAOA_cpm = expr_file_030_cpm[, c(group_name_030_cpm %in% c("BMS_","JYJ_","KYS_","LSJ_"))]
expr_file_030_SAOA_cpm$entrizid = expr_file_030_cpm$entrizid
expr_file_030_SAOA_cpm$Symbol = expr_file_030_cpm$Symbol

save(expr_file_full_cpm, file = "./Data/expr_file_full_cpm.Rdata")
save(expr_file_030_cpm, file = "./Data/expr_file_030_cpm.Rdata")
save(expr_file_090_cpm, file = "./Data/expr_file_090_cpm.Rdata")
save(expr_file_030_SAOA_cpm, file = "./Data/expr_file_030_SAOA_cpm.Rdata")
