# Corestem project - CA Analysis.proj
# 20201102 Sehwan Chun at Corestem, Inc.
# function.R

#### 1. Library Loading ####
packs = c("edgeR","org.Hs.eg.db","data.table","lmerTest", "readxl", "ggpubr","writexl")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. simple DEG ####
anno_DEG = function(DGE_table){
    entrezid = row.names(DGE_table)
    cols = c("ENTREZID", "SYMBOL")
    geneList = biomaRt::select(org.Hs.eg.db,
                               keys = entrezid,
                               columns = cols,
                               keytype = "ENTREZID")
    DGE_table$entrizid = entrezid
    
    #File Cleaning
    for (i in 1:nrow(DGE_table)){
        if(is.na(geneList[i,2]) == F){
            row.names(DGE_table)[i] = geneList[i,2]
        }
    }
    DGE_table$Symbol = row.names(DGE_table)
    return(DGE_table)
}
d030 = function(expr_file, group_count, pairwise){
    expr_file_d030 = expr_file[group_count == "d0" | group_count == "d30"]
    group_count_d030 = group_count[group_count == "d0" | group_count == "d30"]
    
    if(pairwise == T){
        pairwised_number = c()
        for(i in 1:(length(group_count_d030) - 1)){
            if(group_count_d030[i] == "d0" & group_count_d030[i+1] == "d30"){
                pairwised_number = c(pairwised_number,i)
            }
        }
        pairwised_number = sort(c(pairwised_number,pairwised_number+1))
        expr_file_d030 = expr_file_d030[pairwised_number]
        group_count_d030 = group_count_d030[pairwised_number]
    }
    
    group_frame_030 = data.frame(days = group_count_d030)
    design_matrix_030 = model.matrix(~days, data = group_frame_030) #d7 == 1
    
    DGE = DGEList(counts = expr_file_d030, group = group_count_d030)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE, design_matrix_030)
    
    DGETable = exactTest(DGE)$table
    DGETable$FDR = p.adjust(DGETable$PValue, method =  "fdr")
    DGETable = anno_DEG(DGETable)
    return(DGETable)
}
d030_SAOA = function(expr_file, group_count, pairwise){
    
    expr_file_d030 = expr_file[group_count == "d0" | group_count == "d30"]
    group_count_d030 = group_count[group_count == "d0" | group_count == "d30"]
    
    group_name = substr(colnames(expr_file_d030),1,4)
    
    expr_file_d030 = expr_file_d030[, c(group_name %in% c("BMS_","JYJ_","KYS_","LSJ_"))]
    group_count_d030 = group_count_d030[c(group_name %in% c("BMS_","JYJ_","KYS_","LSJ_"))]
    
    if(pairwise == T){
        pairwised_number = c()
        for(i in 1:(length(group_count_d030) - 1)){
            if(group_count_d030[i] == "d0" & group_count_d030[i+1] == "d30"){
                pairwised_number = c(pairwised_number,i)
            }
        }
        pairwised_number = sort(c(pairwised_number,pairwised_number+1))
        expr_file_d030 = expr_file_d030[pairwised_number]
        group_count_d030 = group_count_d030[pairwised_number]
    }
    
    group_frame_030 = data.frame(days = group_count_d030)
    design_matrix_030 = model.matrix(~days, data = group_frame_030) #d7 == 1
    
    DGE = DGEList(counts = expr_file_d030, group = group_count_d030)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE, design_matrix_030)
    
    DGETable = exactTest(DGE)$table
    DGETable$FDR = p.adjust(DGETable$PValue, method =  "fdr")
    DGETable = anno_DEG(DGETable)
    return(DGETable)
}
d090 = function(expr_file, group_count, pairwise){
    expr_file_d090 = expr_file[group_count == "d0" | group_count == "d90"]
    group_count_d090 = group_count[group_count == "d0" | group_count == "d90"]
    
    if(pairwise == T){
        pairwised_number = c()
        for(i in 1:(length(group_count_d090) - 1)){
            if(group_count_d090[i] == "d0" & group_count_d090[i+1] == "d90"){
                pairwised_number = c(pairwised_number,i)
            }
        }
        pairwised_number = sort(c(pairwised_number,pairwised_number+1))
        expr_file_d090 = expr_file_d090[pairwised_number]
        group_count_d090 = group_count_d090[pairwised_number]
    }
    
    group_frame_090 = data.frame(days = group_count_d090)
    design_matrix_090 = model.matrix(~days, data = group_frame_090) #d7 == 1
    
    DGE = DGEList(counts = expr_file_d090, group = group_count_d090)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE, design_matrix_090)
    
    DGETable = exactTest(DGE)$table
    DGETable$FDR = p.adjust(DGETable$PValue, method =  "fdr")
    DGETable = anno_DEG(DGETable)
    return(DGETable)
}
cpm_calculation = function(expr_file){
    DGE = DGEList(counts = expr_file)
    keep = filterByExpr(DGE)
    
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE)
    
    DGECPM = cpm(DGE, log = FALSE)
    
    return(DGECPM)
}
make_group_count = function(expr_file){
    group_count = colnames(expr_file)
    for (i in 1:length(group_count)){
        group_count[i] = gsub("\\.","_", group_count[i])
        group_count[i] = paste0("d",strsplit(group_count[i],"_")[[1]][2])
    }
    return(group_count)
}

#### 3. efficacy index linear mixed model ####
lmm_index_run = function(df, cpm){
    df = df[,c(1,3,4,5)]
    df = reshape2::melt(df)
    df = df[order(df$Samples),]
    row.names(df) = 1:nrow(df)
    
    df$Samples = as.factor(df$Samples)
    df$Days = ifelse(df$variable == "V2", 0, ifelse(df$variable == "V3", 30, 90))
    df_table = data.frame(genes = row.names(cpm))
    
    for (i in 1:nrow(cpm)){
        df$tmp = NA
        for (j in 1:(ncol(cpm)-2)){
            if(length(which(df$Samples == strsplit(colnames(cpm)[j], "_")[[1]][1] &
                            df$Days == strsplit(colnames(cpm)[j], "_")[[1]][2])) != 0){
                k = which(df$Samples == strsplit(colnames(cpm)[j], "_")[[1]][1] &
                              df$Days == strsplit(colnames(cpm)[j], "_")[[1]][2])
                df[k, "tmp"] = cpm[i,j]
            }
        }
        df_table[i,2] = summary(lmer(data = df, formula = value ~ tmp + (1|Samples)))$coefficients[2,1]
        df_table[i,3] = summary(lmer(data = df, formula = value ~ tmp + (1|Samples)))$coefficients[2,5]
        df_table[i,4] = as.numeric(cor.test(df$value, df$tmp)$estimate)
        df_table[i,5] = as.numeric(cor.test(df$value, df$tmp)$p.value)
        
    }
    df_table$FDR_lmm = p.adjust(df_table$V3)
    df_table$FDR_cor = p.adjust(df_table$V5)
    
    colnames(df_table)[2:5] = c("lmm_estimate","lmm_pvalue", "cor_estimate","cor_pvalue")
    
    return(df_table)
}

#### 4. efficacy test ####
efficacy_test_run  = function(efficacy_table){
    tmp_table  = data.frame(test = NA, pvalue = NA)
    tmp_table[1:12,1] = c("wilcox.test_v2_v3",
                      "wilcox.test_v2_v4",
                      "wilcox.test_v2-v1_v3-v2",
                      "var.test_v2_v3",
                      "var.test_v2_v4",
                      "var.test_v2_v1-v3-v2",
                      "t.test_v2_v3",
                      "t.test_v2_v4",
                      "t.test_v2-v1_v3-v2",
                      "pt.test_v2_v3",
                      "pt.test_v2_v4",
                      "pt.test_v2-v1_v3-v2")
    tmp_table[1,2] = wilcox.test(efficacy_table$V2, efficacy_table$V3)$p.value
    tmp_table[2,2] = wilcox.test(efficacy_table$V2, efficacy_table$V4)$p.value
    tmp_table[3,2] = wilcox.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                 (efficacy_table$V3 - efficacy_table$V2))$p.value
    
    tmp_table[4,2] = var.test(efficacy_table$V2, efficacy_table$V3)$p.value
    tmp_table[5,2] = var.test(efficacy_table$V2, efficacy_table$V4)$p.value
    tmp_table[6,2] = var.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                 (efficacy_table$V3 - efficacy_table$V2))$p.value
    
    if(tmp_table[4,2] < 0.05){
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                  var.equal = FALSE)$p.value
    }else{
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                  var.equal = TRUE)$p.value
    }
    if(tmp_table[4,2] < 0.05){
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[4,2] < 0.05){
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                  (efficacy_table$V3 - efficacy_table$V2),
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                  (efficacy_table$V3 - efficacy_table$V2),
                                var.equal = TRUE)$p.value
    }
    
    if(tmp_table[4,2] < 0.05){
        tmp_table[10,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[10,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[4,2] < 0.05){
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[4,2] < 0.05){
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                (efficacy_table$V3 - efficacy_table$V2),
                                paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/2,
                                (efficacy_table$V3 - efficacy_table$V2),
                                paired = T, var.equal = TRUE)$p.value
    }
    

    
    return(tmp_table)
}