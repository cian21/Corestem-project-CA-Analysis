# Corestem project - CA Analysis.proj
# 20201102 Sehwan Chun at Corestem, Inc.
# function.R

#### 1. Library Loading ####
packs = c("edgeR","org.Hs.eg.db","data.table","lmerTest", "readxl", "ggpubr",
          "writexl", "HH","dplyr","variancePartition", "biomaRt", "ggrepel",
          "psych", "ReactomePA", "clusterProfiler", "enrichplot")
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
d030_paired = function(expr_file, group_count){
  expr_file_d030 = expr_file[group_count == "d0" | group_count == "d30"]
  group_count_d030 = group_count[group_count == "d0" | group_count == "d30"]
  
  pairwised_number = c()
  for(i in 1:(length(group_count_d030) - 1)){
    if(group_count_d030[i] == "d0" & group_count_d030[i+1] == "d30"){
      pairwised_number = c(pairwised_number,i)
    }
  }
  pairwised_number = sort(c(pairwised_number,pairwised_number+1))
  expr_file_d030 = expr_file_d030[pairwised_number]
  group_count_d030 = group_count_d030[pairwised_number]
  
  subject_d030 = sort(rep(seq(1,9),2))
  group_frame_030 = data.frame(subject = subject_d030, days = group_count_d030)
  age_d030 = c(58,58,61,61,64,64,39,39,70,70,75,75,57,57,61,61,54,54)
  ageb_d030 = as.factor(c(1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1))
  sex_d030 = as.factor(c(0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1))
  type_d030 = as.factor(c(2,2,2,2,1,1,2,2,3,3,2,2,1,1,3,3,1,1))
  design_matrix_030 = model.matrix(~subject_d030+ageb_d030+type_d030+sex_d030+group_count_d030, data = group_frame_030)
  
  DGE = DGEList(counts = expr_file_d030, group = group_count_d030)
  keep = filterByExpr(DGE)
  DGE = DGE[keep, , keep.lib.sizes = FALSE] 
  DGE = calcNormFactors(DGE)
  DGE = estimateDisp(DGE, design_matrix_030)
  
  DGEfit = glmQLFit(DGE, design_matrix_030)
  DGEqlf = glmQLFTest(DGEfit)
  topTags(DGEqlf)
  return(DGEqlf)
}
d090_paired = function(expr_file, group_count){
  expr_file_d090 = expr_file[group_count == "d0" | group_count == "d90"]
  group_count_d090 = group_count[group_count == "d0" | group_count == "d90"]
  
  pairwised_number = c()
  for(i in 1:(length(group_count_d090) - 1)){
    if(group_count_d090[i] == "d0" & group_count_d090[i+1] == "d90"){
      pairwised_number = c(pairwised_number,i)
    }
  }
  pairwised_number = sort(c(pairwised_number,pairwised_number+1))
  expr_file_d090 = expr_file_d090[pairwised_number]
  group_count_d090 = group_count_d090[pairwised_number]
  
  subject_d090 = sort(rep(seq(1,8),2))
  group_frame_090 = data.frame(subject = subject_d090, days = group_count_d090)
  design_matrix_090 = model.matrix(~subject_d090+group_count_d090, data = group_frame_090)
  
  DGE = DGEList(counts = expr_file_d090, group = group_count_d090)
  keep = filterByExpr(DGE)
  DGE = DGE[keep, , keep.lib.sizes = FALSE] 
  DGE = calcNormFactors(DGE)
  DGE = estimateDisp(DGE, design_matrix_090)
  
  DGEfit = glmQLFit(DGE, design_matrix_090)
  DGEqlf = glmQLFTest(DGEfit)

  return(DGEqlf)
}
dotplot_run = function(genelist, filename){
  tmp = enrichGO(gene = genelist,
                 OrgDb = org.Hs.eg.db, 
                 ont =  "CC", 
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.1,
                 keyType = "SYMBOL")
  
  if(nrow(subset(tmp@result, tmp@result$p.adjust < 0.1))>=1){
    dot_go = clusterProfiler::dotplot(tmp, showCategory = 15, orderBy = "GeneRatio")
    ggsave(plot = dot_go,
           filename = paste0("./Output/",
                             filename,
                             "_CC_enrichment_dotplot.tiff"),
           device = "tiff", width = 12, height = 10, dpi = 300)
  }
  
  tmp = enrichGO(gene = genelist,
                 OrgDb = org.Hs.eg.db, 
                 ont =  "MF", 
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.1,
                 keyType = "SYMBOL")
  
  if(nrow(subset(tmp@result, tmp@result$p.adjust < 0.1))>=1){
    dot_go = clusterProfiler::dotplot(tmp, showCategory = 15, orderBy = "GeneRatio")
    ggsave(plot = dot_go,
           filename = paste0("./Output/",
                             filename,
                             "_MF_enrichment_dotplot.tiff"),
           device = "tiff", width = 12, height = 10, dpi = 300)
  }
  tmp = enrichGO(gene = genelist,
                 OrgDb = org.Hs.eg.db, 
                 ont =  "BP", 
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.1,
                 keyType = "SYMBOL")
  
  if(nrow(subset(tmp@result, tmp@result$p.adjust < 0.1))>=1){
    dot_go = clusterProfiler::dotplot(tmp, showCategory = 15, orderBy = "GeneRatio")
    ggsave(plot = dot_go,
           filename = paste0("./Output/",
                             filename,
                             "_BP_enrichment_dotplot.tiff"),
           device = "tiff", width = 12, height = 10, dpi = 300)
  }
}
basic_plot_run = function(table, filename){
  a = gghistogram(table$PValue, xlab = "p-value", main = "raw p-values", bins = 20)
  b = gghistogram(table$FDR, xlab = "FDR", main = "adjusted p-values(FDR)", bins = 20)
  ggsave(paste0("./Output/",filename,"_hist.tiff"), plot = ggarrange(a,b, ncol = 2), device = "tiff", width = 14, height = 6, dpi = 300)
  
  table$expression = ifelse(table$PValue < 0.01 &
                              abs(table$logFC) >= 0.5, 
                            ifelse(table$logFC> 0.5 ,'Up','Down'),
                            'Stable')
  c =  ggplot(data = table, 
              aes(x = logFC, 
                  y = -log10(table$PValue),
                  color = expression,
                  label = table$Symbol)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("blue", "grey","red"))+
    xlim(c(-1.5, 1.5)) +
    geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (PValue)",
         title="Differential expression")  +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          legend.position="right", 
          legend.title = element_blank()) +
    geom_text_repel(aes(label=ifelse(PValue < 0.01 & abs(logFC) >= 0.5,
                                     row.names(table), '')))
  
  ggsave(paste0("./Output/",filename,"_volcano.tiff"), plot = c, device = "tiff", width = 9, height = 9, dpi = 300)
  
  sign_genes_03 = row.names(subset(table, PValue < 0.01))
  dotplot_run(sign_genes_03, filename)
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
        for (j in 1:(ncol(cpm))){
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
        if (i %% 1000 == 0) {print(i)}
        
    }
    df_table$FDR_lmm = p.adjust(df_table$V3)
    df_table$FDR_cor = p.adjust(df_table$V5)
    
    colnames(df_table)[2:5] = c("lmm_estimate","lmm_pvalue", "cor_estimate","cor_pvalue")
    
    return(df_table)
}
lmm_index_run_short = function(df, cpm){
  df = df[,c(1,3,4)]
  df = reshape2::melt(df)
  df = df[order(df$Samples),]
  row.names(df) = 1:nrow(df)
  
  df$Samples = as.factor(df$Samples)
  df$Days = ifelse(df$variable == "V2", 0,30)
  df_table = data.frame(genes = row.names(cpm))
  
  for (i in 1:nrow(cpm)){
    df$tmp = NA
    for (j in 1:(ncol(cpm))){
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
    if (i %% 1000 == 0) {print(i)}
    
  }
  df_table$FDR_lmm = p.adjust(df_table$V3)
  df_table$FDR_cor = p.adjust(df_table$V5)
  
  colnames(df_table)[2:5] = c("lmm_estimate","lmm_pvalue", "cor_estimate","cor_pvalue")
  
  return(df_table)
}
lmm_index_run_full = function(df, cpm){
  df = df[,c(1,3:9)]
  df = reshape2::melt(df)
  df = df[order(df$Samples),]
  row.names(df) = 1:nrow(df)
  
  df$Samples = as.factor(df$Samples)
  for (i in 1:nrow(df)){
    if (df[i, 2] == "V2") {df[i, 4] = 00}
    if (df[i, 2] == "V3") {df[i, 4] = 30}
    if (df[i, 2] == "V4") {df[i, 4] = 90}
    if (df[i, 2] == "V5") {df[i, 4] = 180}
    if (df[i, 2] == "V6") {df[i, 4] = 210}
    if (df[i, 2] == "V7") {df[i, 4] = 270}
    if (df[i, 2] == "V8") {df[i, 4] = 360}
  }
  colnames(df)[4] = "Days"
  df = df[complete.cases(df),]
  
  df_table = data.frame(genes = row.names(cpm))
  
  for (i in 1:nrow(cpm)){
    df$tmp = NA
    for (j in 1:(ncol(cpm))){
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
    
    if (i %% 1000 == 0) {print(i)}
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
    if(tmp_table[5,2] < 0.05){
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
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
    if(tmp_table[5,2] < 0.05){
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
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
efficacy_test_run_truedate  = function(efficacy_table, type){
    if(type == "SARA" | type == "Gait"){
        visit = visit_date_SARA_gait[1:ncol(visit_date_SARA_gait),]
    }else{
        visit = visit_date_Wearable[1:ncol(visit_date_Wearable),]
    }
    visit = visit[order(visit$Samples),]
    efficacy_table = efficacy_table[order(efficacy_table$Samples),]
    
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
    tmp_table[3,2] = wilcox.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 (efficacy_table$V3 - efficacy_table$V2))$p.value
    
    tmp_table[4,2] = var.test(efficacy_table$V2, efficacy_table$V3)$p.value
    tmp_table[5,2] = var.test(efficacy_table$V2, efficacy_table$V4)$p.value
    tmp_table[6,2] = var.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                              (efficacy_table$V3 - efficacy_table$V2))$p.value
    
    if(tmp_table[4,2] < 0.05){
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[5,2] < 0.05){
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                (efficacy_table$V3 - efficacy_table$V2),
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
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
    if(tmp_table[5,2] < 0.05){
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 (efficacy_table$V3 - efficacy_table$V2),
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 (efficacy_table$V3 - efficacy_table$V2),
                                 paired = T, var.equal = TRUE)$p.value
    }
    
    return(tmp_table)
}
efficacy_test_run_truedate2  = function(efficacy_table, type){
    if(type == "SARA" | type == "Gait"){
        visit = visit_date_SARA_gait[1:ncol(visit_date_SARA_gait),]
    }else{
        visit = visit_date_Wearable[1:ncol(visit_date_Wearable),]
    }
    visit = visit[order(visit$Samples),]
    efficacy_table = efficacy_table[order(efficacy_table$Samples),]
    
    tmp_table  = data.frame(test = NA, pvalue = NA)
    tmp_table[1:12,1] = c("wilcox.test_v2_v3",
                          "wilcox.test_v2_v4",
                          "wilcox.test_v2-v1_v4-v2",
                          "var.test_v2_v3",
                          "var.test_v2_v4",
                          "var.test_v2_v1-v4-v2",
                          "t.test_v2_v3",
                          "t.test_v2_v4",
                          "t.test_v2-v1_v4-v2",
                          "pt.test_v2_v3",
                          "pt.test_v2_v4",
                          "pt.test_v2-v1_v4-v2")
    tmp_table[1,2] = wilcox.test(efficacy_table$V2, efficacy_table$V3)$p.value
    tmp_table[2,2] = wilcox.test(efficacy_table$V2, efficacy_table$V4)$p.value
    tmp_table[3,2] = wilcox.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 ((efficacy_table$V4 - efficacy_table$V2) / 3))$p.value
    
    tmp_table[4,2] = var.test(efficacy_table$V2, efficacy_table$V3)$p.value
    tmp_table[5,2] = var.test(efficacy_table$V2, efficacy_table$V4)$p.value
    tmp_table[6,2] = var.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                              ((efficacy_table$V4 - efficacy_table$V2) / 3))$p.value
    
    if(tmp_table[4,2] < 0.05){
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[7,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[5,2] < 0.05){
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[8,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                ((efficacy_table$V4 - efficacy_table$V2) / 3),
                                var.equal = FALSE)$p.value
    }else{
        tmp_table[9,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                ((efficacy_table$V4 - efficacy_table$V2) / 3),
                                var.equal = TRUE)$p.value
    }
    
    if(tmp_table[4,2] < 0.05){
        tmp_table[10,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[10,2] = t.test(efficacy_table$V2, efficacy_table$V3,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[5,2] < 0.05){
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[11,2] = t.test(efficacy_table$V2, efficacy_table$V4,
                                 paired = T, var.equal = TRUE)$p.value
    }
    if(tmp_table[6,2] < 0.05){
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 ((efficacy_table$V4 - efficacy_table$V2) / 3),
                                 paired = T, var.equal = FALSE)$p.value
    }else{
        tmp_table[12,2] = t.test((efficacy_table$V2 - efficacy_table$V1)/ ((visit$V2 - visit$V1) / 28),
                                 ((efficacy_table$V4 - efficacy_table$V2) / 3),
                                 paired = T, var.equal = TRUE)$p.value
    }
    
    return(tmp_table)
}
efficacy_test_run_truedate_with_plot  = function(efficacy_table, visit_table, index, pair){
  
  
  if(nrow(efficacy_table) != nrow(visit_table)){
    visit_table = visit_table[nrow(efficacy_table),]
  }
  
  if(pair == TRUE){
    tmp = "_paired"
  }else{
    tmp = ""
  }
  
  if(index == "SARA"){
  efficacy_table$V1_V2 = (efficacy_table$V2 - efficacy_table$V1) / 
    ((visit_table$V2 - visit_table$V1) / 28)
  efficacy_table$V2_V3 = (efficacy_table$V3 - efficacy_table$V2) / 
    ((visit_table$V3 - visit_table$V2) / 28)
  efficacy_table$V2_V4 = (efficacy_table$V4 - efficacy_table$V2) / 
    ((visit_table$V4 - visit_table$V2) / 28)
  
  tmp1 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = index,) +
    stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
  
  tmp2 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = index) +
    stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
  
  tmp3 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V4",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = index) +
    stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
  
  tmp4 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V4",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = index) +
    stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
  
  tmp5 = ggpaired(data = efficacy_table[complete.cases(efficacy_table[,c(3,6)]),],
                  cond1 = "V2", cond2 = "V5", line.color = "gray",
                  line.size = 1, palette = "aaas", fill = "condition",
                  xlab = "Visit", ylab = index) +
    stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
  
  tmp6 = ggpaired(data = efficacy_table[complete.cases(efficacy_table[,c(3,6)]),],
                  cond1 = "V2", cond2 = "V5", line.color = "gray",
                  line.size = 1, palette = "aaas", fill = "condition",
                  xlab = "Visit", ylab = index) +
    stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
  
  tmp7 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V3",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
    stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
  
  tmp8 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V3",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
    stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
  
  tmp9 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V4",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
    stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
  
  tmp10 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V4",
                  line.color = "gray", line.size = 1, palette = "aaas",
                  fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
    stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V23_wilcox_boxplot_",index, tmp, ".tiff"),
         tmp1, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V23_t_boxplot_",index, tmp, ".tiff"),
         tmp2, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V24_wilcox_boxplot_",index, tmp, ".tiff"),
         tmp3, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V24_t_boxplot_",index, tmp, ".tiff"),
         tmp4, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V25_wilcox_boxplot_",index, tmp, ".tiff"),
         tmp5, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V25_t_boxplot_",index, tmp, ".tiff"),
         tmp6, width = 8, height = 8)
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V123slope_wilcox_boxplot_",index, tmp, ".tiff"),
         tmp7, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V123slope_t_boxplot_",index, tmp, ".tiff"),
         tmp8, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V124slope_wilcox_boxplot_",index, tmp, ".tiff"),
         tmp9, width = 8, height = 8)
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V124slope_t_boxplot_",index, tmp, ".tiff"),
         tmp10, width = 8, height = 8)
  
  }else{
    efficacy_table$V1_V2 = (efficacy_table$V2 - efficacy_table$V1) / 
      ((visit_table$V2 - visit_table$V1) / 28)
    efficacy_table$V2_V3 = (efficacy_table$V3 - efficacy_table$V2) / 
      ((visit_table$V3 - visit_table$V2) / 28)
    
    tmp1 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3",
                    line.color = "gray", line.size = 1, palette = "aaas",
                    fill = "condition", xlab = "Visit", ylab = index) +
      stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
    
    tmp2 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3",
                    line.color = "gray", line.size = 1, palette = "aaas",
                    fill = "condition", xlab = "Visit", ylab = index) +
      stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
    
    tmp3 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V3",
                    line.color = "gray", line.size = 1, palette = "aaas",
                    fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
      stat_compare_means(paired = pair, method = "wilcox.test", label.x = 0.6)
    
    tmp4 = ggpaired(data = efficacy_table, cond1 = "V1_V2", cond2 = "V2_V3",
                    line.color = "gray", line.size = 1, palette = "aaas",
                    fill = "condition", xlab = "Visit", ylab = paste0(index,"_slope")) +
      stat_compare_means(paired = pair, method = "t.test", label.x = 0.6)
    
    ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V23_wilcox_boxplot_",index, tmp, ".tiff"),
           tmp1, width = 8, height = 8)
    ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V23_t_boxplot_",index, tmp, ".tiff"),
           tmp2, width = 8, height = 8)
    ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V123slope_wilcox_boxplot_",index, tmp, ".tiff"),
           tmp3, width = 8, height = 8)
    ggsave(paste0("../Corestem-project-CA-Analysis/Output/slope/V123slope_t_boxplot_",index, tmp, ".tiff"),
           tmp4, width = 8, height = 8)
    
  }

}
efficacy_test_with_plot = function(efficacy_table, visit_table, index){
  tmp12 = ggpaired(data = efficacy_table, cond1 = "V1", cond2 = "V2", line.color = "gray60", line.size = 1, palette = "aaas",fill = "condition", xlab = "Visit", ylab = paste0(index, " score"), title = paste0("Comparison of ", index, " scores by visit (V1 vs V2)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V12_wilcox_",index,"_paired_boxplot.tiff"), tmp12, width = 8, height = 8)
  
  tmp23 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V3", line.color = "gray60", line.size = 1, palette = "aaas",fill = "condition", xlab = "Visit", ylab = paste0(index, " score"), title = paste0("Comparison of ", index, " scores by visit (V2 vs V3)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V23_wilcox_",index,"_paired_boxplot.tiff"), tmp23, width = 8, height = 8)
  
  tmp24 = ggpaired(data = efficacy_table, cond1 = "V2", cond2 = "V4", line.color = "gray60", line.size = 1, palette = "aaas",fill = "condition", xlab = "Visit", ylab = paste0(index, " score"), title = paste0("Comparison of ", index, " scores by visit (V2 vs V4)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V24_wilcox_", index, "_paired_boxplot.tiff"), tmp24, width = 8, height = 8)
  
  tmp25 = ggpaired(data = efficacy_table[-2,], cond1 = "V2", cond2 = "V5", line.color = "gray60", line.size = 1, palette = "aaas",fill = "condition", xlab = "Visit", ylab = paste0(index, " score"), title = paste0("Comparison of ", index, " scores by visit (V2 vs V5)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=5,label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V25_wilcox_", index, "_paired_boxplot.tiff"), tmp25, width = 8, height = 8)
  
  efficacy_table$V12 = (efficacy_table$V2 - efficacy_table$V1) / ((visit_table$V2 - visit_table$V1) / 28)
  efficacy_table$V23 = (efficacy_table$V3 - efficacy_table$V2) / ((visit_table$V3 - visit_table$V2) / 28)
  efficacy_table$V24 = (efficacy_table$V4 - efficacy_table$V2) / ((visit_table$V4 - visit_table$V2) / 28)
  efficacy_table$V25 = (efficacy_table$V5 - efficacy_table$V2) / ((visit_table$V5 - visit_table$V2) / 28)
  
  print(efficacy_table)
  
  tmp1223 = ggpaired(data = efficacy_table, cond1 = "V12", cond2 = "V23", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V3)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=min((c(efficacy_table$V12,efficacy_table$V23)) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V1223_wilcox_", index, "_paired_boxplot.tiff"), tmp1223, width = 8, height = 8)
  
  tmp1224 = ggpaired(data = efficacy_table, cond1 = "V12", cond2 = "V24", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V4)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=min((c(efficacy_table$V12,efficacy_table$V24)) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V1224_wilcox_", index, "_paired_boxplot.tiff"), tmp1224, width = 8, height = 8)
  
  tmp1225 = ggpaired(data = efficacy_table[-2,], cond1 = "V12", cond2 = "V25", line.color = "gray60", line.size = 1, palette = "aaas", fill = "condition", xlab = "Period", ylab = paste0(index, " slope"), title = paste0("Comparison of ", index, " slope (V1-V2 vs V2-V5)")) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 0.6) + 
    stat_summary(fun.data = function(x) data.frame(y=min(na.omit((c(efficacy_table$V12,efficacy_table$V25))) - 1),label = paste(round(mean(x), digits = 2),"\U00B1",round(sd(x), digits = 2))), geom="text")
  
  ggsave(paste0("../Corestem-project-CA-Analysis/Output/index/V1225_wilcox_", index, "_paired_boxplot.tiff"), tmp1225, width = 8, height = 8)
  
}

