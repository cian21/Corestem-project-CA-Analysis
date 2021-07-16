# Corestem project - CA Analysis.proj
# 20210712 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for CA data analysis

#### 1. Source Loading ####
load("./Data/CA analysis files.image")
source("./Script/Functions/CA analysis - functions.R")

#### 2. Cpm cleaning ####
# Cleaning gene expression data ---------------------------
row.names(expr_file) = expr_file$gene_id
expr_file$gene_id = NULL
group_name = substr(colnames(expr_file),1,4)
group_count = make_group_count(expr_file)

# DEG analysis with selected data (edgeR) ---------------------------
expr_file_030 = expr_file[,which(substr(group_name,4,4) == "_" & group_count %in% c("d0","d30"))] #select samples
expr_file_030_gc = make_group_count(expr_file_030)
expr_file_030_DEG = d030(expr_file_030, expr_file_030_gc, pairwise = T)

basic_plot_run(expr_file_030_DEG, "CA_0030")
sign_gene_030 = subset(expr_file_030_DEG, PValue < 0.001)$Symbol

# DEG analysis with selected data (DREAM) ---------------------------
expr_file_030 = expr_file[,which(group_count %in% c("d0","d30"))] #select samples
expr_file_030_cpm = cpm_calculation(expr_file_030)
expr_file_030_cpm = anno_DEG(as.data.frame(expr_file_030_cpm))
expr_file_030_cpm = expr_file_030_cpm[,-c((ncol(expr_file_030_cpm)-1),(ncol(expr_file_030_cpm)))] #got CPM 

# Setting DREAM parameter and design matrix ---------------------------
design_df = data.frame(row.names = colnames(expr_file_030_cpm),
                       Individual = substr(colnames(expr_file_030_cpm),1,3),
                       Treatment = ifelse(substr(colnames(expr_file_030_cpm),5,6) %in% c(0, "_0"), "N", "Y"))

for (i in 1:nrow(design_df)){
    design_df[i,3] = gsub("[1-9]","",efficacy_file_SARA[which(efficacy_file_SARA$Samples == substr(colnames(expr_file_030_cpm)[i],1,3)), 11])
    design_df[i,4] = efficacy_file_SARA[which(efficacy_file_SARA$Samples == substr(colnames(expr_file_030_cpm)[i],1,3)), 13]
}

colnames(design_df)[c(3,4)] = c("Subgroup", "Age")
form = ~ Treatment + (1|Individual) + (1|Age) + + (1|Subgroup) 
vobjDream = voomWithDreamWeights(expr_file_030_cpm, form, design_df)

fitmm = dream(vobjDream, form, design_df)
table = topTable(fitmm, adjust.method = "fdr", number = nrow(fitmm$p.value))

# Check p-value distribution ---------------------------
a = gghistogram(subset(topTable(fitmm, adjust.method = "fdr", number = 13094))$P.Value, xlab = "p-value", main = "raw p-values", bins = 20)
b = gghistogram(subset(topTable(fitmm, adjust.method = "fdr", number = 13094))$adj.P.Val, xlab = "FDR", main = "adjusted p-values(FDR)", bins = 20)
ggsave(paste0("./Output/CA_IIT_DEG_0030D_hist.tiff"), plot = ggarrange(a,b, ncol = 2), device = "tiff", width = 14, height = 6, dpi = 300)

# Check volcano plot ---------------------------
table$expression = ifelse(table$P.Value < 0.01 &
                              abs(table$logFC) >= 0.5, 
                          ifelse(table$logFC> 0.5 ,'Up','Down'),
                          'Stable')
c =  ggplot(data = table, 
            aes(x = logFC, 
                y = -log10(table$P.Value),
                color = expression,
                label = table$Symbol)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (PValue)",
         title="Differential expression")  +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          legend.position="right", 
          legend.title = element_blank()) +
    geom_text_repel(aes(label=ifelse(P.Value < 0.01 & abs(logFC) >= 0.5,
                                     row.names(table), '')))

ggsave(paste0("./Output/CA_IIT_DEG_0030D_volcano.tiff"), plot = c, device = "tiff", width = 9, height = 9, dpi = 300)

#### 3. Extract full cpm ####
expr_file_full_cpm = cpm_calculation(expr_file)
expr_file_full_cpm = anno_DEG(as.data.frame(expr_file_full_cpm))
expr_file_full_cpm = expr_file_full_cpm[,-c((ncol(expr_file_full_cpm)-1),(ncol(expr_file_full_cpm)))]
save(expr_file_full_cpm, file = "./Data/expr_file_full_cpm.Rdata")
