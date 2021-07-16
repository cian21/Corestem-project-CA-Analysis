# Corestem project - FACS analysis.proj
# 20201028 Sehwan Chun at Corestem, Inc.
# 3.2. CA PCR analysis

#### 1. source Loading ####
load("./Data/CA analysis files.image")
source("./Script/Functions/CA analysis - functions.R")

plot_run = function(sub, total, genes){
    sub = sub[,2:4]
    sub_melted = melt(sub)
    sub_melted$variable = factor(sub_melted$variable, level = 
                                             c("HC", "SAOA","MSA-C"))
    
    sub_melted$value = log(sub_melted$value)
    my_comparisons =  list(c("HC", "SAOA"),
                            c("SAOA", "MSA-C"),
                            c("HC", "MSA-C"))
    
    my_palette = c("#0000AC","#FC4E07", "#BB3099",  "#EE0099", "#00AFBB", "#E7B800")
    
    sub_plot = ggboxplot(data = sub_melted,
              x = "variable",
              y = "value", 
              color = "variable",
              add = "jitter",
              add.params = list(color = "variable"),
              palette = my_palette) +
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                           label = "p.format")
    
    sub_plot = ggpar(sub_plot, ylab = paste0("log(relative gene expression of ", genes, ")"),
                     xlab = "group", ggtheme = theme_classic(base_size = 13),
                     main = paste0("Relative gene expression of ", genes),
                     legend = "none")
    
    sub_plot
    
    total = total[,2:3]
    total_melted = melt(total)
    total_melted$variable = factor(total_melted$variable, level = 
                                     c("HC", "CA"))
    
    total_melted$value = log(total_melted$value)
    my_comparisons =  list(c("HC", "CA"))
    
    total_plot = ggboxplot(data = total_melted,
                         x = "variable",
                         y = "value", 
                         color = "variable",
                         add = "jitter",
                         add.params = list(color = "variable"),
                         palette = my_palette) +
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                           label = "p.format")
    
    total_plot = ggpar(total_plot, ylab = paste0("log(relative gene expression of ", genes, ")"),
                     xlab = "group", ggtheme = theme_classic(base_size = 13),
                     main = paste0("Relative gene expression of ", genes),
                     legend = "none")
    
    total_plot
    
    double_plot = ggarrange(sub_plot, total_plot, ncol = 2)
    my_file = paste0("./Output/Relative gene expression of ",genes," plot.tiff")
    
    ggsave(my_file,
           double_plot,
           device = "tiff",
           width = 15,
           height = 7.5,
           units = "in",
           dpi = 300)
}

plot_run(CA_PCR_file_ADIPOQ_CA_sub, CA_PCR_file_ADIPOQ_CA_total, "ADIPOQ")
plot_run(CA_PCR_file_FABP4_CA_sub, CA_PCR_file_FABP4_CA_total, "FABP4")
plot_run(CA_PCR_file_LPL_CA_sub, CA_PCR_file_LPL_CA_total, "LPL")
plot_run(CA_PCR_file_PPARG2_CA_sub, CA_PCR_file_PPARG2_CA_total, "PPARG2")

