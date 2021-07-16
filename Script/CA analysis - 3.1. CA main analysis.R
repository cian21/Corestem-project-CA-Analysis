# Corestem project - FACS analysis.proj
# 20210712 Sehwan Chun at Corestem, Inc.
# 3.1. CA main analysis

#### 1. source Loading ####
load("./Data/expr_file_full_cpm.Rdata")
load("./Data/CA analysis expr.image")
source("./Script/Functions/CA analysis - functions.R")

#### 2. revise column names in full cpm file ####
expr_file_full_cpm_rev = expr_file_full_cpm 
for (i in 1:ncol(expr_file_full_cpm_rev)){
    tmp = colnames(expr_file_full_cpm_rev)[i]
    if (substr(tmp,4,4) == 2) {
        tmp_name = substr(tmp,1,3)
        tmp_day = as.numeric(strsplit(tmp, "_")[[1]][2]) + 180
        colnames(expr_file_full_cpm_rev)[i] = paste0(tmp_name,"_",tmp_day)
    }
}

#### 3. lmm running with SARA ####
SARA_lmm_0030 = lmm_index_run_short(efficacy_file_SARA, expr_file_full_cpm_rev)
SARA_lmm_003090 = lmm_index_run(efficacy_file_SARA, expr_file_full_cpm_rev)
SARA_lmm_00_360 = lmm_index_run_full(efficacy_file_SARA, expr_file_full_cpm_rev)
write.table(SARA_lmm_00_360, "./Output/CA_SARA_lmm_00_360.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE)
dotplot_run(subset(SARA_lmm_00_360, FDR_cor < 0.1)$genes, "SARA_full")
#Swayneck_lmm_030 = lmm_index_run(efficacy_file_Swayneck, expr_file_030_cpm)
#Swaywaist_lmm_030 = lmm_index_run(efficacy_file_Swaywaist, expr_file_030_cpm)
#Gaitvelocity_lmm_030 = lmm_index_run(efficacy_file_Gaitvelocity, expr_file_030_cpm)
#Gaitcadence_lmm_030 = lmm_index_run(efficacy_file_Gaitcadence, expr_file_030_cpm)
#Gaitsupport_lmm_030 = lmm_index_run(efficacy_file_Gaitsupport, expr_file_030_cpm)

#### 4. efficacy test with SARA ####
SARA_efficacy_test = efficacy_test_with_plot(efficacy_file_SARA, visit_date_SARA_gait, "SARA")
#Swayneck_efficacy_test = efficacy_test_with_plot(efficacy_file_Swayneck, "sway area (neck)")
#Swaywaist_efficacy_test = efficacy_test_with_plot(efficacy_file_Swaywaist, "sway area (waist)")
#Gaitcadence_efficacy_test = efficacy_test_with_plot(efficacy_file_Gaitcadence, "gait rite (cadence)")
#Gaitsupport_efficacy_test = efficacy_test_with_plot(efficacy_file_Gaitsupport, "gait rite (mean base of support)")
#Gaitvelocity_efficacy_test = efficacy_test_with_plot(efficacy_file_Gaitvelocity, "gait rite (velocity)")
SARA_efficacy_test = efficacy_test_with_plot(subset(efficacy_file_SARA, Samples %in% c("KYS","LSJ","JYJ","BMS")),"SARA")
SARA_efficacy_test = efficacy_test_with_plot(subset(efficacy_file_SARA, Samples %in% c("BBH","LJO","OSH")),"SARA")
SARA_efficacy_test = efficacy_test_with_plot(subset(efficacy_file_SARA, Samples %in% c("KKH","KDB")),"SARA")
