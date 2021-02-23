# Corestem project - FACS analysis.proj
# 20201028 Sehwan Chun at Corestem, Inc.
# 3.1. CA main analysis

#### 1. source Loading ####
load("./Data/CA analysis expr.image")
load("./Data/expr_file_full_cpm.Rdata")
load("./Data/expr_file_030_cpm.Rdata")
load("./Data/expr_file_090_cpm.Rdata")
load("./Data/expr_file_030_SAOA_cpm.Rdata")
source("./Script/Functions/CA analysis - functions.R")

#### 2. simple DEG running ####
DEG_d030 = d030(expr_file, group_count, pairwise = TRUE)
DEG_d090 = d090(expr_file, group_count, pairwise = TRUE)
DEG_d030_SAOA = d030_SAOA(expr_file, group_count, pairwise = TRUE)

#### 3. lmm running with 0,30 ####
SARA_lmm_030 = lmm_index_run(efficacy_file_SARA, expr_file_030_cpm)
Swayneck_lmm_030 = lmm_index_run(efficacy_file_Swayneck, expr_file_030_cpm)
Swaywaist_lmm_030 = lmm_index_run(efficacy_file_Swaywaist, expr_file_030_cpm)
Gaitvelocity_lmm_030 = lmm_index_run(efficacy_file_Gaitvelocity, expr_file_030_cpm)
Gaitcadence_lmm_030 = lmm_index_run(efficacy_file_Gaitcadence, expr_file_030_cpm)
Gaitsupport_lmm_030 = lmm_index_run(efficacy_file_Gaitsupport, expr_file_030_cpm)

#### 4. efficacy test1 ####
SARA_efficacy_test = efficacy_test_run(efficacy_file_SARA)
Swayneck_efficacy_test = efficacy_test_run(efficacy_file_Swayneck)
Swaywaist_efficacy_test = efficacy_test_run(efficacy_file_Swaywaist)
Gaitcadence_efficacy_test = efficacy_test_run(efficacy_file_Gaitcadence)
Gaitsupport_efficacy_test = efficacy_test_run(efficacy_file_Gaitsupport)
Gaitvelocity_efficacy_test = efficacy_test_run(efficacy_file_Gaitvelocity)
total_efficacy_test = cbind(SARA_efficacy_test,
                            Swayneck_efficacy_test$pvalue,
                            Swaywaist_efficacy_test$pvalue,
                            Gaitcadence_efficacy_test$pvalue,
                            Gaitsupport_efficacy_test$pvalue,
                            Gaitvelocity_efficacy_test$pvalue)

#### 5. efficacy test2 ####
SARA_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_SARA, "SARA")
Swayneck_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_Swayneck, "Wear")
Swaywaist_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_Swaywaist, "Wear")
Gaitcadence_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_Gaitcadence, "Gait")
Gaitsupport_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_Gaitsupport ,"Gait")
Gaitvelocity_efficacy_test_td = efficacy_test_run_truedate(efficacy_file_Gaitvelocity ,"Gait")
total_efficacy_test_td = cbind(SARA_efficacy_test_td,
                            Swayneck_efficacy_test_td$pvalue,
                            Swaywaist_efficacy_test_td$pvalue,
                            Gaitcadence_efficacy_test_td$pvalue,
                            Gaitsupport_efficacy_test_td$pvalue,
                            Gaitvelocity_efficacy_test_td$pvalue)

SARA_efficacy_test_td2 = efficacy_test_run_truedate2(efficacy_file_SARA, "SARA")


save.image("./Output/CA_analysis.Rdata")
