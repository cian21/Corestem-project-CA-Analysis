# Corestem project - CA Analysis.proj
# 20210712 Sehwan Chun at Corestem, Inc.
# 1.1. load environment for analysis

#### 1. Library Loading ####
packs = c("edgeR","org.Hs.eg.db","data.table", "lmerTest",
          "readxl", "ggpubr","writexl")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. Files Loading ####
# Load Gene expression data ---------------------------
expr_file = "./Data/Rawdata/gene_count_matrix_210709.csv"
expr_file = read.csv(expr_file, stringsAsFactors = F, header = T)
expr_file = expr_file[complete.cases(expr_file),]

# Load efficacy indices data ---------------------------
efficacy_file = "./Data/Rawdata/CA_efficacy.xlsx"
efficacy_file_SARA = read_xlsx(efficacy_file, sheet = 1)

# Non-used 
#efficacy_file_Swayneck = read_xlsx(efficacy_file, sheet = 2)
#efficacy_file_Swaywaist = read_xlsx(efficacy_file, sheet = 3)
#efficacy_file_Gaitvelocity = read_xlsx(efficacy_file, sheet = 4)
#efficacy_file_Gaitcadence = read_xlsx(efficacy_file, sheet = 5)
#efficacy_file_Gaitsupport = read_xlsx(efficacy_file, sheet = 6)

# Load visit date data ---------------------------
visit_date_file = "./Data/Rawdata/유효성검사 평가일.xlsx"
visit_date_SARA_gait = read_xlsx(visit_date_file, sheet = 1)
#visit_date_Wearable = read_xlsx(visit_date_file, sheet = 2)

#### 3. Saving image file ####
save.image(file = "./Data/CA analysis files.image")

