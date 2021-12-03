dir_path <- paste0(getwd(), "/data-raw/")

# 1) Start by loading the necessary files.

# a) DEET files.
# DEET_metadata: DEET's metadata of studies.
load(paste0(dir_path, "DEET_metadata.rda"))
rownames(DEET_metadata) <- DEET_metadata$DEET.ID

# DEET_DE: DEET's DE of studies.
load(paste0(dir_path, "DEET_DE_final.rda"))


# b) ActivePathway (AP) files.
gmt_BP <- ActivePathways::read.GMT(paste0(dir_path, "Human_GO_AllPathways_with_GO_iea_June_01_2021_symbol.gmt"))
gmt_TF <- ActivePathways::read.GMT(paste0(dir_path, "Human_TranscriptionFactors_MSigdb_June_01_2021_symbol.gmt"))
DEET_gmt_DE <- ActivePathways::read.GMT(paste0(dir_path, "DEET_DE.gmt"))
DEET_gmt_BP <- ActivePathways::read.GMT(paste0(dir_path, "DEET_BP.gmt"))
DEET_gmt_TF <- ActivePathways::read.GMT(paste0(dir_path, "DEET_TF.gmt"))

usethis::use_data(DEET_metadata, DEET_DE, gmt_BP, gmt_TF, DEET_gmt_DE,
                  DEET_gmt_BP, DEET_gmt_TF, internal = TRUE, overwrite = TRUE)
