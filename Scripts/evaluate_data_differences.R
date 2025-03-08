# PURPOSE: to calculate differences for number of chela measured between survey data output by crabpack and 
# chela data provided by Jon Richar

# NOTES: 

# - For tanner and snow crab, it seems like data < 2008 and < 2020, respectively, use only shell condition 2 crab
# - Other discrepancies in recent years could be do to special project data being included in Jon's data but not in
# data output by crabpack
# - ***Only shell 2 crab chela are used in the stock assessment for maturity estimates***


# LOAD LIBS/PARAMS ----------
source("./Scripts/load_libs_params.R")


# TANNER CRAB ---------------
cp.dat <- readRDS("./Data/tanner_survey_specimenEBS.rda")$specimen
jon.dat <- read.csv("./Data/ebscrab_bairdi_chela_TS_combined.csv")

cp.dat %>%
  filter(YEAR %in% 1989:2024, HAUL_TYPE !=17,
         case_when((YEAR <2008) ~ SHELL_CONDITION == 2,
                   TRUE ~ SHELL_CONDITION == SHELL_CONDITION)) %>%
  group_by(YEAR) %>%
  reframe(tot_ms = sum(is.na(CHELA_HEIGHT) == FALSE)) -> pp

jon.dat %>%
  mutate(YEAR = as.numeric(substr(CRUISE, 1, 4))) %>%
  group_by(YEAR) %>%
  reframe(tot_msrd_jon = sum(is.na(CHELA_HEIGHT) == FALSE)) -> qq

pp %>%
  filter(YEAR %in% qq$YEAR) %>%
  rename(tot_msrd_survey = tot_ms) -> tt

right_join(tt, qq) %>%
  mutate(diff = tot_msrd_survey - tot_msrd_jon) %>%
  dplyr::select(YEAR, tot_msrd_survey, tot_msrd_jon, diff) %>%
  as.data.frame() -> diff.bairdi

write.csv(diff.bairdi, "./Output/diff_chelamsrd_bairdi.csv")

# SNOW CRAB ---------------
cp.dat <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen
jon.dat <- read.csv("./Data/opilio_chela_height_TS.csv")


cp.dat %>%
  filter(YEAR %in% 1989:2024, HAUL_TYPE !=17,
         case_when((YEAR <2020) ~ SHELL_CONDITION == 2,
                   TRUE ~ SHELL_CONDITION == SHELL_CONDITION)) %>%
  group_by(YEAR) %>%
  reframe(tot_ms = sum(is.na(CHELA_HEIGHT) == FALSE)) -> pp

jon.dat %>%
  mutate(YEAR = as.numeric(substr(CRUISE, 1, 4))) %>%
  group_by(YEAR) %>%
  reframe(tot_msrd_jon = sum(is.na(CHELA_HEIGHT) == FALSE)) -> qq

pp %>%
  filter(YEAR %in% qq$YEAR) %>%
  rename(tot_msrd_survey = tot_ms) -> tt

right_join(tt, qq) %>%
  mutate(diff = tot_msrd_survey - tot_msrd_jon) %>%
  dplyr::select(YEAR, tot_msrd_survey, tot_msrd_jon, diff) %>%
  as.data.frame() -> diff.snow

write.csv(diff.snow, "./Output/diff_chelamsrd_snow.csv")


# Comparing maturity ratios between Emily and Jon's estimates ----
# read in Emily's data
agg.prop.dat <- read.csv("./Output/opilio_propmat_agg.csv")

# Read in maturity output from crabpack (Jon's data)
ratio <- readRDS("./Data/snow_survey_maturityEBS.rda")$male_mat_ratio
pars <- readRDS("./Data/snow_survey_maturityEBS.rda")$model_parameters

# calculate differences for male maturity ratios
agg.prop.dat %>%
  dplyr::select(YEAR, MIDPOINT, NUM_IMMATURE, NUM_MATURE, TOTAL_CRAB, PROP_MATURE) %>%
  rename(SIZE_BIN = MIDPOINT, TOTAL_CRAB_EM = TOTAL_CRAB, PROP_MATURE_EM = PROP_MATURE, 
         NUM_IMMATURE_EM = NUM_IMMATURE, NUM_MATURE_EM = NUM_MATURE) -> em.dat

ratio %>%
  dplyr::select(YEAR, SIZE_BIN, NUM_IMMATURE, NUM_MATURE, TOTAL_CRAB, PROP_MATURE) %>%
  right_join(., em.dat) %>%
  filter(is.na(TOTAL_CRAB) == FALSE) %>%
  mutate(diff_imm = NUM_IMMATURE - NUM_IMMATURE_EM,
         diff_mat = NUM_MATURE - NUM_MATURE_EM,
         diff_total = TOTAL_CRAB - TOTAL_CRAB_EM,
         diff_prop = PROP_MATURE - PROP_MATURE_EM) -> ll

ll %>%
  filter(diff_total > 0.0001) -> probs

# calculate differences for model parameters
jon.pars <- pars
diff <- read.csv("./Output/maturity_model_params.csv") %>%
  rename(A_EST_EM = A_EST, A_SE_EM = A_SE, B_EST_EM = B_EST, B_SE_EM = B_SE) %>%
  right_join(., jon.pars) %>%
  mutate(diff_AEST = A_EST_EM - A_EST,
         diff_ASE = A_SE_EM - A_SE,
         diff_BEST = B_EST_EM - B_EST,
         diff_BSE = B_SE_EM - B_SE)




