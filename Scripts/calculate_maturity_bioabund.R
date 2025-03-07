# PURPOSE: to calculate biomass and abundance use morphometric maturity cutlines to define immature vs. mature for 
# snow crab shell 2 males

# Author: Emily Ryznar (with help from Jon Richar's SQL scripts)

## LOAD LIBS/PARAMS ----------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

## LOAD/FILTER DATA ----------------------------------------------------------------------------------------
spec <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>% # survey specimen data from crabpack
        filter(SEX == 1,
               SHELL_CONDITION == 2,
               HAUL_TYPE != 17) %>%
        filter(YEAR > 1988)

haul <- readRDS("./Data/snow_survey_specimenEBS.rda")$haul

params <- read.csv("./Output/maturity_model_params.csv") # maturity model parameters

params <- readRDS("./Data/snow_survey_maturityEBS.rda")$model_parameters

## CALCULATE MORPHOMETRIC MATURITY BIOMASS/ABUNDANCE -------------------------------------------------------
# aggregate data into 1mm bins, apply mat model params 
morph.dat <- spec %>% # this data already has 1MM size bins
              right_join(params, .) %>% # join in model params
              rename(COUNT = SAMPLING_FACTOR) %>%
              filter(is.na(A_EST) == FALSE) %>%
              group_by(YEAR, STATION_ID, HAUL_TYPE, AREA_SWEPT, LATITUDE, LONGITUDE, REGION, DISTRICT, STRATUM,
                       TOTAL_AREA, SPECIES, SEX, SIZE_1MM) %>%
              reframe(COUNT = sum(COUNT), # aggregate by 1mm bins
                      CALCULATED_WEIGHT_1MM = sum(CALCULATED_WEIGHT_1MM)) %>%
              mutate(MATURE_COUNT = ((1/(1 + exp(-A_EST * (SIZE_1MM - B_EST)))) * COUNT), # apply ogive params to get proportion mature
                     MATURE_WEIGHT = ((1/(1 + exp(-A_EST * (SIZE_1MM - B_EST)))) * CALCULATED_WEIGHT_1MM),
                     MATURE_COUNT = case_when((SIZE_1MM <= 50) ~  0, # set mat males <= 50mm to zero
                                              TRUE ~ MATURE_COUNT),
                     MATURE_WEIGHT = case_when((SIZE_1MM <= 50) ~  0, # set mat males <= 50mm (60 for bairdi) to zero
                                           TRUE ~ MATURE_WEIGHT),
                     IMMATURE_COUNT = COUNT - MATURE_COUNT, # calculate immature count and KG
                     IMMATURE_WEIGHT = CALCULATED_WEIGHT_1MM - MATURE_WEIGHT) %>%
              rename(TOTAL_COUNT = COUNT, TOTAL_CALCULATED_WEIGHT_1MM = CALCULATED_WEIGHT_1MM) %>%
              pivot_longer(., cols = c("MATURE_COUNT", "IMMATURE_COUNT", "MATURE_WEIGHT", "IMMATURE_WEIGHT"), names_to = "MAT_SEX", values_to = "VALUE") %>%
              mutate(TYPE = case_when((grepl("_COUNT", MAT_SEX) == "TRUE") ~ "MORPH_COUNT",
                                      TRUE ~ "MORPH_WEIGHT"), # morphometric counts and weights
                      MAT_SEX = case_when((grepl("IMMATURE_", MAT_SEX) == "TRUE") ~ "Immature Male",
                                         TRUE ~ "Mature Male")) %>%
              pivot_wider(., names_from = "TYPE", values_from = "VALUE") %>%
              mutate(CPUE = MORPH_COUNT/AREA_SWEPT,
                     CPUE_KG = (MORPH_COUNT * MORPH_WEIGHT)/AREA_SWEPT/1000) %>%
              group_by(REGION, YEAR, STATION_ID, HAUL_TYPE, LATITUDE, LONGITUDE, DISTRICT, STRATUM,
                       TOTAL_AREA, SPECIES, SIZE_1MM, MAT_SEX) %>%
              reframe(MORPH_COUNT = sum(MORPH_COUNT),
                      MORPH_WEIGHT = sum(MORPH_WEIGHT),
                      CPUE = sum(CPUE),
                      CPUE_KG = sum(CPUE_KG))

years <- c(1989:2019, 2021:2024) 

station_haul_cpue = data.frame()

for(ii in 1:length(years)){
  morph.dat %>%
    filter(YEAR %in% years[ii]) -> data.crab
  
  # Specify unique haul types by year
  HT = unique(data.crab$HAUL_TYPE)
  
  # Specify unique matsex combos
  mat_sex_combos <- c("Mature Male", "Immature Male")
  
  # Join in zero-catch stations
  cpue.dat <- morph.dat %>%
                right_join(expand_grid(REGION = unique(morph.dat$REGION),
                                       SPECIES = unique(morph.dat$SPECIES),
                                       DISTRICT = unique(morph.dat$DISTRICT),
                                       MAT_SEX = mat_sex_combos,
                                       SIZE_1MM = 1:max(morph.dat$SIZE_1MM, na.rm = T),
                                       HAUL_TYPE = HT,
                                       haul %>%
                                         filter(YEAR %in% years[ii]) %>%
                                         distinct(YEAR, STATION_ID, MID_LATITUDE, MID_LONGITUDE, STRATUM, TOTAL_AREA) %>%
                                         rename(LATITUDE = MID_LATITUDE, LONGITUDE = MID_LONGITUDE))) %>%
                replace_na(list(MORPH_COUNT = 0, MORPH_WEIGHT = 0, CPUE = 0, CPUE_KG = 0)) 
  
  # Add data each year
  station_haul_cpue = rbind(station_haul_cpue, cpue.dat)
}

# Calculate abundance and biomass
bio_abund_df <- station_haul_cpue %>%
                group_by(YEAR, STATION_ID, STRATUM, DISTRICT, SPECIES, MAT_SEX, TOTAL_AREA) %>% # summarize across haul types
                dplyr::mutate(MORPH_COUNT = sum(MORPH_COUNT), MORPH_WEIGHT = sum(MORPH_WEIGHT),
                                 CPUE = sum(CPUE), CPUE_KG = sum(CPUE_KG)) %>%
                ungroup() %>%
                #Scale to abundance by strata
                group_by(YEAR, STRATUM, SPECIES, MAT_SEX) %>%
                dplyr::reframe(AREA = TOTAL_AREA,
                               MEAN_CPUE = mean(CPUE),
                               N_CPUE = n(),
                               VAR_CPUE = (var(CPUE)*(AREA^2))/N_CPUE,
                               MEAN_CPUE_KG = mean(CPUE_KG),
                               N_CPUE_KG = n(),
                               VAR_CPUE_KG = (var(CPUE_KG)*(AREA^2))/N_CPUE_KG,
                               ABUNDANCE = (MEAN_CPUE * AREA),
                               BIOMASS = (MEAN_CPUE_KG * AREA),
                               N_STATIONS = length(unique(STATION_ID))) %>%
                distinct() %>%
                #Sum across strata
                group_by(YEAR, SPECIES, MAT_SEX) %>%
                dplyr::reframe(AREA = sum(AREA),
                               MEAN_CPUE = sum(MEAN_CPUE),
                               SD_CPUE = sqrt(sum(VAR_CPUE)),
                               N_CPUE = sum(N_CPUE),
                               MEAN_CPUE_KG = sum(MEAN_CPUE_KG),
                               SD_CPUE_KG = sqrt(sum(VAR_CPUE_KG)),
                               N_CPUE_KG = sum(N_CPUE_KG),
                               ABUNDANCE = sum(ABUNDANCE),
                               ABUNDANCE_CI = 1.96*(SD_CPUE),
                               BIOMASS = sum(BIOMASS),
                               BIOMASS_CI = 1.96*(SD_CPUE_KG),
                               N_STATIONS = sum(N_STATIONS))
          
# Convert abundance to millions, biomass to tons
bio_abund_df2 <- bio_abund_df %>%
                    mutate(ABUNDANCE = ABUNDANCE/1e6,
                           ABUNDANCE_CI = ABUNDANCE_CI/1e6,
                           BIOMASS = BIOMASS/1000,
                           BIOMASS_CI = BIOMASS/1000) 

# Compare to Jon's output
mat.male <- bio_abund_df2 %>%
              dplyr::select(YEAR, MAT_SEX, ABUNDANCE, ABUNDANCE_CI, BIOMASS, BIOMASS_CI) %>%
              filter(MAT_SEX == "Mature Male") 
