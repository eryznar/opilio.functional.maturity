# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
df.dat <- read.csv("./Data/opilio_directedfishery_catch.csv") %>%
            dplyr::select(Year, Retained_kt) %>%
            rename(directedfish_biomass = Retained_kt)# directed fishery data

survey.dat <- read.csv("./Data/opilio_survey_biomass.csv") %>% # survey data
                mutate(small_male_biomass = gsub("\\s*\\([^\\)]+\\)", "", small_male_biomass),
                       large_male_biomass = gsub("\\s*\\([^\\)]+\\)", "", large_male_biomass),
                       small_male_biomass = as.numeric(as.character(gsub(",", "", small_male_biomass))),
                       large_male_biomass = as.numeric(as.character(gsub(",", "", large_male_biomass))),
                       small_male_biomass = small_male_biomass/1000,
                       large_male_biomass= large_male_biomass/1000) # convert to kilotons 
                       #Year = c(1988:2019, 2021:2024) - 1) # lagging survey data back one year 

morph.dat <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE) %>%
  mutate(CUTOFF = BETA0 + BETA1*(log(SIZE_1MM)),
         MATURE = case_when((log(CHELA_HEIGHT) > CUTOFF) ~ 1,
                            TRUE ~ 0)) %>%
  dplyr::select(YEAR, SIZE_1MM, MATURE) %>%
  rename(Year = YEAR) %>%
  group_by(Year) %>%
  mutate(size_at_mat = ifelse(MATURE == 1, mean(SIZE_1MM[MATURE == 1], na.rm = TRUE), NA))%>%
  ungroup()


right_join(df.dat, survey.dat) %>%
  right_join(., morph.dat) -> model.dat

# FIT MODELS ---------------------------------------------------------------------------------------------
model.dat %>%
  filter(is.na(size_at_mat) == FALSE) %>%
  dplyr::select(!SIZE_1MM, MATURE) %>%
  distinct() -> model.dat2

mod.1 <- gam(size_at_mat ~ s(directedfish_biomass) + s(small_male_biomass) + s(large_male_biomass), data = model.dat2)
mod.2 <- gam(MATURE ~ s(directedfish_biomass) + s(small_male_biomass) + s(large_male_biomass) + s(SIZE_1MM), data = model.dat)

