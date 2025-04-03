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

bin.dat <- morph.dat %>%
  mutate(bin = case_when((SIZE_1MM>=55 & SIZE_1MM<=65) ~ "Small (55-65mm)", # 
                         (SIZE_1MM>=95 & SIZE_1MM<=105) ~ "Large (95-105mm)")) %>%
  filter(is.na(bin) == FALSE) %>%
  group_by(Year, bin) %>%
  mutate(total_crab = n(),
         total_mature = sum(MATURE == 1),
         prop_mature = total_mature/total_crab) %>%
  ungroup()



right_join(df.dat, survey.dat) %>%
  right_join(., bin.dat) -> model.dat

# FIT MODELS ---------------------------------------------------------------------------------------------
model.dat %>%
  filter(is.na(size_at_mat) == FALSE) %>%
  dplyr::select(Year, directedfish_biomass, small_male_biomass, large_male_biomass, size_at_mat, bin) %>%
  distinct() -> model.dat2

# small bin (55-65)
SaM.dat <- model.dat2 %>% filter(bin == "Small (55-65mm)") %>%
                distinct()

prop.dat <- model.dat %>% filter(bin == "Small (55-65mm)") %>%
              dplyr::select(!c(SIZE_1MM, MATURE, size_at_mat)) %>%
              distinct()

small.SaM <- gam(size_at_mat ~ s(directedfish_biomass, k = 3) + s(small_male_biomass, k = 3) + s(large_male_biomass,  k = 3), data = SaM.dat)
small.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 3) + s(small_male_biomass, k = 3) + s(large_male_biomass, k = 3), data = prop.dat)

# large bin (95-105mm)
SaM.dat <- model.dat2 %>% filter(bin == "Large (95-105mm)") %>%
  distinct()

prop.dat <- model.dat %>% filter(bin == "Large (95-105mm)") %>%
  dplyr::select(!c(SIZE_1MM, MATURE, size_at_mat)) %>%
  distinct()

large.SaM <- gam(size_at_mat ~ s(directedfish_biomass, k = 3) + s(small_male_biomass, k = 3) + s(large_male_biomass,  k = 3), data = SaM.dat)
large.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 3) + s(small_male_biomass, k = 3) + s(large_male_biomass, k = 3), data = prop.dat)

