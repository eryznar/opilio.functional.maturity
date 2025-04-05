# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
df.dat <- read.csv("./Data/opilio_directedfishery_catch.csv") %>%
  dplyr::select(Year, Retained_kt) %>%
  rename(directedfish_biomass = Retained_kt)# directed fishery data

survey.dat <- read.csv("./Data/opilio_survey_malebiomass.csv") %>% # survey data
                mutate(small_male_biomass = gsub("\\s*\\([^\\)]+\\)", "", small_male_biomass),
                       large_male_biomass = gsub("\\s*\\([^\\)]+\\)", "", large_male_biomass),
                       large_male_biomass_sh2 = gsub("\\s*\\([^\\)]+\\)", "", large_male_biomass_sh2),
                       small_male_biomass = as.numeric(as.character(gsub(",", "", small_male_biomass))),
                       large_male_biomass = as.numeric(as.character(gsub(",", "", large_male_biomass))),
                       large_male_biomass_sh2 = as.numeric(as.character(gsub(",", "", large_male_biomass_sh2))),
                       small_male_biomass = small_male_biomass/1000,
                       large_male_biomass= large_male_biomass/1000,
                       large_male_biomass_sh2= large_male_biomass_sh2/1000) %>%
  filter(Year > 2000) # convert to kilotons 
                       #Year = c(1988:2019, 2021:2024) - 1) # lagging survey data back one year 

morph.dat <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE) %>%
  mutate(CUTOFF = BETA0 + BETA1*(log(SIZE_1MM)),
         MATURE = case_when((log(CHELA_HEIGHT) > CUTOFF) ~ 1,
                            TRUE ~ 0)) %>%
  dplyr::select(YEAR, SIZE_1MM, MATURE) %>%
  rename(Year = YEAR) #%>%
  # group_by(Year) %>%
  # mutate(size_at_mat = ifelse(MATURE == 1, mean(SIZE_1MM[MATURE == 1], na.rm = TRUE), NA))%>%
  # ungroup()

params <- read.csv("./Output/maturity_model_params.csv") %>%
  dplyr::select(YEAR, B_EST) %>%
  rename(Year = YEAR, size_at_mat = B_EST)

# Load March-April ice data
ice <- read.csv("./Data/ERA5ice_1972.2024.csv") %>%
  filter(name == "Mar-Apr ice cover") %>%
  dplyr::select(year, value) %>%
  rename(Year = year, MarApr_ice = value)


# Bin specimen data
bin.dat <- morph.dat %>%
  mutate(bin = case_when((SIZE_1MM>=55 & SIZE_1MM<=65) ~ "Small (55-65mm)", # 
                         (SIZE_1MM>=95 & SIZE_1MM<=105) ~ "Large (95-105mm)")) %>%
  #filter(is.na(bin) == FALSE) %>%
  group_by(Year, bin) %>%
  mutate(total_crab = n(),
         total_mature = sum(MATURE == 1),
         prop_mature = total_mature/total_crab) %>%
  ungroup()


right_join(df.dat, survey.dat) %>%
  right_join(., bin.dat) %>%
  right_join(., ice) %>%
  #filter(is.na(size_at_mat) == FALSE) %>%
  dplyr::select(!c(SIZE_1MM, MATURE)) %>%
  distinct() %>%
  right_join(., params) -> model.dat

# FIT MODELS ---------------------------------------------------------------------------------------------
# proportion mature in small bin (55-65) ----
small.dat <- model.dat %>% 
                filter(bin == "Small (55-65mm)") %>%
                distinct() %>%
            na.omit()

small.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 7) + s(small_male_biomass, k = 4) + s(large_male_biomass_sh2, k = 5) + s(MarApr_ice, k = 4), data = small.dat)

# diagnostics
summary(small.prop)
gam.check(small.prop)
appraise(small.prop)
draw(small.prop)

# proportion mature in large bin (95-105mm) ----
large.dat <- model.dat %>% 
              filter(bin == "Large (95-105mm)") %>%
              distinct() %>%
            na.omit()

large.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + s(small_male_biomass, k = 4) + s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), data = large.dat)

# diagnostics
summary(large.prop)
gam.check(large.prop)
appraise(large.prop)
draw(large.prop)

# size at maturity ----
SaM.dat <- model.dat %>%
              dplyr::select(!c(bin, total_crab, total_mature, prop_mature)) %>%
              distinct() %>%
           na.omit()

SaM.mod <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + s(small_male_biomass, k = 4) + s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), data = SaM.dat)

# diagnostics
summary(SaM.mod)
gam.check(SaM.mod)
appraise(SaM.mod)
draw(SaM.mod)

