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
                       large_male_biomass_sh2= large_male_biomass_sh2/1000, 
                       Year = c(1988:2019, 2021:2024) - 0) # lagging survey data back one year 

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

small.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 7) + s(small_male_biomass, k = 4) 
                  + s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), method =  "REML", data = small.dat)

small.prop.beta <- gam(prop_mature ~ s(directedfish_biomass, k = 7) + s(small_male_biomass, k = 4) 
                  + s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), family = betar(link = "logit"), method =  "REML", data = small.dat)


small.prop.beta.int <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + s(small_male_biomass, large_male_biomass_sh2, k = 4)
                         + s(MarApr_ice, k = 4), family = betar(link = "logit"), method = "REML", data = small.dat)

AICc(small.prop, small.prop.beta, small.prop.beta.int)


# diagnostics
summary(small.prop)
gam.check(small.prop)
appraise(small.prop)
draw(small.prop)

summary(small.prop.beta)
gam.check(small.prop.beta)
appraise(small.prop.beta)
draw(small.prop.beta)

summary(small.prop.beta.int)
gam.check(small.prop.beta.int)
appraise(small.prop.beta.int)
draw(small.prop.beta.int)


# proportion mature in large bin (95-105mm) ----
large.dat <- model.dat %>% 
              filter(bin == "Large (95-105mm)") %>%
              distinct() %>%
            na.omit()

large.prop <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + 
                    s(small_male_biomass, k = 4) + s(large_male_biomass, k = 4) + s(MarApr_ice, k = 4), method = "REML", data = large.dat)

large.prop.int <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + 
                    s(small_male_biomass, large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), method = "REML", data = large.dat)

large.prop.beta <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + 
                    s(small_male_biomass, k = 4) + s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), family = betar(link = "logit"), method = "REML",
                    data = large.dat)


large.prop.beta.int <- gam(prop_mature ~ s(directedfish_biomass, k = 4) + 
                       s(small_male_biomass, large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), family = betar(link = "logit"), method = "REML",
                     data = large.dat)

AICc(large.prop, large.prop.int, large.prop.beta, large.prop.beta.int)

# diagnostics
summary(large.prop)
gam.check(large.prop)
appraise(large.prop)
draw(large.prop)

summary(large.prop.int)
gam.check(large.prop.int)
appraise(large.prop.int)
draw(large.prop.int)

summary(large.prop.beta)
gam.check(large.prop.beta)
appraise(large.prop.beta)
draw(large.prop.beta)

summary(large.prop.beta.int)
gam.check(large.prop.beta.int)
appraise(large.prop.beta.int)
draw(large.prop.beta.int)

# size at maturity ----
SaM.dat <- model.dat %>%
              dplyr::select(!c(bin, total_crab, total_mature, prop_mature)) %>%
              distinct() %>%
           na.omit()

SaM.mod <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + s(small_male_biomass, k = 4) + 
                 s(large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), method = "REML", data = SaM.dat)
SaM.mod.int <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                     s(small_male_biomass, large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), method = "REML", data = SaM.dat)

SaM.mod.tw.int <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                     s(small_male_biomass, large_male_biomass_sh2, k = 4) + s(MarApr_ice, k = 4), family = tw(link = "log"), method = "REML", data = SaM.dat)

AICc(SaM.mod, SaM.mod.int, SaM.mod.tw.int)

# diagnostics
summary(SaM.mod)
gam.check(SaM.mod)
appraise(SaM.mod)
draw(SaM.mod)

summary(SaM.mod.int)
gam.check(SaM.mod.int)
appraise(SaM.mod.int)
draw(SaM.mod.int)

summary(SaM.mod.tw.int)
gam.check(SaM.mod.tw.int)
appraise(SaM.mod.tw.int)
draw(SaM.mod.tw.int)

