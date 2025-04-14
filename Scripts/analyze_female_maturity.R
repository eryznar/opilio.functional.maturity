# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------

# Survey specimen data
bin.dat <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  #dplyr::select(colnames(chela_10.13)) %>%
  filter(HAUL_TYPE !=17, SEX == 2, SHELL_CONDITION == 2) %>% # filter for females, sh2, not HT17
  mutate(MATURE = case_when((CLUTCH_SIZE > 0) ~ 1,
                            TRUE ~ 0)) %>%
  dplyr::select(YEAR, SIZE_1MM, MATURE, SAMPLING_FACTOR, CALCULATED_WEIGHT_1MM) %>%
  rename(Year = YEAR) %>%
  mutate(bin = case_when((SIZE_1MM>40 & SIZE_1MM<=50) ~ "Small (40-50mm)", # apply small, medium, and large bins for modeling below
                         (SIZE_1MM>50 & SIZE_1MM<=60) ~ "Medium (50-60mm)",
                         (SIZE_1MM>60 & SIZE_1MM<=70) ~ "Large (60-70mm)")) %>%
  filter(is.na(bin) == FALSE) %>%
  group_by(Year, bin) %>%
  reframe(total_crab = sum(SAMPLING_FACTOR),
          total_kg = sum(SAMPLING_FACTOR*CALCULATED_WEIGHT_1MM)/1000,
          total_mature = sum(SAMPLING_FACTOR[MATURE == 1]), # summing the sampling factor of only crab that are mature in each bin
          total_mature_kg = sum(SAMPLING_FACTOR[MATURE == 1] * CALCULATED_WEIGHT_1MM[MATURE == 1])/1000,
          prop_mature = total_mature/total_crab,
          prop_mature_kg = total_mature_kg/total_kg)  # calculate proportion mature in each bin

# Join binned survey specimen data with abundance from survey for same size bins
abund.dat <- read.csv("./Data/surveyabund_snowfemales40.70.csv") %>%
  dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT, bin) %>%
  rename(Year = YEAR, abundance = ABUNDANCE, biomass = BIOMASS_MT) %>%
  pivot_wider(., values_from = c("abundance", "biomass"), names_from = "bin") %>%
  rename(abundance_small = `abundance_Small (40-50mm)`, abundance_medium = `abundance_Medium (50-60mm)`, abundance_large = `abundance_Large (60-70mm)`,
         biomass_small = `biomass_Small (40-50mm)`, biomass_medium = `biomass_Medium (50-60mm)`, biomass_large = `biomass_Large (60-70mm)`)

abund.bin.dat <- bin.dat %>%
  dplyr::select(Year, bin, prop_mature, prop_mature_kg) %>%
  pivot_wider(., names_from = "bin", values_from = c(prop_mature, prop_mature_kg)) %>%
  rename(propmature_small = `prop_mature_Small (40-50mm)`, propmature_medium = `prop_mature_Medium (50-60mm)`, propmature_large = `prop_mature_Large (60-70mm)`,
         propmature_small_kg = `prop_mature_kg_Small (40-50mm)`, propmature_medium_kg = `prop_mature_kg_Medium (50-60mm)`, propmature_large_kg = `prop_mature_kg_Large (60-70mm)`) %>%
  right_join(., abund.dat) %>%
  mutate(mature_abundance_small = (abundance_small * propmature_small)/1e6,
         mature_abundance_medium = (abundance_medium * propmature_medium)/1e6,
         mature_abundance_large = (abundance_large * propmature_large)/1e6,
         mature_biomass_small = (biomass_small * propmature_small)/1000,
         mature_biomass_medium = (biomass_medium * propmature_medium)/1000,
         mature_biomass_large = (biomass_large * propmature_small)/1000,
         abundance_small = abundance_small/1e6,
         abundance_medium = abundance_medium/1e6,
         abundance_large = abundance_large/1e6,
         biomass_small = biomass_small/1000,
         biomass_medium = biomass_medium/1000,
         biomass_large= biomass_large/1000,
         total_matfem_abundance = mature_abundance_small+mature_abundance_medium+mature_abundance_large,
         total_matfem_biomass = mature_biomass_small+mature_biomass_medium+mature_biomass_large) %>%
  dplyr::select(Year, propmature_large, propmature_medium, propmature_small, propmature_large_kg, propmature_medium_kg, propmature_small_kg, 
                total_matfem_abundance, total_matfem_biomass)


# Size at 50% maturity timeseries from maturity model
params <- read.csv("./Output/female_maturity_model_params.csv") %>%
  dplyr::select(YEAR, B_EST) %>%
  rename(Year = YEAR, size_at_mat = B_EST)

# Load March-April ice data
ice <- read.csv("./Data/ERA5ice_1972.2024.csv") %>%
  filter(name == "Mar-Apr ice cover") %>%
  dplyr::select(year, value) %>%
  rename(Year = year, MarApr_ice = value) %>%
  mutate(Year = Year - 0) # ice lag

# Load male dat
male.dat <- read.csv("./Data/snow_male_GAM_modeldat.csv") %>%
  mutate(matmale_abundance_55.105 = mature_abundance_small+mature_abundance_large, 
         matmale_biomass_55.105 = mature_biomass_small+mature_biomass_large) %>%
  dplyr::select(Year, matmale_abundance_55.105, matmale_biomass_55.105)

male.survey.bioabund <- read.csv("./Data/allmalesnow_surveybioabund.csv") %>%
  dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT) %>%
  rename(Year = YEAR, all_male_abundance = ABUNDANCE, all_male_biomass = BIOMASS_MT) %>%
  mutate(all_male_abundance = all_male_abundance/1e6,
         all_male_biomass= all_male_biomass/1000)

male.dat <- right_join(male.dat, male.survey.bioabund)

# Bind all dataframes into df for modeling
model.dat <- right_join(abund.bin.dat, ice) %>%
  #filter(is.na(size_at_mat) == FALSE) %>%
  right_join(., male.dat) %>%
  right_join(., params) %>%
  mutate(total_matfem_abundance = log(total_matfem_abundance + 1),
         total_matfem_biomass = log(total_matfem_biomass + 1),
         lag1_matfem_abundance = lag(total_matfem_abundance, n = 1),
         lag1_matfem_biomass = lag(total_matfem_biomass, n = 1),
         lag1_matmale_abundance = lag(matmale_abundance_55.105, n = 1),
         lag1_matmale_biomass = lag(matmale_biomass_55.105, n = 1),
         lag1_allmale_abundance= lag(all_male_abundance, n = 1),
         lag1_allmale_biomass= lag(all_male_biomass, n = 1),
         lag3_MarApr_ice = lag(MarApr_ice, n = 3))



# FIT MODELS ---------------------------------------------------------------------------------------------
# proportion mature in small bin (40-50) ----
  # abundance covariates ----
  small.main.nolag.aa <- gam(propmature_small ~
                              s(total_matfem_abundance, k = 4) +
                              s(matmale_abundance_55.105, k = 4)+
                              s(all_male_abundance, k = 4)+
                              s(MarApr_ice, k = 4),
                            family = betar(link = "logit"),
                            method =  "REML",
                            data = model.dat)
  
  
  small.main.lagabund.aa <- gam(propmature_small ~
                                  s(lag1_matfem_abundance, k = 4) +
                                  s(lag1_matmale_abundance, k = 4)+
                                  s(lag1_allmale_abundance, k = 4)+
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"),
                                method =  "REML",
                                data = model.dat)
  
  small.main.icelag.aa <- gam(propmature_small ~
                                s(total_matfem_abundance, k = 4) +
                                s(matmale_abundance_55.105, k = 4)+
                                s(all_male_abundance, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  small.main.iceabundlag.aa <- gam(propmature_small ~
                                     s(lag1_matfem_abundance, k = 4) +
                                     s(lag1_matmale_abundance, k = 4)+
                                     s(lag1_allmale_abundance, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"),
                                   method =  "REML",
                                   data = model.dat)
  
  # biomass covariates ----
  small.main.nolag.bb <- gam(propmature_small_kg ~
                               s(total_matfem_biomass, k = 4) +
                               s(matmale_biomass_55.105, k = 4)+
                               s(all_male_biomass, k = 4)+
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"),
                             method =  "REML",
                             data = model.dat)
  
  
  small.main.lagabund.bb <- gam(propmature_small_kg ~
                                  s(lag1_matfem_biomass, k = 4) +
                                  s(lag1_matmale_biomass, k = 4)+
                                  s(lag1_allmale_biomass, k = 4)+
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"),
                                method =  "REML",
                                data = model.dat)
  
  small.main.icelag.bb <- gam(propmature_small_kg ~
                                s(total_matfem_biomass, k = 4) +
                                s(matmale_biomass_55.105, k = 4)+
                                s(all_male_biomass, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  small.main.iceabundlag.bb <- gam(propmature_small_kg ~
                                     s(lag1_matfem_biomass, k = 4) +
                                     s(lag1_matmale_biomass, k = 4)+
                                     s(lag1_allmale_biomass, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"),
                                   method =  "REML",
                                   data = model.dat)
  
  
  # model comparison ----
  AICc(small.main.nolag.aa, small.main.lagabund.aa, small.main.icelag.aa, small.main.iceabundlag.aa,
       small.main.nolag.bb, small.main.lagabund.bb, small.main.icelag.bb, small.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  
  # diagnostics ----
  summary(small.main.lagabund.aa)
  gam.check(small.main.lagabund.aa)
  appraise(small.main.lagabund.aa)
  draw(small.main.lagabund.aa)
  
  summary(small.main.icelag.aa)
  gam.check(small.main.icelag.aa)
  appraise(small.main.icelag.aa)
  draw(small.main.icelag.aa)
  
# proportion mature in medium bin (50-60) ----
  # abundance covariates ----
  medium.main.nolag.aa <- gam(propmature_medium ~
                               s(total_matfem_abundance, k = 4) +
                               s(matmale_abundance_55.105, k = 4)+
                               s(all_male_abundance, k = 4)+
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"),
                             method =  "REML",
                             data = model.dat)
  
  
  medium.main.lagabund.aa <- gam(propmature_medium ~
                                  s(lag1_matfem_abundance, k = 4) +
                                  s(lag1_matmale_abundance, k = 4)+
                                  s(lag1_allmale_abundance, k = 4)+
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"),
                                method =  "REML",
                                data = model.dat)
  
  medium.main.icelag.aa <- gam(propmature_medium ~
                                s(total_matfem_abundance, k = 4) +
                                s(matmale_abundance_55.105, k = 4)+
                                s(all_male_abundance, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  medium.main.iceabundlag.aa <- gam(propmature_medium ~
                                     s(lag1_matfem_abundance, k = 4) +
                                     s(lag1_matmale_abundance, k = 4)+
                                     s(lag1_allmale_abundance, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"),
                                   method =  "REML",
                                   data = model.dat)
  
  # biomass covariates ----
  medium.main.nolag.bb <- gam(propmature_medium_kg ~
                               s(total_matfem_biomass, k = 4) +
                               s(matmale_biomass_55.105, k = 4)+
                               s(all_male_biomass, k = 4)+
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"),
                             method =  "REML",
                             data = model.dat)
  
  
  medium.main.lagabund.bb <- gam(propmature_medium_kg ~
                                  s(lag1_matfem_biomass, k = 4) +
                                  s(lag1_matmale_biomass, k = 4)+
                                  s(lag1_allmale_biomass, k = 4)+
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"),
                                method =  "REML",
                                data = model.dat)
  
  medium.main.icelag.bb <- gam(propmature_medium_kg ~
                                s(total_matfem_biomass, k = 4) +
                                s(matmale_biomass_55.105, k = 4)+
                                s(all_male_biomass, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  medium.main.iceabundlag.bb <- gam(propmature_medium_kg ~
                                     s(lag1_matfem_biomass, k = 4) +
                                     s(lag1_matmale_biomass, k = 4)+
                                     s(lag1_allmale_biomass, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"),
                                   method =  "REML",
                                   data = model.dat)
  
  
  # model comparison ----
  AICc(medium.main.nolag.aa, medium.main.lagabund.aa, medium.main.icelag.aa, medium.main.iceabundlag.aa,
       medium.main.nolag.bb, medium.main.lagabund.bb, medium.main.icelag.bb, medium.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  
  # diagnostics ----
  summary(medium.main.icelag.aa)
  gam.check(medium.main.icelag.aa)
  appraise(medium.main.icelag.aa)
  draw(medium.main.icelag.aa)
  
  summary(medium.main.lagabund.aa)
  gam.check(medium.main.lagabund.aa)
  appraise(medium.main.lagabund.aa)
  draw(medium.main.lagabund.aa)
  
  
# proportion mature in large bin (60-70) ----
  # abundance covariates ----
  large.main.nolag.aa <- gam(propmature_large ~
                                s(total_matfem_abundance, k = 4) +
                                s(matmale_abundance_55.105, k = 4)+
                                s(all_male_abundance, k = 4)+
                                s(MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  large.main.lagabund.aa <- gam(propmature_large ~
                                   s(lag1_matfem_abundance, k = 4) +
                                   s(lag1_matmale_abundance, k = 4)+
                                   s(lag1_allmale_abundance, k = 4)+
                                   s(MarApr_ice, k = 4),
                                 family = betar(link = "logit"),
                                 method =  "REML",
                                 data = model.dat)
  
  large.main.icelag.aa <- gam(propmature_large ~
                                 s(total_matfem_abundance, k = 4) +
                                 s(matmale_abundance_55.105, k = 4)+
                                 s(all_male_abundance, k = 4)+
                                 s(lag3_MarApr_ice, k = 4),
                               family = betar(link = "logit"),
                               method =  "REML",
                               data = model.dat)
  
  
  large.main.iceabundlag.aa <- gam(propmature_large ~
                                      s(lag1_matfem_abundance, k = 4) +
                                      s(lag1_matmale_abundance, k = 4)+
                                      s(lag1_allmale_abundance, k = 4)+
                                      s(lag3_MarApr_ice, k = 4),
                                    family = betar(link = "logit"),
                                    method =  "REML",
                                    data = model.dat)
  
  # biomass covariates ----
  large.main.nolag.bb <- gam(propmature_large_kg ~
                                s(total_matfem_biomass, k = 4) +
                                s(matmale_biomass_55.105, k = 4)+
                                s(all_male_biomass, k = 4)+
                                s(MarApr_ice, k = 4),
                              family = betar(link = "logit"),
                              method =  "REML",
                              data = model.dat)
  
  
  large.main.lagabund.bb <- gam(propmature_large_kg ~
                                   s(lag1_matfem_biomass, k = 4) +
                                   s(lag1_matmale_biomass, k = 4)+
                                   s(lag1_allmale_biomass, k = 4)+
                                   s(MarApr_ice, k = 4),
                                 family = betar(link = "logit"),
                                 method =  "REML",
                                 data = model.dat)
  
  large.main.icelag.bb <- gam(propmature_large_kg ~
                                 s(total_matfem_biomass, k = 4) +
                                 s(matmale_biomass_55.105, k = 4)+
                                 s(all_male_biomass, k = 4)+
                                 s(lag3_MarApr_ice, k = 4),
                               family = betar(link = "logit"),
                               method =  "REML",
                               data = model.dat)
  
  
  large.main.iceabundlag.bb <- gam(propmature_large_kg ~
                                      s(lag1_matfem_biomass, k = 4) +
                                      s(lag1_matmale_biomass, k = 4)+
                                      s(lag1_allmale_biomass, k = 4)+
                                      s(lag3_MarApr_ice, k = 4),
                                    family = betar(link = "logit"),
                                    method =  "REML",
                                    data = model.dat)
  
  
  # model comparison ----
  AICc(large.main.nolag.aa, large.main.lagabund.aa, large.main.icelag.aa, large.main.iceabundlag.aa,
       large.main.nolag.bb, large.main.lagabund.bb, large.main.icelag.bb, large.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  
  # diagnostics ----
  summary(large.main.lagabund.aa)
  gam.check(large.main.lagabund.aa)
  appraise(large.main.lagabund.aa)
  draw(large.main.lagabund.aa)
  
  
  summary(large.main.lagabund.bb)
  gam.check(large.main.lagabund.bb)
  appraise(large.main.lagabund.bb)
  draw(large.main.lagabund.bb)
  
  
# size at maturity ----
  # abundance covariates ----
  SAM.main.nolag.aa <- gam(size_at_mat ~
                               s(total_matfem_abundance, k = 4) +
                               s(matmale_abundance_55.105, k = 4)+
                               s(all_male_abundance, k = 4)+
                               s(MarApr_ice, k = 4),
                             method =  "REML",
                             data = model.dat)
  
  
  SAM.main.lagabund.aa1 <- gam(size_at_mat ~
                                  s(lag1_matfem_abundance, k = 4) +
                                  s(lag1_matmale_abundance, k = 4)+
                                  s(lag1_allmale_abundance, k = 4)+
                                  s(MarApr_ice, k = 4),
                                method =  "REML",
                                data = model.dat)
  
  SAM.main.lagabund.aa2 <- gam(size_at_mat ~
                                 s(lag1_matfem_abundance, k = 4) +
                                 s(lag1_matmale_abundance, k = 4)+
                                 s(lag1_allmale_abundance, k = 4)+
                                 s(MarApr_ice, k = 4),
                               family = tw(link = "log"),
                               method =  "REML",
                               data = model.dat)
  
  SAM.main.icelag.aa <- gam(size_at_mat ~
                                s(total_matfem_abundance, k = 4) +
                                s(matmale_abundance_55.105, k = 4)+
                                s(all_male_abundance, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              method =  "REML",
                              data = model.dat)
  
  
  SAM.main.iceabundlag.aa <- gam(size_at_mat ~
                                     s(lag1_matfem_abundance, k = 4) +
                                     s(lag1_matmale_abundance, k = 4)+
                                     s(lag1_allmale_abundance, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   method =  "REML",
                                   data = model.dat)
  
  # biomass covariates ----
  SAM.main.nolag.bb <- gam(size_at_mat ~
                               s(total_matfem_biomass, k = 4) +
                               s(matmale_biomass_55.105, k = 4)+
                               s(all_male_biomass, k = 4)+
                               s(MarApr_ice, k = 4),
                             method =  "REML",
                             data = model.dat)
  
  
  SAM.main.lagabund.bb <- gam(size_at_mat ~
                                  s(lag1_matfem_biomass, k = 4) +
                                  s(lag1_matmale_biomass, k = 4)+
                                  s(lag1_allmale_biomass, k = 4)+
                                  s(MarApr_ice, k = 4),
                                method =  "REML",
                                data = model.dat)
  
  SAM.main.icelag.bb <- gam(size_at_mat ~
                                s(total_matfem_biomass, k = 4) +
                                s(matmale_biomass_55.105, k = 4)+
                                s(all_male_biomass, k = 4)+
                                s(lag3_MarApr_ice, k = 4),
                              method =  "REML",
                              data = model.dat)
  
  
  SAM.main.iceabundlag.bb <- gam(size_at_mat ~
                                     s(lag1_matfem_biomass, k = 4) +
                                     s(lag1_matmale_biomass, k = 4)+
                                     s(lag1_allmale_biomass, k = 4)+
                                     s(lag3_MarApr_ice, k = 4),
                                   method =  "REML",
                                   data = model.dat)
  
  
  # model comparison ----
  AICc(SAM.main.nolag.aa, SAM.main.lagabund.aa1, SAM.main.lagabund.aa2, SAM.main.icelag.aa, SAM.main.iceabundlag.aa,
       SAM.main.nolag.bb, SAM.main.lagabund.bb, SAM.main.icelag.bb, SAM.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  # diagnostics ----
  summary(SAM.main.lagabund.aa)
  gam.check(SAM.main.lagabund.aa)
  appraise(SAM.main.lagabund.aa)
  draw(SAM.main.lagabund.aa)
  
  summary(SAM.int.icelag.aa)
  gam.check(SAM.int.icelag.aa)
  appraise(SAM.int.icelag.aa)
  draw(SAM.int.icelag.aa)
  
  summary(SAM.main.icelag.bb)
  gam.check(SAM.main.icelag.bb)
  appraise(SAM.main.icelag.bb)
  draw(SAM.main.icelag.bb)
  
  summary(SAM.int.icelag.bb)
  gam.check(SAM.int.icelag.bb)
  appraise(SAM.int.icelag.bb)
  draw(SAM.int.icelag.bb)


