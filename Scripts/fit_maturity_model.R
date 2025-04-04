# PURPOSE: to use morphometric maturity cutlines to define immature vs. mature for snow crab shell 2 males

# Author: Emily Ryznar, Jon Richar, Shannon Hennessey

# GENERAL WORKFLOW:
# 1) Calculate distribution minima between ln(chela height) density clouds by ln(carapace width) interval
# 2) Run linear model between minima and chela height midpoints to create cutlines (extract model params)
    # could maybe try fitting a nonlinear model to determine cutline, wouldn't necessitate log trans
# 3) Use cutline params to define mature/immature on chela-measured crab from survey data
# 4) Fit nls model between proportion mature and carapace width midpoint by year, extract parameters
# 5) Apply nls model parameters (ogives) to all male survey crab to determine maturity

# LOAD LIBS/PARAMS ----
source("./Scripts/load_libs_params.R")

# Read in maturity output from crabpack
ratio <- readRDS("./Data/snow_survey_maturityEBS.rda")$male_mat_ratio
pars <- readRDS("./Data/snow_survey_maturityEBS.rda")$model_parameters

# Read in minima data and run model to create cutlines
minima <- read.csv("./Output/opilio_cutline_minima.csv") %>%
          mutate(BETA0 = coef(lm(minima ~ midpoint))[1],
                 BETA1 = coef(lm(minima ~ midpoint))[2])
#   
# minima <- read.csv("./Data/Jon_minima.csv") %>%
#   mutate(BETA0 = coef(lm(y ~ x))[1],
#          BETA1 = coef(lm(y ~ x))[2])

#mod <- lm(y ~ x, dat = minima)

# Read and filter survey data
haul <- readRDS("./Data/snow_survey_specimenEBS.rda")$haul %>%
  filter(#YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?,
          HAUL_TYPE !=17)

haul %>%
  dplyr::select(HAULJOIN, YEAR, MID_LATITUDE, MID_LONGITUDE) %>%
  rename(LATITUDE = MID_LATITUDE, LONGITUDE = MID_LONGITUDE)-> hauljoin

chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SEX == 1,
         SHELL_CONDITION == 2, 
         HAUL_TYPE !=17, 
         #YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?
         is.na(CHELA_HEIGHT) == FALSE) %>%
  dplyr::select(HAULJOIN, SPECIES, REGION, DISTRICT, YEAR, SIZE, CHELA_HEIGHT, SAMPLING_FACTOR)

# Read in Jon data for 2010 and 2013
jon.dat <- read.csv("./Data/Jon_opilio_chela_height_TS.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE, 1, 4)),
         SPECIES = "SNOW",
         REGION = "EBS",
         DISTRICT = "ALL") %>%
  dplyr::filter(YEAR %in% c(2010, 2013), 
                #HAUL_TYPE !=17, need to make sure this was filtered by Jon
                SEX == 1, 
                SHELL_CONDITION == 2) %>%
  dplyr::select(HAULJOIN, SPECIES, REGION, DISTRICT, YEAR, WIDTH, CHELA_HEIGHT, SAMPLING_FACTOR) %>%
  rename(SIZE = WIDTH)

# Join survey data and Jon's data from 2010 and 2013
chela <- rbind(jon.dat, chela)

# Pool and transform data using natural log
chela %>%
  mutate(LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT)) -> chela2

# Subset data into size intervals at ln(CW) of 10mm and add in cutline params
chela2 %>%
  #filter(LN_CW >= 3.9 & LN_CW <= 4.6) %>% # these were filtered out by Jon
  mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),
         BIN2 = BIN,
         BETA0 = unique(minima$BETA0), # add in cutline params
         BETA1 = unique(minima$BETA1)) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2) -> bin.dat

# Calculate proportion mature by size bin
prop.dat <- bin.dat %>%
            mutate(MATURITY = case_when(LN_CH >= (BETA0 + BETA1*LN_CW) ~ "MATURE",
                                        LN_CH < (BETA0 + BETA1*LN_CW) ~ "IMMATURE"),
                   WIDTH_1MM = floor(SIZE)) %>%
            unite("AREA", c(SPECIES, REGION, DISTRICT), sep = "_") %>%
            group_by(AREA, YEAR, HAULJOIN, BIN, MATURITY) %>%
            # Sum n_crabs per size bin
            reframe(COUNT = sum(SAMPLING_FACTOR))  %>%
            #dplyr::select(AREA, YEAR, HAULJOIN, BIN, MATURITY, COUNT) %>%
            # Expand/add 0-count size bins
            right_join(., right_join(hauljoin, expand_grid(BIN = unique(cut_width(c(0:200), width = 10, center = 5, closed = "left")),
                                             AREA = c("SNOW_EBS_ALL"), 
                                             MATURITY = c("MATURE", "IMMATURE"),
                                             YEAR = c(1989:2019, 2021:2024)), # also joining to unique hauljoins by year
                                  relationship = "many-to-many")) %>%
            replace_na(list(COUNT = 0)) %>%
            pivot_wider(names_from = MATURITY, values_from = COUNT) %>%
            mutate(TOTAL_CRAB= IMMATURE + MATURE,
                   PROP_MATURE = MATURE/TOTAL_CRAB) %>%
            separate(BIN, sep = ",", into = c("LOWER", "UPPER")) %>%
            mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                   UPPER = as.numeric(gsub('.$', '', UPPER)),
                   MIDPOINT = (UPPER + LOWER)/2,
                   PROP_MATURE = case_when((LOWER <= 40) ~ 0, # set all crab in 40-50 bin or smaller to immature; previous code assigns 0 and 1 to crab if no crab were caught
                                           (LOWER >= 130) ~ 1, # set all crab in 130-140 bin or larger to mature
                                           TRUE ~ PROP_MATURE),
                   SIZE_BIN = paste0(LOWER, "-", UPPER)) %>%
            separate(AREA, sep = "_", into = c("SPECIES", "REGION", "DISTRICT")) %>%
            # Remove specific years: EBS snow - 2008, 2012, 2014, 2016, 2020; tanner - 1989, 2013, 2015, 2020, tannerE - 2011 too
            # filter(!((SPECIES == "TANNER" & YEAR %in% c(1989, 2013, 2015, 2020)) | (DISTRICT == "E166" & YEAR == 2011) |
            #            (SPECIES == "SNOW" & REGION == "EBS" & YEAR %in% c(2008, 2012, 2014, 2016, 2020)))) %>%
            # Format final output file
            rename(NUM_IMMATURE = IMMATURE,
                   NUM_MATURE = MATURE) %>%
            #replace_na(list(PROP_MATURE = 0)) %>% # SHOULD NOT REPLACE NAs with ZERO IF THERE WERE NO CRAB CAUGHT
            dplyr::select(SPECIES, REGION, DISTRICT, YEAR, SIZE_BIN, MIDPOINT, LATITUDE, LONGITUDE, NUM_IMMATURE, NUM_MATURE, TOTAL_CRAB, PROP_MATURE)
          
write.csv(prop.dat, "./Output/opilio_propmat_latlon.csv")

# Aggregate data by year
agg.prop.dat <- prop.dat %>%
                group_by(SPECIES, REGION, DISTRICT, YEAR, SIZE_BIN, MIDPOINT) %>%
                reframe(NUM_IMMATURE = sum(NUM_IMMATURE),
                        NUM_MATURE = sum(NUM_MATURE)) %>%
                mutate(TOTAL_CRAB = NUM_IMMATURE + NUM_MATURE,
                       PROP_MATURE = NUM_MATURE/TOTAL_CRAB,
                       LOWER = as.numeric(sub("-.*", "", SIZE_BIN)),
                       UPPER =  as.numeric(sub(".*-", "", SIZE_BIN)),
                       PROP_MATURE = case_when((LOWER < 40) ~ 0, # set all crab in 40-50 bin smaller to immature; previous code assigns 0 and 1 to crab if no crab were caught
                                               (LOWER >= 130) ~ 1, # set all crab in 130-140 bin or larger to mature
                                               TRUE ~ PROP_MATURE)) %>%
                #replace_na(list(PROP_MATURE = 0)) %>% # SHOULD NOT REPLACE NAs with ZERO IF THERE WERE NO CRAB CAUGHT
                dplyr::select(!c(LOWER, UPPER)) %>%
                group_by(YEAR) %>%
                arrange(., MIDPOINT, .by_group = TRUE) %>%
                ungroup()

write.csv(agg.prop.dat, "./Output/opilio_propmat_agg.csv")


## Fit size at maturity model --------------------------------------------------
# Jon filters out 2008, 2012, 2014, 2016 due to non-convergence?
# 2020 is omitted bc we didn't have a survey
# 2010 and 2013 are problematic for me I don't have any data (only via special projects?)

yrs <- c(1989:2007,2009:2011, 2013, 2015, 2017:2019, 2021:2024)
params <- data.frame()
preds <- data.frame()

for(ii in 1:length(yrs)){
  print(paste("Fitting year", yrs[ii]))
  
  #filter by year
  agg.prop.dat %>%
    filter(YEAR == yrs[ii]) -> mod.dat
  
  #fit nls model
  mod <- nls(PROP_MATURE ~ (1/(1 + exp(-a*(as.numeric(MIDPOINT) - b)))),
              data = mod.dat,
              start = list(a = 0.10, b = 60.0),
              na.action = na.omit, 
              nls.control(maxiter = 5000))
  
  #pull out params
  A_EST <- summary(mod)$coefficients[1,1]
  A_SE <- summary(mod)$coefficients[1,2]
  B_EST <- summary(mod)$coefficients[2,1]
  B_SE <- summary(mod)$coefficients[2,2]
  
  #create summary df
  out <- data.frame(SPECIES = unique(mod.dat$SPECIES),
                    REGION = unique(mod.dat$REGION),
                    DISTRICT = unique(mod.dat$DISTRICT),
                    YEAR = yrs[ii],
                    A_EST = A_EST,
                    A_SE = A_SE,
                    B_EST = B_EST,
                    B_SE = B_SE)
  
  params <- rbind(params, out)
  
  #create prediction df
  out2 <- data.frame(SPECIES = unique(mod.dat$SPECIES),
                    REGION = unique(mod.dat$REGION),
                    DISTRICT = unique(mod.dat$DISTRICT),
                    YEAR = yrs[ii],
                    MIDPOINT = mod.dat$MIDPOINT,
                    PROP_MATURE = predict(mod, mod.dat))
  
  preds <- rbind(preds, out2)
  
}

write.csv(params, "./Output/maturity_model_params.csv")

##Plot maturity ogives
miss_yr <- c(2008, 2012, 2014, 2016, 2020)

dummy <- data.frame(SPECIES = "SNOW", REGION = "EBS", DISTRICT = "ALL", YEAR = miss_yr, A_EST = NA, A_SE = NA, B_EST = NA, B_SE = NA)

params <- rbind(params, dummy)

ggplot(preds, aes(MIDPOINT, PROP_MATURE, group = YEAR, color = "YEAR"))+
  geom_line(aes(color = YEAR))+
  theme_bw()+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")+
  ggtitle("Proportion mature at size")+
  ylab("Proportion mature")+
  xlab("Survey year")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

## Plot average size at 50% maturity
ggplot(params, aes(YEAR, B_EST, group = 1))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin = B_EST - B_SE, ymax = B_EST + B_SE)) +
  theme_bw()+
  ylab("Carapace width (mm)")+
  xlab("Survey year")+
  ggtitle("Snow crab size at 50% maturity")+
  scale_x_continuous(breaks = seq(min(params$YEAR), max(params$YEAR), by = 5))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

miss_yr <- c(2013, 2015, 2020)

dummy <- data.frame(SPECIES = "TANNER", REGION = "EBS", DISTRICT = rep(c("ALL", "E166", "W166"), each =3), YEAR = rep(miss_yr, 3), A_EST = NA, A_SE = NA, B_EST = NA, B_SE = NA)

t.params <- rbind(t.params, dummy)

ggplot()+
  geom_line(t.params %>% filter(DISTRICT == "W166", YEAR > 1989), mapping = aes(YEAR, scale(B_EST)[,1], group = 1), color = "blue")+
  geom_point(size = 2)+
  #geom_errorbar(t.params %>% filter(DISTRICT == "ALL", YEAR > 1989), mapping = aes(x = YEAR, ymin = B_EST - B_SE, ymax = B_EST + B_SE), color = "blue") +
  geom_line(params, mapping = aes(YEAR, scale(B_EST)[,1], group = 1))+
  #geom_errorbar(params, mapping = aes(x = YEAR, ymin = B_EST - B_SE, ymax = B_EST + B_SE)) +
  theme_bw()+
  ylab("Carapace width (mm)")+
  xlab("Survey year")+
  ggtitle("Tanner crab size at 50% maturity")+
  scale_x_continuous(breaks = seq(min(params$YEAR), max(params$YEAR), by = 5))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggplot(t.params %>% filter(DISTRICT == "W166", YEAR > 1989), aes(YEAR, B_EST, group = 1))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin = B_EST - B_SE, ymax = B_EST + B_SE)) +
  theme_bw()+
  ylab("Carapace width (mm)")+
  xlab("Survey year")+
  ggtitle("Tanner crab size at 50% maturity")+
  scale_x_continuous(breaks = seq(min(params$YEAR), max(params$YEAR), by = 5))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

var.1 <- t.params %>% filter(DISTRICT == "ALL", YEAR > 1989)
var.2 <- params %>% filter(YEAR > 1989)

cor.test(scale(var.1$B_EST)[,1], scale(var.2$B_EST)[,1])
