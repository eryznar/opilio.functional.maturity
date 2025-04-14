# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
spec.dat <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(HAUL_TYPE !=17, SEX == 2, SHELL_CONDITION == 2) # filter for females, sh2, not HT17

plot.dat <- spec.dat %>% filter(YEAR>1988, CLUTCH_SIZE >0)
ggplot()+
  geom_histogram(plot.dat %>% filter(YEAR>1988), mapping = aes(SIZE_1MM), fill = "lightgrey", color = "black")+
  scale_x_continuous(breaks = seq(min(plot.dat$SIZE_1MM), max(plot.dat$SIZE_1MM), by = 10))+
  theme_bw()+
  facet_wrap(~YEAR, scales = "free")
  


haul <- readRDS("./Data/snow_survey_specimenEBS.rda")$haul %>%
  filter(#YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?,
    HAUL_TYPE !=17)

haul %>%
  dplyr::select(HAULJOIN, YEAR, MID_LATITUDE, MID_LONGITUDE) %>%
  rename(LATITUDE = MID_LATITUDE, LONGITUDE = MID_LONGITUDE)-> hauljoin


# Survey specimen data
prop.dat <- spec.dat %>%
  #filter(SIZE >= 15 & SIZE <= 85) %>%
  #dplyr::select(colnames(chela_10.13)) %>%
  mutate(MATURE = case_when((CLUTCH_SIZE > 0) ~ 1,
                            TRUE ~ 0),
         BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2,
         WIDTH_1MM = floor(SIZE)) %>%
  group_by(YEAR, WIDTH_1MM, SAMPLING_FACTOR, MATURE, HAULJOIN, MIDPOINT) %>%
  reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
  # right_join(., right_join(hauljoin, expand_grid(BIN = unique(cut_width(c(15:85), width = 10, center = 5, closed = "left")),
  #                                                MATURE = 0:1,
  #                                                YEAR = c(1975:2019, 2021:2024)), # also joining to unique hauljoins by year
  #                          relationship = "many-to-many")) %>%
  # replace_na(list(COUNT = 0))  %>%
  distinct() %>%
  group_by(YEAR, WIDTH_1MM, MATURE, MIDPOINT) %>%
  reframe(TOTAL_MATURE = sum(COUNT[MATURE == 1]), # calculate proportion mature by 1mm size bin
          TOTAL_CRAB = sum(COUNT),
          PROP_MATURE = TOTAL_MATURE/TOTAL_CRAB)

# Fit size at maturity model
yrs <- c(1989:2019, 2021:2024)
params <- data.frame()
preds <- data.frame()

for(ii in 1:length(yrs)){
  print(paste("Fitting year", yrs[ii]))
  
  #filter by year
  prop.dat %>%
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
  out <- data.frame(SPECIES = "SNOW",
                    REGION = "EBS",
                    DISTRICT = "ALL",
                    YEAR = yrs[ii],
                    A_EST = A_EST,
                    A_SE = A_SE,
                    B_EST = B_EST,
                    B_SE = B_SE)
  
  params <- rbind(params, out)
  
  #create prediction df
  out2 <- data.frame(SPECIES = "SNOW",
                     REGION = "EBS",
                     DISTRICT = "ALL",
                     YEAR = yrs[ii],
                     MIDPOINT = mod.dat$MIDPOINT,
                     PROP_MATURE = predict(mod, mod.dat))
  
  preds <- rbind(preds, out2)
  
}


write.csv(params, "./Output/female_maturity_model_params.csv")

##Plot maturity ogives
miss_yr <- c(2020)

dummy <- data.frame(SPECIES = "SNOW", REGION = "EBS", DISTRICT = "ALL", YEAR = miss_yr, A_EST = NA, A_SE = NA, B_EST = NA, B_SE = NA)

params <- rbind(params, dummy)

ggplot(preds, aes(MIDPOINT, PROP_MATURE, group = YEAR, color = YEAR))+
  geom_line(aes(color = YEAR), linewidth = 1)+
  theme_bw()+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")+
  ggtitle("Proportion mature at size")+
  ylab("Proportion mature (female)")+
  xlab("Carapace width (mm)")+
  scale_color_gradient2(low = "cadetblue", midpoint = 2006, high = "darkred")+
  #scale_color_viridis_d()+
  scale_x_continuous(breaks = seq(min(preds$MIDPOINT), max(preds$MIDPOINT), by = 5))+
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


