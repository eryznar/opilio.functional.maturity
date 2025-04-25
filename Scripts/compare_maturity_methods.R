# Set params ----
current_year <- 2024
years <- c(1989:current_year)

# Load data ----
# chela <- read.csv("./Data/snow_chela_UPDATED.csv") %>% # already != HT 17, only shell 2, no special projects
#   filter(YEAR %in% years) %>%
#   mutate(RATIO = SIZE/CHELA_HEIGHT) %>% 
#   filter(RATIO > 2 & RATIO < 35) %>% # filter extreme measurements
#   mutate(LN_CW = log(SIZE),
#          LN_CH = log(CHELA_HEIGHT)) %>%
#   #filter(LN_CW >= 3.9 & LN_CW <= 4.6) %>% # these were filtered out by Jon
#   mutate(BIN = cut_width(LN_CW, width = 0.025, center = 0.0125, closed = "left", dig.lab = 4),
#          BIN2 = BIN) %>%
#   separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
#   mutate(LOWER = as.numeric(sub('.', '', LOWER)),
#          UPPER = as.numeric(gsub('.$', '', UPPER)),
#          MIDPOINT = (UPPER + LOWER)/2) %>%
#   mutate(cc = case_when((YEAR == 2024)~"recent",
#                         TRUE ~ "not")) %>%
#   mutate(cutoff = BETA0 + BETA1*(LN_CW), # apply cutline model
#          mature = case_when((LN_CH > cutoff) ~ 1,
#                             TRUE ~ 0))

minima <- read.csv("./Output/opilio_cutline_minima.csv")


# crabpack specimen
cp_chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
          #dplyr::select(colnames(chela_10.13)) %>%
          filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
                 YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
          mutate(ratio = SIZE/CHELA_HEIGHT,
                 LN_CH = log(CHELA_HEIGHT),
                 LN_CW = log(SIZE)) %>%
          filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
          dplyr::select(!c(ratio)) %>%
          mutate(cutoff = BETA0 + BETA1*LN_CW, # apply cutline model
                 MATURE = case_when((LN_CH > cutoff) ~ 1,
                                    TRUE ~ 0))

# Plot 1: cutline with 2024 highlighted ----
# Fit linear model to evaluate ln(CH) minima (imm/mat division) and ln(CW) bin midpoint
mod <- lm(minima ~ midpoint, minima)

prd <- predict(mod, data = min.dat, se.fit=TRUE, interval="confidence", level=0.95)

cbind(min.dat %>% dplyr::select(!c(lwr, upr)), prd$fit) -> min.dat2

labs <- data.frame(x = 4, y = 3, lab = paste0("R-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001"))


ggplot()+
  #geom_point(min.dat2, mapping = aes(midpoint, minima), color = "red")+
  geom_point(cp_chela, mapping = aes(LN_CW, LN_CH, color = as.factor(MATURE)), alpha = 0.5)+
  #geom_density(plot.dat, mapping = aes(LN_CH), color = "red", linetype = "dashed", linewidth = 1)+
  theme_bw()+
  facet_wrap(~YEAR)+
  #geom_line(min.dat2, mapping = aes(midpoint, fit), linewidth = 1.5, color = "black", linetype = "dashed")+
  #geom_ribbon(min.dat2, mapping = aes(x = midpoint, ymin = lwr, ymax = upr), alpha = 0.4, fill = "black")+
  #annotate("text", x = 4, y = 3, label = paste0("R-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001"))+
  ylab("ln(chela height)")+
  xlab("ln(carapace width)")+
  geom_abline(slope = BETA1, intercept = BETA0, linewidth= 1, linetype = "dashed")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("./Figures/year_facetted_cutline.png", width = 11, height = 8.5)

# Plot 2: legacy weighting of 0s and 1s by sampling factor per 10mm bin and 5mm bin
# Subset data into size intervals at ln(CW) of 10mm and add in cutline params
cp_chela %>%
  mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2,
         SIZE_1MM = floor(SIZE)) %>%
  group_by(YEAR, BIN, MIDPOINT, MATURE) %>%
  # Sum n_crabs per size bin
  reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
  mutate(bin.size = "10mm") ->  bin.10mm

cp_chela %>%
  mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 10mm
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2,
         SIZE_1MM = floor(SIZE)) %>%
  group_by(YEAR, BIN, MIDPOINT, MATURE) %>%
  # Sum n_crabs per size bin
  reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
  mutate(bin.size = "5mm") ->  bin.5mm


rbind(bin.10mm, bin.5mm) -> plot.dat


  ggplot(plot.dat %>% filter(bin.size == "10mm"), aes(MIDPOINT, MATURE, size = log(COUNT + 10))) +
    geom_point(alpha = 0.5)+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = sampling factor, 10mm bins")+
    scale_size_continuous(breaks = seq(2,8, by = 1), name = "log(SF + 10)")+
    theme_bw()+
     xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_SFweight_10mm.png", height = 8.5, width = 11, units = "in")
  
  ggplot(plot.dat %>% filter(bin.size == "5mm"), aes(MIDPOINT, MATURE, size = log(COUNT + 10))) +
    geom_point(alpha = 0.5, color = "cadetblue")+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = sampling factor, 5mm bins")+
    scale_size_continuous(breaks = seq(2,8, by = 1), name = "log(SF + 10)")+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_SFweight_5mm.png", height = 8.5, width = 11, units = "in")
  

# Plot 3: Proposed weighting of 0s and 1s by CPUE by 1mm bin
  # Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd
  snow_cpue <- calc_cpue(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                         species = "SNOW",
                         years = years,
                         sex = "male",
                         shell_condition = "new_hardshell",
                         bin_1mm = TRUE) %>%
               dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)
   
  
  cp_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
     mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
            BIN2 = BIN) %>%
       separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
       mutate(LOWER = as.numeric(sub('.', '', LOWER)),
              UPPER = as.numeric(gsub('.$', '', UPPER)),
              MIDPOINT = (UPPER + LOWER)/2) %>%
       group_by(YEAR, BIN, MIDPOINT, MATURE, MAT_TEXT) %>%
       # Sum n_crabs per size bin
       reframe(CPUE = sum(CPUE)) %>%
       mutate(bin.size = "10mm") ->  bin.10mm
  
  
  cp_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
    mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 5mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2) %>%
    group_by(YEAR, BIN, MIDPOINT, MATURE, MAT_TEXT) %>%
    # Sum n_crabs per size bin
    reframe(CPUE = sum(CPUE)) %>%
    mutate(bin.size = "5mm") ->  bin.5mm
  
  
  
  
  rbind(bin.10mm, bin.5mm) %>%
    na.omit()-> plot.dat
  
  
  ggplot(plot.dat %>% filter(!YEAR %in% c(2008, 2010, 2013, 2014), bin.size == "10mm"), aes(MIDPOINT, MATURE, size = log(CPUE + 10))) +
    geom_point(alpha = 0.5)+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = CPUE, 10mm bins")+
    scale_size_continuous(breaks = seq(2,8, by = 1), name = "log(CPUE + 10)")+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_CPUEweight_10mm.png", height = 8.5, width = 11, units = "in")
  
  ggplot(plot.dat %>% filter(!YEAR %in% c(2008, 2010, 2013, 2014), bin.size == "5mm"), aes(MIDPOINT, MATURE, size = log(CPUE + 10))) +
    geom_point(alpha = 0.5,  color = "cadetblue")+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = CPUE, 5mm bins")+
    scale_size_continuous(breaks = seq(2,8, by = 1), name = "log(CPUE + 10)")+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_CPUEweight_5mm.png", height = 8.5, width = 11, units = "in")
  
  