# PURPOSE: to fit binomial GAM to mature (0,1) weighted by CPUE or unweighted

# Author: Emily Ryznar

# NOTES:
# Make sure you're using updated chela database compiled by shannon, though can't use it for everything bc
# it doesn't have sampling factor...

# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd
snow_cpue <- calc_cpue(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                       species = "SNOW",
                       years = years,
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)


# Join 1mm cpue to chela crab
mod.dat <- sh_chela %>%
              mutate(SIZE_1MM = floor(SIZE),
                     MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                          TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
              right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin
              mutate(YEAR = as.factor(YEAR)) %>%
              mutate(MATURE = case_when((SIZE <=35) ~ 0,
                                           (SIZE >= 135) ~ 1,
                                        TRUE ~ MATURE)) %>%
          na.omit()


# Fit unweighted binomial gam ----
mod.1 <- bam(MATURE ~ YEAR + s(SIZE), family = binomial(link = "logit"), data = mod.dat)

prd <- predict(mod.1, data.frame(YEAR = mod.dat$YEAR, SIZE = mod.dat$SIZE),
               type = "response", se.fit = TRUE)$fit

plot_dat<- data.frame(YEAR = mod.dat$YEAR, MATURE = mod.dat$MATURE, SIZE = mod.dat$SIZE, prd = prd)

ggplot(plot_dat, aes(SIZE, prd, color = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()

# # Predict on one year
# newdat<-data.frame(YEAR=2021,
#                    SIZE=seq(40,140,5))
# 
# prd <- predict(mod.1, newdat,
#                type = "response", se.fit = TRUE)$fit
# 
# plot_dat<- cbind(newdat, prd)
# 
# ggplot(plot_dat, aes(SIZE, prd))+
#   geom_line(linewidth = 1)+
#   theme_bw()

# Fit cpue-weighted binomial gam ----
mod.2 <- bam(MATURE ~ YEAR + s(SIZE), family = binomial(link = "logit"), data = mod.dat, weights = CPUE)

prd <- predict(mod.2, data.frame(YEAR = mod.dat$YEAR, SIZE = mod.dat$SIZE),
               type = "response", se.fit = TRUE)$fit

plot_dat2<- data.frame(YEAR = mod.dat$YEAR, MATURE = mod.dat$MATURE, SIZE = mod.dat$SIZE, prd = prd)

ggplot(plot_dat2, aes(SIZE, prd, color = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()

# # Predict on one year
# newdat<-data.frame(YEAR=2021,
#                    SIZE=seq(40,140,5))
# 
# prd <- predict(mod.2, newdat,
#                type = "response", se.fit = TRUE)$fit
# 
# plot_dat2<- cbind(newdat, prd)
# 
# ggplot(plot_dat2, aes(SIZE, prd))+
#   geom_line(linewidth = 1)+
#   theme_bw()

# 

# Combine predictions from both models
all.dat <- rbind(plot_dat %>% mutate(type = "Unweighted"), plot_dat2 %>% mutate(type = "CPUE-weighted"))


ggplot()+
  geom_line(all.dat %>% filter(type == "Unweighted"), mapping = aes(SIZE, prd, color = as.factor(1), group = YEAR), alpha = 0.5, linewidth = 1)+
  geom_line(all.dat %>% filter(type == "CPUE-weighted"), mapping = aes(SIZE, prd, color = as.factor(2), group = YEAR), alpha = 0.5, linewidth = 1)+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Unweighted", "CPUE-weighted"), name = "")+
  theme_bw()+
  ylab("Probability of terminal molt")+
  xlab("Carapace width (mm)")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.direction = "horizontal",
        legend.position = "bottom")

ggsave("./Figures/binomialGAMs_pmolt_predictions.png", width = 10, height = 8)
  
