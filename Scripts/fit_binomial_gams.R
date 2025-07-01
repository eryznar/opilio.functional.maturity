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
        mutate(MATURE = case_when((SIZE <40) ~ 0,
                                  (SIZE >= 115) ~ 1,
                                  TRUE ~ MATURE),
               log.CPUE = as.integer(round(log(CPUE+10))),
               CPUE = as.integer(round(CPUE))) %>%
        na.omit() 

# Binomial GAM for loop to fit by year
yrs <- unique(mod.dat$YEAR)
plot.dat <- data.frame()

for(ii in 1:length(yrs)){
  mod.dat %>%
    filter(YEAR %in% yrs[ii]) -> mod.dat2
  
  print(paste0("Fitting ", yrs[ii]))
  
  new.dat <- data.frame(YEAR = unique(mod.dat2$YEAR), SIZE = seq(min(mod.dat2$SIZE), max(mod.dat2$SIZE), by = 1))
  
  mod.1 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2)
  prd.1 <- predict(mod.1, new.dat, type = "response", se.fit = TRUE)
  
  mod.2 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = CPUE)
  prd.2 <- predict(mod.2, new.dat, type = "response", se.fit = TRUE)
  
  
  mod.3 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = log.CPUE)
  prd.3 <- predict(mod.3, new.dat, type = "response", se.fit = TRUE)
  
  
  out <- rbind(data.frame(new.dat, prd = prd.1, type = "Unweighted"),
               data.frame(new.dat, prd = prd.2, type = "CPUE-weighted"),
               data.frame(new.dat, prd = prd.3, type = "log(CPUE)-weighted")) %>%
    mutate(prd.fit = case_when((SIZE <=44) ~ 0,
                               (SIZE >= 115) ~ 1,
                               TRUE ~ prd.fit),
           N = nrow(mod.dat2))
  
  plot.dat <- rbind(plot.dat, out)
}

labelz <- paste0(unique(plot.dat$YEAR), " \n(N=", unique(plot.dat$N), ")")
names(labelz) <- unique(plot.dat$YEAR)

ggplot()+
  geom_line(plot.dat %>% filter(type == "Unweighted"), mapping = aes(SIZE, prd.fit, color = as.factor(1), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "CPUE-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(2), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "log(CPUE)-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(3), group = YEAR), linewidth = 1)+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "Log(CPUE)-weighted"), name = "")+
  scale_fill_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "Log(CPUE)-weighted"), name = "")+
  theme_bw()+
  facet_wrap(~YEAR, labeller = labeller(YEAR = labelz))+
  # geom_ribbon(plot.dat %>% filter(type == "Unweighted"),
  #           mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(1)), alpha = 0.25)+
  # geom_ribbon(plot.dat %>% filter(type == "CPUE-weighted"),
  #             mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(2)), alpha = 0.25)+
  # ylab("Probability of terminal molt")+
  xlab("Carapace width (mm)")+
  scale_x_continuous(limits = c(25, 138), breaks = seq(25, 140, by = 25))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 11))

ggsave("./Figures/log_cpue_untwd_ogives.png", width = 11, height = 8.5, units = "in")

 
## Solve for 50% maturity ----
# Specify function
find_l50 <- function(year_val, model) {
  f <- function(size_val) {
    predict(model, newdata = data.frame(SIZE = size_val, YEAR = year_val), type = "response") - 0.5
  }
  l50 <- uniroot(f, lower = min(mod.dat$SIZE), upper = max(mod.dat$SIZE))$root
  
  df <- data.frame(year = year_val, L50 = l50)
  
  return(df)
}

# Run functions
unique(mod.dat$YEAR) %>%
purrr::map_df(~find_l50(., mod.1)) %>%
  mutate(type = "Unweighted") -> uw.out


unique(mod.dat$YEAR) %>%
  purrr::map_df(~find_l50(., mod.2)) %>%
  mutate(type = "CPUE-weighted") -> w.out

# Combine and plot
l50.df <- rbind(uw.out, w.out) %>%
  mutate(year = as.numeric(as.character(year)))

all.years <- data.frame(year = seq(min(l50.df$year), max(l50.df$year), by = 1))

missing.yrs <- right_join(l50.df, all.years) %>%
                filter(is.na(L50) == TRUE) %>%
                pull(year)

l50.df <- rbind(l50.df, data.frame(year = missing.yrs, L50 = rep(NA, 2*length(missing.yrs)), 
                         type = rep(c("CPUE-weighted", "Unweighted"), each = length(missing.yrs))) )
 
 
ggplot(l50.df, aes(year, L50, group = type, color = type))+
  geom_point()+
  geom_line()+
  xlab("Year")+
  theme_bw()

## SAME WORKFLOW AS ABOVE BUT WITH CRABPACK DATA TO COMPARE WITH SAMPLING FACTOR----

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
mod.dat <- cp_chela %>%
  mutate(SIZE_1MM = floor(SIZE),
         MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                              TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
  right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin
  mutate(YEAR = as.factor(YEAR)) %>%
  mutate(MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE),
         CPUE = as.integer(round(CPUE)),
         SAMPLING_FACTOR = as.integer(round(SAMPLING_FACTOR))) %>%
  filter(is.na(MATURE) == FALSE)


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

# Fit SF-weighted binomial gam ----
mod.3 <- bam(MATURE ~ YEAR + s(SIZE), family = binomial(link = "logit"), data = mod.dat, weights = SAMPLING_FACTOR)

prd <- predict(mod.3, data.frame(YEAR = mod.dat$YEAR, SIZE = mod.dat$SIZE),
               type = "response", se.fit = TRUE)$fit

plot_dat3<- data.frame(YEAR = mod.dat$YEAR, MATURE = mod.dat$MATURE, SIZE = mod.dat$SIZE, prd = prd)

ggplot(plot_dat3, aes(SIZE, prd, color = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()


# Combine predictions from both models
all.dat <- rbind(plot_dat %>% mutate(type = "Unweighted"), plot_dat2 %>% mutate(type = "CPUE-weighted"),
                 plot_dat3 %>% mutate(type = "SF-weighted"))


ggplot()+
  geom_line(all.dat %>% filter(type == "Unweighted"), mapping = aes(SIZE, prd, color = as.factor(1), group = YEAR), linewidth = 1)+
  geom_line(all.dat %>% filter(type == "CPUE-weighted"), mapping = aes(SIZE, prd, color = as.factor(2), group = YEAR),  linewidth = 1)+
  geom_line(all.dat %>% filter(type == "SF-weighted"), mapping = aes(SIZE, prd, color = as.factor(3), group = YEAR), linewidth = 1)+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "SF-weighted"), name = "")+
  theme_bw()+
  facet_wrap(~YEAR)+
  ylab("Probability of terminal molt")+
  xlab("Carapace width (mm)")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14))

#ggsave("./Figures/binomialGAMs_pmolt_predictions.png", width = 11, height = 8.5)

## Solve for 50% maturity ----
# Specify function
find_l50 <- function(year_val, model) {
  f <- function(size_val) {
    predict(model, newdata = data.frame(SIZE = size_val, YEAR = year_val), type = "response") - 0.5
  }
  l50 <- uniroot(f, lower = min(mod.dat$SIZE), upper = max(mod.dat$SIZE))$root
  
  df <- data.frame(year = year_val, L50 = l50)
  
  return(df)
}

# Run functions
unique(mod.dat$YEAR) %>%
  purrr::map_df(~find_l50(., mod.1)) %>%
  mutate(type = "Unweighted, binomial GAM") -> uw.out


unique(mod.dat$YEAR) %>%
  purrr::map_df(~find_l50(., mod.2)) %>%
  mutate(type = "CPUE-weighted, binomial GAM") -> w.out

unique(mod.dat$YEAR) %>%
  purrr::map_df(~find_l50(., mod.3)) %>%
  mutate(type = "SF-weighted, binomial GAM") -> sf.out

# get crabpack data for original size at mat ts using sigmoid model
mat <- readRDS("./Data/snow_survey_maturityEBS.rda")$model_parameters %>%
      dplyr::select(YEAR, B_EST, B_SE) %>%
      rename(year = YEAR, L50 = B_EST, SE = B_SE) %>%
  mutate(type = "SF-weighted, legacy")



# Combine and plot
l50.df <- rbind(uw.out, w.out, sf.out, mat) %>%
  mutate(year = as.numeric(as.character(year))) %>%
  filter(!year %in% c(2010, 2012))

all.years <- data.frame(year = seq(min(l50.df$year), max(l50.df$year), by = 1))

missing.yrs <- right_join(l50.df, all.years) %>%
  filter(is.na(L50) == TRUE) %>%
  pull(year)

missing.yrs <- c(2008, 2010, 2012, 2013, 2016, 2020)


l50.df <- rbind(l50.df, data.frame(year = missing.yrs, L50 = rep(NA, 4*length(missing.yrs)), 
                                   type = rep(c("CPUE-weighted, binomial GAM", 
                                                "Unweighted, binomial GAM", "SF-weighted, binomial GAM", 
                                                "SF-weighted, legacy"), each = length(missing.yrs))))




ggplot(l50.df, aes(year, L50, group = type, color = type))+
  geom_point()+
  geom_line()+
  xlab("Year")+
  scale_x_continuous(breaks = seq(min(l50.df$year), max(l50.df$year), by = 2))+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon", "darkgreen"), labels = c("Unweighted, binomial GAM", 
                                                                                        "CPUE-weighted, binomial GAM",
                                                                                        "SF-weighted, binomial GAM",
                                                                                        "SF-weighted, legacy"), name = "")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave("./Figures/size_at_50maturity_comparison.png", width = 11, height = 8.5)
ll <- data.frame()

sizemat <- function(model, type, year){
  # Example code
    new_data <- data.frame(YEAR = year,  SIZE = seq(min(mod.dat$SIZE), max(mod.dat$SIZE), length.out = 100))
    
    preds <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
    
    log_odds <- preds$fit
    prob <- plogis(log_odds)  # Convert log-odds to probability
    l50 <- new_data$SIZE[which.min(abs(prob - 0.5))]
    
    ll <-  data.frame(l50 = l50, yy = year,
                               tt = type)
  
  
  return(ll)
}

unique(mod.dat$YEAR) %>%
purrr::map_df(~sizemat(mod.1, "Unweighted, binomial GAM", .)) -> out.1

unique(mod.dat$YEAR) %>%
  purrr::map_df(~sizemat(mod.2, "CPUE-weighted, binomial GAM", .)) -> out.2

unique(mod.dat$YEAR) %>%
  purrr::map_df(~sizemat(mod.3, "SF-weighted, binomial GAM", .)) -> out.3


rbind(out.1, out.2, out.3) %>%
  rename(L50 = l50, year = yy, type = tt) %>%
  mutate(SE = NA) %>%
  rbind(mat) %>%
  filter(year != 2012) %>%
  mutate(year = as.numeric(as.character(year))) -> bin.gam

all.years <- data.frame(year = seq(min(bin.gam$year), max(bin.gam$year), by = 1))

missing.yrs <- right_join(bin.gam, all.years) %>%
  filter(is.na(L50) == TRUE) %>%
  pull(year)


bin.gam2<- rbind(bin.gam, data.frame(year = missing.yrs, L50 = rep(NA, 4*length(missing.yrs)), 
                                     SE = rep(NA, 4*length(missing.yrs)),
                                   type = rep(c("CPUE-weighted, binomial GAM", 
                                                "Unweighted, binomial GAM", "SF-weighted, binomial GAM", 
                                                "SF-weighted, legacy"), each = length(missing.yrs))))




ggplot(bin.gam2, aes(year, L50, group = type, color = type))+
  geom_point()+
  geom_line()+
  geom_errorbar(bin.gam2, mapping = aes(ymin = L50-SE, ymax = L50+SE))+
  xlab("Year")+
  scale_x_continuous(breaks = seq(min(l50.df$year), max(l50.df$year), by = 3))+
  # scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon", "darkgreen"),
  #                    labels = c("Unweighted, binomial GAM", 
  #                               "CPUE-weighted, binomial GAM",
  #                               "SF-weighted, binomial GAM",
  #                               "SF-weighted, legacy"), name = "")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggsave("./Figures/size_at_50maturity_comparison.png", width = 11, height = 8.5)

