# FluVax_prior_infection_study
# Data analysis code for each study figure
# Figure 1a Demographic beeswarm plot
LS_wide <- read.csv("HI_timecourse.csv",header = T, stringsAsFactors = F)
LS_wide$prior_H3_2 [LS_wide$prior_H3==0] <- "No"
LS_wide$prior_H3_2 [LS_wide$prior_H3==1] <- "Yes"

LS_wide$Sex <- LS_wide$SexS

p <- ggplot(data = subset(LS_wide, virus==5), aes(prior_H3_2, Age, color = Sex)) +
  geom_quasirandom(width=0.3, size=1) +
  scale_color_manual(values = c("#d6604d", "#4393c3")) +
  scale_y_continuous(breaks = seq(10, 90,10)) +
  geom_abline(intercept=49, slope=0, linetype="dashed", col="gray30", lwd=0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6, margin = margin(2,0,0,0)),
        axis.title.y = element_text(("Age, years"),size = 7, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=5, margin = margin(0,0,0,0)),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=7),
        legend.position = "top")
p
--
# Figure 2a Antibody titres against vaccine A(H3N2) strain (HK14e) by prior infection ande seroconversion status
library(tidyverse)

# read and format data

data <- read.csv("HI_long_diff.csv",header = T, stringsAsFactors = F)

data <- subset(data, virus  %in% c(5))

data <- subset(data, time  <= 6)

data$Conv[data$t1_otherDiff >=2] <- 1

data$Conv[data$t1_otherDiff <2 |data$t1_otherDiff == "NA"] <- 0

data_renamed <- data %>%
  select(
    pid = Subject_ID, timepoint = time, sex = Sex, age = Age,
    prior_H3 = prior_H3,conv = Conv, l2hi = L2titre, titre = Titer
      )  

data_extra <- data_renamed %>%
  mutate(
    exposure_group = case_when(
      prior_H3 == 0  ~ "no",
      prior_H3 == 1 ~ "yes",
    ) %>%
      factor(c(
        "no", "yes"
      )),
    timepoint_lbl = factor(
      timepoint, 1:6, c("pre", "d4", "d7", "d14", "d21","d280")
    ),
    conv_lbl = factor(
     conv, 0:1, c("no", "yes")
    ),
  ) 

# calculate gmts
summarise_logmean <- function(arr) {
  logarr <- log(arr)
  logmean <- mean(logarr)
  logse <- sd(logarr) / sqrt(length(arr))
  logerr_margin <- qnorm(0.975) * logse
  tibble(
    mean = exp(logmean),
    low = exp(logmean - logerr_margin),
    high = exp(logmean + logerr_margin)
  )
}

gmts <- data_extra %>%
  group_by(exposure_group,timepoint_lbl) %>%
  summarise(summarise_logmean(titre), .groups = "drop")

# plot data
yticks <- seq(2.32, 14.32, 1)
ylabs <- c(5,10,20,40,80,160,320,640,1280,2560,5120,10240,20480)  

# plot each value
p <- ggplot(data_extra, aes(exposure_group, l2hi,shape = conv_lbl, color = exposure_group )) +
  facet_grid(~ timepoint_lbl) +
  geom_jitter(width=0.3, height=0.3, size=1) +
  scale_color_manual(values = c("#808080","#0099FF"))+
  scale_shape_manual(values = c(1,3)) +
  scale_y_continuous( "HI titre", breaks = seq(2.32,14.32,1),labels=ylabs) +
  geom_abline(intercept=5.32, slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position=c(0.2,1))
p

# plot gmts to be overlaid
gmt_plot <- ggplot(gmts,aes(exposure_group,  mean, ymin = low, ymax = high, color = exposure_group)) +
  facet_grid(~ timepoint_lbl) +
  geom_crossbar()+
  geom_linerange()+
  scale_color_manual(values = c("#333333","#0033FF"))+
  scale_y_continuous("GMT",trans = "log2", breaks = 5 * 2^(0:15)) +
  coord_cartesian(ylim=c(5, 20000)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        legend.text = element_text(size=8),
        legend.position=c(0.12,0.9))
gmt_plot
