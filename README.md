# FluVax_prior_infection_study
# Data analysis code for each study figure
# Figure 1a Demographic beeswarm plot
library(ggbeeswarm)

# read and format data
LS_wide <- read.csv("HI_timecourse.csv",header = T, stringsAsFactors = F)

LS_wide$prior_H3_2 [LS_wide$prior_H3==0] <- "No"

LS_wide$prior_H3_2 [LS_wide$prior_H3==1] <- "Yes"

LS_wide$Sex <- LS_wide$SexS

# plot individual values with median line
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
summarise_logmean <- function(arr) 
{
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

# plot labels, log2 to absolute
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
--
# Figure 2b-c
library(mgcv)

library(tidyverse)

library(ggpubr)

# read and format data
LS_long <- read.csv("HI_long.csv",header = T, stringsAsFactors = F)
# group by birth decade - 4 different birth decades. First extract year from DoB and create new column YoB
LS_long$DoBS <- as.Date(LS_long$DoBS)

LS_long$YoB <- substring(LS_long$DoBS,1,4)

LS_long$YoB2 [LS_long$YoB < 1960] <- "1935-59 (n=24)"

LS_long$YoB2 [LS_long$YoB>=1960 & LS_long$YoB < 1970] <- "1960-69 (n=36)"

LS_long$YoB2 [LS_long$YoB>=1970 & LS_long$YoB < 1980] <- "1970-79 (n=27)"

LS_long$YoB2 [LS_long$YoB>=1980] <- "1980-96 (n=13)"

# convert characters to factors
LS_long$YoB2 <- factor(LS_long$YoB2)

LS_long$Subject_ID <- factor(LS_long$Subject_ID)

# Exclude Townsville 99 from all analysis
LS_long <- subset(LS_long, !Short_Name %in% c("Townsville/2/99"))
                                              
# plot labels, log2 to absolute
# y axis labels
yticks <- seq(2.32, 13.32, 1)

ylabs <- c("nd",10,20,40,80,160,320,640,1280,2560,5120,10240)
# x ticks clade and year all viruses reduced
LS_long$YearClCode2[LS_long$YearClCode2==2018] <- 2018.5

LS_long$YearClCode2[LS_long$YearClCode2==2017.5] <- 2017.75

xticks <- c(1968,1972,1975,1977,1979.5,1981.5,1987,1989.5,1992.5,1995.5,1996.5,1999,2002,2002.5,2004,2005,2007,2008,2009,2010,2011.5,2012,2013,2013.25, 2014,2014.25,2016,2017,2017.75,2018.5)

xlabels <- c("1968","1972","1975","1977","1979","1982","1987","1989", "1993","1995","1997","1999","2002","","2004","2005","2007","2008","2009","2010","2011","","2013","", "2014", "", "2016", "2017", "", "2018")

# run gam for fig 2b, pre-vaccine landscape by age group
model1.bd = gam(L2titre ~ s(YearClCode, by = YoB2) + YoB2 + s(Subject_ID, bs="re"), 
                data = subset(LS_long, time==1), method = "REML")
                
# plot gam for fig 2b
fig2b <- ggplot(data = LS_long, aes(x = YearClCode, y = L2titre)) + 
  geom_jitter(data = subset(LS_long, time==1), width=0.3, height=0.4, alpha=0.2, aes(colour=YoB2), size=0.5) +
  stat_smooth(data = subset(LS_long, time==1), method="gam", formula=formula(model1.bd)) +
  geom_smooth(data = subset(LS_long, time==1), aes(colour=YoB2, fill=YoB2)) +
  scale_color_manual(values = c("#547BD3","#C77CFF","#CC6677","#EBA85F")) +
  scale_fill_manual(values = c("#547BD3","#C77CFF","#CC6677","#EBA85F")) +
  geom_abline(intercept=5.32, slope=0, linetype="dotted", col="gray30", lwd=1.15) +
  geom_vline(xintercept = 2014, linetype="dotted", col="gray30", lwd=1.15) +
  xlab("A(H3N2) virus isolation year") + 
  scale_x_continuous(breaks=xticks, labels = xlabels) +
  ylab(expression("HI titer")) + 
  coord_cartesian(ylim=c(2.32, 12.5)) + 
  scale_y_continuous(breaks=yticks,labels=ylabs) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position="top")
--
 fig2b 
 
 ggsave("gam titre by time Fig 2c.pdf", fig2c, unit = "cm", width = 18, height = 12 )
-
# run gams for fig 2c, landscape by vaccination time-point, all vaccinees
model1.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==1), method = "REML")
                    
model2.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==2), method = "REML")
                    
model3.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==3), method = "REML")
                    
model4.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==4), method = "REML")
                    
model5.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==5), method = "REML")
                    
model6.cohort = gam(L2titre ~ s(YearClCode) + s(Subject_ID, bs="re"), 
                    data = subset(LS_long, time==6), method = "REML")
        
# plot gams for fig 2c
fig2c <- ggplot(data = LS_long, aes(x = YearClCode, y = L2titre)) + 
  geom_jitter(data = subset(LS_long, time==1), alpha=0.2, width=0.4, height=0.4, size=0.5, colour = "#6e6e6e") + 
  geom_jitter(data = subset(LS_long, time==3), alpha=0.2, width=0.4, height=0.4, size=0.5, colour = "#EBA85F") +
  geom_jitter(data = subset(LS_long, time==4), alpha=0.2, width=0.4, height=0.4, size=0.5, colour = "#CC6677") + 
  geom_jitter(data = subset(LS_long, time==5), alpha=0.2, width=0.4, height=0.4, size=0.5, colour = "#547BD3") +
  geom_jitter(data = subset(LS_long, time==6), alpha=0.2, width=0.4, height=0.4, size=0.5, colour = "#C77CFF") + 
  stat_smooth(data = subset(LS_long, time==1), method="gam", formula=formula(model1.cohort)) +
  geom_smooth(data = subset(LS_long, time==1), colour = "#6e6e6e", fill = "#6e6e6e", alpha=0.3) + 
  stat_smooth(data = subset(LS_long, time==3), method="gam", formula=formula(model3.cohort)) +
  geom_smooth(data = subset(LS_long, time==3), colour = "#EBA85F", fill = "#EBA85F",alpha=0.65) +
  stat_smooth(data = subset(LS_long, time==4), method="gam", formula=formula(model4.cohort)) +
  geom_smooth(data = subset(LS_long, time==4), colour = "#CC6677", fill = "#CC6677", alpha=0.5) + 
  stat_smooth(data = subset(LS_long, time==5), method="gam", formula=formula(model5.cohort)) +
  geom_smooth(data = subset(LS_long, time==5), colour = "#547BD3", fill = "#547BD3",alpha=0.6) +
  stat_smooth(data = subset(LS_long, time==6), method="gam", formula=formula(model6.cohort)) +
  geom_smooth(data = subset(LS_long, time==6), colour = "#C77CFF", fill = "#C77CFF") + 
  geom_abline(intercept=5.32, slope=0, linetype="dotted", col="gray30", lwd=1.15) +
  geom_vline(xintercept = 2014, linetype="dotted", col="gray30", lwd=1.15) +
  xlab("A(H3N2) virus isolation year") + 
  scale_x_continuous(breaks=xticks, labels = xlabels) +
  ylab(expression("HI titer")) + 
  coord_cartesian(ylim=c(2.32, 12.5)) + 
  scale_y_continuous(breaks=yticks,labels=ylabs) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, margin = margin(2,0,0,0)),
        axis.title.x=element_text(size = 8, margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position="top") 

fig2c

