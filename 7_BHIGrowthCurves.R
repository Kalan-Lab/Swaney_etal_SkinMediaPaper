
# Script 7 ----------------------------------------------------------------
# Supplemental Figure 1

# Skin Media Manuscript
# BHI growth curves
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(tidyverse)
library(plater)
library(ggplot2)

# set working directory to data folder
setwd("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/BHI_GrowthCurves/")

# Strain Set 1

set1 <- read_excel("2022_9_2_SkinMediaStrains_BHI_1.xlsx") # read in growth curve data
gathered1 <- set1 %>% gather(key = Wells, value = OD600, 2:ncol(set1)) # convert to long format

layout1 <- read_plate("2022_9_2_SkinMediaStrains_BHI_1_setup.csv") # read in palte layout

gathered1 <- left_join(gathered1,layout1)
gathered1 <- gathered1 %>% filter(!is.na(Strain))
gathered1$Replicate <- 1

# Strain Set 2
set2 <- read_excel("2022_9_5_SkinMediaStrains_BHI_2.xlsx")
gathered2 <- set2 %>% gather(key = Wells, value = OD600, 2:ncol(set2))

layout2 <- read_plate("2022_9_5_SkinMediaStrains_BHI_2_setup.csv")

gathered2 <- left_join(gathered2,layout2)
gathered2 <- gathered2 %>% filter(!is.na(Strain))
gathered2$Replicate <- 2

gathered_all <- rbind(gathered1, gathered2) # Join data from strain sets 1 and 2

# Strain Set 3

set3 <- read_excel("2022_9_7_SkinMediaStrains_BHI_3.xlsx")
gathered3 <- set3 %>% gather(key = Wells, value = OD600, 2:ncol(set2))

layout3 <- read_plate("2022_9_7_SkinMediaStrains_BHI_3_setup.csv")

gathered3 <- left_join(gathered3,layout3)
gathered3 <- gathered3 %>% filter(!is.na(Strain))
gathered3$Replicate <- 3

gathered_all <- rbind(gathered_all, gathered3) # Join data from strain sets 1-3

# Strain Set 4

set4 <- read_excel("2022_9_9_SkinMediaStrains_BHI_4.xlsx")
gathered4 <- set4 %>% gather(key = Wells, value = OD600, 2:ncol(set2))

layout4 <- read_plate("2022_9_9_SkinMediaStrains_BHI_4_setup.csv")

gathered4 <- left_join(gathered4,layout4)
gathered4 <- gathered4 %>% filter(!is.na(Strain))
gathered4$Replicate <- 4

gathered_all <- rbind(gathered_all, gathered4) # Join data from strain sets 1-4

# Average all replicates and plot  --------------------------------------------

# re-name LK369
gathered_all$Strain <- ifelse(gathered_all$Strain == "Microbacterium sp. 002456035 LK369","Microbacterium sp. LK369", gathered_all$Strain)

gathered_all$Time <- round(gathered_all$Time, digits = 2) # round time to 2 decimal places

gathered_all <- gathered_all %>% group_by(Strain, Time, Replicate) %>% dplyr::summarise(mean=mean(OD600)) # average technical replicates

# Pull out media control data
blank <- gathered_all %>% filter(Strain == "Blank")
colnames(blank) <- c("Strain","Time","Replicate","blank_mean")

# Join together strain data with their respective media control OD600s
samples <- gathered_all %>% filter(Strain !="Blank")
samples <- left_join(samples,blank[,c(2:4)])

# Subtract the media control OD600 from the sample OD600
samples$subtracted_OD600 <- samples$mean-samples$blank_mean

# Take the average and standard deviation of the biological replicates
samples_avg_BHI <- samples %>% group_by(Strain, Time) %>% dplyr::summarise(avg_OD600 = mean(subtracted_OD600),
                                                                                      avg_SD = sd(subtracted_OD600))

samples_avg_BHI$Strain <- factor(samples_avg_BHI$Strain, levels = c("Corynebacterium amycolatum LK19","Kocuria rhizophila LK221","Dermacoccus nishinomiyaensis LK1128",
                                                  "Staphylococcus epidermidis LK593","Corynebacterium glucuronolyticum LK488","Kocuria marina_A LK478",
                                                  "Dermabacter hominis LK522","Staphylococcus epidermidis LK717","Corynebacterium kefirresidentii LK1134",
                                                  "Microbacterium sp. LK369","Citrobacter freundii LK704","Staphylococcus aureus USA300",
                                                  "Dietzia cinnamea LK439","Micrococcus luteus LK410","Klebsiella pneumoniae LK469","Sphingobacterium hotanense LK485"))

# Final plot - Supplemental Figure 1A
BHI_plot <- samples_avg_BHI %>% filter(Strain != "Blank") %>%
  ggplot(aes(x=Time,y=avg_OD600)) + 
  geom_ribbon(aes(x=Time,ymin=avg_OD600-avg_SD, ymax=avg_OD600+avg_SD) ,alpha=0.15) + 
  geom_smooth(method = 'loess',color="black", se = FALSE, size=1.25) + 
  facet_wrap(~Strain) + 
  theme_bw() + ylab("OD600") +
  xlab("Time (hours)")

#ggsave(plot=BHI_plot, filename="/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/BHI_Curves.pdf",
      #width = 12, height=10, units = "in", device = "pdf")

# # plot individual averaged points
# samples_avg_BHI %>% filter(Strain != "Blank") %>%
#   ggplot(aes(x=Time,y=avg_OD600)) + 
#   geom_point() + 
#   facet_wrap(~Strain) + 
#   theme_bw() + ylab("OD600") +
#   xlab("Time (hours)")

