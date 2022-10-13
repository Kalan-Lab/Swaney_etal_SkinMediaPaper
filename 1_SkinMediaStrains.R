
# Script 1 ----------------------------------------------------------------
# Figure 1 

# Skin Media Manuscript
# Strain information
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(tidyverse)
library(ggplot2)

strains <- read_excel("/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/SupplementalMaterial.xlsx", sheet = "S2 Strain Information", skip = 1)

# color palette
cols <- c("#103a66","#4e80b5","#a1c2e6","#107c91", "#10adb0","#94d2bd", "#cae3db",
          "#c2a242", "#ebd8a0",
          "#b35539", 
          "#821b5b")

# Plot Figure 1B ----------------------------------------------------------

strains$Microenvironment <- factor(strains$Microenvironment, levels = c("Sebaceous","Moist","Dry","Foot"))

strains$Genus <- factor(strains$Genus, levels = c("Corynebacterium","Dermabacter","Dermacoccus","Dietzia","Kocuria","Microbacterium",
                                                  "Micrococcus","Citrobacter","Klebsiella","Staphylococcus","Sphingobacterium"))

strains %>% filter(Site != "NA") %>% 
  ggplot(aes(x=Microenvironment, fill=Genus)) + geom_bar(color="black") + theme_bw() + 
  scale_fill_manual(values=cols) + ylab("Number of strains") +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size=0.1, color="black" ) ,
    panel.grid.minor.y = element_line(size=0.1, color="black")
  )

