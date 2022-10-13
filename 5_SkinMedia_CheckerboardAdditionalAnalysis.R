
# Script 5 ----------------------------------------------------------------
# Figures 5B-D

# Skin Media Manuscript
# Checkerboard assay - Sweat/sebum preference, principal components analysis, upset plot
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(ggfortify)
library(ComplexUpset)

setwd("/Users/mhswaney/GitHub/scripts/mhswaney/Skin_media_project/SkinMediaManuscript/")

source("4_SkinMedia_Checkerboards.R") # load data for checkerboard analysis

# Read in strain information
strains <- read_excel("/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/SupplementalMaterial.xlsx", sheet = "S2 Strain Information", skip = 1)
strains[11,1] <- "Microbacterium sp. LK369" # re-name

# color palette for genera
cols <- c("#103a66","#4e80b5","#a1c2e6","#107c91", "#10adb0","#94d2bd", "#cae3db", #Actinobacteria
          "#c2a242", "#ebd8a0", #Proteobacteria
          "#b35539", #Firmicutes
          "#821b5b") #Bacteroidetes

#### Analysis ----------------------------------------------------------------

# Sweat and sebum preference based on Multiple Linear Regression -------------

mlr_plot <- left_join(mlr_plot,strains) # join mlr data with strain info

mlr_plot$Genus <- factor(mlr_plot$Genus, levels = c("Corynebacterium","Dermabacter","Dermacoccus","Dietzia","Kocuria","Microbacterium",
                                          "Micrococcus","Citrobacter","Klebsiella","Staphylococcus","Sphingobacterium"))
mlr_plot$Microenvironment <- factor(mlr_plot$Microenvironment, levels = c("Sebaceous","Moist","Dry","Foot","NA"))

# plot sweat and sebum preference based on MLR - genus (Figure 5C)
ggplot(mlr_plot, aes(x_start, y_start, xend = x_end, yend = y_end, color=Genus), size=)+geom_segment(arrow = arrow(length = unit(0.1,"inches")),size=1.5, lineend = "round", linejoin="round") + 
  theme_bw() +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0) +  scale_color_manual(values=cols) + 
  ylim(-3,3) + xlim(-3,3) + ylab("Sebum preference") + xlab("Sweat preference")

# plot sweat and sebum preference based on MLR - microenvironment (Supplemental Figure 4)
ggplot(mlr_plot, aes(x_start, y_start, xend = x_end, yend = y_end, color=Microenvironment), size=)+geom_segment(arrow = arrow(length = unit(0.1,"inches")),size=1.5, lineend = "round", linejoin="round") + 
  theme_bw() +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0) +  scale_color_manual(values=c("#C32446","#F2B342","#3E7274","#2B3B4D","#B5B5B5")) + 
  ylim(-3,3) + xlim(-3,3) + ylab("Sebum preference") + xlab("Sweat preference")


# Classify sweat and sebum preferences ------------------------------------

mlr_plot$SebumPref <- ifelse(mlr_plot$sebum_pval < 0.05 & mlr_plot$y_end < 0, "Low sebum",
                           ifelse(mlr_plot$sebum_pval < 0.05 & mlr_plot$y_end > 0, "High sebum", "No preference"))

mlr_plot$SweatPref <- ifelse(mlr_plot$sweat_pval < 0.05 & mlr_plot$x_end < 0, "Low sweat",
                             ifelse(mlr_plot$sweat_pval < 0.05 & mlr_plot$x_end > 0, "High sweat", "No preference"))


# Principal Component Analysis --------------------------------------------

# prepare data
t_allStrains <- as.data.frame(t(allStrains[,2:ncol(allStrains)]))
colnames(t_allStrains) <- allStrains$SweatSebum
t_allStrains$Strain <- sub("_", " ",rownames(t_allStrains))
t_allStrains$Strain <- sub("_", " ",t_allStrains$Strain)
t_allStrains[which(t_allStrains$Strain == "Kocuria marina A_LK478"),"Strain"] <- "Kocuria marina_A LK478"

# join data with strain info
t_allStrains <- left_join(t_allStrains, strains)
rownames(t_allStrains) <- t_allStrains$Strain

# perform PCA
pca <- prcomp(t_allStrains[,c(1:64)], scale. = FALSE)
summary(pca)

t_allStrains$Genus <- factor(t_allStrains$Genus, levels = c("Corynebacterium","Dermabacter","Dermacoccus","Dietzia","Kocuria","Microbacterium",
                                                            "Micrococcus","Citrobacter","Klebsiella","Staphylococcus","Sphingobacterium"))
t_allStrains$Microenvironment <- factor(t_allStrains$Microenvironment, levels=c("Sebaceous","Moist","Dry","Foot","NA"))

# plot PCA - Figure 5B
pca_plot <- ggplot2::autoplot(pca, data=t_allStrains, colour="Genus", size=6, shape = "Microenvironment") + 
  theme_bw() + scale_color_manual(values=cols) +
  scale_shape_manual(values=c("Sebaceous" = 19, "Moist" = 17, "Dry" = 15, "Foot" = 18,"NA" = 4)) 

#ggsave(plot=pca_plot, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure5/CheckerboardPCA.pdf",
 # device= "pdf",width=5.75, height=4, units="in")

# Hierarchical clustering - For Script 6 --------------------------------------

# compute Euclidean distance matrix
dist_mat <- dist(t_allStrains[,c(1:64)], method = "euclidean")

# perform hierarchical clustering - this hclust object is used in Script 6 for the tanglegram analysis
hclust_avg <- hclust(dist_mat, method="average")

#save(hclust_avg,file = "/Users/mhswaney/GitHub/scripts/mhswaney/Skin_media_project/SkinMediaManuscript/CheckerboardDendrogram.RData")

#plot(hclust_avg)

# Classify strain sweat/sebum preferences ---------------------------------

forupset <- mlr_plot

forupset$HighSweat <- 0
forupset$LowSweat <- 0
forupset$NoSwPref <- 0
forupset$HighSebum <- 0
forupset$LowSebum <- 0
forupset$NoSebPref <- 0

# assign binary for sweat/sebum preference groups
forupset$HighSweat <- ifelse(forupset$SweatPref == "High sweat", 1, forupset$HighSweat)
forupset$LowSweat <- ifelse(forupset$SweatPref == "Low sweat", 1, forupset$LowSweat)
forupset$NoSwPref <- ifelse(forupset$SweatPref == "No sweat preference", 1, forupset$NoSwPref)
forupset$HighSebum <- ifelse(forupset$SebumPref == "High sebum", 1, forupset$HighSebum)
forupset$LowSebum <- ifelse(forupset$SebumPref == "Low sebum", 1, forupset$LowSebum)
forupset$NoSebPref <- ifelse(forupset$SebumPref == "No sebum preference", 1, forupset$NoSebPref)

# join upset plot data with strain information
forupset <- forupset[,c(7,14:19)]
forupset <- left_join(forupset,strains)

forupset$Genus <- factor(forupset$Genus, levels = c("Corynebacterium","Dermabacter","Dermacoccus","Dietzia","Kocuria","Microbacterium",
                                                    "Micrococcus","Citrobacter","Klebsiella","Staphylococcus","Sphingobacterium"))

groups <- c("HighSweat","LowSweat","NoSwPref","HighSebum","LowSebum","NoSebPref")

# Upset plot - Figure 5D
upset(forupset,groups, 
      base_annotations=list(
        'Intersection size'=intersection_size(aes(fill=Genus)) + scale_fill_manual(values = cols)))

