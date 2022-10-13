# Script 3 ----------------------------------------------------------------
# Figure 3 and Supplemental Figure 2B

# Skin Media Manuscript
# Sebum growth curves
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(tidyverse)
library(data.table)
library(growthrates)

# set working directory to data folder
setwd("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/")

map <- read_excel("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/File_Mapping_GrowthCurves.xlsx", sheet = "Sebum")

# Read in strain information
strains <- read_excel("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/File_Mapping_GrowthCurves.xlsx", sheet = "Strains")
strains[11,1] <- "Microbacterium sp. LK369" # re-name

# color palette for genera
cols <- c("#103a66","#4e80b5","#a1c2e6","#107c91", "#10adb0","#94d2bd", "#cae3db", #Actinobacteria
          "#c2a242", "#ebd8a0", #Proteobacteria
          "#b35539", #Firmicutes
          "#821b5b") #Bacteroidetes

# color palette for sebum concentration
rainbow <- c("#C92C31","#F16704","#FFB01F","#DECB36","#BCE64C","#63BD71","#0A9396","#005F73")

# Functions ---------------------------------------------------------------

# Function to combine growth curve data from individual experiments
ReadSebum <- function(mapping_file, type){
  
  File <- mapping_file[1] # file location 
  Wells1_3 <- mapping_file[2] # strain 1
  Wells4_6 <- mapping_file[3] # strain 2
  Wells7_9 <- mapping_file[4] # strain 3
  Date <- mapping_file[5] # experiment date
  Rep <- mapping_file[6] # if any strains had multiple biological replicates on the plate
  
  data <- read_excel(File)
  data.new <- data[seq(1, nrow(data), 4), ] # get every 4th time point
  
  gathered <- data.new %>% gather(key = Well, value = OD600, 2:ncol(data.new))# gather data into long format
  
 # assign sebum concentrations
   gathered$Concentration <- ifelse(startsWith(gathered$Well,"A"), 0, 
                                   ifelse(startsWith(gathered$Well,"B"), 0.0075,
                                          ifelse(startsWith(gathered$Well,"C"), 0.01,
                                                 ifelse(startsWith(gathered$Well,"D"), 0.025,
                                                        ifelse(startsWith(gathered$Well,"E"), 0.05,
                                                               ifelse(startsWith(gathered$Well,"F"), 0.075,
                                                                      ifelse(startsWith(gathered$Well,"G"), 0.1,
                                                                             ifelse(startsWith(gathered$Well,"H"), 0.25, NA))))))))
   # assign strain information
  gathered$Strain <- ifelse(endsWith(gathered$Well, "10") | endsWith(gathered$Well, "11") | endsWith(gathered$Well, "12"), "Blank",
                            ifelse(endsWith(gathered$Well, "4") | endsWith(gathered$Well, "5") | endsWith(gathered$Well, "6"), Wells4_6,
                                   ifelse(endsWith(gathered$Well, "7") | endsWith(gathered$Well, "8") | endsWith(gathered$Well, "9"), Wells7_9,
                                          ifelse(endsWith(gathered$Well, "1") | endsWith(gathered$Well, "2") | endsWith(gathered$Well, "3"), Wells1_3, NA))))
  
  gathered$Type <- "Sebum"
  gathered$Date <- Date
  gathered$Time <- round(gathered$Time,2)
  
  gathered$Replicate <- Rep
 
  # assign replicate values depending on the number of biological replicates in a plate
  rep1_2 <- c(rep(c(1),165), rep(c(2),165), rep(c(1),330))
  rep1_2 <- c(rep(rep1_2,8))
  rep2_3 <- c(rep(c(1),165), rep(c(1),165), rep(c(2),165),rep(c(1),165))
  rep2_3 <- c(rep(rep2_3,8))
  rep1_2_3 <- c(rep(c(1),165), rep(c(2),165), rep(c(3),165),rep(c(1),165))
  rep1_2_3 <- c(rep(rep1_2_3,8))
  
  gathered$Replicate <- ifelse(gathered$Replicate == "1_2", rep1_2, 
                               ifelse(gathered$Replicate == "2_3",rep2_3,
                                      ifelse(gathered$Replicate == "1_2_3", rep1_2_3, 1)))
  
  # inoculated an additional column for C. amycolatum
  gathered$Strain <- ifelse(gathered$Strain == "Blank" & gathered$Date == "2022_2_22" & (gathered$Well == "A10" | gathered$Well == "B10" | gathered$Well == "C10" | 
                                                                                           gathered$Well == "D10" | gathered$Well == "E10" | gathered$Well =="F10" | 
                                                                                           gathered$Well == "G10" | gathered$Well == "H10"), 
                            "Corynebacterium amycolatum LK19", 
                            gathered$Strain)
  
  
  #exclude well - contamination
  gathered <- gathered %>% filter(!(Well == "D7" & Strain == "Staphylococcus epidermidis LK593" & Date == "2022_4_20"))
  
  
  #exclude column - used to inoculate a strain
  gathered <- gathered %>% filter(!(Strain == "Blank" & Date == "2022_4_11" & (Well == "A10" | Well == "B10" | Well == "C10" | 
                                                                                 Well == "D10" | Well == "E10" | Well =="F10" | 
                                                                                 Well == "G10" | Well == "H10")))
  return(gathered)       
  
}

# Growth Curves ----------------------------------------------------------------

# Read in sebum growth curves and bind data together
DataList <- apply(map, 1, FUN = ReadSebum)
DataAllStrains <- rbindlist(DataList)

# Calculations
# take the mean of the OD600 for technical replicates - both strains and blanks
# for each individual experiment, subtract the mean OD600 of the strain minus the mean OD600 of the blank, for each sweat/sebum concentration
# then, take the mean and standard deviation of all biological replicates (across multiple dates) for each sweat/sebum concentration and time point

# Take average of technical replicates
grouped <- DataAllStrains %>% group_by(Strain, Concentration, Time, Date, Replicate) %>% dplyr::summarise(mean=mean(OD600))

# Pull out media control data
blank <- grouped %>% filter(Strain == "Blank")
colnames(blank) <- c("Strain","Concentration","Time","Date","Replicate","blank_mean")

# Join together strain data with their respective media control OD600s
samples <- grouped %>% filter(Strain !="Blank")
samples <- left_join(samples,blank[,c(2:4,6)])

# Subtract the media control OD600 from the sample OD600
samples$subtracted_OD600 <- samples$mean-samples$blank_mean

# Take the average and standard deviation of the biological replicates
samples_avg <- samples %>% group_by(Strain, Concentration, Time) %>% dplyr::summarise(avg_OD600 = mean(subtracted_OD600),
                                                                               avg_SD = sd(subtracted_OD600))
samples_avg <- samples_avg %>% filter(Strain != "NA")

samples_avg$Strain <- factor(samples_avg$Strain, levels = c("Corynebacterium amycolatum LK19","Kocuria rhizophila LK221","Dermacoccus nishinomiyaensis LK1128",
                                                            "Staphylococcus epidermidis LK593","Corynebacterium glucuronolyticum LK488","Kocuria marina_A LK478",
                                                            "Dermabacter hominis LK522","Staphylococcus epidermidis LK717","Corynebacterium kefirresidentii LK1134",
                                                            "Microbacterium sp. LK369","Citrobacter freundii LK704","Staphylococcus aureus USA300",
                                                            "Dietzia cinnamea LK439","Micrococcus luteus LK410","Klebsiella pneumoniae LK469","Sphingobacterium hotanense LK485"))

# Plot sebum growth curves - Figure 3A
sebum_big <- samples_avg %>%
  ggplot(aes(x=Time,y=avg_OD600,color=as.factor(Concentration))) +  
  theme_bw() +
  facet_wrap(~Strain) + 
  geom_ribbon(aes(x=Time,ymin=avg_OD600-avg_SD, ymax=avg_OD600+avg_SD, fill=as.factor(Concentration)),alpha=0.08, color=NA) + 
  geom_smooth(se = FALSE,aes(linetype=as.factor(Concentration)), size=1.25) +
  scale_color_manual(values=rainbow) + 
  labs(color="Sebum\nconcentration (%)",linetype="Sebum\nconcentration (%)") + 
  xlab("Time (hours)") + ylab("OD600") + scale_linetype_manual(values=c("0"=1,"0.0075"=2,"0.01"=2,"0.025"=2,"0.05"=1,"0.075"=2,"0.1"=1,"0.25"=2)) + 
  scale_fill_manual(values=rainbow, guide="none")+
  ylim(-0.1,1.3)

#ggsave(plot=sebum_big, filename="/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/SebumCurves.pdf",
       #width = 12, height=10, units = "in", device = "pdf")

# Area Under the Curve ----------------------------------------------------

conc <- c("0","0.0075","0.01","0.025","0.05","0.075","0.1","0.25")
strain <- unique(DataAllStrains$Strain)
strain <- strain[which(strain != "NA")]
strain <- strain[which(strain != "Blank")]

# Calculate area under the curve for each sebum concentration for each strain
AUC <- lapply(strain, FUN=function(X){
  strain_filter <- samples_avg %>% filter(Strain == X)
  areas <- lapply(conc, FUN=function(Y){
    conc_filter <- strain_filter %>% filter(Concentration == Y)
    fit <- loess(conc_filter$avg_OD600~conc_filter$Time) # Fit LOESS model
    df <- data.frame(x = conc_filter$Time)
    df <- transform(df, y.pred = predict(fit, df))
    auc <- pracma::trapz(df$x, df$y.pred)
  })
  names(areas) <- conc
  return(areas)
})


# combine all data for each strain and sebum concentration
aucauc <- rbindlist(AUC)
aucauc$Strain <- strain
auc_df <- melt(aucauc)
colnames(auc_df) <- c("Strain","Concentration","AUC")
auc_df$Strain <- factor(auc_df$Strain, levels = c("Corynebacterium amycolatum LK19","Kocuria rhizophila LK221","Dermacoccus nishinomiyaensis LK1128",
                                                            "Staphylococcus epidermidis LK593","Corynebacterium glucuronolyticum LK488","Kocuria marina_A LK478",
                                                            "Dermabacter hominis LK522","Staphylococcus epidermidis LK717","Corynebacterium kefirresidentii LK1134",
                                                            "Microbacterium sp. LK369","Citrobacter freundii LK704","Staphylococcus aureus USA300",
                                                            "Dietzia cinnamea LK439","Micrococcus luteus LK410","Klebsiella pneumoniae LK469","Sphingobacterium hotanense LK485"))

# Plot sebum AUC - Figure 3B
auc_plot <- ggplot(auc_df, aes(x=Concentration, y=AUC, fill=Concentration)) + 
  theme_bw() + 
  geom_col(color="black") + 
  facet_wrap(~Strain) + scale_fill_manual(values = rainbow) + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        text=element_text(size=5)) + xlab("Sebum concentration (%)") + 
  ylim(0,32)

#ggsave(plot=auc_plot, filename="/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/SebumAUC.pdf",
       #width = 5.45, height=5, units = "in", device = "pdf")

# Linear Regression -------------------------------------------------------

# normalize AUC data out of 1
auc_df_min_max <- auc_df %>% group_by(Strain) %>% summarise(minAUC = min(AUC), maxAUC = max(AUC))
auc_df <- left_join(auc_df, auc_df_min_max)
auc_df$NormAUC <- (auc_df$AUC - auc_df$minAUC) / (auc_df$maxAUC - auc_df$minAUC)

# linear regression with AUC data for each strain
fitted_models <- auc_df %>% 
  group_by(Strain) %>% 
  do(model=lm(NormAUC~as.numeric(Concentration), data = .))

# extract p-values from fitted models
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# extract coefficients from fitted models
coef <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$coefficients
  c <- f[2]
  return(c)
}

# create a data frame with extracted coefficients and p-values for each strain
auc_lm <- apply(fitted_models, 1, FUN=function(X) {
  strain <- X[[1]]
  coef <- coef(X[[2]])
  pval <- lmp(X[[2]])
  return(data.frame(Strain = strain, Coef = coef, Pval= pval))
})

auc_lm <- rbindlist(auc_lm)
auc_lm <- left_join(auc_lm, strains)
auc_lm$Genus <- factor(auc_lm$Genus, levels = c("Corynebacterium","Dermabacter","Dermacoccus","Dietzia","Kocuria","Microbacterium",
                                                "Micrococcus","Citrobacter","Klebsiella","Staphylococcus","Sphingobacterium"))

# plot for Figure 3C
auc_point <- ggplot(auc_lm, aes(x=Coef, y=-log10(Pval), color = Genus)) + 
  theme_classic() + 
  geom_point(size=7) +
  scale_color_manual(values=cols) +
  xlim(-0.3, 0.3) + 
  ylab("-log10(p-value)") +
  xlab("Slope coefficient") + 
  geom_hline(yintercept = -log10(0.05), color = "grey") + 
  geom_vline(xintercept = 0, color = "grey") +
  annotate("text", x=-0.16, y=4.5, label = "Low sebum") +
  annotate("text", x=0.16, y=4.5, label = "High sebum") + 
  annotate("text",x=-0.2, y=0.5, label = "No preference") + 
  theme(legend.position = "bottom")

# classify sebum preferences
auc_lm$SebumPref <- ifelse(auc_lm$Pval < 0.05 & auc_lm$Coef < 0, "Low sebum",
                           ifelse(auc_lm$Coef > 0 & auc_lm$Pval < 0.05, "High sebum", "No preference"))
auc_lm$SebumPref <- factor(auc_lm$SebumPref, levels = c("No preference","High sebum","Low sebum"))

# preference plot for Figure 3C
pref_auc <- ggplot(auc_lm, aes(x=SebumPref, fill=Genus)) + geom_bar(color="black") + 
  theme_bw() +
  scale_fill_manual(values=cols) +
  scale_x_discrete(drop=FALSE) + 
  xlab("Sebum preference predicted by AUC") + 
  ylab("Number of strains") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major.x = element_blank() , 
        panel.grid.major.y = element_line(size=0.1, color="black",) ,
        panel.grid.minor.y = element_blank()) + 
  scale_y_continuous(breaks = seq(0,14, by=2)) + 
  theme(legend.position = "none")

legend <- as_ggplot(get_legend(auc_point))

#ggsave(legend, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/Legend.pdf",
      # width = 6, height = 1, units = "in", device = "pdf")

#ggsave(pref_auc, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/SebumPrefLM.pdf",
      # width = 3, height=4, units = "in", device = "pdf")

#ggsave(auc_point, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/SebumLMPoint.pdf",
      # width = 3.5, height=4, units = "in", device = "pdf")

# Maximum Growth Rate -----------------------------------------------------

# calculate max growth rate
L <- all_easylinear(avg_OD600 ~ Time | Strain + Concentration, data=samples_avg)
maxgrowth <- results(L)

maxgrowth$Concentration <- factor(maxgrowth$Concentration, levels = c("0","0.0075","0.01","0.025","0.05","0.075","0.1","0.25"))
maxgrowth$Strain <- factor(maxgrowth$Strain, levels = c("Corynebacterium amycolatum LK19","Kocuria rhizophila LK221","Dermacoccus nishinomiyaensis LK1128",
                                                  "Staphylococcus epidermidis LK593","Corynebacterium glucuronolyticum LK488","Kocuria marina_A LK478",
                                                  "Dermabacter hominis LK522","Staphylococcus epidermidis LK717","Corynebacterium kefirresidentii LK1134",
                                                  "Microbacterium sp. LK369","Citrobacter freundii LK704","Staphylococcus aureus USA300",
                                                  "Dietzia cinnamea LK439","Micrococcus luteus LK410","Klebsiella pneumoniae LK469","Sphingobacterium hotanense LK485"))

# Supplemental Figure 3A
growthrate <- ggplot(maxgrowth, aes(x=as.factor(Concentration), y=mumax, color=as.factor(Concentration))) + geom_point(size=3) + 
  facet_wrap(~Strain) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text=element_text(size=5),
        legend.position = "none") +
  xlab("Sebum Concentration") + ylab("Maximum growth rate") +
  scale_color_manual(values=rainbow)

#ggsave(plot=growthrate, filename="/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure3/SebumGrowthRate.pdf",
       #width = 5.45, height=5, units = "in", device = "pdf")
