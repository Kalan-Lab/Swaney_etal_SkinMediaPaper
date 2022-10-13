
# Script 4 ----------------------------------------------------------------
# Figure 4

# Skin Media Manuscript
# Checkerboard assay - Heat maps, multiple linear regression analysis, and prepare data frame
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(purrr)
library(tidyverse)
library(reshape2)
library(relaimpo)
library(pheatmap)

#### Create heat maps for each experiment -------------------------------------

# Code is the same for all strains, apart from reading in results from file
# Commented code is included for analyzing results for the first two strains

# Blanks ------------------------------------------------------------------

# read in experiments
Blank_1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_17_blank.txt",
                          skip=32, n_max = 8, delim = "\t")[,2:9])
Blank_2 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_24_blank_checkerboard.txt",
                          skip=32, n_max = 8, delim = "\t")[,2:9])
Blank_3 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_5_blank_checkerboard.txt",
                             skip=34, n_max = 8, delim = "\t")[,2:9])

# create list with replicates
list_reps <- list(Blank_1, Blank_2, Blank_3)

# take averages and round to second decimal place
blank_avg <- reduce(list_reps, `+`) / length(list_reps)
blank_avg <- round(as.matrix(blank_avg),3)

# Corynebacterium amycolatum LK19 -----------------------------------------

# read in experiments - 3 biological replicates and 3 OD600 technical replicates (because of cell clumping)
Rep1_1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                            skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep1_2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep1_3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 122, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2_1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                            skip = 166, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2_2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 210, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2_3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 254, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3_1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                            skip = 298, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3_2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 342, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3_3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_8_LK19_checkerboard.txt",
                             skip = 386, n_max = 8, delim="\t")[,2:9]) - blank_avg

# create list with replicates
list_reps <- list(Rep1_1, Rep1_2, Rep1_3, Rep2_1, Rep2_2, Rep2_3, Rep3_1, Rep3_2, Rep3_3)

# take averages and round to second decimal place
avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),3)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 3 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,9)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

# plot heat map
#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Corynebacterium amycolatum LK19",
                  #    display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30,filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK19_Checkerboard.pdf")

# data in long format
data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

# z-score normalization
data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

# create df for binding all data together
allStrains <- data
allStrains$Corynebacterium_amycolatum_LK19 <- allStrains$OD600

# Join sweat and sebum concentrations into one variable
allStrains <- allStrains %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Corynebacterium_amycolatum_LK19)

# z-score normalization
data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

# multiple linear regression
bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

# relative importance of sweat and sebum concentration
calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

# plot sweat or sebum concentration vs OD600 
#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point(aes(color=as.factor(SebumConc))) + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point(aes(color=as.factor(SweatConc)))  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
mlr_plot <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
mlr_plot$Strain <- "Corynebacterium amycolatum LK19"


# Corynebacterium glucuronolyticum LK488 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_19_LK488_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK488_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK488_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg

# create list with replicates
list_reps <- list(Rep1, Rep2, Rep3)

# take averages and round to second decimal place
avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 3 decimal points
with_sd <- avg

sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

# plot heat map
#plot_LK488 <- pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Corynebacterium sp. LK488",display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK488_Checkerboard.pdf")

# data in long format
data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

# z-score normalization
data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

# Join sweat and sebum concentrations into one variable
data_long <- data
data_long$Corynebacterium_glucuronolyticum_LK488 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Corynebacterium_glucuronolyticum_LK488)

# bind data
allStrains <- left_join(allStrains,data_long)

# z-score normalization
data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

# multiple linear regression
bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

# relative importance of sweat and sebum concentration
calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

# plot sweat or sebum concentration vs OD600 
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Corynebacterium glucuronolyticum LK488"
mlr_plot <- rbind(mlr_plot,new)

# Corynebacterium keferresidentii LK1134 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK1134_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK1134_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK1134_checkerboard.txt",
                          skip = 116, n_max = 8, delim="\t")[,2:9]) - blank_avg

# create list with replicates
list_reps <- list(Rep1, Rep2, Rep3)

# take averages and round to second decimal place
avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 3 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Corynebacterium keferresidentii LK1134",
        # display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK1134_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Corynebacterium_kefirresidentii_LK1134 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Corynebacterium_kefirresidentii_LK1134)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)   
new$Strain <- "Corynebacterium kefirresidentii LK1134"
mlr_plot <- rbind(mlr_plot,new)

# Kocuria rhizophila LK221 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK221_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK221_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK221_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Kocuria sp. LK221",
                       #display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK221_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Kocuria_rhizophila_LK221 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Kocuria_rhizophila_LK221)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))

data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))


bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Kocuria rhizophila LK221"
mlr_plot <- rbind(mlr_plot,new)

# Kocuria marina_A LK478 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK478_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK478_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK478_checkerboard.txt",
                           skip = 116, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Kocuria sp. LK478",
                     #  display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK478_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Kocuria_marina_A_LK478 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Kocuria_marina_A_LK478)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Kocuria marina_A LK478"
mlr_plot <- rbind(mlr_plot,new)

# Dietzia cinnamea LK439 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_24_LK439_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_24_LK439_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9])- blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_24_LK439_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9])- blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Dietzia sp. LK439",
                       # display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK439_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Dietzia_cinnamea_LK439 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Dietzia_cinnamea_LK439)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Dietzia cinnamea LK439"
mlr_plot <- rbind(mlr_plot,new)

# Citrobacter freundii LK704 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK704_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_24_LK704_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_24_LK704_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Citrobacter sp. LK704",
                     #  display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK704_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Citrobacter_freundii_LK704 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Citrobacter_freundii_LK704)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Citrobacter freundii LK704"
mlr_plot <- rbind(mlr_plot,new)

# Klebsiella pneumoniae LK469 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK469_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK469_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_3_31_LK469_checkerboard.txt",
                           skip = 116, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Klebsiella sp. LK469",
                     #  display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK469_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Klebsiella_pneumoniae_LK469 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Klebsiella_pneumoniae_LK469)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  
new$Strain <- "Klebsiella pneumoniae LK469"
mlr_plot <- rbind(mlr_plot,new)

# Sphingobacterium hotanense LK485 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK485_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK485_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK485_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Sphingobacterium sp. LK485",
                      # display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK485_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Sphingobacterium_hotanense_LK485 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Sphingobacterium_hotanense_LK485)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Sphingobacterium hotanense LK485"
mlr_plot <- rbind(mlr_plot,new)

# Staphylococcus epidermidis LK593 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_3_LK593_Checkerboard.txt",
                          skip = 30, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_19_LK593_checkerboard.txt",
                           skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_19_LK593_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points

with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Staphylococcus epidermidis LK593",
                       #display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK593_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Staphylococcus_epidermidis_LK593 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Staphylococcus_epidermidis_LK593)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Staphylococcus epidermidis LK593"
mlr_plot <- rbind(mlr_plot,new)

# Staphylococcus epidermidis LK717 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK717_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK717_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_19_LK717_checkerboard.txt",
                           skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Staphylococcus epidermidis LK717",
                  #     display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK717_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Staphylococcus_epidermidis_LK717 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Staphylococcus_epidermidis_LK717)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Staphylococcus epidermidis LK717"
mlr_plot <- rbind(mlr_plot,new)

# Staphylococcus aureus LK376 (USA300) -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK376_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK376_checkerboard.txt",
                           skip = 116, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK376_checkerboard.txt",
                           skip = 158, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Staphylococcus aureus USA300",
        # display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK376_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Staphylococcus_aureus_USA300 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Staphylococcus_aureus_USA300)
allStrains <- left_join(allStrains,data_long)


data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))


bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Staphylococcus aureus USA300"
mlr_plot <- rbind(mlr_plot,new)

# Microbacterium sp. LK369 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_5_LK369_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9])
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_5_LK369_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9])
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_5_LK369_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9])

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±",
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Microbacterium sp. LK369",
        # display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK369_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Microbacterium_sp._LK369 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Microbacterium_sp._LK369)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Microbacterium sp. LK369"
mlr_plot <- rbind(mlr_plot,new)

# Dermabacter hominis LK522 -----------------------------------------

# read in experiments - 6 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_3_LK522_Checkerboard.txt",
                          skip = 30, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK522_checkerboard.txt",
                           skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_10_LK522_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep4 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK522_checkerboard.txt",
                          skip = 32, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep5 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK522_checkerboard.txt",
                           skip = 74, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep6 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_2_LK522_checkerboard.txt",
                           skip = 116, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3, Rep4, Rep5, Rep6)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),3)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,6)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Dermabacter sp. LK522",
                     #  display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK522_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Dermabacter_hominis_LK522 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Dermabacter_hominis_LK522)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
ggplot(data, aes(x=SweatConc, y=OD600, color = as.factor(SebumConc))) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Dermabacter hominis LK522"
mlr_plot <- rbind(mlr_plot,new)

# Dermacoccus nishinomiyaensis LK1128 -----------------------------------------

Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK1128_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK1128_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_2_22_LK1128_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,3)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Dermacoccus sp. LK1128",
                    #   display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/CheckerboardFigures/LK1128_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Dermacoccus_nishinomiyaensis_LK1128 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Dermacoccus_nishinomiyaensis_LK1128)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Dermacoccus nishinomiyaensis LK1128"
mlr_plot <- rbind(mlr_plot,new)

# Micrococcus luteus LK410 -----------------------------------------

# read in experiments - 3 biological replicates
Rep1 <- suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_22_LK410_checkerboard.txt",
                          skip = 34, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep2 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_22_LK410_checkerboard.txt",
                           skip = 78, n_max = 8, delim="\t")[,2:9]) - blank_avg
Rep3 <-  suppressWarnings(readr::read_delim("/Users/mhswaney/Lab_data_permanent/Skin_media_project/DataAnalysis/Checkerboard/2022_4_22_LK410_checkerboard.txt",
                           skip = 122, n_max = 8, delim="\t")[,2:9]) - blank_avg

list_reps <- list(Rep1, Rep2, Rep3)

avg <- reduce(list_reps, `+`) / length(list_reps)
avg <- round(as.matrix(avg),2)

# unlist the replicates, put them into an array with dimensions 8 x 8 and generate 9 arrays from the list
# then apply through the arrays and take the standard deviation of each cell
# round the value to 2 decimal points
with_sd <- avg
sd <- round(apply(array(unlist(list_reps), c(8,8,4)), c(1,2), sd),3)

rownames(avg) <- c("2X","1X","0.5X","0.25X","0.125X","0.0625X","0.031X","0X")
colnames(avg) <- c("0.125%","0.0625%","0.031%","0.016%","0.008%","0.004%","0.002%","0%")

# add the avg OD value to the standard deviation in each cell
with_sd[-1] <- paste0(as.matrix(with_sd[-1]), "±", 
                      as.matrix(sd[-1]))

#pheatmap(avg,cluster_rows = FALSE, cluster_cols = FALSE, main = "Micrococcus sp. LK410",
       #  display_numbers = data.matrix(with_sd), cellwidth = 53, cellheight = 30, filename = "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/InProgress/CheckerboardFigures/LK410_Checkerboard.pdf")

data <- melt(avg) %>% rename(SweatConc = Var1, SebumConc = Var2, OD600 = value)

data$OD600 <- (data$OD600  - mean(data$OD600 )) / sd(data$OD600)

data_long <- data
data_long$Micrococcus_luteus_LK410 <- data_long$OD600
data_long <- data_long %>% unite("SweatSebum", SweatConc:SebumConc, remove=TRUE) %>% dplyr::select(SweatSebum,Micrococcus_luteus_LK410)
allStrains <- left_join(allStrains,data_long)

data$SweatConc <- as.numeric(substr(as.character(data$SweatConc), 1, nchar(as.character(data$SweatConc))-1))
data$SebumConc <- as.numeric(substr(as.character(data$SebumConc), 1, nchar(as.character(data$SebumConc))-1))
data$SweatConc <- (data$SweatConc-min(data$SweatConc))/(max(data$SweatConc)-min(data$SweatConc))
data$SebumConc <- (data$SebumConc-min(data$SebumConc))/(max(data$SebumConc)-min(data$SebumConc))

bothfit <- lm(data$OD600~data$SebumConc + data$SweatConc)
summary(bothfit)

calc.relimp(bothfit, type=c("lmg","last","first","pratt"), rela=TRUE)

#dev.off()
#ggplot(data, aes(x=SweatConc, y=OD600)) + geom_point() + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SebumConc)), size=0.5,method="lm")
#ggplot(data, aes(x=SebumConc, y=OD600)) + geom_point()  + theme_classic() + geom_smooth(se = FALSE,aes(color=as.factor(SweatConc)), size=0.5,method="lm")

#extract MLR coefficients
sebumcoef <- bothfit$coefficients[2]
sweatcoef <- bothfit$coefficients[3]
sebumpval <- summary(bothfit)$coefficients[11]
sweatpval <- summary(bothfit)$coefficients[12]

# create df for binding all regression data together
new <- data.frame(x_start=0, y_start=0, x_end=sweatcoef, y_end=sebumcoef, 
                       sweat_pval = sweatpval,sebum_pval = sebumpval)  

new$Strain <- "Micrococcus luteus LK410"
mlr_plot <- rbind(mlr_plot,new)
