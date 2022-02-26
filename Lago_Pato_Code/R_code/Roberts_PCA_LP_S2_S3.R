# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Set working directory ---------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021")
#check working directory
getwd()


# Load libraries ----------------------------------------------------------

##load libraries
library(tidyverse)
library(tidypaleo)
library(ggplot2)
library(readr)
library(vegan) # for paleo analysis
library(rioja) # for paleo analysis
library(repr)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(ellipse)  # for PCA and cluster
library(dplyr)  # for transforming data
library(factoextra) # for PCA and cluster
library(reshape2)
library(GGally)
library(ggsci)
library(ggdendro)
library(dendextend)
library(dynamicTreeCut)
library(colorspace)
library(cluster)
library(ggpubr)
library(cowplot) # for plotting
library(ggfortify) # for time series
library(magrittr) # for piping %>%
library(mgcv)
library(gridExtra)
library(gtable)
library(repr)
library(patchwork)
library(cowplot)
library(gridGraphics)
#ITRAX-specific 
library(itraxR) # not used 
library(PeriodicTable)
library(egg)
library(devtools)
library(tidyr)
library(tidypaleo)
library(rbacon)
library(corrplot)
library(chemometrics)
library(stringr)
library(egg)
library(grid)
library(data.table)
library(bestNormalize)
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)
#RColorBrewer
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

# Customise BrBG colour palette from 11 to 18 colours & get colour codes to copy  --------
nb.cols <- 11
BrBG1 <- colorRampPalette(brewer.pal(11, "BrBG"))(nb.cols)
BrBG1
# BrBG1 11 colour codes 
c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30")

# Custom RdYlBlu colour palettes from 11 to 18 colours & get col --------
nb.cols <- 11
RdYlBu <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
RdYlBu
c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

nb.cols <- 5
RdYlBu1 <- colorRampPalette(brewer.pal(5, "RdYlBu"))(nb.cols)
RdYlBu1
# RdYlBu1 colour codes 
c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

nb.cols <- 7
RdYlBu2 <- colorRampPalette(brewer.pal(7, "RdYlBu"))(nb.cols)
RdYlBu2
# RdYlBu1 colour codes 
c("D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")

nb.cols <- 10
RdYlBu3 <- colorRampPalette(brewer.pal(10, "RdYlBu"))(nb.cols)
RdYlBu3
c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

nb.cols <- 4
RdYlBu4 <- colorRampPalette(brewer.pal(4, "RdYlBu"))(nb.cols)
RdYlBu4
c("#D7191C", "#FDAE61", "#ABD9E9", "#2C7BB6")

nb.cols <- 11
RdBu1 <- colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols)
RdBu1
c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")


# ITRAX PCA LP08 & LP16 colour scheme and elements list -------------------------------------

# LP08 and LP16 Unit colours for Unit 1, 5B, 5C, 6 based on RdYlBu & RdBu
unit_col_LP08 <- c("#4575B4","#FDAE61", "#F46D43", "#B2182B") #Unit 1, 5B, 5C, 6
unit_col_LP16 <- c("#4575B4", "#74ADD1", "#ABD9E9", "#ABD9E9", "#FDAE61", "#F46D43", "#B2182B") #Unit 1, 2, 3, 4, 5A, 5B, 5C, 6
cluster_col5 <-  c("#A50026","#F46D43","#FEE090","#ABD9E9","#4575B4")

# Elements in LP08 PCA filtered as >0.5% mean and >0.1 max TSN% value and well defined signal (i.e., not noise)

#Different elements
elements_LP08 <- c("Si","S","K","Ca","Ti","Mn","Fe","Zn","Br","Rb","Sr","Zr","inc","coh","inc_coh", "coh_inc") #n=12 elements

elements_LP16 <- c("Si","S","K","Ca","Ti","V","Cr","Mn","Fe","Ni","Zn","As","Br","Rb","Sr","Zr","Ba","inc","coh","inc_coh","coh_inc") #n=17 elements

#Combined match
elements_LP08_16 <- c("Si","S","K","Ca","Ti","Mn","Fe","Zn","Br","Rb","Sr","Zr","inc","coh","inc_coh", "coh_inc") #n=12 elements as LP08

# Load period table data to create a list of element symbols in the $symb column 
data("periodicTable")
# symb(1:116)
# TO DO - write code to map elements onto periodic table output - so code will work with any combination automatically 


# Import and select, filter data  ------------------------------------------------------------

# Link to folder with datafiles in: https://www.dropbox.com/sh/u4d9jwujokw8qb0/AACyalT2Mpn4Pr6Daj_bqGR9a?dl=0
# First, choose with interval dataset to import - i.1, 200um, 2mm or 1cm for each record
PCA_db_LP08 <- read_csv("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Data/LP08_200um_TSN.csv")
PCA_db_LP08

PCA_db_LP16 <- read_csv("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Data/LP16_200um_TSN.csv")
PCA_db_LP16

# Dataframe to take forward for analysis and plotting  
PCA_df_LP08 <- select(PCA_db_LP08, depth, elements_LP08_16, Subunit, Group) %>% 
  replace(is.na(.), 0) %>% #convert NA to 0
  mutate_if(is.logical,as.numeric) #to convert SO2 FALSE to 0 - overcomes parsing problem with SO2
# To convert NA to 0 in base R: EPMA_df[is.na(EPMA_df)] <- 0
PCA_df_LP08

# OR

PCA_df_LP16 <- select(PCA_db_LP16, depth, elements_LP08_16, Subunit, Group) %>% 
  replace(is.na(.), 0) %>% #convert NA to 0
  mutate_if(is.logical,as.numeric) #to convert SO2 FALSE to 0 - overcomes parsing problem with SO2
# To convert NA to 0 in base R: EPMA_df[is.na(EPMA_df)] <- 0
PCA_df_LP16

# Convert data to long format for ggplot and tidypaleo facet plotting --------------------------

# select out the parameters to plot, and convert to long format - i.e., one value per row per variable - for later plotting .... perhaps ... 
PCA_df_long_LP08 <- select(PCA_df_LP08,  depth, elements_LP08_16, Subunit, Group) %>%
  pivot_longer(elements_LP08_16, names_to = "param", values_to = "value")
  #relocate(param, .before = Type)
PCA_df_long_LP08

PCA_df_long_LP16 <- select(PCA_df_LP16,  depth, elements_LP08_16, Subunit, Group) %>%
  pivot_longer(elements_LP08_16, names_to = "param", values_to = "value")
#relocate(param, .before = Type)
PCA_df_long_LP16

# Figure S2 & S3 - Correlation & PCA -----------------------------------------------------

# Choose dataset to input for processing for PCA in base R-----------------------------------------------

# Set up and select columns to use - allows to make subset at this point if different from above
var_sub_n <- PCA_df_LP08
ME_n <- elements_LP08_16

# OR

var_sub_n <- PCA_df_LP16
ME_n <- elements_LP08_16


# # Import into PCA -------------------------------------------------------

kta <- PCA_df_LP08

# OR

kta <- PCA_df_LP16

# remove rows if necessary
#kta <- kta[,(var_sub_LP08),drop=FALSE]
# kta <- kta[!kta$Group=='ANT',]
# kta <- kta[!kta$Group=='DI',]

# Transform data ----------------------------------------------------------

head(kta)
#plot(kta[, ME_n], pch=19, cex = 0.05)

library(bestNormalize)
kta.bc <- select(kta, Si:coh_inc) %>%
  replace(is.na(.), 0) %>% #convert NA to 0
  mutate_if(is.logical,as.numeric) #to convert any FALSE to 0 - overcomes parsing problem
# To convert NA to 0 in base R: EPMA_df[is.na(EPMA_df)] <- 0
# Apply bestNormalise to see which transformation makes the data "most Normal'  - using Ti as example vector
kta.bc
hist(kta.bc$Br)
(kta.bc.best <- bestNormalize(kta.bc$Br, r = 1, k = 5))
kta.bc.best
# write.csv(kta.bc.best,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta_best_norm.csv", row.names = TRUE)

(arcsinh_obj <- arcsinh_x(kta.bc$Ti))
(yeojohnson_obj <- yeojohnson(kta.bc$Ti))
(orderNorm_obj <- orderNorm(kta.bc$Ti))

# Standardise and centre variables -  Z-scores using scale() function for PCA analysis - values are mean of 0 and +/-1 of 1 std dev
# copy the original filename to a new filename, apply a z-score transform, look at it,plot it
kta.Z <- kta
kta.Z[, ME_n] <- scale(kta[ME_n], center = TRUE, scale = TRUE)
kta.Z[is.na(kta.Z)] <- 0
kta.Z
#plot(kta.Z[, ME_LP08], pch=19, cex = 0.05)

# Sqrt transform kta = %kta sq root-transformed, cuberoot = %TSBN cuberoot transformed
kta.sqrt <- kta
kta.sqrt[, ME_n] <- sqrt (kta.sqrt[ME_n])
kta.sqrt[is.na(kta.sqrt)] <- 0
kta.sqrt
#plot(kta.sqrt[, LP08], pch=19, cex = 0.05)

# Standardise and centre the sqrt dataset - calculate Z-scores of sqrt transformed data 
kta.sqrt.Z <- kta.sqrt
kta.sqrt.Z[, ME_n] <- scale(kta.sqrt[ME_n], center = TRUE, scale = TRUE)
kta.sqrt.Z
#plot(kta.sqrt.Z[, ME_LP08], pch=19, cex = 0.05)

# Log transform kta (ln = natural log) using ln(x+1) to set 0 from inf. back to 0
kta.ln <- kta
kta.ln[, ME_n] <- log (kta.ln[ME_n]+1)
kta.sqrt[is.na(kta.ln)] <- 0
kta.ln
#plot(kta.sqrt[, LP08], pch=19, cex = 0.05)

# Standardise and centre the ln dataset - calculate Z-scores of naturallog transformed data 
kta.ln.Z <- kta.ln
kta.ln.Z[, ME_n] <- scale(kta.ln[ME_n], center = TRUE, scale = TRUE)
kta.ln.Z
#plot(kta.sqrt.Z[, ME_LP08], pch=19, cex = 0.05)

#  Write transformed data to file --------------------------------------------------
write.csv(kta.Z,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta_Z.csv", row.names = TRUE)
write.csv(kta.sqrt,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta_sqrt.csv", row.names = TRUE)
write.csv(kta.sqrt.Z,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta.sqrt.Z.csv", row.names = TRUE)
write.csv(kta.ln,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta_ln.csv", row.names = TRUE)
write.csv(kta.ln.Z,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/kta.ln.Z.csv", row.names = TRUE)

# Correlation plots  --------------------------------------------

#clear plot window
dev.off()
#reset parameters to default
.pardefault <- par()
par(.pardefault)

library(GGally)


# Basic correlation matrix -------------------------------------------------------

# Correlation plot - untransformed - all elements, density plots and correlation matrxi summary
theme_set(theme_bw(base_size=6))
# use this to see where positive/significant correlations as an overview overall
ggcorr(kta[, ME_n], method = c("everything", "pearson"), size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 

## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr_matrix_untrans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr_matrix_untrans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation plot - log transformed, scaled and centered - all elements, density plots and correlation matrxi summary
theme_set(theme_bw(base_size=6))
# use this to see where positive/significant correlations as an overview overall
ggcorr(kta.ln[, ME_n], method = c("everything", "pearson"), size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 

## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr_matrix_ln_trans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr_matrix_ln_trans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation plot - log transformed, scaled and centered - all elements, density plots and correlation matrxi summary
theme_set(theme_bw(base_size=6))
# use this to see where positive/significant correlations as an overview overall
ggcorr(kta.sqrt[, ME_n], method = c("everything", "pearson"), size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 

## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr_matrix_sqrt_trans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr_matrix_sqrt_trans.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation-density matrix ----------------------------------------------

# Correlation density matrix plot on untransformed data
ggpairs(kta, columns = ME_n, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr-den_matrix_untrans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr-den_matrix_untrans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Correlation density matrix plot on log-transformed data
ggpairs(kta.ln, columns = ME_n, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr-den_matrix_ln_trans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr-den_matrix__ln_trans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Correlation density matrix plot on sqrt-transformed data
ggpairs(kta.ln, columns = ME_n, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
## LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr-den_matrix_sqrt_trans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
## LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr-den_matrix__sqrt_trans.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Correlation-density plot - grouped by Unit ------------------------------

# Correlation density Unit plot - untransformed data - scatterplot, density dist and stats for each unit
# limit this to 5 x 5 matrix as text difficult to read 
# and assign the matrix plot to a variable and then iterating through each plot with a for-loop
theme_set(theme_bw(base_size=8))
plot.matrix <- ggpairs(kta, columns = ME_n, upper = list(continuous = wrap("cor", size = 1)),
                       lower = list(continuous = wrap("points", alpha = 1, size=0.5)),
                       ggplot2::aes(colour = Group, title="Correlation plot by Group"))
plot.matrix

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr-den-unit_matrix_untrans.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr-den-unit_matrix_untrans.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Correlation density Unit plot - log transformed, scaled and centered data - scatterplot, density dist and stats for each unit
# limit this to 5 x 5 matrix as text difficult to read 
# and assign the matrix plot to a variable and then iterating through each plot with a for-loop
theme_set(theme_bw(base_size=8))
plot.matrix <- ggpairs(kta.ln, columns = ME_n, upper = list(continuous = wrap("cor", size = 1)),
                       lower = list(continuous = wrap("points", alpha = 1, size=0.5)),
                       ggplot2::aes(colour = Group, title="Correlation plot by Group"))
plot.matrix

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Corr-den-unit_matrix_trans.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Corr-den-unit_matrix_trans.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")


# new version from Nov 2021 update - need to chnage size of the density plot lines 
tmp.plot <- ggpairs(
  data = kta, 
  mapping = ggplot2::aes(colour = Eruption_ID, size = 0.5, labelSize = 5),
  columns = ME_n,
  diag = list(continuous = "densityDiag", size = 0.5, combo = "box_no_facet"), 
  upper = list(continuous = wrap("cor", size = 2.5)),
  lower = list(continuous = wrap("points", alpha = 1, size=0.5)),
  axisLabels = "show",
) 

tmp.plot 


# PCA  --------------------------------------------------------------------

## Visualising PCA with the "factoextra" package and ggplot2

#clear plot window
dev.off()

# select dataset to take forward to PCA analysis

# ln transformed, scaled and centered data
pca.kta <- kta.ln.Z [ME_n]
head(pca.kta)
#plot(pca.kta, pch=19, cex = 0.05)

# OR 

# sqrt transformed scaled and centered data
pca.kta <- kta.sqrt.Z [ME_n]
head(pca.kta)
#plot(pca.kta, pch=19, cex = 0.05)

# Examine correlation matrix for PCA df
cor_kta <- cor(pca.kta)
round(cor_kta, 2)
# p-values
library("Hmisc")
cor1_kta <- rcorr(as.matrix(pca.kta))
cor1_kta

# Examine co-variance in the PCA df
res.cov <- cov(pca.kta)
round(res.cov,2)
# p-values
cor1_kta <- rcov(as.matrix(pca.kta))
cor1_kta


# Make PCA not including Unit columns (1 in csv file) - depth is --------

p <- prcomp(pca.kta,  scale = TRUE)
# The coordinates for PCA are contained within the prcomp object as a matrix:res.pca$x

# Make PC covariance matrix 
r = eigen(cov(pca.kta))
p$rotation

# Scree plots above and below to assess how many PC to carry forward --------

# Percentage of variance explained by the first four PC dimensions
ev <- p$sdev^2
var_percent <-ev[1:4] / sum(ev)
var_percent
# Plot as a scree plot 
pca_scree <- fviz_eig(p) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

pca_scree

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_scree.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_scree.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# # Screeplot of PC Inertia & writing stats to file -----------------------------------------------

if(!is.null(dev.list())) dev.off()
layout(matrix(1:1, ncol=1))
p_inertia <- screeplot(p)

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_inertia.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_inertia.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Show data
PC14 <- p$x[,1:4]
head(PC14)
#Show summary stats for variables using standardize package - and standardize data (to finish)
x.stats <- summary(PC14)
x.stats

# write LP08 & LP16 PC data and stats to csv
write.csv(pca.kta,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/pca.kta.csv", row.names = TRUE)
write.csv(PC14,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/PC1-4.csv", row.names = TRUE)
write.csv(var_percent,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/PC1-4_var_percent.csv", row.names = TRUE)
write.csv(x.stats,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Outputs/ITRAX/PC1-4_stats.csv", row.names = TRUE)

# Plot PC axes ------------------------------------------------------------

# Plot the first 2-4 principal components, so use p$x[,1] (PCA 1) and p$x[,2] (PCA 2) etc. to select the relevant data. 
# Plot datapoints to check - use change axes to = 1,3 for PC1 vs PC3
# cos2 = quality of individual points on the factor map
# contib = relative contribution of individual datapoints to PCA score
fviz_pca_ind(p, geom = c("point", "text"), axes = c(1,2), labelsize = 1, col.ind = "cos2", 
             gradient.cols = c("grey90","#2E9FDF", "#FC4E07"),
             repel = FALSE, title = "PCA - Biplot with depths labelled") + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_contrib.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_contrib.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# PC axes loadings --------------------------------------------------------

if(!is.null(dev.list())) dev.off()
#Create a barplot for loadings on PC1 axis
barplot(p$rotation[,1], main="PC 1 Loadings Plot", las=2) # the loadings for the variables for first principal component (PC) are p$rotation[,1], 
# create a numeric object ("n.pc1") containing values to position text underneath the bars 
# p$rotation[,1] > 0 sets the condition or test for ifelse() to meet. 
# yes argument sets the value specified if the condition is met,no argument returns value if not
# c.pc1 adds colours and replot
n.pc1 <- ifelse(p$rotation[,1] > 0, yes=-0.01, no=p$rotation[,1]-0.01)
c.pc1 <- ifelse(p$rotation[,1] > 0, yes="darkblue", no="darkred")
# plot diagram
par(mar=c(8,3,2,1)) # Set margins
b1 <- barplot(p$rotation[,1], main="PC1 Loadings Plot", col=c.pc1, las=2, axisnames=FALSE)
abline(h=0) # Add horizontal line
text(x=b1, y=n.pc1, labels=names(p$rotation[,1]), adj=1, srt=90, xpd=FALSE) # Add variable names
# axisnames=FALSE stops barplot() plotting automatically
# creates the bar chart as an object to extract the "midpoints" of each bar to correctly position the variable names when using the text() function.
# adj=1 sets the alignment of the variable names to right align.
# srt=90 changes the direction of the text to a 90 degree angle (vertical).
# xpd=TRUE tells R to plot the text outside the plot region, and within the figure region

# Plot all PC axes loadings together --------------------------------------

if(!is.null(dev.list())) dev.off()
# Change colour of bar plot
c.pc1 <- ifelse(p$rotation[,1] > 0, "darkblue","darkred")
c.pc2 <- ifelse(p$rotation[,2] > 0, "darkblue","darkred")
c.pc3 <- ifelse(p$rotation[,3] > 0, "darkblue","darkred")
c.pc4 <- ifelse(p$rotation[,4] > 0, "darkblue","darkred")
# Get position for variable names
n.pc1 <- ifelse(p$rotation[,1] > 0, -0.01, p$rotation[,1]-0.01)
n.pc2 <- ifelse(p$rotation[,2] > 0, -0.01, p$rotation[,2]-0.01)
n.pc3 <- ifelse(p$rotation[,3] > 0, -0.01, p$rotation[,3]-0.01)
n.pc4 <- ifelse(p$rotation[,4] > 0, -0.01, p$rotation[,4]-0.01)

# Plot PC1-4 loadings - 8 cm wide plot ------------------------------------
layout(matrix(1:4, ncol=1)) # Set up layout
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149608 inches but R converts to 3.149606 - 
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.2,0.2), pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 
#dev.off() #need to include this to write to pdf file fully - will delete screen plot
#par(mar=c(2,4,4,2), oma=c(8,0,0,0)) # Set up margins

# Plot PC 1, 2, 3 & 4 together add variable names
b1 <- barplot(p$rotation[,1], main="PC1 Loadings Plot", col=c.pc1, las=3, axisnames=FALSE)
abline(h=0)
text(x=b1, y=ifelse(p$rotation[,1] > 0, -0.01, p$rotation[,1]-0.01), labels=names(p$rotation[,1]), adj=1, srt=90, xpd=NA)
b2 <- barplot(p$rotation[,2], main="PC2 Loadings Plot", col=c.pc2, las=3, axisnames=FALSE)
abline(h=0)
text(x=b2, y=ifelse(p$rotation[,2] > 0, -0.01, p$rotation[,2]-0.01), labels=names(p$rotation[,2]), adj=1, srt=90, xpd=NA)
b3 <- barplot(p$rotation[,3], main="PC3 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,3] > 0, -0.01, p$rotation[,3]-0.01), labels=names(p$rotation[,3]), adj=1, srt=90, xpd=NA)
b4 <- barplot(p$rotation[,4], main="PC4 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,4] > 0, -0.01, p$rotation[,4]-0.01), labels=names(p$rotation[,4]), adj=1, srt=90, xpd=NA)

# print plot screen to save - will be same dimensions are screen plot
# LP08 save
dev.print(pdf, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_load.pdf")
# LP16 save
dev.print(pdf, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_load.pdf")

# Plot variable lines and data on chosen PC axes biplots ----------------------------------

if(!is.null(dev.list())) dev.off()

theme_set(theme_minimal(base_size=16))
var1 <- fviz_pca_var(p, axes = c(1,2), col.var = "contrib",
             gradient.cols = c("brown", "blue"), #"Blues",
             title = "PCA - Variable Plot: PC1/PC2")

var2 <- fviz_pca_var(p, axes = c(1,3), col.var = "contrib",
             gradient.cols = c("brown", "blue"), #"Blues",
             title = "PCA - Variable Plot: PC1/PC3")

var3 <- fviz_pca_var(p, axes = c(2,3), col.var = "contrib",
             gradient.cols = c("brown", "blue"), #"Blues",
             title = "PCA - Variable Plot: PC2/PC3")

var4 <- fviz_pca_var(p, axes = c(2,4), col.var = "contrib",
                     gradient.cols = c("brown", "blue"), #"Blues",
                     title = "PCA - Variable Plot: PC2/PC4")

library(cowplot)
pca_var <- plot_grid(var1, var2, var3, var4,
  labels = "AUTO", ncol = 2
)
pca_var

file4 <- tempfile("PCA_var", fileext = ".pdf")
save_plot(file4, pca_var, ncol = 2, base_asp = 1)

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_var.pdf", 
       height = c(30), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_var.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Final PCA graph -------------------------------------------------------------

# define custom colour scheme based on Brewer RdYlBu
# LP08 and LP16 Uni colours for Unit 1, 5B, 5C, 6 based on RdYlBu & RdBu
unit_col_LP08 <- c("#4575B4","#FDAE61", "#F46D43", "#B2182B") #Unit 1, 5B, 5C, 6
unit_col_LP16 <- c("#4575B4", "#74ADD1", "#ABD9E9", "#ABD9E9", "#FDAE61", "#F46D43", "#B2182B") #Unit 1, 2, 3, 4, 5A, 5B, 5C, 6

# choose unit colour first

unit_col <- unit_col_LP08

# OR

unit_col <- unit_col_LP16


# Create PCA PC1 vs PC2 plot with legend for later on --------------------------------

# Biplot of individuals and variables - Keep only the labels for variables - Color by groups, add ellipses
# build four plots without a legend - but use pca0 to check working ok and add legend back in later
# legend omitted for pca1, pca2, pca3, pca4 - and added back in in cowplot - comment out for pca1 to check working ok
theme_set(theme_classic(base_size=12))
pca0 <- fviz_pca_biplot(p, axes = c(1,2), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC1/PC2 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="bottom")
pca0
# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PC1_vs_PC2_legendbot.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PC1_vs_PC2_legendbot.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Make 4 part plot with most commonly useful PC axes combinations -------

# PC1 vs PC2
theme_set(theme_classic(base_size=12))
pca1 <- fviz_pca_biplot(p, axes = c(1,2), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC1/PC2 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="none")
pca1
# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PC1_vs_PC2_no_legend.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PC1_vs_PC2_no_legend.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# PC1 vs PC3
pca2 <- fviz_pca_biplot(p, axes = c(1,3), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC1/PC3 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="none")
pca2
# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PC1_vs_PC3.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PC1_vs_PC3.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# PC2 vs PC3
pca3 <- fviz_pca_biplot(p, axes = c(2,3), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC2/PC3 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="none")
pca3
# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PC2_vs_PC3.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP016 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PC2_vs_PC3.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# PC2 vs PC4
pca4 <- fviz_pca_biplot(p, axes = c(2,4), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC2/PC4 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="none")
pca4
# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PC2_vs_PC4.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PC2_vs_PC4.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Final 4 part plot -----------------------------------------------------------------------
# Build grid plot without legends and save
pca_grid <- plot_grid(pca1, pca2, pca3, pca4, 
                      labels = "AUTO", label_size = 16, ncol = 2)
pca_grid

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_grid_final.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_grid_final.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# use PC0 for bottom horizontal legend
pca0

# for right vertical legend do this ... 
pca_legend <- fviz_pca_biplot(p, axes = c(1,2), label = "var", col.var = "blue",
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Group,
                        addEllipses=TRUE, ellipse.level=0.95, 
                        palette = unit_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC1/PC2 and Variables (95% CI ellipses)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position="right")

# hide legend and replot - then save to temp file to transfer to ggsave
grobs <- ggplotGrob(pca_legend)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
pca_grid_legend <- plot_grid(pca_grid, legend, labels = "AUTO",
                       ncol = 2, rel_widths = c(1, .1))
pca_grid_legend

file4 <- tempfile("PCA_final", fileext = ".pdf")
save_plot(file4,pca_grid_legend , ncol = 2, base_asp = 1)

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/PCA_final_legend.pdf", 
       height = c(30), width = c(35), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/PCA_final_legend.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Hierarchical clustering --------------------------------------------------------

#clear plot window
dev.off()
#reset parameters to default
.pardefault <- par()
par(.pardefault)

# Choose dataset to cluster
kta.cluster <- kta.sqrt.Z[ME_n]
# OR
kta.cluster <- kta.ln.Z[ME_n]

head(kta.cluster)

#  -----------------------------

# Determine the optimal number of clusters to use -----------------------------------------------------------------------

# USE THIS ONE - optimal number of clusters using Gap Statistic Method - bootstrapping MCMC approach (Tibshirani et al., 2001)
gap_stat <- clusGap(kta.cluster, FUN = kmeans, nstart = 25,
                    K.max = 15, B = 50)
p9 <- fviz_gap_stat(gap_stat) 
p9 +theme(text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)
)

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Cluster_optimal.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Cluster_optimal.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# optimal number of clusters usng average silhouette method (https://uc-r.github.io/kmeans_clustering#silo) 
fviz_nbclust(kta.cluster, FUN = hcut, method = "silhouette")

# k-means clustering - Dim 1 vs Dim 2 ------------------------------------------------------

# Visualize kmeans clustering using repel = TRUE to avoid overplotting - nstart = no. of groups defined by optimal plot
theme_set(theme_classic(base_size=16))
km.res <- kmeans(kta.cluster, 5, nstart = 5)
fviz_cluster(km.res, geom = c("point"), kta.cluster, ellipse.type = "euclid", ellipse.level=0.95, 
             palette = BrBG1, ggtheme = theme_minimal()+
               theme(text = element_text(size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     legend.position="right"))

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Cluster_dim.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Cluster_dim.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Hierachical clustering - Vegan - Dimension plot ------------------------------------------

theme_set(theme_classic(base_size=12))
# Use hcut() which computes hclust and cut the tree - 
hc.cut <- hcut(kta.cluster, k = 5, hc_method = "ward.D2") #can be complete, single, average or wards
# Visualize dendrogram - 
fviz_dend(hc.cut, show_labels = TRUE, labelsize = 1, rect = TRUE) + 
  theme(text = element_text(size = 7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

fviz_cluster(hc.cut, ellipse.type = "convex", ellipse.level=0.95, 
             palette = unit_col, ggtheme = theme_minimal()+
               theme(text = element_text(size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     legend.position="right"))

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Cluster_convex.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Cluster_convex.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# K means clustering dendrogram ---------------------------------------------------------------

#clear plot window
dev.off()

# Color as a function of the cluster - using hc as defined above - need to define colours to match number of groups
par(mar=c(4,1,1,7))
# custom colour for branching
cluster_col5 <-  c("#A50026","#F46D43","#FEE090","#ABD9E9","#4575B4")
# Run clustering with just data matrix
kta.matrix <- as.matrix(kta.sqrt.Z[, ME_n])
head(kta.matrix)
# Link row names from original matrix
rownames(kta.matrix) <- kta.sqrt.Z$depth
dend <- as.dendrogram(hclust(d = dist(x = kta.matrix), "ward.D2"))

#par(mfrow = c(1,2), mar = c(2,2,1,0))
dend <- dend %>%
  color_branches(k = 5) %>%
  set("branches_lwd", 0.75) %>%
  set("branches_lty", 1) %>%
  set("labels_cex", 0.01) %>%
  set("labels_col", value = cluster_col5, k=5) %>%
  set("branches_k_color", value = cluster_col5, k=5)
# use this if havent set the colours: dend <- color_labels(dend, k = 10)
# code is the same as:labels_colors(dend)  <- get_leaves_branches_col(dend)
p10 <- plot(dend, horiz=TRUE, axes=TRUE)

# LP08 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP08/Cluster_final.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# LP16 save
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/ITRAX/LP16/Cluster_final.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

#get a list of cluster group members
dend_list <- get_subdendrograms(dend, 5)
cluster_list <- lapply(dend_list, labels)
library(data.table) # install if not installed already
fwrite(list(cluster_list), file = "cluster_list.csv")

#10 colour scheme 
#c("skyblue", "steelblue3", "blue3", "burlywood3", "chocolate4", "salmon2", "red", "darkred", "darkorange",  "darkmagenta")

# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Figures/FigureS12B.pdf", height = c(29.7), width = c(21), dpi = 600, units = "cm")

# Dendroplot in ggplot ----------------------------------------------------

# Color as a function of the cluster - using hc as defined above - need to define colours to match number of groups
par(mar=c(4,1,1,7))
# Run clustering with Long ID names 
kta.matrix <- as.matrix(kta.sqrt.Z[, ME_n])
head(kta.matrix)
# Link row names from original matrix
rownames(kta.matrix) <- kta.sqrt.Z$LongID
kta_dendro <- as.dendrogram(hclust(d = dist(x = kta.matrix), "ward.D2"))
# Create dendro
dendro.plot <- ggdendrogram(data = kta_dendro, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 1))
# Preview the plot
print(dendro.plot)


# END ---------------------------------------------------------------------

