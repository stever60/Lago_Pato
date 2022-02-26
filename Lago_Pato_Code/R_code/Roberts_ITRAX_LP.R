# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Load libraries & colours ----------------------------------------------------------

##load libraries
library(tidyverse)
library(tidypaleo)
library(readr)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(cowplot) # for plotting
library(vegan)
library(rioja)
library(ellipse)  # for PCA and cluster
library(factoextra) # for PCA and cluster
library(reshape2)
library(GGally)
library(ggsci)
library(ggdendro)
library(dendextend)
library(dynamicTreeCut)
library(colorspace)
library(cluster)
library(magrittr) # for piping %>%
library(mgcv)
library(gtable)
library(repr)
library(bestNormalize)
library(sjmisc)
library(chemometrics)
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)
#RColorBrewer
display.brewer.all()
display.brewer.pal(11,"BrBG")
display.brewer.all(colorblindFriendly = TRUE)
# Show BrBG colour palette with 11 colours & get colour codes to copy  --------
nb.cols <- 11
PiYG1 <- colorRampPalette(brewer.pal(11, "PiYG"))(nb.cols)
PiYG1
nb.cols <- 9
Greys1 <- colorRampPalette(brewer.pal(9, "Greys"))(nb.cols)
Greys1
nb.cols <- 9
Greens1 <- colorRampPalette(brewer.pal(9, "Greens"))(nb.cols)
Greens1
nb.cols <- 9
Blues1 <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)
Blues1

# -------------------------------------------------------------------------

# PARTS -------------------------------------------------------------------

# PART 1 - LP08 Composite SH20-M1H
# PART 2 - Lp16 Composite Sh20-M1H

# -------------------------------------------------------------------------

# TO DO
# - import and composite ITRAX datasets done in excel into R / Tidyverse
# - write code to add age-depth model columns from BACON output

# -------------------------------------------------------------------------

# Set working directory ---------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021")
#check working directory
getwd()

# PART 1 - LP08-ITRAX-SH20  ------------------------------------------------------------------

# LP08-ITRAX-SH20 - Import cps data ---------------------------------------------
library(dplyr)

# Location of LP08 and LP16-ITRAX-COMP-SH20 composite datasets in dropbox - add figure doi when done 
# https://www.dropbox.com/s/6ovlxp4bgx9v9p6/LP08_ITRAX_COMP_2mm_SH20.csv?dl=0
# https://www.dropbox.com/s/fqtewao1yqnfvvu/LP016_ITRAX_COMP_2mm_SH20.csv?dl=0
# https://www.dropbox.com/s/dk2o9b9497zirco/LP16_ITRAX_COMP_200um_SH20.csv?dl=0

# Composite LP08 - ITRAX cps data & Sh20 BACON age-depth model (added in Excel)
# Find . replace all LP08 with LP16 touse LP16 dataset
LP08_x.raw <- read_csv("Data/LP08_ITRAX_COMP_2mm_SH20.csv")
LP08_x.raw

# OR

LP16_x.raw <- read_csv("Data/LP16_ITRAX_COMP_2mm_SH20.csv")
LP16_x.raw

# OR

LP16_x.raw <- read_csv("Data/LP16_ITRAX_COMP_200um_SH20.csv")
LP16_x.raw

# Filter data to keep only valid data and remove kcps <-2sd and MSE >2sd ----------------------------
LP08_kcps.mean <- mean(LP08_x.raw$kcps)
LP08_kcps.sd <- 2*sd(LP08_x.raw$kcps)
LP08_kcps.thres <- LP08_kcps.mean - LP08_kcps.sd 
LP08_kcps.thres

LP08_MSE.mean <- mean(LP08_x.raw$MSE)
LP08_MSE.sd <- 2*sd(LP08_x.raw$MSE)
LP08_MSE.thres <- LP08_MSE.mean + LP08_MSE.sd 
LP08_MSE.thres

# filter by validity and then remove analyses above MSE threshold and below kcps - 2s & 

# Leave as is (if already filtered)
LP08_COMP <- LP08_x.raw

# OR filter based on MSE and kcps thresholds
LP08_COMP <- filter(LP08_x.raw, validity == "1") %>% 
  filter(MSE < LP08_MSE.thres) %>%
  filter(kcps > LP08_kcps.thres)

# Remove unwanted columns - might need to change these - and write to file
LP08_COMP_filter <- select(LP08_COMP, 
                           -c(filename, 
                              min, max,
                              median, mean_2s_range:Foffset,
                              D1, S1, S2, S3)) %>% 
  rename(SH20_age = mean)

write.csv(LP08_COMP_filter,"Output/LP08/1_cps_filter.csv", row.names = FALSE)

# Define elements lists -------------------------------------------------
library(PeriodicTable)

# list of all possible elements
elements <- c(symb(1:117), "Mo_inc", "Mo_coh", "TS_sum", "cps_sum")

# list of all elements from Q-spec matching
LP08_elements <- select(LP08_COMP_filter, c(Mg:Mo_coh)) %>% 
  names()
LP08_elements

# Ar (tube gas), Ta & W (splutter)
machine_elements <- select(LP08_COMP_filter, c(Ar, Ta, W)) %>% 
  names()

# REE 
# REE <- select(LP08_COMP_filter, c(La:Ho)) %>% 
#  names()
# REE

# Machine generated elements Ar (tube gas), Ta & W (splutter), REE removed
# LP08_elements1 <- select(LP08_COMP_filter, 
#                        -c(Core:MSE, all_of(machine_elements), 
#                           Nb:Cs, REE, Ir, Pt)) %>% names()
# LP08_elements1

# Calculate as % of normalising factors TS (Total Scatter), cps_sum & inc/coh -------------------------------------------------------

# Add to columns - find a way to replace [9:67] with column headings as "Mg":"Mo_coh"
LP08_COMP_filter

# Type in row numbers of elements and scatter to calculate  across 
LP08_rowsums <- 9:43

LP08_COMP_filter1 <- LP08_COMP_filter %>% 
  replace(is.na(.), 0) %>%
  mutate(TS_sum = Mo_inc + Mo_coh) %>% 
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  mutate(coh_inc = Mo_coh / Mo_inc) %>%
  mutate(cps_sum = rowSums(.[LP08_rowsums]))
#mutate(cps_sum = row_sums(LP08_elements))

LP08_COMP_filter1

# create clr file from ratio data using chemometrics package - this doesnt seem right - only works for elements where not a lot of zeroes
library(chemometrics)
LP08_COMP_filter1_text <- select(LP08_COMP_filter1, Core:MSE)
LP08_COMP_filter1_text

LP08_COMP_filter1_clr1 <- select(LP08_COMP_filter1, Fe) %>% 
  clr%>% 
  as_tibble()
LP08_COMP_filter1_clr1

LP08_COMP_filter1_clr_Fe <- bind_cols(LP08_COMP_filter1_text, LP08_COMP_filter1_clr1)
LP08_COMP_filter1_clr_fe

# Standardise and centre cps data  - using column names --------------------------------------

# list of element column names for plotting
LP08_col <- select(LP08_COMP_filter1, c(Mg:inc_coh)) %>% names()
LP08_col

# list of all element and other parameter column names 
LP08_col1 <- select(LP08_COMP_filter1, c(kcps:cps_sum)) %>% names()
LP08_col1

LP08_COMP_filter1.Z <- LP08_COMP_filter1
LP08_COMP_filter1.Z[, LP08_col] <- scale(LP08_COMP_filter1[, LP08_col], center = TRUE, scale = TRUE)
LP08_COMP_filter1.Z

# Convert cps and cps Z-scores to long format for plotting -------------------------------------------------

LP08_COMP_filter1_long <- select(LP08_COMP_filter1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
LP08_COMP_filter1_long 

LP08_COMP_filter1_long.Z <- select(LP08_COMP_filter1.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
LP08_COMP_filter1_long.Z 

# Generate cps summary stats  --------------------------------------------
library(psych)

LP08_summary <- LP08_COMP_filter1 %>%
  select(all_of(LP08_col1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter cps dataset and stats to file --------------------------------------------------
write.csv(LP08_COMP_filter1,"Output/LP08/2_cps_filter1.csv", row.names = FALSE)
write.csv(LP08_summary,"Output/LP08/2.1_cps_filter1_stats.csv", row.names = FALSE)
write.csv(LP08_COMP_filter1.Z,"Output/LP08/3_cps_filter1_Z.csv", row.names = FALSE)

# Filter cps element column list based on mean and max cps values -------------------

LP08_mean50 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 50)
}
LP08_mean200 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 200)
}
LP08_max50 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 50)
}

LP08_m50 <- select(LP08_COMP_filter1, 
                   Core, depth_cm, SH20_age, kcps, MSE, 
                   Mg:Mo_coh & where(LP08_mean50)) %>% print()
LP08_m200 <- select(LP08_COMP_filter1, 
                    Core, depth_cm, SH20_age, kcps, MSE, 
                    Mg:Mo_coh & where(LP08_mean200)) %>% print()
LP08_mx50 <- select(LP08_COMP_filter1, 
                    Core, depth_cm, SH20_age, kcps, MSE, 
                    Mg:Mo_coh & where(LP08_max50), TS_sum:cps_sum) %>% print()

# Create lists of column names to take forwLP08
LP08_col_m50 <- select(LP08_COMP_filter1, Mg:Mo_coh 
                       & -machine_elements 
                       & where(LP08_mean50)) %>% names() %>% print()
LP08_col_m200 <- select(LP08_COMP_filter1, Mg:Mo_coh 
                        & -machine_elements 
                        & where(LP08_mean200)) %>% names() %>% print()
LP08_col_m50_mx50 <- select(LP08_COMP_filter1, Mg:Mo_coh 
                            & -machine_elements 
                            & where(LP08_max50) 
                            & where(LP08_mean50)) %>% names() %>% print()

# Correlation matrices -------------------------------------------------------

library(GGally)
library(dplyr)
# Correlation plot with max200 elements - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=2))

ggcorr(LP08_COMP_filter1[, LP08_col1], method = c("everything", "pearson"), 
       size = 2, label = FALSE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/ITRAX/LP08/Fig 0_Corr_matrix_cps_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Manual list of elements with significant correlations - determined visually from correlation matrix
plot_elements1_LP08 <-c("Si", "P", "S", "K", "Ca", "Ti", "Cr", "Mn", "Fe", "Co", "Cu", "Zn", "Ga", "As", "Se",
                        "Br", "Rb", "Sr", "Zr", "Ba", "Pb", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")

# Stats generated list of elements to take forwLP08 based on m50
plot_elements2_LP08 <- c(LP08_col_m50, "inc_coh", "coh_inc")

ggcorr(LP08_COMP_filter1[,  plot_elements2_LP08], method = c("everything", "pearson"), 
       size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/ITRAX/LP08/Fig 0_Corr_matrix_cps_m50.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
theme_set(theme_bw(base_size=8))
ggpairs(LP08_COMP_filter1, columns = plot_elements2_LP08, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Figures/ITRAX/LP08/Fig 0__Corr-den_matrix_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary cps plots vs depth --------------------------------------------------------

library(tidypaleo)

# Figure 1 - cps elements1
theme_set(theme_bw(8))
LP08_Fig1 <- LP08_COMP_filter1_long  %>%
  filter(param %in% plot_elements1_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements1_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Signifcantly correlated elements cps summary")
LP08_Fig1
ggsave("Figures/ITRAX/LP08/Fig 1_cps_elements1.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 2 - cps elements2 based on mean>50 cps stats
theme_set(theme_bw(8))
LP08_Fig2 <- LP08_COMP_filter1_long  %>%
  filter(param %in% plot_elements2_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements2_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps based on mean>50 cps stats")
LP08_Fig2
ggsave("Figures/ITRAX/LP08/Fig 2_cps_elements2.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 3 - cps elements2 Z-scores
theme_set(theme_bw(8))
LP08_Fig3 <- LP08_COMP_filter1_long.Z  %>%
  filter(param %in% plot_elements2_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements2_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps Z-scores based on mean>50 cps stats")
LP08_Fig3
ggsave("Figures/ITRAX/LP08/Fig 3_cps_elements2_Z.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Normalization (elements and scatter parameters) ---------------------------------------------------------

# Inc normalised
LP08_inc_norm <- LP08_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Mo_inc))
# coh/inc normalised - Boyle 2015
LP08_coh_inc_norm <- LP08_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ (Mo_coh/Mo_inc)))
# Total scatter (TS) normalised - Kylander et al 2011/2012
LP08_TS_norm <- LP08_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ TS_sum))
# cps_sum normalised
LP08_cps_sum_norm <- LP08_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ cps_sum))
# Ti normalised
LP08_Ti_norm <- LP08_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Ti))
# Natural log Ti normalized & replace -Inf with NA --------
LP08_Ln_Ti_norm <- LP08_Ti_norm  %>% 
  mutate(across(c("Mg":"Mo_coh"),log)) 
is.na(LP08_Ln_Ti_norm)<-sapply(LP08_Ln_Ti_norm, is.infinite)
# Replace NA with 0 - if needed
# LP08_Ln_Ti_norm[is.na(LP08_Ln_Ti_norm)]<-0

#Check = 1
LP08_inc_norm$Mo_inc
LP08_Ti_norm$Ti
#Check = 0
LP08_Ln_Ti_norm$Ti

# StandLP08ize and center Ti-normalized data  --------------------------------------
LP08_Ln_Ti_norm.Z <- LP08_Ln_Ti_norm
LP08_Ln_Ti_norm.Z[, LP08_col1] <- scale(LP08_Ln_Ti_norm[, LP08_col1], center = TRUE, scale = TRUE)
LP08_Ln_Ti_norm.Z
# can replace Ti normalized with other normalisation parameters above

# write LP08_Ln_Ti_norm.Z to file for Part 3 & 4
write.csv(LP08_Ln_Ti_norm.Z,"Output/LP08/LP08_Ln_Ti_norm.Z.csv", row.names = FALSE)

# Elements and scatter as %cps sum - Bertrand et al 2021 -------- produces the same result as %TSN sum

LP08_cps_sum_norm_pc <- LP08_cps_sum_norm %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(cps_pc_sum = rowSums(.[LP08_rowsums]))

#Check everything = 100
head(LP08_cps_sum_norm_pc$cps_pc_sum)

# define new element list to include %cps sum (check=100 in output)
LP08_col3 <- select(LP08_cps_sum_norm_pc, c(Mg:cps_pc_sum)) %>% 
  names()
LP08_col3

# Create %cps_sum summary stats table
LP08_cps_sum_pc_summary <- LP08_cps_sum_norm_pc %>%
  select(all_of(LP08_col3)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Elements and scatter as percentage of TSN sum - Roberts et al 2017 --------- produces the same result as %cps sum
LP08_TSN_sum <- LP08_TS_norm %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_sum = rowSums(.[LP08_rowsums]))
LP08_TSN_sum

LP08_TSN_pc <- LP08_TSN_sum %>% mutate(across(c("Mg":"Mo_coh"),
                                              .fns = ~./TSN_sum)) %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum = rowSums(.[LP08_rowsums]))
LP08_TSN_pc

#Check everything = 100
head(LP08_TSN_pc$TSN_pc_sum)

# define new element list to include %TSN sum (check=100 in output)
LP08_col4 <- select(LP08_TSN_pc, c(Mg:TSN_pc_sum)) %>% 
  names()
LP08_col4

# Create TSN stats table
LP08_TSN_pc_summary <- LP08_TSN_pc %>%
  select(all_of(LP08_col4)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Filter elements based on %cps_sum (%TSN) max > 0.5 and mean >0.1 - *** replace %TSN with %cps sum here *** ----------------
LP08_max0.5 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 0.5)
}
LP08_mean0.1 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 0.1)
}

LP08_mx0.5 <- select(LP08_TSN_pc, Mg:Mo_coh & where(LP08_max0.5)) %>% print()
LP08_m0.1 <- select(LP08_TSN_pc, Mg:Mo_coh & where(LP08_mean0.1)) %>% print()

# define new element list of filtered elements and their %TSN sum (check=100 in output)
LP08_TSN_mx0.5 <- select(LP08_TSN_pc, Mg:Mo_coh 
                         & -machine_elements 
                         & where(LP08_max0.5)) %>% names()

LP08_TSN_m0.1 <- select(LP08_TSN_pc, Mg:Mo_coh 
                        & -machine_elements 
                        & where(LP08_mean0.1)) %>% names()
LP08_TSN_mx0.5
LP08_TSN_m0.1

# New element list based on  LP08_col_m0.1 + Si, P , S
LP08_TSN_m0.1_list <- c("Si", "P", "S", LP08_TSN_m0.1)
LP08_TSN_m0.1_list

# Add filtered element TSN sum to TSN_pc dataset & rename
LP08_TSN_pc1 <- LP08_TSN_pc %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum1 = rowSums(.[LP08_TSN_m0.1_list]))
LP08_TSN_pc1
# check filtered element total <100 but >95%
head(LP08_TSN_pc1$TSN_pc_sum1)

# Create TSN stats table with filtered element %TSN Sum added to end 
LP08_col5 <- select(LP08_TSN_pc1, Mg:TSN_pc_sum1) %>% names()

LP08_TSN_pc1_summary <- LP08_TSN_pc1 %>%
  select(all_of(LP08_col5)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write normalised datasets, %cps & %TSN and stats summaries to file --------------------------------------------------
write.csv(LP08_Ti_norm,"Output/LP08/4_Ti_norm.csv", row.names = FALSE)
write.csv(LP08_inc_norm,"Output/LP08/5_inc_norm.csv", row.names = FALSE)
write.csv(LP08_TS_norm,"Output/LP08/6_TS_norm.csv", row.names = FALSE)
write.csv(LP08_cps_sum_norm,"Output/LP08/7_cps_sum_norm.csv", row.names = FALSE)
write.csv(LP08_TSN_pc,"Output/LP08/8_TSN_pc.csv", row.names = FALSE)
write.csv(LP08_TSN_pc_summary,"Output/LP08/8.1_TSN_pc_stats.csv", row.names = FALSE)
write.csv(LP08_TSN_pc1,"Output/LP08/8.2_TSN_pc1.csv", row.names = FALSE)
write.csv(LP08_TSN_pc1_summary,"Output/LP08/8.3_TSN_pc1_stats.csv", row.names = FALSE)
write.csv(LP08_cps_sum_norm_pc,"Output/LP08/9_cps_sum_pc.csv", row.names = FALSE)
write.csv(LP08_cps_sum_pc_summary,"Output/LP08/9.1_cps_sum_pc_stats.csv", row.names = FALSE)
write.csv(LP08_Ln_Ti_norm,"Output/LP08/10_Ln_Ti_norm.csv", row.names = FALSE)
write.csv(LP08_Ln_Ti_norm.Z,"Output/LP08/11_Ln_Ti_norm_Z.csv", row.names = FALSE)

# Convert to long format ------------------------------------------------
LP08_inc_norm_long <- select(LP08_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
LP08_inc_norm_long

LP08_coh_inc_norm_long <- select(LP08_coh_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
LP08_coh_inc_norm_long

LP08_Ln_Ti_norm_long <- select(LP08_Ln_Ti_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
LP08_Ln_Ti_norm_long

LP08_Ln_Ti_norm.Z_long <- select(LP08_Ln_Ti_norm.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col1)) %>%
  pivot_longer(all_of(LP08_col1), names_to = "param", values_to = "value")
LP08_Ln_Ti_norm.Z_long

LP08_cps_sum_norm_pc_long <- select(LP08_cps_sum_norm_pc,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col3)) %>%
  pivot_longer(all_of(LP08_col3), names_to = "param", values_to = "value")
LP08_cps_sum_norm_pc_long

LP08_TSN_pc1_long <- select(LP08_TSN_pc1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(LP08_col5)) %>%
  pivot_longer(all_of(LP08_col5), names_to = "param", values_to = "value")
LP08_TSN_pc1_long

# Summary normalised plots--------------------------------------

# # Figure 4 cps/inc. -----------------------------------------------------

# manually user defined
plot_elements3_LP08 <- c("S", "K", "Ca", "Ti", "Mn", "Fe", "Zn", 
                         "Br", "Sr", "Zr", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")
# OR 

# stats defined based on mean %TSN >0.1% - %TSN output is the same as %cps sum output
plot_elements3_LP08 <- c(LP08_TSN_m0.1_list, "inc_coh", "coh_inc")

# can replace %TSN_M0.1 element list with cps LP08_col_m50 element list
# plot_elements3_LP08 <- c(LP08_col_m50, "inc_coh", "coh_inc")

# Unit depths LP08 & LP16
unit_depths_LP08 <- c(45, 175, 324, 378, 470)
unit_depths_LP16 <- c(35, 47, 67, 74, 83, 110, 295) 

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
theme_set(theme_bw(8))
LP08_Fig4a <- LP08_inc_norm_long  %>%
  filter(param %in% plot_elements3_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements3_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc.", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Incoherent (inc.) scatter normalised (filtered elements): CONISS XRF")

LP08_Fig4b <- LP08_coh_inc_norm_long  %>%
  filter(param %in% plot_elements3_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements3_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/(coh./inc.)", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Coherent/Incoherent (coh./inc.) scatter ratio normalised (filtered elements): CONISS XRF")
ggarrange(LP08_Fig4a, LP08_Fig4b, nrow = 2)
ggsave("Figures/ITRAX/LP08/Fig 4_cps_inc&coh_inc.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 5 %TSN----------------------------------------------------------------

plot_elements4_LP08 <- c("S", "K", "Ca", "Ti", "Mn", "Fe", "Zn", 
                         "Br", "Sr", "Zr", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")
# OR

# stats defined 
plot_elements4_LP08 <- c(LP08_TSN_m0.1_list, "inc_coh", "coh_inc")

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
# Roberts et al 2017 - CONISS with broken stick Hellingers dist defined zones - red dashed lines
theme_set(theme_bw(8))
LP08_Fig5 <- LP08_TSN_pc1_long  %>%
  filter(param %in% plot_elements4_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements4_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%TSN", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%TSN (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss1_LP08 <- LP08_TSN_pc1_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 6 - %TSN with CONISS
LP08_Fig6 <- LP08_Fig5 +
  layer_dendrogram(coniss1_LP08, aes(y = depth_cm), param = "CONISS") +
  layer_zone_boundaries(coniss1_LP08, aes(y = depth_cm, col = "blue", lty = 2, alpha = 0.7)) +
  ggtitle("%TSN (filtered elements): CONISS XRF(red) ITRAX (black)")

ggarrange(LP08_Fig5, LP08_Fig6, nrow = 2)
ggsave("Figures/ITRAX/LP08/Fig 5_6_TSN_pc_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

#Figure 7 - TSN with CONISS and Roberts et al (2017) zone boundaries
#LP08_Fig6 +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%TSN (filtered elements): Coniss zone comparison")
#ggsave("Figures/ITRAX/LP08/Fig 7_TSN_pc_CONISS_comparison.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 8  %cps_sum
#define colour scheme for 6 groups (5 guano and 1 non-guano)
guano_zone_colours <- c("#BDBDBD", "#00441B", "#00441B", "#006D2C", "#238B45", "#74C476")
theme_set(theme_bw(8))
LP08_Fig8 <- LP08_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements3_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements3_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm, aes(colour = guano_zone_colours))) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%cps_sum", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%cps sum (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss2_LP08 <- LP08_cps_sum_norm_pc_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 9 %cps_sum with CONISS
LP08_Fig9 <- LP08_Fig8 +
  layer_dendrogram(coniss2_LP08, aes(y = depth_cm), param = "CONISS") +
  layer_zone_boundaries(coniss2_LP08, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
  ggtitle("%cps sum (filtered elements): CONISS XRF(red) ITRAX (black)")
ggarrange(LP08_Fig8, LP08_Fig9, nrow = 2)
ggsave("Figures/ITRAX/LP08/Fig 8_9_cps_sum_pc_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 10 - compare to Roberts et al. 2017 CONISS zoning 
#LP08_Fig9 +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%cps sum (filtered elements): CONISS comparison")
#ggsave("Figures/LP08/Fig 10_cps_sum_pc_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")


# create a new element list without Ti
plot_elements5_LP08 <- purrr::discLP08(plot_elements4_LP08,.p = ~stringr::str_detect(.x,"Ti"))
plot_elements5_LP08

# Figure 11 Ti normalised
theme_set(theme_bw(8))
LP08_Fig11 <- LP08_Ln_Ti_norm_long  %>%
  filter(param %in% plot_elements5_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements5_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Ln(cps/Ti)", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised (filtered elements): CONISS XRF")
LP08_Fig11
ggsave("Figures/ITRAX/LP08/Fig 11_Ti_norm.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 12 - Ti normalised as Z-scores 
theme_set(theme_bw(8))
LP08_Fig12 <- LP08_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements5_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements5_LP08)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "LN(cps/Ti) Z-scores", y = "Depth (cm)") +
  geom_hline(yintercept = unit_depths_LP08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss3_LP08 <- LP08_Ln_Ti_norm.Z_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 13
LP08_Fig13 <- LP08_Fig12 + labs(x = "Ln(cps/Ti Z-scores)", y = "Depth (cm)") +
  layer_dendrogram(coniss3_LP08, aes(y = depth_cm), param = "CONISS ITRAX") +
  layer_zone_boundaries(coniss3_LP08, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF(red) ITRAX (black)")
ggarrange(LP08_Fig12, LP08_Fig13, nrow = 2)
ggsave("Figures/ITRAX/LP08/Fig 12_13_Ti_Z_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 14
#LP08_Fig12 + labs(x = "cps/Ti Z-scores", y = "Depth (cm)") +
#  layer_dendrogram(coniss3_LP08, aes(y = depth_cm), param = "CONISS") +
#  layer_zone_boundaries(coniss3_LP08, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS comparison")
#ggsave("Figures/LP08/Fig 14_Ti_Z_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Plots vs age ------------------------------------------------------------

unit_ages_LP08 <- c(1380, 3310, 5600, 7370, 10160, 21230, 29800)


# Figure 15 %cps sum vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
LP08_Fig15 <- LP08_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements4_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements4_LP08)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "%cps sum") +
  geom_vline(xintercept = unit_ages_Lp08, colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log %cps sum (filtered elements): CONISS XRF")
ggsave("Figures/ITRAX/LP08/Fig 15_cps_sum_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Figure 16 cps/Ti as Z-scores vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
LP08_Fig16 <- LP08_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements4_LP08) %>%
  mutate(param = fct_relevel(param, plot_elements4_LP08)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps/Ti Z-score") +
  geom_vline(xintercept = unit_ages_LP08, colour = "red", lty = 2, alpha = 0.7)+
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")
ggarrange(LP08_Fig15, LP08_Fig16, nrow = 2)
ggsave("Figures/LP08/Fig 16_Ti_Z_CONISS_2017_age.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")


# END ---------------------------------------------------------------------


