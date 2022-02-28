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
library(GGally)
library(bestNormalize)
library(sjmisc)
library(chemometrics)
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)

# Set working directory ---------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021")
#check working directory
getwd()

# PART 1 Matched 4 cm interval subsample and ITRAX 1cm  dataset ------------------------------------------------------------------


# Set up variable lists  ---------------------------------------------

# elements are %TSN

all_variables <- c("GRD", "DMAR", "WC", "DW", "LOI550", "LOI950_c", "TOC_N", "TOC_B", 
                   "Ca", "Ti","Mn", "Fe","Br", "Sr", "inc", "coh",
                   "LnBr_Ti",	"LnSr_Ti",	"LnMn_Ti",	"LnCa_Ti",	"LnFe_Mn", "inc_coh", "coh_inc", "TSN_RS")

key_variables <- c("GRD", "DW", "TOC_N", "LOI950_c", "Ti",  "Br", "Ca", "inc_coh", "coh_inc", "TSN_RS")

TSN_variables <- c( "Ca", "Ti","Mn", "Fe","Br", "Sr", "inc", "coh")

Ln_variables<- c("LnBr_Ti",	"LnSr_Ti",	"LnMn_Ti",	"LnCa_Ti",	"LnFe_Mn")


# Import data  ----------------------------------------------------------

# data import
db <- read_csv("Data/LP08_Sedimentology_4cm.csv") %>% 
  filter(!if_any(everything(), is.na))
db
  
db1 <-  select(db, c(Strat_depth,Unit, Subunit, SH20_mean_age, all_variables))
db1

# Create log ratio, centered and centred, standardised (z-scores) log ratios for all data --------------------------

db_ln <- db1
db_ln[, TSN_variables]<- log(db_ln[TSN_variables]) %>% 
  as_tibble()
db_ln

db_cln <- db1
db_cln[, TSN_variables] <- scale(db_cln[TSN_variables], center = TRUE) %>% 
  as_tibble()
db_cln

db_cln.Z <- db_ln
db_cln.Z[, TSN_variables] <- scale(db_cln.Z[TSN_variables], center = TRUE, scale = TRUE) %>% 
  as_tibble()
db_cln.Z

# Note1: centred log ratio oxide and elemental data should be the same when centered
# Note2: centered and standardised (scaled) produces the same output as centering alone for log ratio data

#  Write data without blanks to Output folder --------------------------------------------------
write.csv(db,"Output/Correlation/LP08_4cm/db.csv", row.names = FALSE)
write.csv(db1,"Output/Correlation/LP08_4cm/db1.csv", row.names = FALSE)
write.csv(db_ln,"Output/Correlation/LP08_4cm/ln.csv", row.names = FALSE)
write.csv(db_cln,"Output/Correlation/LP08_4cm/cln.csv", row.names = FALSE)
write.csv(db_cln.Z,"Output/Correlation/LP08_4cm/cln_Z.csv", row.names = FALSE)


# Generate stats & write to file  --------------------------------------------
library(psych)

db1_summary <- db1 %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

dbcln_summary <- db_cln %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

db1_summary_key <- db1 %>%
  select(all_of(key_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

dbcln_summary_key <- db_cln %>%
  select(all_of(key_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter to file --------------------------------------------------
write.csv(db1_summary,"Output/Correlation/LP08_4cm/stats_all.csv", row.names = FALSE)
write.csv(dbcln_summary,"Output/Correlation/LP08_4cm/stats_cln_all.csv", row.names = FALSE)
write.csv(db1_summary_key,"Output/Correlation/LP08_4cm/stats_key.csv", row.names = FALSE)
write.csv(dbcln_summary_key,"Output/Correlation/LP08_4cm/stats_cln_key.csv", row.names = FALSE)


# Convert to long format --------------------------------

# Convert to long format - i.e., one value per row per variable - for later
db1_long <- select(db1, Strat_depth, Unit, Subunit, SH20_mean_age, all_of(all_variables)) %>%
  pivot_longer(all_of(all_variables), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
db1_long

# Convert to long format - i.e., one value per row per variable - for later
db1_long_key <- select(db1, Strat_depth, Unit, Subunit, SH20_mean_age, all_of(key_variables)) %>%
  pivot_longer(all_of(key_variables), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
db1_long_key


# Correlation and covariance matrices -------------------------------------

# Examine correlation matrix & p-values for db1 all variables 
db1_var <- select(db1, all_variables)
db1_var
cor <- cor(db1_var)
round(cor, 2)

# p-values
library("Hmisc")
cor_p <- rcorr(as.matrix(db1_var))
cor_p

# Examine correlation matrix & p-values for db1 key variables 
db1_key <- select(db1, key_variables)
db1_key
cor_key <- cor(db1_key)
round(cor_key, 2)

# p-values
library("Hmisc")
cor_p_key <- rcorr(as.matrix(db1_key))
cor_p_key

# Examine co-variance of %TSN elements only (only data measured in same way, with the same units)
db1_TSN <- select(db1, TSN_variables)
cov_key <- cov(db1_TSN)
round(cov_key,2)


# Correlation summary for all data ---------------
library(GGally)

# Correlation plot - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=8))
ggcorr(db1[, all_variables], method = c("everything", "pearson"),
       size = 3, label = TRUE, label_size = 3, label_alpha = TRUE, label_round=2) 
ggsave("Output/Correlation/LP08_4cm/Fig 1A_Corr_matrix_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
# Stars meaning = *** p<0.001; ** p<0.001; * p<0.05; . p<0.1; nothing = p>0.1
theme_set(theme_bw(base_size=8))
ggpairs(db1, columns = all_variables, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Output/Correlation/LP08_4cm/Fig 1B_Corr-den_matrix_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix - by Unit
theme_set(theme_bw(base_size=8))
ggpairs(db1, columns = all_variables, upper = list(continuous = wrap("cor", size = 1.5)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.5)),
        ggplot2::aes(colour = Unit, title="Correlation plot by Unit", alpha = 0.5))
ggsave("Output/Correlation/LP08_4cm/Fig 1C_Corr-den-unit_matrix_all.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation summary for key correlations ---------------

# Correlation plot - overview
theme_set(theme_bw(base_size=12))
ggcorr(db1[, key_variables], method = c("everything", "pearson"),
       size = 6, label = TRUE, label_size = 6, label_alpha = TRUE, label_round=2) 
ggsave("Output/Correlation/LP08_4cm/Fig 1A_Corr_matrix_key.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
# Stars meaning = *** p<0.001; ** p<0.001; * p<0.05; . p<0.1; nothing = p>0.1
theme_set(theme_bw(base_size=12))
ggpairs(db1, columns = key_variables, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1)),
        title="Correlation-density plot")
ggsave("Output/Correlation/LP08_4cm/Fig 1B_Corr-den_matrix_key.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix - by Unit
theme_set(theme_bw(base_size=12))
ggpairs(db1, columns = key_variables, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1)),
        ggplot2::aes(colour = Unit, title="Correlation plot by Unit", alpha = 0.5))
ggsave("Output/Correlation/LP08_4cm/Fig 1C_Corr-den-unit_matrix_key.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation stats for key correlations ---------------

# Get full stats for some key correlation tests (if needed)
res_1 <- cor.test(db1$TOC_N, db1$inc_coh, method = "pearson")
print(res_1)
res_2 <- cor.test(db1$TOC_N, db1$Br, method = "pearson")
print(res_2)
res_3 <- cor.test(db1$TSN_RS, db1$coh_inc, method = "pearson")
print(res_3)
res_4 <- cor.test(db1$DW, db1$coh_inc, method = "pearson")
print(res_4)
res_5 <- cor.test(db1$LOI950_c, db1$Ca, method = "pearson")
print(res_5)

# Save console outputs above to file
sink("Output/Correlation/LP08_4cm/LP08_4cm_stats_key.txt")
print(res_1)
print(res_2)
print(res_3)
print(res_4)
print(res_5)
sink(file = NULL)
# reset back to console output
sink(file = NULL)
res_1
