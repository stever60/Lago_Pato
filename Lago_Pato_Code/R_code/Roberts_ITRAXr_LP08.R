
# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# setup workspace ----
library(itraxR)
library(tidyverse) # all core tidyverse packages
library(tidypaleo) # Dewey Dunnington's ggplot extensions for palaeo-style plots
library(readr)
library(ggpubr)

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021/Data/ITRAX/LP08/")
#check working directory
getwd() 

# START NOTES ----------------------------------------------------------

# contains extra lines D1, S1, S2, S3; D1 = Fe a*2 = 12.8 keV used in itrax.R example file

# Renaming columns in  old filenames to match ITRAX.r formats before using itraximport  ----------------------------------

# result.txt = original txt output file 
# Results = kcps renamed as cps
# Results1 = kcps & cps retained; cps calculated as element and scatter sum by code below - use this one 
# Results2 = cps calculated as element and scatter sum by code below & kcps removed - don't need this one
# need to run it for each core and change folder name using find and replace 
# Reanalysis data contains extra column callled revalidity; validity is the original run validity

# import and rename cps as kcps - saved as Results.tsv
df<- read_tsv("LP081A/result.txt", col_names = TRUE, skip = 2) %>% 
  rename(cps = kcps, `Fe a*2` = D1) %>% 
  write_tsv("LP081A/Results.tsv")

# import, calculate cps_sum and rename as cps, retain kcps column - saved as Results1.tsv
df1_rowsums <- c(12:51, 55, 56)
df1 <- read_tsv("LP081A/result.txt", col_names = TRUE, skip = 2) %>%
  mutate(cps_sum = rowSums(.[df1_rowsums])) %>% 
  rename(cps = cps_sum, `Fe a*2` = D1) %>% 
  relocate(cps, .before = MSE) %>% 
  write_tsv("LP081A/Results1.tsv")

# list of all elements in LP08 and LP16 files
cps_elementsList_LP08 <- select(df1, c(Mg:Pb, `Mo inc`, `Mo coh`)) %>% 
  names()
cps_elementsList_LP08

# Now do this in Finder -----------------------------------------------

# open tsv Results and Results1 in excel 
# add first two row back into all new txt files using copy/paste
# rename .tsv files to .txt in Finder
# now ready to import into itrax.R and runs as normal
# choose whether to use Results.txt or Results1.txt
# below uses Results1.txt
# repeat for all sections 

# SECTION 1 - DATA STRUCTURE & TESTING FILE IMPORTS WORK OK ----------------------------------------------

# import the data - Results file has cps instead of kcps so that ITRAX.r will run - original file is in 2008 folder 
itrax_import("LP081A/Results.txt", depth_top = 0) %>% # import xrf data
  ggplot(aes(x=`Mo inc`/`Mo coh`, y=depth)) + # setup plot area
  geom_point() + # plots points             
  geom_lineh() + # plots the line
  scale_y_reverse() + # reverses the scale
  theme_bw()

# look at the metadata - need to add extra data for earlier Aber run cores to match Manchester metadata - see template 
itrax_meta("LP081A/document.txt")
itrax_meta("LP081B/document.txt")

# plot the optical image vs position - 
# optical1 is Photoshopenhanced image in Photoshop - HDR - photorealistic/low contrast, brightness+100, contrast+50
itrax_image(file = "LP081B/optical1.tif", # define location of image file
            meta = "LP081B/document.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
)

# plot the radiograph vs position - radiograph1 is Photshop enhanced image - HDR - photorealistic, brightness+100, contrast+50
itrax_radiograph(file = "LP081B/radiograph1.tif", # define location of radiograph image
                 meta = "LP081B/document.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure


# plot the optical image vs position - check LP081B
# LP081B was compressed and then realigned using LOI, MS and XRF subsample data
# LP081B position adjusted intervals of 3.07692307692307 i.e., 65 cm should be 100 cm long core (100/65 x 2 mm)
# document start and end cood adjusted to 0 and 1000 - optical step size adjusted upwLP08s by *100/65 multiple too
# optical1 is Photoshop enhanced image in Photoshop - HDR - photorealistic/low contrast, brightness+100, contrast+50
itrax_image(file = "LP081B/optical1.tif", # define location of image file
            meta = "LP081B/document.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
)

# plot the radiograph vs position - radiograph1 is Photoshop enhanced image - HDR - photorealistic, brightness+100, contrast+50
itrax_radiograph(file = "LP081B/radiograph1.tif", # define location of radiograph image
                 meta = "LP081B/document.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

itrax_qspecsettings("LP081A/settings.dfl") # parse some q-spec settings

itrax_spectra(filename = "LP081A/sumspectra.spe", # define a raw spectra file or sumspectra file
              parameters = "LP081A/settings.dfl", # define an associated settings file
              plot = TRUE # suppress the plot
) 


# SECTION 2 - IMPORTING DATA ------------------------------------------------

LP081A_S1 <- list(metadata   = itrax_meta("LP081A/document.txt"),
                 xrf        = itrax_import("LP081A/Results1.txt", 
                                           depth = 86, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP081A/optical1.tif",
                                          meta = "LP081A/document.txt"),
                 radiograph = itrax_radiograph(file = "LP081A//radiograph1.tif",
                                               meta = "LP081A/document.txt",
                                               trim = as.numeric(itrax_meta("LP081A/document.txt")[6:7,2])))

# LP081B using document and Results files - as measured, start depth adjusted to 900 mm 
# LP081B_1 using document1 and Results1 - start depth at 650 mm retained - decompressed and then realigned based on LOI, MS and XRF subsample data
# LP081B Results1 position data has been adjusted to intervals of 3.07692307692307 i.e., 65 cm section should be 100 cm long core (100/65 x 2 mm)

# LP081B using Results.txt and document.txt
LP081B_S2 <- list(metadata   = itrax_meta("LP081B/document.txt"),
                 xrf        = itrax_import("LP081B/Results1.txt", 
                                           depth = 600, #650 - field depth
                                           parameters = "all"),
                 image      = itrax_image(file = "LP081B/optical1.tif",
                                          meta = "LP081B/document.txt"),
                 radiograph = itrax_radiograph(file = "LP081B//radiograph1.tif",
                                               meta = "LP081B/document.txt",
                                               trim = as.numeric(itrax_meta("LP081B/document.txt")[6:7,2])))

# LP081B using Results1.txt and document1.txt - as in Roberts et al. (2017) based on age/XRF info
# core decompresssed/stretched from 65 cm as measured length to 100 cm 
#LP081B_S2 <- list(metadata   = itrax_meta("LP081B/document1.txt"),
#              xrf        = itrax_import("LP081B/Results1.txt", 
#                                        depth = 544, 
#                                        parameters = "all"),
#              image      = itrax_image(file = "LP081B/optical1.tif",
#                                       meta = "LP081B/document.txt"),
#              radiograph = itrax_radiograph(file = "LP081B//radiograph1.tif",
#                                            meta = "LP081B/document.txt",
#                                            trim = as.numeric(itrax_meta("LP081B/document1.txt")[6:7,2])))


LP081C_S3 <- list(metadata   = itrax_meta("LP081C/document.txt"),
                 xrf        = itrax_import("LP081C/Results1.txt", 
                                           depth = 1098, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP081C/optical1.tif",
                                          meta = "LP081C/document.txt"),
                 radiograph = itrax_radiograph(file = "LP081C//radiograph1.tif",
                                               meta = "LP081C/document.txt",
                                               trim = as.numeric(itrax_meta("LP081C/document.txt")[6:7,2])))

LP081D_S4 <- list(metadata   = itrax_meta("LP081D/document.txt"),
                 xrf        = itrax_import("LP081D/Results1.txt", 
                                           depth = 1888, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP081D/optical1.tif",
                                          meta = "LP081D/document.txt"),
                 radiograph = itrax_radiograph(file = "LP081D//radiograph1.tif",
                                               meta = "LP081D/document.txt",
                                               trim = as.numeric(itrax_meta("LP081D/document.txt")[6:7,2])))

LP081E_S5 <- list(metadata   = itrax_meta("LP081E/document.txt"),
                 xrf        = itrax_import("LP081E/Results1.txt", 
                                           depth = 2748, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP081E/optical1.tif",
                                          meta = "LP081E/document.txt"),
                 radiograph = itrax_radiograph(file = "LP081E//radiograph1.tif",
                                               meta = "LP081E/document.txt",
                                               trim = as.numeric(itrax_meta("LP081E/document.txt")[6:7,2])))

LP081F_S5 <- list(metadata   = itrax_meta("LP081F/document.txt"),
                  xrf        = itrax_import("LP081F/Results1.txt", 
                                            depth = 2748, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP081F/optical1.tif",
                                           meta = "LP081F/document.txt"),
                  radiograph = itrax_radiograph(file = "LP081F//radiograph1.tif",
                                                meta = "LP081F/document.txt",
                                                trim = as.numeric(itrax_meta("LP081F/document.txt")[6:7,2])))

LP081G_S5 <- list(metadata   = itrax_meta("LP081G/document.txt"),
                  xrf        = itrax_import("LP081G/Results1.txt", 
                                            depth = 2748, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP081G/optical1.tif",
                                           meta = "LP081G/document.txt"),
                  radiograph = itrax_radiograph(file = "LP081G//radiograph1.tif",
                                                meta = "LP081G/document.txt",
                                                trim = as.numeric(itrax_meta("LP081G/document.txt")[6:7,2])))

LP081H_S5 <- list(metadata   = itrax_meta("LP081H/document.txt"),
                  xrf        = itrax_import("LP081H/Results1.txt", 
                                            depth = 2748, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP081H/optical1.tif",
                                           meta = "LP081H/document.txt"),
                  radiograph = itrax_radiograph(file = "LP081H//radiograph1.tif",
                                                meta = "LP081H/document.txt",
                                                trim = as.numeric(itrax_meta("LP081H/document.txt")[6:7,2])))

LP081I_S5 <- list(metadata   = itrax_meta("LP081I/document.txt"),
                  xrf        = itrax_import("LP081I/Results1.txt", 
                                            depth = 2748, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP081I/optical1.tif",
                                           meta = "LP081I/document.txt"),
                  radiograph = itrax_radiograph(file = "LP081I//radiograph1.tif",
                                                meta = "LP081I/document.txt",
                                                trim = as.numeric(itrax_meta("LP081I/document.txt")[6:7,2])))


# join the xrf data for the sections together ---
LP08_xrf <- itrax_join(list(S1 = LP081A_S1$xrf, # S1 will be the "label" given to the core section
                           S2 = LP081B_S2$xrf, 
                           S3 = LP081C_S3$xrf,
                           S4 = LP081D_S4$xrf, 
                           S5 = LP081E_S5$xrf,
                           S6 = LP081F_S6$xrf, 
                           S7 = LP081G_S7$xrf,
                           S8 = LP081H_S8$xrf, 
                           S9 = LP081I_S9$xrf,)
)

LP08_xrf
write.csv(OUTPUT/LP08/LP08_xrf,"LP08_xrf.csv", row.names = FALSE)

# Figure 1 - Summary overlaps plot - all sections -------------------------

Fig1.1 <- ggplot(data = na.omit(LP08_xrf), mapping = aes(x = depth, y = Ti)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.2 <- ggplot(data = na.omit(LP08_xrf), mapping = aes(x = depth, y = Br)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.3 <- ggplot(data = na.omit(LP08_xrf), mapping = aes(x = depth, y = `Mo coh`/`Mo inc`)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
ggarrange(Fig1.1, Fig1.2, Fig1.3, ncol = 3)
ggsave("Figures/Fig1_overlaps.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")


# SECTION 3 - TIDYING DATA ----------------------------------------------------------

# cps filtering using Fe rather than  Fe a*2 & adjust cps to between 10,000-20,000
Fig2 <- ggplot(data = LP08_xrf, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
Fig2
ggsave("Figures/Fig2_tolerance_Fe_count.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# cps tolerance filter - could use >mean+/-2s cps (based on kcps or cps_sum)
cps.min.thres <- 30000
cps.max.thres <- 70000

#  OR 

cps.mean <- mean(LP08_xrf$cps)
cps.sd <- 2*sd(LP08_xrf$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

LP08_xrf  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw()
ggsave("Figures/Fig3_tolerance_cps.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# MSE tolerance filter - 2 used here but could use >mean+2s - which is 1.764766 for LP08 record 

MSE.thres <- 2 # use this for LP08

#  OR

MSE.mean <- mean(LP08_xrf$MSE)
MSE.sd <- 2*sd(LP08_xrf$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

LP08_xrf %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw()
ggsave("Figures/Fig4_tolerance_MSE.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Surface slope tolerance filter
slope.min.thres = -0.2
slope.max.thres = 0.2

#  OR - used for LP08
slope1 <-  LP08_xrf$`sample surface` - lag(LP08_xrf$`sample surface`)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 2*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

LP08_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_tolerance)) +
  scale_y_continuous(limits = c(-0.55, 0.55), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()
ggsave("Figures/Fig5_tolerance_sur_slope_.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Combining 'validity' flags   
LP08_xrf <- LP08_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) # %>% filter(qc == TRUE) #to remove from LP08_xrf rows that dont pass QC
# plot summary
theme_set(theme_bw(8))
Fig6.1 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.2 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.3 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.4 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig6.1, Fig6.2, Fig6.3, Fig6.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig6_tolerance_combined.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

LP08_xrf 

# Combining 'validity' flags - removing data that doesn't pass - replotting as line only
LP08_xrf1 <- LP08_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=-0.3 | slope >=0.3 | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=20000 | cps >=60000 | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=2, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) %>%  
  filter(qc == TRUE) #to remove from LP08_xrf rows that dont pass QC

# plot summary
theme_set(theme_bw(8))
Fig7.1 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.2 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.3 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.4 <- ggplot(data = LP08_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig7.1, Fig7.2, Fig7.3, Fig7.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig7_tolerance_filtered.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correcting for dead time & Uncertainties ---------------------------------------------- 

# This doesn't work on cores scanned at Aber - no Dt (Dwell time) column in results.txt file

# ggplot(data = LP08_xrf, aes(x = depth, y = Dt)) +
#  scale_x_reverse() + 
#  scale_y_continuous(sec.axis = sec_axis( trans=~(.+(1-mean(LP08_xrf$Dt, na.rm = TRUE))), name="Correction Factor")) +
#  geom_line() +
#  geom_hline(yintercept = mean(LP08_xrf$Dt, na.rm = TRUE), linetype = "dotted")


# Noisy data ------------------------------------------------------------

# see also section: # Calculate as % of normalising factors TS (Total Scatter), cps_sum - in Bertrand et al.R

# Using autocorrelation to detect noisy signals - non-noisy data should be AC, higher, outside 95% limits, showing some order/pattern
library(forecast)
library(ggpubr)

# Individual elements
Fig8 <- ggarrange(
  ggAcf(LP08_xrf$Ca) + ylim(c(NA,1)), ggAcf(LP08_xrf$Ti) + ylim(c(NA,1)), 
  ggAcf(LP08_xrf$Fe) + ylim(c(NA,1)),  ggAcf(LP08_xrf$Sr) + ylim(c(NA,1)),
  ggAcf(LP08_xrf$P) + ylim(c(NA,1)),ggAcf(LP08_xrf$Cu) + ylim(c(NA,1)), 
  ggAcf(LP08_xrf$Zn) + ylim(c(NA,1)), ggAcf(LP08_xrf$S) + ylim(c(NA,1)),
  ggAcf(LP08_xrf$Ni) + ylim(c(NA,1)), ggAcf(LP08_xrf$Cs) + ylim(c(NA,1)),
  ggAcf(LP08_xrf$Ba) + ylim(c(NA,1)), ggAcf(LP08_xrf$Pb) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
Fig8
ggsave("Figures/Fig8_ACF.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# all elements in summary plot
elementsList <- select(LP08_xrf, c(Mg:`Mo coh`)) %>% 
  names()
elementsList

apply(LP08_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == 5) %>% arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line()
ggsave("Figures/Fig9_ACF_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# identify acceptable variables
# use coh and inc as stop points for well measured: 0.5 takes down to Mo inc, 0.23 goes down to Mo coh 
apply(LP08_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == 5) %>% arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == 5) %>%
  filter(value >= 0.23) %>%
  pull(elements) %>% 
  ordered() -> myElements
myElements

# get acceptable rows and variables and make into long format and then plot
# add P for LP08 regression analysis as not in selected
# remove Ar, Ta, W (if selected above) - these are detector generated elements
LP08_xrf %>% 
  filter(qc == TRUE) %>% # pivot long
  select(P, any_of(myElements), depth, label) %>% 
  select(-c(Ar)) %>% 
  tidyr::pivot_longer(!c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc"))) %>%
  # plot
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig10_filtered_elements.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


## Visualising Raw Data



# SECTION 4 - PLOTTING ----------------------------------------------------

allelements <- 5:45

theme_set(theme_paleo(8))
xrfStrat <- LP08_xrf %>% 
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc/coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method
  # select(Fe, Ti, Mn,`coh/inc`, `inc/coh`,TS_sum, cps_sum, depth, label) %>% # a smaller set of elements defined manually to test.
  select(P, S, any_of(myElements), `coh/inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar)) %>% 
  tidyr::pivot_longer(!c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))  
# note that the levels controls the order

ggplot(xrfStrat, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  #tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig11_filtered_elements_plot.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Plots with tolerance filtered row removed but no colour
ggplot(xrfStrat, aes(x = peakarea, y = depth)) +
  geom_lineh() + #aes(color = label)
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  #tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig12_multi_tolerance_filtered_elements.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

xrfStrat

# ADDING CORE IMAGES - TO DO / FINISH  ------------------------------------------------------------------

Fig2 <- ggplot() +
  scale_y_continuous(limits = rev(range(LP08_xrf$depth))) +
  scale_x_continuous(breaks = round(c(0, max(as.numeric(colnames(LP081B_S2$image$image)))*2)),
                     limits = c(0, max(as.numeric(colnames(LP081B_S2$image$image)))*2)) +
  coord_fixed(ratio = 1) +
  labs(y = "Depth [mm]", x = "[mm]") +
  annotation_custom(rasterGrob(LP081A_S1$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP081A_S1$xrf$depth),
                    ymin = min(LP081A_S1$xrf$depth),
                    xmin = min(as.numeric(colnames(LP081A_S1$image$image))),
                    xmax = max(as.numeric(colnames(LP081A_S1$image$image)))
  ) +
  annotation_custom(rasterGrob(LP081B_S2$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP081B_S2$xrf$depth),
                    ymin = min(LP081B_S2$xrf$depth),
                    xmin = max(as.numeric(colnames(LP081B_S2$image$image))),
                    xmax = max(as.numeric(colnames(LP081B_S2$image$image)))*2
  ) +
  annotation_custom(rasterGrob(LP081C_S3$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP081C_S3$xrf$depth),
                    ymin = min(LP081C_S3$xrf$depth),
                    xmin = min(as.numeric(colnames(LP081C_S3$image$image))),
                    xmax = max(as.numeric(colnames(LP081C_S3$image$image)))
  ) +
  annotation_custom(rasterGrob(LP081D_S4$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP081D_S4$xrf$depth),
                    ymin = min(LP081D_S4$xrf$depth),
                    xmin = min(as.numeric(colnames(LP081D_S4$image$image))),
                    xmax = max(as.numeric(colnames(LP081D_S4$image$image)))*2
  ) +
  annotation_custom(rasterGrob(LP081E_S5$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP081E_S5$xrf$depth),
                    ymin = min(LP081E_S5$xrf$depth),
                    xmin = min(as.numeric(colnames(LP081E_S5$image$image))),
                    xmax = max(as.numeric(colnames(LP081E_S5$image$image)))
  )

Fig2
ggsave("Figures/Fig2_core_image_overlap.pdf", 
       height = c(15), width = c(5), dpi = 600, units = "cm")

egg::ggarrange(Fig2 + theme_paleo(), 
               Fig1 + theme(axis.title.y = element_blank(),
                            axis.text.y  = element_blank(),
                            axis.ticks.y = element_blank()), 
               ncol = 2, 
               widths = c(1, 4) # these are relative. For c(1, 5), the first plot will be 1/5th the width of the second.
)


# SECTION 5 - TRANSFORMING DATA -------------------------------------------

# validity filtered
LP08_xrfNorm <- LP08_xrf %>% # n = 2005
  mutate(`coh_inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc_coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method
  
  # transform
  mutate(across(any_of(elementsList)) /`Mo inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all Inf to NA
  
  # identify acceptable observations - validity and/or qc
  # comment in/out to choose either/or, or both
  filter(validity == TRUE) %>% # n = 1996
  #filter(qc == TRUE) %>% # n = 1707
  
  # identify acceptable variables
  select(P, S, any_of(myElements), `coh_inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar))

# pivot
LP08_xrfNorm_long <-  tidyr::pivot_longer(LP08_xrfNorm, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))

# plot
ggplot(LP08_xrfNorm_long, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo. inc.", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig13_Mo inc_normalised.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Validity and qc filtered
LP08_xrfNorm_qc <- LP08_xrf %>% # n = 2005
  mutate(`coh_inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc_coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method
  
  # transform
  mutate(across(any_of(elementsList)) /`Mo inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all Inf to NA
  
  # identify acceptable observations - validity and/or qc
  # comment in/out to choose either/or, or both
  filter(validity == TRUE) %>% # n = 1996
  filter(qc == TRUE) %>% # n = 1707
  
  # identify acceptable variables
  select(P, S, any_of(myElements), `coh_inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar))

# pivot
LP08_xrfNorm_qc_long <-  tidyr::pivot_longer(LP08_xrfNorm_qc, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))

# plot
ggplot(LP08_xrfNorm_long_qc, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo. inc.", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig13A_Mo inc_normalised_qc.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


# SMOOTHING  ---------------------------------------------------------------

LP08_xrfSmooth <- LP08_xrf %>%
  # uses a 10 point running mean (2 cm for this data); 5 before, 5 after - 1 cm i.e., 5 point RM 2.5 before/after doesnt work
  mutate(across(any_of(elementsList), 
                function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
  )
  ) 

ggplot(LP08_xrfSmooth, mapping = aes(x = depth, y = Ca)) + 
  geom_line(data = LP08_xrf, col = "grey80") + 
  geom_line() + 
  scale_x_reverse() +
  theme_paleo()
ggsave("Figures/Fig14_Ca_smoothed.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed stratigraphic diagram -------------------------------------------------
# smoothed data has to be labelled and combined with the original data so it can be faceted.
# make the xrf plot with running means

# make new dataset LP08_xrf1 with coh/inc and cps_sum included
LP08_xrf1 <- LP08_xrf %>%
  mutate(coh_inc = `Mo coh`/`Mo inc`) %>%
  mutate(inc_coh = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements]))# added by sjro to match Aber normalisation method
LP08_xrf1

# make new element list
elementsList1 <- select(LP08_xrf1, c(Mg:`Mo coh`, coh_inc, inc_coh, cps_sum)) %>% names()
elementsList1

# Smoothed cps plot - final join, smooth and plot with elements of most interest
full_join(y = LP08_xrf1 %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP08_xrf1 %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  filter(validity == TRUE) %>%
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  
  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) cps; validity filtered")
#axis.text.x = element_blank(),
#axis.ticks.x = element_blank()) 
ggsave("Figures/Fig15_smoothed_cps.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Smoothed Mo inc normalised plot - validity filtered - elements of most interest
full_join(y = LP08_xrfNorm %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP08_xrfNorm %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because LP08_xrfNorm has already been filtered
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  
  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo inc.", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) inc. normalised; validity filtered")
#axis.text.x = element_blank(),
#axis.ticks.x = element_blank())
ggsave("Figures/Fig16_smoothed_incNorm.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed Mo inc normalised plot - validity and qc filtered - elements of most interest
full_join(y = LP08_xrfNorm_qc %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP08_xrfNorm_qc %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because LP08_xrfNorm has already been filtered
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  
  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo inc.", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) inc. normalised; validity and qc filtered")
#axis.text.x = element_blank(),
#axis.ticks.x = element_blank())
ggsave("Figures/Fig17_smoothed_incNorm_qc.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed %cps_sum and Ti-log normalised plots - validity filtered - elements of most interest ***** TO DO **** 
# get code from Bertrand et al.R section 
# log normalised

# Write to file  --------------------------------------------------------
write.csv(LP08_xrf,"Output/LP08_xrf.csv", row.names = FALSE)
write.csv(LP08_xrf1,"Output/LP08_xrf1.csv", row.names = FALSE)
write.csv(LP08_xrfNorm,"Output/LP08_xrfNorm.csv", row.names = FALSE)
write.csv(LP08_xrfNorm_qc,"Output/LP08_xrfNorm_qc.csv", row.names = FALSE)


# SECTION 6 - MULTIVARIATE METHODS -------------------------------------------



# SECTION 7 - CALIBRATING DATA -------------------------------------------

