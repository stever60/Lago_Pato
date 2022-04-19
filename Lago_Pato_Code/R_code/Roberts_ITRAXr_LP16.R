
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
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2022/Data/ITRAX/LP16/LP16_16_5_Mo/")
#check working directory
getwd() 

# Renaming result.txt file - only run this once at start ----------------------------------------------------------

# contains extra lines D1, S1, S2, S3; D1 = Fe a*2 = 12.8 keV used in itrax.R example file

# Renaming columns in  old filenames to match ITRAX.r formats before using itraximport  ----------------------------------

# result.txt = original txt output file 
# Results = kcps renamed as cps
# Results1 = kcps & cps retained; cps calculated as element and scatter sum by code below - use this one 
# Results2 = cps calculated as element and scatter sum by code below & kcps removed - don't need this one
# need to run it for each core and change folder name using find and replace 
# Reanalysis data sometimes contains extra column called revalidity; validity is the original run validity

# import and rename cps as kcps - saved as Results.tsv
df<- read_tsv("LP161C/result.txt", col_names = TRUE, skip = 2) %>% 
  rename(cps = kcps, `Fe a*2` = D1) %>% 
  write_tsv("LP161C/Results.tsv")
df

# import, calculate cps_sum and rename as cps, retain kcps column - saved as Results1.tsv
df1_rowsums <- c(11:48, 53, 54)
df1 <- read_tsv("LP161C/result.txt", col_names = TRUE, skip = 2) %>%
  mutate(cps_sum = rowSums(.[df1_rowsums])) %>% 
  rename(cps = cps_sum, `Fe a*2` = D1) %>% 
  relocate(cps, .before = MSE) %>% 
  write_tsv("LP161C/Results1.tsv")
df1

# list of all elements in LP16 and LP16 files
cps_elementsList_LP16 <- select(df1, c(Mg:U, `Mo inc`, `Mo coh`)) %>% 
  names()
cps_elementsList_LP16

# Need to do this in Finder -----------------------------------------------

# open tsv Results and Results1 in excel 
# add first two row back into all new txt files using copy/paste
# rename .tsv files to .txt in Finder
# now ready to import into itrax.R and runs as normal
# choose whether to use Results.txt or Results1.txt
# below uses Results1.txt
# repeat for all sections 

# SECTION 1 - DATA STRUCTURE & TESTING FILE IMPORTS WORK OK ----------------------------------------------

# import the data - Results file has cps instead of kcps so that ITRAX.r will run - original file is in 2008 folder 
itrax_import("LP161A/Results1.txt", depth_top = 0) %>% # import xrf data
  ggplot(aes(x=`Mo inc`/`Mo coh`, y=depth)) + # setup plot area
  geom_point() + # plots points             
  geom_lineh() + # plots the line
  scale_y_reverse() + # reverses the scale
  theme_bw()

# Check metadata in document.txt and create new document1.txt --------------

itrax_meta("LP161A/document.txt")

# Make sure extra data for earlier Aber run cores are added to match Manchester metadata
# Check start/end/Step size etc. are correct
# Use ARD1A document.txt in ITRAX.r file as a template or copy/paste missing info from below

#Optical Start	285.0	Optical End	785.0
#Step size	500	microns	
#Start temperature		∞C	
#Start humidity		%	
#Start vacuum	-95	kPa	
#Stop temperature		∞C	
#Stop humidity		%	
#Stop vacuum	-94.8	kPa

# 1 - Make a copy of document.txt and 
 
# 2 - For each core, need to make sure that the XRF start/stop position (mm) coordinates and Optical start/end coordinates 
# are the same if using optical1 image that has been cropped to fit the XRF scan limits
# otherwise the data will be misaligned

# 3 - Check adiograph start/stop has the same start/stop coordinates as the XRF scan

# 4 - Save changes as document1.txt

itrax_meta("LP161A/document1.txt") #

# Adjust optical image and radiograph in Photoshop ------------------------

# optica1 = autotone/autocontrast/autocolour then +100 brightness 
# radiograph1 = autotone/autocontrast/autocolour then +150 brightness x2 for LP2A-LP4B; LP1A, 1C +150 x1 brightness 

# LP161A images -----------------------------------------------------------

# plot the optical image vs position - check LP161A - optical1 is Photoshop enhanced and cropped image
opt <- itrax_image(file = "LP161A/optical1.tif", # define location of image file
            meta = "LP161A/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image
radx <- itrax_radiograph(file = "LP161A/radiograph1.tif", # define location of radiograph image
                 meta = "LP161A/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure


# LP161B images -----------------------------------------------------------

# plot the optical image vs position - check LP161A - optical1 is Photoshop enhanced and cropped image
opt <- itrax_image(file = "LP161B/optical1.tif", # define location of image file
                   meta = "LP161B/document1.txt", # define location of associated metadata
                   plot = TRUE,
                   trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image
radx <- itrax_radiograph(file = "LP161B/radiograph1.tif", # define location of radiograph image
                         meta = "LP161B/document1.txt", # define location of associated metadata
                         plot = TRUE,
                         trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure


# LP161C Images -----------------------------------------------------------

# plot the optical image vs position - check LP161C - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP161C/optical1.tif", # define location of image file
            meta = "LP161C/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP161C/radiograph1.tif", # define location of radiograph image
                 meta = "LP161C/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

# LP162A Images -----------------------------------------------------------

# plot the optical image vs position - check LP162A - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP162A/optical1.tif", # define location of image file
            meta = "LP162A/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP162A/radiograph1.tif", # define location of radiograph image
                 meta = "LP162A/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

# LP162B Images -----------------------------------------------------------

# plot the optical image vs position - check LP162B - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP162B/optical1.tif", # define location of image file
            meta = "LP162B/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP162B/radiograph1.tif", # define location of radiograph image
                 meta = "LP162B/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

# LP163A_200um Images -----------------------------------------------------------

# plot the optical image vs position - check LP163A_200um - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP163A_200um/optical1.tif", # define location of image file
            meta = "LP163A_200um/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP163A_200um/radiograph1.tif", # define location of radiograph image
                 meta = "LP163A_200um/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure


# LP163B_200um Images -----------------------------------------------------------

# plot the optical image vs position - check LP163B_200um - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP163B_200um/optical1.tif", # define location of image file
            meta = "LP163B_200um/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP163B_200um/radiograph1.tif", # define location of radiograph image
                 meta = "LP163B_200um/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

# LP164A Images -----------------------------------------------------------

# plot the optical image vs position - check LP164A - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP164A/optical1.tif", # define location of image file
            meta = "LP164A/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP164A/radiograph1.tif", # define location of radiograph image
                 meta = "LP164A/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure

# LP164B Images -----------------------------------------------------------

# plot the optical image vs position - check LP164B - optical1 is Photoshop enhanced and cropped image
itrax_image(file = "LP164B/optical1.tif", # define location of image file
            meta = "LP164B/document1.txt", # define location of associated metadata
            plot = TRUE,
            trim = FALSE
) %>%  #pipe (send) the output of that to...
  str() # summarise the structure

# plot the radiograph vs position - radiograph1 is Photoshop enhanced and cropped image 
itrax_radiograph(file = "LP164B/radiograph1.tif", # define location of radiograph image
                 meta = "LP164B/document1.txt", # define location of associated metadata
                 plot = TRUE,
                 trim = FALSE
) %>% # pipe (send) the output of that to...
  str() # summarise the structure


# Check Qpec settings and energy spectra ----------------------------------

itrax_qspecsettings("LP161C/settings.dfl") # parse some q-spec settings

itrax_spectra(filename = "LP161C/sumspectra.spe", # define a raw spectra file or sumspectra file
              parameters = "LP161C/settings.dfl", # define an associated settings file
              plot = TRUE # suppress the plot
) 

# SECTION 2 - IMPORTING DATA ------------------------------------------------

# LP161A using document1 and Results1 files - 1A is below 1C

# LP161C using Results1.txt and document1.txt
LP161C_S1 <- list(metadata   = itrax_meta("LP161C/document1.txt"),
                  xrf        = itrax_import("LP161C/Results1.txt", 
                                            depth = 67,
                                            parameters = "all"),
                  image      = itrax_image(file = "LP161C/optical1.tif",
                                           meta = "LP161C/document1.txt"),
                  radiograph = itrax_radiograph(file = "LP161C//radiograph1.tif",
                                                meta = "LP161C/document1.txt",
                                                trim = as.numeric(itrax_meta("LP161C/document1.txt")[6:7,2])))

LP161A_S2 <- list(metadata   = itrax_meta("LP161A/document1.txt"),
                 xrf        = itrax_import("LP161A/Results1.txt", 
                                           depth = 280, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP161A/optical1.tif",
                                          meta = "LP161A/document1.txt"),
                 radiograph = itrax_radiograph(file = "LP161A//radiograph1.tif",
                                               meta = "LP161A/document1.txt",
                                               trim = as.numeric(itrax_meta("LP161A/document1.txt")[6:7,2])))

LP161B_S3 <- list(metadata   = itrax_meta("LP161B/document1.txt"),
                  xrf        = itrax_import("LP161B/Results1.txt", 
                                            depth = 400, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP161B/optical1.tif",
                                           meta = "LP161B/document1.txt"),
                  radiograph = itrax_radiograph(file = "LP161B//radiograph1.tif",
                                                meta = "LP161B/document1.txt",
                                                trim = as.numeric(itrax_meta("LP161B/document1.txt")[6:7,2])))

LP162A_S4 <- list(metadata   = itrax_meta("LP162A/document1.txt"),
                 xrf        = itrax_import("LP162A/Results1.txt", 
                                           depth = 560, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP162A/optical1.tif",
                                          meta = "LP162A/document1.txt"),
                 radiograph = itrax_radiograph(file = "LP162A//radiograph1.tif",
                                               meta = "LP162A/document1.txt",
                                               trim = as.numeric(itrax_meta("LP162A/document1.txt")[6:7,2])))

LP162B_S5 <- list(metadata   = itrax_meta("LP162B/document1.txt"),
                 xrf        = itrax_import("LP162B/Results1.txt", 
                                           depth = 935, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP162B/optical1.tif",
                                          meta = "LP162B/document1.txt"),
                 radiograph = itrax_radiograph(file = "LP162B//radiograph1.tif",
                                               meta = "LP162B/document1.txt",
                                               trim = as.numeric(itrax_meta("LP162B/document1.txt")[6:7,2])))

LP163A_S6 <- list(metadata   = itrax_meta("LP163A_200um/document1.txt"),
                 xrf        = itrax_import("LP163A_200um/Results1.txt", 
                                           depth = 1295, 
                                           parameters = "all"),
                 image      = itrax_image(file = "LP163A_200um/optical1.tif",
                                          meta = "LP163A_200um/document1.txt"),
                 radiograph = itrax_radiograph(file = "LP163A_200um//radiograph1.tif",
                                               meta = "LP163A_200um/document1.txt",
                                               trim = as.numeric(itrax_meta("LP163A_200um/document1.txt")[6:7,2])))

LP163B_S7 <- list(metadata   = itrax_meta("LP163B_200um/document1.txt"),
                  xrf        = itrax_import("LP163B_200um/Results1.txt", 
                                            depth = 1691, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP163B_200um/optical1.tif",
                                           meta = "LP163B_200um/document1.txt"),
                  radiograph = itrax_radiograph(file = "LP163B_200um//radiograph1.tif",
                                                meta = "LP163B_200um/document1.txt",
                                                trim = as.numeric(itrax_meta("LP163B_200um/document1.txt")[6:7,2])))

LP164A_S8 <- list(metadata   = itrax_meta("LP164A/document1.txt"),
                  xrf        = itrax_import("LP164A/Results1.txt", 
                                            depth = 2058, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP164A/optical1.tif",
                                           meta = "LP164A/document1.txt"),
                  radiograph = itrax_radiograph(file = "LP164A//radiograph1.tif",
                                                meta = "LP164A/document1.txt",
                                                trim = as.numeric(itrax_meta("LP164A/document1.txt")[6:7,2])))

LP164B_S9 <- list(metadata   = itrax_meta("LP164B/document1.txt"),
                  xrf        = itrax_import("LP164B/Results1.txt", 
                                            depth = 2475.2, 
                                            parameters = "all"),
                  image      = itrax_image(file = "LP164B/optical1.tif",
                                           meta = "LP164B/document1.txt"),
                  radiograph = itrax_radiograph(file = "LP164B//radiograph1.tif",
                                                meta = "LP164B/document1.txt",
                                                trim = as.numeric(itrax_meta("LP164B/document1.txt")[6:7,2])))


# join the xrf data for the sections together ---
LP16_xrf <- itrax_join(list(S1 = LP161C_S1$xrf, # S1 will be the "label" given to the core section
                           S2 = LP161A_S2$xrf,
                           S3 = LP161B_S3$xrf,
                           S4 = LP162A_S4$xrf,
                           S5 = LP162B_S5$xrf, 
                           S6 = LP163A_S6$xrf,
                           S7 = LP163B_S7$xrf, 
                           S8 = LP164A_S8$xrf,
                           S9 = LP164B_S9$xrf)
)

LP16_xrf
write.csv(LP16_xrf,"LP16_xrf.csv", row.names = FALSE)

# Figure 1 - Summary overlaps plot - all sections -------------------------

Fig1.1 <- ggplot(data = na.omit(LP16_xrf), mapping = aes(x = depth, y = Ti)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.2 <- ggplot(data = na.omit(LP16_xrf), mapping = aes(x = depth, y = Br)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.3 <- ggplot(data = na.omit(LP16_xrf), mapping = aes(x = depth, y = `Mo coh`/`Mo inc`)) +
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
Fig2 <- ggplot(data = LP16_xrf, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
Fig2
ggsave("Figures/Fig2_tolerance_Fe_count.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# cps tolerance filter - could use >mean+/-2s cps (based on kcps or cps_sum)
cps.min.thres <- 30000
cps.max.thres <- 70000

#  OR 

cps.mean <- mean(LP16_xrf$cps)
cps.sd <- 3*sd(LP16_xrf$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

LP16_xrf  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw()
ggsave("Figures/Fig3_tolerance_cps.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# MSE tolerance filter - 2 used here but could use >mean+2s - which is 1.764766 for LP16 record 

MSE.thres <- 2 # use this for LP16

#  OR

MSE.mean <- mean(LP16_xrf$MSE)
MSE.sd <- 3*sd(LP16_xrf$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

LP16_xrf %>%
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

#  OR - used for LP16
slope1 <-  LP16_xrf$`sample surface` - lag(LP16_xrf$`sample surface`)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 3*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

LP16_xrf %>%
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
LP16_xrf <- LP16_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) # %>% filter(qc == TRUE) #to remove from LP16_xrf rows that dont pass QC
# plot summary
theme_set(theme_bw(8))
Fig6.1 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  #geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.2 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  #geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.3 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  #geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.4 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  #geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig6.1, Fig6.2, Fig6.3, Fig6.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig6_tolerance_combined.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

LP16_xrf 

# Combining 'validity' flags - removing data that doesn't pass - replotting as line only
LP16_xrf1 <- LP16_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=-0.3 | slope >=0.3 | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=20000 | cps >=60000 | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=2, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) %>%  
  filter(qc == TRUE) #to remove from LP16_xrf rows that dont pass QC

# plot summary
theme_set(theme_bw(8))
Fig7.1 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.2 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.3 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.4 <- ggplot(data = LP16_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig7.1, Fig7.2, Fig7.3, Fig7.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig7_tolerance_filtered.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correcting for dead time & Uncertainties ---------------------------------------------- 

# This doesn't work on cores scanned at Aber - no Dt (Dwell time) column in results.txt file

# ggplot(data = LP16_xrf, aes(x = depth, y = Dt)) +
#  scale_x_reverse() + 
#  scale_y_continuous(sec.axis = sec_axis( trans=~(.+(1-mean(LP16_xrf$Dt, na.rm = TRUE))), name="Correction Factor")) +
#  geom_line() +
#  geom_hline(yintercept = mean(LP16_xrf$Dt, na.rm = TRUE), linetype = "dotted")


# Noisy data ------------------------------------------------------------

# see also section: # Calculate as % of normalising factors TS (Total Scatter), cps_sum - in Bertrand et al.R

# Using autocorrelation to detect noisy signals - non-noisy data should be AC, higher, outside 95% limits, showing some order/pattern
library(forecast)
library(ggpubr)

# Individual elements
Fig8 <- ggarrange(
  ggAcf(LP16_xrf$Ca) + ylim(c(NA,1)), ggAcf(LP16_xrf$Ti) + ylim(c(NA,1)), 
  ggAcf(LP16_xrf$Fe) + ylim(c(NA,1)),  ggAcf(LP16_xrf$Sr) + ylim(c(NA,1)),
  ggAcf(LP16_xrf$P) + ylim(c(NA,1)),ggAcf(LP16_xrf$Cu) + ylim(c(NA,1)), 
  ggAcf(LP16_xrf$Zn) + ylim(c(NA,1)), ggAcf(LP16_xrf$S) + ylim(c(NA,1)),
  ggAcf(LP16_xrf$Ni) + ylim(c(NA,1)), ggAcf(LP16_xrf$Cs) + ylim(c(NA,1)),
  ggAcf(LP16_xrf$Ba) + ylim(c(NA,1)), ggAcf(LP16_xrf$Pb) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
Fig8
ggsave("Figures/Fig8_ACF.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# all elements in summary plot
elementsList <- select(LP16_xrf, c(Mg:`Mo coh`)) %>% 
  names()
elementsList

apply(LP16_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
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
apply(LP16_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
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
# add P for LP16 regression analysis as not in selected
# remove Ar, Ta, W (if selected above) - these are detector generated elements
LP16_xrf %>% 
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
xrfStrat <- LP16_xrf %>% 
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
  scale_y_continuous(limits = rev(range(LP16_xrf$depth))) +
  scale_x_continuous(breaks = round(c(0, max(as.numeric(colnames(LP161C_S2$image$image)))*2)),
                     limits = c(0, max(as.numeric(colnames(LP161C_S2$image$image)))*2)) +
  coord_fixed(ratio = 1) +
  labs(y = "Depth [mm]", x = "[mm]") +
  annotation_custom(rasterGrob(LP161A_S1$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP161A_S1$xrf$depth),
                    ymin = min(LP161A_S1$xrf$depth),
                    xmin = min(as.numeric(colnames(LP161A_S1$image$image))),
                    xmax = max(as.numeric(colnames(LP161A_S1$image$image)))
  ) +
  annotation_custom(rasterGrob(LP161C_S2$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP161C_S2$xrf$depth),
                    ymin = min(LP161C_S2$xrf$depth),
                    xmin = max(as.numeric(colnames(LP161C_S2$image$image))),
                    xmax = max(as.numeric(colnames(LP161C_S2$image$image)))*2
  ) +
  annotation_custom(rasterGrob(LP162A_S3$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP161C_S3$xrf$depth),
                    ymin = min(LP161C_S3$xrf$depth),
                    xmin = min(as.numeric(colnames(LP161C_S3$image$image))),
                    xmax = max(as.numeric(colnames(LP161C_S3$image$image)))
  ) +
  annotation_custom(rasterGrob(LP162B_S4$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP162B_S4$xrf$depth),
                    ymin = min(LP162B_S4$xrf$depth),
                    xmin = min(as.numeric(colnames(LP162B_S4$image$image))),
                    xmax = max(as.numeric(colnames(LP162B_S4$image$image)))*2
  ) +
  annotation_custom(rasterGrob(LP163A_200um_S5$image$image,
                               width = unit(1, "npc"),
                               height = unit(1, "npc")),
                    ymax = max(LP163A_200um_S5$xrf$depth),
                    ymin = min(LP163A_200um_S5$xrf$depth),
                    xmin = min(as.numeric(colnames(LP163A_200um_S5$image$image))),
                    xmax = max(as.numeric(colnames(LP163A_200um_S5$image$image)))
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
LP16_xrfNorm <- LP16_xrf %>% # n = 2005
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
LP16_xrfNorm_long <-  tidyr::pivot_longer(LP16_xrfNorm, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))

# plot
ggplot(LP16_xrfNorm_long, aes(x = peakarea, y = depth)) +
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
LP16_xrfNorm_qc <- LP16_xrf %>% # n = 2005
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
LP16_xrfNorm_qc_long <-  tidyr::pivot_longer(LP16_xrfNorm_qc, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))

# plot
ggplot(LP16_xrfNorm_long_qc, aes(x = peakarea, y = depth)) +
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

LP16_xrfSmooth <- LP16_xrf %>%
  # uses a 10 point running mean (2 cm for this data); 5 before, 5 after - 1 cm i.e., 5 point RM 2.5 before/after doesnt work
  mutate(across(any_of(elementsList), 
                function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
  )
  ) 

ggplot(LP16_xrfSmooth, mapping = aes(x = depth, y = Ca)) + 
  geom_line(data = LP16_xrf, col = "grey80") + 
  geom_line() + 
  scale_x_reverse() +
  theme_paleo()
ggsave("Figures/Fig14_Ca_smoothed.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed stratigraphic diagram -------------------------------------------------
# smoothed data has to be labelled and combined with the original data so it can be faceted.
# make the xrf plot with running means

# make new dataset LP16_xrf1 with coh/inc and cps_sum included
LP16_xrf1 <- LP16_xrf %>%
  mutate(coh_inc = `Mo coh`/`Mo inc`) %>%
  mutate(inc_coh = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements]))# added by sjro to match Aber normalisation method
LP16_xrf1

# make new element list
elementsList1 <- select(LP16_xrf1, c(Mg:`Mo coh`, coh_inc, inc_coh, cps_sum)) %>% names()
elementsList1

# Smoothed cps plot - final join, smooth and plot with elements of most interest
full_join(y = LP16_xrf1 %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP16_xrf1 %>% 
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
full_join(y = LP16_xrfNorm %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP16_xrfNorm %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because LP16_xrfNorm has already been filtered
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
full_join(y = LP16_xrfNorm_qc %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = LP16_xrfNorm_qc %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because LP16_xrfNorm has already been filtered
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
write.csv(LP16_xrf,"Output/LP16_xrf.csv", row.names = FALSE)
write.csv(LP16_xrf1,"Output/LP16_xrf1.csv", row.names = FALSE)
write.csv(LP16_xrfNorm,"Output/LP16_xrfNorm.csv", row.names = FALSE)
write.csv(LP16_xrfNorm_qc,"Output/LP16_xrfNorm_qc.csv", row.names = FALSE)


# SECTION 6 - MULTIVARIATE METHODS -------------------------------------------



# SECTION 7 - CALIBRATING DATA -------------------------------------------

