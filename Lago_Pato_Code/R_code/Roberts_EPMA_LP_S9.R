# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021")
#check working directory
getwd()

# Load libraries ----------------------------------------------------------

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
library(tidyverse)
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
library(itraxR)
library(PeriodicTable)
library(tidypaleo)
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


# Expand the BrBG colour palette from 11 to 18 colours & get colour codes to copy  --------
nb.cols <- 11
BrBG1 <- colorRampPalette(brewer.pal(11, "BrBG"))(nb.cols)
BrBG1
# BrBG1 11 colour codes 
# "#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30"

nb.cols <- 10
BrBG2 <- colorRampPalette(brewer.pal(10, "BrBG"))(nb.cols)
BrBG2

nb.cols <- 5
RdYlBu1 <- colorRampPalette(brewer.pal(5, "RdYlBu"))(nb.cols)
RdYlBu1
# RdYlBu1 colour codes 
#"#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"

nb.cols <- 7
RdYlBu1 <- colorRampPalette(brewer.pal(7, "RdYlBu"))(nb.cols)
RdYlBu1
# RdYlBu1 colour codes 
#"#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"

nb.cols <- 10
RdYlBu2 <- colorRampPalette(brewer.pal(10, "RdYlBu"))(nb.cols)
RdYlBu2
# "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"

# # Gabriel et al specific colour schemes ---------------------------------

#6: SAM, ANT, NAP, DI, K33, K58 based on RdYlBu direction -1
SAM <- c("#4575B4", "#ABD9E9", "#FEE090", "#DFC27D", "#F46D43", "#A50026")

#5: SAM, NAP, DI, K33, K58 based on RdYlBu direction -1
SAM1 <- c("#4575B4", "#FEE090", "#DFC27D", "#F46D43", "#A50026")

#4: SAM, NAP, K33, K58 based on RdYlBu direction -1
SAM2 <- c("#4575B4", "#FEE090",  "#F46D43", "#A50026")

##7: SAM, ANT, NAP, DI, LP16_12.5, LP16_32, LP16_57 based on RdYlBu direction -1
SAM3 <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")


# Roberts et al colour and symbol schemes ---------------------------------

# Set up Eruption_ID grouping / colours 
##19: A1	MB2	So_UK	H1	UK	Ax	R1	MB1	MB2	H2	Hb	Ha	H1991	A1	MBX	MB1910	MB1910? LP16_12.5, LP16_32, LP16_57 based on RdYlBu & BrBG
eruptions <- c("A1","Ax","H1991","H1","H2","Ha","Hb","MB1","MB2","MBX", "MB1910", "MB1910?",
               "R1","R2","So_UK", "UK", "LP16_12.5_UK", "LP16_32_UK", "LP16_57_UK")

erup_col <- c("#FDAE61", "#DFC27D", "#BF812D","#FDAE61", "#FEE090", "#ABD9E9", "#74ADD1", "#4575B4", "#313695",
          "#003C30","#01665E","#35978F", "#ABD9E9","#4575B4", "#DFC27D", "#BF812D", "#A50026", "#D73027", "#F46D43")

erup_symbols <- c(3, 3, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 22, 22, 21, 4, 23, 25, 24)

erup_size <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2)


# TO DO - NORMALISE RAW DATA TO 100%  - DONT RUN YET ----------------------

# Import log transformed data file .raw or .inter from above 
x <- read.csv("x_inter.log.csv", 
              header=TRUE, row.names=1, sep=",", check.names=FALSE)
head(x)

##create % sum dataset if not as % data already or sums to 100% if not already 
totals<-apply(x, 1, sum)
x.pc<-x/totals*100
mx<-apply(x.pc, 2, max)
x.pc.5<-x.pc[, mx >0.5]


# START - Import and filter data to remove totals <95% for ACID and <97% for BS & INT data  ------------------------------------------------------------

# Whole EMPA_db.csv file: https://www.dropbox.com/s/v2mkattn7k4m9gs/EPMA_db.csv?dl=0 
EPMA_db <- read_csv("/Users/Steve/Dropbox/BAS/Data/R/Papers/Gabriel_et_al_2021/Data/EPMA_db.csv")


# remove Blockley and Langdon unpublished data and any other data based on the reference
EPMA_db1 <- EPMA_db[!EPMA_db$Geochem_Ref1=='Blockley_Langdon_unpublished',]
EPMA_db1.1 <- EPMA_db1[!EPMA_db1$Geochem_Ref1=='Gabriel et al_submitted',]
EPMA_db2.1 <- EPMA_db1.1[!EPMA_db1.1$Group=='NAP',]
EPMA_db2.2 <- EPMA_db2.1[!EPMA_db2.1$Group=='ANT',]
EPMA_db2 <- EPMA_db2.2[!EPMA_db2.2$Group=='DI',]

# Database to take forward for analysis and plotting  
EPMA_df <- select(EPMA_db2, ID:Total_n, Type, TAS_Class, Group, Eruption_ID:Volcano) %>% 
  replace(is.na(.), 0) %>% #convert NA to 0
  mutate_if(is.logical,as.numeric) #to convert SO2 FALSE to 0 - overcomes parsing problem with SO2
# To convert NA to 0 in base R: EPMA_df[is.na(EPMA_df)] <- 0
EPMA_df

# filter by Type and then remove analyses with low totals
# filtering based on Hunt and Hill (1995) <95% and Rutledal et al (2020): https://www.sciencedirect.com/science/article/pii/S0277379119309965#bib40 
EPMA_BS <- filter(EPMA_df, Type == "Basic") %>% 
  filter(Total > 97)
EPMA_BS

EPMA_INT <- filter(EPMA_df, Type == "Intermediate") %>% 
  filter(Total > 97)
EPMA_INT

EPMA_ACID <- filter(EPMA_df, Type == "Acid") %>% 
  filter(Total > 95)
EPMA_ACID

EPMA_ACID <- filter(EPMA_df, Type == "Acid") %>% 
  filter(Total > 95)
EPMA_ACID

EPMA_RHY <- filter(EPMA_df, TAS_Class == "Rhyolite") %>% 
  filter(Total > 95)
EPMA_RHY

# bind function - assumes order of first one etc. so then need to reorder each as SAM, ANT etc. within Basic, Int, Acid
EPMA_ALL_filter <- bind_rows(EPMA_BS, EPMA_INT,EPMA_ACID)
EPMA_ALL_filter$Eruption_ID <- factor(EPMA_ALL_filter$Eruption_ID, 
                                levels = eruptions)
EPMA_ALL_filter %>% 
  arrange(Eruption_ID)

#BS+INT data- arrnage by groups for plotting and PCA
EPMA_BSINT_filter <- bind_rows(EPMA_BS, EPMA_INT)
EPMA_BSINT_filter$Eruption_ID <- factor(EPMA_BSINT_filter$Eruption_ID, 
                                levels = eruptions)
EPMA_BSINT_filter %>% 
  arrange(Eruption_ID)

#ACID data- arrange by groups for plotting and PCA
EPMA_ACID_filter <- EPMA_ACID
EPMA_ACID_filter$Eruption_ID <- factor(EPMA_ACID_filter$Eruption_ID, 
                                 levels = eruptions)
EPMA_ACID_filter %>% 
  arrange(Eruption_ID)

#RHY data- arrange by groups for plotting and PCA
EPMA_RHY_filter <- EPMA_RHY
EPMA_RHY_filter$Eruption_ID <- factor(EPMA_RHY_filter$Eruption_ID, 
                                 levels = eruptions)
EPMA_RHY_filter %>% 
  arrange(Eruption_ID)

# export to csv file 
# all filtered data to file
write_csv(EPMA_ALL_filter, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/EPMA_ALL_filter.csv")
# BS+INT filtered data to file
write_csv(EPMA_BSINT_filter, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/EPMA_BS+INT_filter.csv")
# ACID filtered data to file
write_csv(EPMA_ACID_filter, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/EPMA_ACID_filter.csv")
# RHY filtered data to file
write_csv(EPMA_RHY_filter, "/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/EPMA_RHY_filter.csv")

# Convert data to long format for ggplot and tidypaleo facet plotting --------------------------

# select out the parameters to plot, and convert to long format - i.e., one value per row per variable - for later plotting .... perhaps ... 
EPMA_long <- select(EPMA_ALL_filter, ID, LabID, LongID, SiO2:K2O, Type, Group, Eruption_ID, Volcano) %>%
  pivot_longer(c(`SiO2`:`K2O`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Type)
EPMA_long

# convert normalised data only to long format for plotting
EPMA_long_n <- select(EPMA_ALL_filter, ID, LabID, LongID, SiO2_n:K2O_n, Type, Group, Eruption_ID, Volcano) %>%
  pivot_longer(c(`SiO2_n`:`K2O_n`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Type)
EPMA_long_n

# convert normalised RHYOLITIC DATA data only to long format for plotting
EPMA_long_RHY <- select(EPMA_RHY_filter , ID, LabID, LongID, SiO2_n:K2O_n, Type, Group, Eruption_ID, Volcano) %>%
  pivot_longer(c(`SiO2_n`:`K2O_n`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Type)
EPMA_long_RHY

#  Create TAS diagram background plot  ------------------------------------

# create a data frame
d = data.frame(x = c(40, 80), y = c(0,15))
theme_set(theme_bw(base_size=12))

##make the TAS template
p_ALL <- ggplot(data=d, mapping=aes(x=x, y=y)) +
  geom_blank() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(0,15), breaks=seq(0,15,2), expand = c(0, 0)) + 
  scale_x_continuous(limits=c(40,80), breaks=seq(40,80,5), expand = c(0, 0)) +
  labs(y=expression(Na[2]*O + K[2]*O*~ wt~'%'), x=expression(SiO[2]*~ wt~'%'))+
  annotate("segment", x=45, xend=45, y=0, yend=5)+
  annotate("segment", x=45, xend=52, y=5, yend=5)+
  annotate("segment", x=52, xend=69, y=5, yend=8)+
  annotate("segment", x=76.5, xend=69, y=0, yend=8)+
  annotate("segment", x=69, xend=69, y=8, yend=13)+
  annotate("segment", x=45, xend=61.32, y=5, yend=13.7)+
  annotate("segment", x=52, xend=52, y=0, yend=5)+
  annotate("segment", x=57, xend=57, y=0, yend=5.9)+
  annotate("segment", x=63, xend=63, y=0, yend=6.9)+
  annotate("segment", x=52, xend=49.4, y=5, yend=7.3)+
  annotate("segment", x=57, xend=53.05, y=5.9, yend=9.25)+
  annotate("segment", x=63, xend=57.6, y=6.9, yend=11.7)+
  annotate("segment", x=41, xend=45, y=3, yend=3)+
  annotate("segment", x=41, xend=41, y=0, yend=3)+
  annotate("segment", x=41, xend=41, y=3, yend=7, linetype="dashed")+
  annotate("segment", x=41, xend=45, y=7, yend=9.4, linetype="dashed")+
  annotate("segment", x=45, xend=52.5, y=9.4, yend=14)+
  annotate("segment", x=49.4, xend=45, y=7.3, yend=9.4)+
  annotate("segment", x=53, xend=48.4, y=9.3, yend=11.5)+
  annotate("segment", x=57.6, xend=50.3, y=11.7, yend=15) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
        )

textsize <- 4
c = "gray40"
## Add the text annotations
tas <- p_ALL + annotate("text", label = "Basalt", x = 48.5, y = 2, size=textsize, colour = c) +
  annotate("text", label = "Basaltic\n andesite", x = 54.8, y = 2.7, size=textsize, colour = c) +
  annotate("text", label = "Andesite", x = 60, y = 3.5, size=textsize, colour = c) +
  annotate("text", label = "Dacite", x = 67.5, y = 4.2, size=textsize, colour = c) +
  annotate("text", label = "Rhyolite", x = 75, y = 8, size=textsize, colour = c) +
  annotate("text", label = "Trachy- \n basalt", x = 48.8, y = 5.7, size=textsize, colour = c) +
  annotate("text", label = "Basaltic \n trachy- \n andesite", x = 52.5, y = 7, size=textsize, colour = c) +
  annotate("text", label = "Trachy- \n andesite", x = 57.8, y = 8.2, size=textsize, colour = c) +
  annotate("text", label = "Trachydacite", x = 65, y = 9, size=textsize, colour = c) +
  annotate("text", label = "Trachyte", x = 62.5, y = 11.5, size=textsize, colour = c) +
  annotate("text", label = "Picro- \n basalt", x = 43, y = 1.5, size=textsize, colour = c) +
  annotate("text", label = "Basanite", x = 44, y = 6, size=textsize, colour = c) +
  annotate("text", label = "Tephrite", x = 43.5, y = 7, size=textsize, colour = c) +
  annotate("text", label = "Phono- \n tephrite", x = 48.5, y = 9.5, size=textsize, colour = c) +
  annotate("text", label = "Tephri- \n phonolite", x = 52.5, y = 11.5, size=textsize, colour = c) +
  annotate("text", label = "Phonolite", x = 57, y = 14, size=textsize, colour = c) +
  annotate("text", label = "Foidite", x = 42, y = 12, size=textsize, colour = c)

tas

# Figure S9A - adding geochem data to TAS plot  ----------------------------

# Manual colour scheme - SAM defined at start 
p1 <- tas + 
  geom_point(data = EPMA_ALL_filter, aes(x=SiO2_n, y=Na2O_n+K2O_n, colour = Eruption_ID, fill = Eruption_ID, 
                                         shape = Eruption_ID, size = Eruption_ID), size = 1, alpha = 1)  + 
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  scale_size_manual(values=erup_size) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10), 
        legend.justification = c(1,1),
        legend.position = "bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))

# wait a few seconds then run this to avoid timing out  
p1

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS9A.pdf", height = c(20), width = c(20), dpi = 600, units = "cm")

# Figure S9B Make a zoom in plot of rhyolite -----------------------------------------

# create a data frame zoom in on rhyloite 
e = data.frame(x = c(60, 80), y = c(2,12))
theme_set(theme_bw(base_size=12))

##make the TAS template
p_RHY <- ggplot(data=e, mapping=aes(x=x, y=y)) +
  geom_blank() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(0,13), breaks=seq(0,13,2), expand = c(0, 0)) + 
  scale_x_continuous(limits=c(63,80), breaks=seq(63,80,2), expand = c(0, 0)) +
  labs(y=expression(Na[2]*O + K[2]*O*~ wt~'%'), x=expression(SiO[2]*~ wt~'%')) +
  annotate("segment", x=76.5, xend=69, y=0, yend=8)+
  annotate("segment", x=69, xend=69, y=8, yend=13)+
  annotate("segment", x=63, xend=69, y=6.9, yend=8) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

textsize <- 4
c = "gray40"
## Add the text annotations
tas_RHY <- p_RHY +
  annotate("text", label = "Dacite", x = 67.5, y = 4.2, size=textsize, colour = c) +
  annotate("text", label = "Rhyolite", x = 75, y = 8, size=textsize, colour = c) +
  annotate("text", label = "Trachydacite", x = 66, y = 10, size=textsize, colour = c)
  
tas_RHY

# Add the data Manual colour scheme - SAM defined at start 
p2 <- tas_RHY + 
  geom_point(data = EPMA_ALL_filter, aes(x=SiO2_n, y=Na2O_n+K2O_n, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1, alpha = 1)  + 
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10), 
        legend.justification = c(1,1),
        legend.position = "bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))

# wait a few seconds then run this to avoid timing out  
p2

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS9B.pdf", height = c(16), width = c(16), dpi = 600, units = "cm")

#  Figure S9C - selected bi-plots ---------------------------------------------------------------

theme_set(theme_bw(base_size=12))
#clear plot window
dev.off()

# K vs Si
p3g <- ggplot(data = EPMA_ALL_filter, aes(x = SiO2_n, y = K2O_n))
p3 <- p3g + geom_point(data = EPMA_ALL_filter, aes(x=SiO2_n, y=K2O_n, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1, alpha = 1)  +
  labs(x=expression(SiO[2]*~ '(wt %)'), y=expression(K[2]*O*~ '(wt %)'))+
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12), 
        legend.justification = c(1,1),
        legend.position = "bottom", 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
p3

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS9C.pdf", height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure S10A-F -------------------------------------------------------------

# Ca vs Fe
p5g <- ggplot(data = EPMA_ALL_filter, aes(x = FeO_n, y = CaO_n))
p5 <- p5g + geom_point(data = EPMA_ALL_filter, aes(x = FeO_n, y = CaO_n, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1.5, alpha = 1)  +
  labs(y=expression('CaO (wt %)'), x=expression(FeO*~ '(wt %)')) +
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12), 
        legend.justification = c(1,1),
        legend.position = "bottom", 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=16, colour = "black"), axis.title=element_text(size=16, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
p5

# Ca vs K
p6g <- ggplot(data = EPMA_ALL_filter, aes(x = K2O_n, y = CaO_n))
p6 <- p6g + geom_point(data = EPMA_ALL_filter, aes(x = K2O, y = CaO, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1.5, alpha = 1)  +
  labs(y=expression('CaO (wt %)'), x= expression(K[2]*O*~ '(wt %)')) +
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12), 
        legend.justification = c(1,1),
        legend.position = "bottom", 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=16, colour = "black"), axis.title=element_text(size=16, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
p6

# Mg vs Fe
p7g <- ggplot(data = EPMA_ALL_filter, aes(x = FeO_n, y = MgO_n))
p7 <- p7g + geom_point(data = EPMA_ALL_filter, aes(x = FeO_n, y = MgO_n, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1.5, alpha = 1)  +
  labs(y=expression('MgO (wt %)'), x= expression('FeO (wt %)')) +
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12), 
        legend.justification = c(1,1),
        legend.position = "bottom", 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=16, colour = "black"), axis.title=element_text(size=16, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
p7

# K vs Na
p8g <- ggplot(data = EPMA_ALL_filter, aes(x = Na2O_n, y = K2O_n))
p8 <- p8g + geom_point(data = EPMA_ALL_filter, aes(x = Na2O_n, y = K2O_n, colour = Eruption_ID, fill = Eruption_ID, shape = Eruption_ID), size = 1.5, alpha = 1)  +
  labs(y=expression(K[2]*O ~ ('wt%')), x=expression(Na[2]*O ~ ('wt%'))) +
  scale_shape_manual(values=erup_symbols) + # 21:25 are outline = colour and fill = fill shapes
  scale_fill_manual(values = erup_col) +
  scale_color_manual(values = erup_col) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12), 
        legend.justification = c(1,1),
        legend.position = "bottom", 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=16, colour = "black"), axis.title=element_text(size=16, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
p8


# Arrange plots into a grid  ----------------------------------------------

# Figure S9 --------------------------------------------------------------

ggarrange(p1, # First row with line plot
          # Second row with labels
          ggarrange(p3, p2, ncol = 2, labels = c("B", "C")), 
          nrow = 2, common.legend = TRUE,  legend = "top", 
          labels = "A" # Label of the first plot
) 

# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS9.pdf", height = c(30), width = c(20), dpi = 600, units = "cm")

# alternative plotting - no zoom RHY plot included 
#ggarrange(p1, p4, common.legend = TRUE, align = c("hv"),legend = "bottom", labels = c("A", "B", "C","D", "E", "F"), font.label = list(size = 16, color = "black"),ncol = 1, nrow = 2) 

# Figure S10 --------------------------------------------------------------
ggarrange(p5, p6, p7, p8, 
          common.legend = TRUE, 
          align = c("hv"),
          legend = "top", 
          labels = c("A", "B","C", "D"),
          font.label = list(size = 16, color = "black"),
          ncol = 2, nrow = 3) 

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS10.pdf", height = c(16), width = c(12), dpi = 600) #, units = "cm")


# Figure S11 - Correlation & PCA -----------------------------------------------------

# Data input and processing for PCA in base R-----------------------------------------------
#clear plot window
dev.off()
#reset parameters to default
.pardefault <- par()
par(.pardefault)

# Set up and select columns to use - as measured data
var_sub <- c("ID", "LongID", "SiO2",	"TiO2",	"Al2O3", "FeO",	"MnO", "MgO", "CaO", "Na2O",
             "K2O",	"Total", "Type", "Group",	"Eruption_ID",	"Volcano")
ME <- c("SiO2",	"TiO2",	"Al2O3", "FeO",	"MnO", "MgO", "CaO", "Na2O","K2O")

# or normalised data
var_sub_n <- c("ID", "LongID", "SiO2_n",	"TiO2_n",	"Al2O3_n", "FeO_n",	"MnO_n", "MgO_n", "CaO_n", "Na2O_n",
               "K2O_n",	"Total_n", "Type", "Group",	"Eruption_ID",	"Volcano")
ME_n <- c("SiO2_n",	"TiO2_n",	"Al2O3_n", "FeO_n",	"MnO_n", "MgO_n", "CaO_n", "Na2O_n", "K2O_n")

# Choose which dataset to import from TAS/bi-plotsfor PCA analysis on 
# 1) ALL DATA
kta <- EPMA_ALL_filter 
kta <- kta[,(var_sub_n),drop=FALSE]
kta

# OR

# 2) BS+INT filtered data to file
kta <-  EPMA_BSINT_filter.csv
kta <- kta[,(var_sub_n),drop=FALSE]
kta

# OR

# 3) ACID filtered data to file
kta <- EPMA_ACID_filter
kta <- kta[,(var_sub_n),drop=FALSE]
# remove ANT data - only 2 rows - not enough data for correlation and PCA in ACID dataset
kta <- kta[!kta$Group=='ANT',]
kta

# OR

# 4) RHY filtered data to file - choose this one for the 
kta <- EPMA_RHY_filter
kta <- kta[,(var_sub_n),drop=FALSE]
# remove ANT and DI data - only 2 rows of each- not enough data for correlation and PCA in ACID dataset
kta <- kta[!kta$Group=='ANT',]
kta <- kta[!kta$Group=='DI',]
kta

head(kta)
plot(kta[, ME_n], pch=19, cex = 0.05)

# Transform data ----------------------------------------------------------

#centre variables and  Z-scores using scale() function - for standardised (scaled) and centred PCA analysis - values are mean of 0 and +/-1 of 1 std dev
#copy the original filename to a new filename and then apply a z-score transform to it, plot it
kta.Z <- kta
kta.Z[, ME_n] <- scale(kta[ME_n], center = TRUE, scale = TRUE)
kta.Z[is.na(kta.Z)] <- 0
kta.Z
head(kta.Z)
plot(kta.Z[, ME_n], pch=19, cex = 0.05)

# Transform the %kta data by sqrt = %kta sq root-transformed, cuberoot = %TSBN cuberoot transformed
kta.sqrt <- kta
kta.sqrt[, ME_n] <- sqrt (kta.sqrt[ME_n])
kta.sqrt[is.na(kta.sqrt)] <- 0
head(kta.sqrt)
plot(kta.sqrt[, ME_n], pch=19, cex = 0.05)

# Calculate Z-scores of sqrt transformed data 
#copy the original sqrt filename to a new filename and then apply a z-score transform to it, plot it
kta.sqrt.Z <- kta.sqrt
kta.sqrt.Z[, ME_n] <- scale(kta.sqrt[ME_n], center = TRUE, scale = TRUE)
head(kta.sqrt.Z)
plot(kta.sqrt.Z[, ME_n], pch=19, cex = 0.05)

#  Write transformed data to file --------------------------------------------------
write.csv(kta.sqrt,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/kta_sqrt.csv", row.names = TRUE)
write.csv(kta.Z,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/kta_Z.csv", row.names = TRUE)
write.csv(kta.sqrt.Z,"/Users/Steve/Dropbox/BAS/Data/R/Roberts_et_al_2020_Pato/Outputs/kta_sqrt_Z.csv", row.names = TRUE)

# Correlation plots  --------------------------------------------

#Correlation plot - all elements, density plots and correlation matrxi summary
library(GGally)
# use this to see where positive/significant correlations as an overview overall
ggcorr(kta[, ME_n], method = c("everything", "pearson"), label = TRUE, label_alpha = TRUE, label_round=2) 

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/EPMA/Figure_corr_matrix.pdf", height = c(20), width = c(20), dpi = 600, units = "cm")

#then run this with +correlations only 
ggpairs(kta, columns = ME_n, upper = list(continuous = wrap("cor", size = 5)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation plot") 

## save plot in working directory
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/EPMA/Figure_corr-den_matrix.pdf", height = c(20), width = c(20), dpi = 600, units = "cm")

# Then run this with +correlations only - will plot scatterplot, density dist and stats for each unit
# limit this to 5 x 5 matrix as text difficult to read 
# and assign the matrix plot to a variable and then iterating through each plot with a for-loop
theme_set(theme_bw(base_size=8))
plot.matrix <- ggpairs(kta, columns = ME_n, upper = list(continuous = wrap("cor", size = 2.5)),
                       lower = list(continuous = wrap("points", alpha = 1, size=0.5)),
                       ggplot2::aes(colour = Eruption_ID, title="Correlation plot by Group"))

plot.matrix

# STOP & WAIT FOR 5 secs - set up loop to change colours to RdYlBu
for(i in 1:9) {
  for(j in 1:9){
    plot.matrix[i,j] <- plot.matrix[i,j] +
      scale_fill_brewer(palette = "RdYlBu", direction=-1) +
      scale_color_brewer(palette = "RdYlBu", direction=-1)
  }
}

plot.matrix

# save plot in working directory - A4 - square dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS5.pdf", height = c(21), width = c(21), dpi = 600, units = "cm")


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
pca.kta <- kta.sqrt.Z [ME_n]
head(pca.kta)
plot(pca.kta, pch=19, cex = 0.05)

# Make PCA not including Unit columns (1 in csv file) - depth is ID
p <- prcomp(pca.kta,  scale = TRUE)
# The coordinates for PCA are contained within the prcomp object as a matrix:res.pca$x

# Make PC variable matrix 
r = eigen(cov(pca.kta))
p$rotation

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

# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS11A.pdf", height = c(16), width = c(16), dpi = 600, units = "cm")

# Clear plot window and plot a Screeplot of PC Inertia
if(!is.null(dev.list())) dev.off()
layout(matrix(1:1, ncol=1))
p_inertia <- screeplot(p)

# Use the scree plots above and below to assess how many PC to keep - usually PC1, 2, 3. 

# Show data being used run
PC14 <- p$x[,1:4]
head(PC14)
#Show summary stats for variables using standardize package - and standardize data (to finish)
x.stats <- summary(PC14)
x.stats

#write PC data and stats to csv to use later if needed
write.csv(pca.kta,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/kta_sqrt_pca.csv", row.names = TRUE)
write.csv(PC14,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/PC1-4.csv", row.names = TRUE)
write.csv(var_percent,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/PC1-4_var_percent.csv", row.names = TRUE)
write.csv(x.stats,"/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Outputs/PC1-4_stats.csv", row.names = TRUE)

# Plot PC axes ------------------------------------------------------------

# Plot the first 2-4 principal components, so use p$x[,1] (PCA 1) and p$x[,2] (PCA 2) etc. to select the relevant data. 
# Plot datapoints to check - use change axes to = 1,3 for PC1 vs PC3
# cos2 = quality of individual points on the factor map
# contib = relative contribution of individual datapoints to PCA score
fviz_pca_ind(p, geom = c("point", "text"), axes = c(1,2), labelsize = 2, col.ind = "cos2", 
             gradient.cols = c("grey90","#2E9FDF", "#FC4E07"),
             repel = FALSE, title = "PCA - Biplot with depths labelled") + 
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# PC axes loadings --------------------------------------------------------

#Create a barplot for loadings on PC1 axis - can change this 
if(!is.null(dev.list())) dev.off()
#Loadings plots for PC1 only
barplot(p$rotation[,1], main="PC 1 Loadings Plot", las=2)
#the loadings for the variables for our first principal component (PC) are p$rotation[,1], while the second principal component loadings are p$rotation[,2] and so on
n.pc1 <- ifelse(p$rotation[,1] > 0, yes=-0.01, no=p$rotation[,1]-0.01)
#fix the labelling - create a numeric object ("n.pc1") containing values which will position text underneath the bars in the bar chart using the ifelse() function. The first argument p$rotation[,1] > 0 sets the condition or test for ifelse() to meet. The yes argument will return the value specified if the condition is met, while the no argument will return the value specified if it fails to meet the condition.
c.pc1 <- ifelse(p$rotation[,1] > 0, yes="blue", no="red2")
#add colours and replot
par(mar=c(8,3,2,1)) # Set margins
b1 <- barplot(p$rotation[,1], main="PC 1 Loadings Plot", col=c.pc1, las=2, axisnames=FALSE)
abline(h=0) # Add horizontal line
text(x=b1, y=n.pc1, labels=names(p$rotation[,1]), adj=1, srt=90, xpd=FALSE) # Add variable names
#use axisnames=FALSE to stop barplot() plotting them automatically
#this creates the bar chart as an object so that we can extract the "midpoints" of each bar to correctly position the variable names when using the text() function.
#adj=1 sets the alignment of the variable names to right align.
#srt=90 changes the direction of the text to a 90 degree angle (vertical).
#xpd=TRUE tells R that it can plot the text outside the plot region, and within the figure region


# Plot all PC axes loadings together --------------------------------------

if(!is.null(dev.list())) dev.off()
# Change colour of bar plot
c.pc1 <- ifelse(p$rotation[,1] > 0, yes="blue", no="red2")
c.pc2 <- ifelse(p$rotation[,2] > 0, "blue", "red2")
c.pc3 <- ifelse(p$rotation[,3] > 0, "blue", "red2")
c.pc4 <- ifelse(p$rotation[,4] > 0, "blue", "red2")
# Get position for variable names
n.pc1 <- ifelse(p$rotation[,1] > 0, -0.01, p$rotation[,1]-0.01)
n.pc2 <- ifelse(p$rotation[,2] > 0, -0.01, p$rotation[,2]-0.01)
n.pc3 <- ifelse(p$rotation[,3] > 0, -0.01, p$rotation[,3]-0.01)
n.pc4 <- ifelse(p$rotation[,4] > 0, -0.01, p$rotation[,4]-0.01)

# plot PC1 and PC2 together

# Plot PC1-4 loadings - 8 cm wide plot ------------------------------------
layout(matrix(1:4, ncol=1)) # Set up layout
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149608 inches but R converts to 3.149606 - 
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.2,0.2), pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 
#dev.off() #need to include this to write to pdf file fully - will delete screen plot
#par(mar=c(2,4,4,2), oma=c(8,0,0,0)) # Set up margins

# Plot PC 1, 2, 3 & 4 together add variable names
b1 <- barplot(p$rotation[,1], main="PC 1 Loadings Plot", col=c.pc1, las=3, axisnames=FALSE)
abline(h=0)
text(x=b1, y=ifelse(p$rotation[,1] > 0, -0.01, p$rotation[,1]-0.01), labels=names(p$rotation[,1]), adj=1, srt=90, xpd=NA)
b2 <- barplot(p$rotation[,2], main="PC 2 Loadings Plot", col=c.pc2, las=3, axisnames=FALSE)
abline(h=0)
text(x=b2, y=ifelse(p$rotation[,2] > 0, -0.01, p$rotation[,2]-0.01), labels=names(p$rotation[,2]), adj=1, srt=90, xpd=NA)
b3 <- barplot(p$rotation[,3], main="PC 3 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,3] > 0, -0.01, p$rotation[,3]-0.01), labels=names(p$rotation[,3]), adj=1, srt=90, xpd=NA)
b4 <- barplot(p$rotation[,4], main="PC 4 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,4] > 0, -0.01, p$rotation[,4]-0.01), labels=names(p$rotation[,4]), adj=1, srt=90, xpd=NA)

# Open a pdf file to save to file - sort out the margins on this plot
pdf("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS11B.pdf") 
# 2. Create a plot
# Plot PC1-4 loadings - 8 cm wide plot ------------------------------------
layout(matrix(1:4, ncol=1)) # Set up layout
## Translate 8 cm graph plot size *from* cm *to* inches:
plotinch <- 8 / cm(1) # -> 8 cm  is  3.149608 inches but R converts to 3.149606 - 
#par(mai=c(4,4,3,4), pin=c(plotinch, plotinch/2), xaxs='i') #plot size in inches
par(mai=c(0.1,0.1,0.2,0.2), pin=c(plotinch, plotinch/2),  mgp=c(2,0.5,0), xaxs='i', omi=c(0,0,0,0)) 
#dev.off() #need to include this to write to pdf file fully - will delete screen plot
#par(mar=c(2,4,4,2), oma=c(8,0,0,0)) # Set up margins

# Plot PC 1, 2, 3 & 4 together add variable names
b1 <- barplot(p$rotation[,1], main="PC 1 Loadings Plot", col=c.pc1, las=3, axisnames=FALSE)
abline(h=0)
text(x=b1, y=ifelse(p$rotation[,1] > 0, -0.01, p$rotation[,1]-0.01), labels=names(p$rotation[,1]), adj=1, srt=90, xpd=NA)
b2 <- barplot(p$rotation[,2], main="PC 2 Loadings Plot", col=c.pc2, las=3, axisnames=FALSE)
abline(h=0)
text(x=b2, y=ifelse(p$rotation[,2] > 0, -0.01, p$rotation[,2]-0.01), labels=names(p$rotation[,2]), adj=1, srt=90, xpd=NA)
b3 <- barplot(p$rotation[,3], main="PC 3 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,3] > 0, -0.01, p$rotation[,3]-0.01), labels=names(p$rotation[,3]), adj=1, srt=90, xpd=NA)
b4 <- barplot(p$rotation[,4], main="PC 4 Loadings Plot", col=c.pc3, las=3, axisnames=FALSE)
abline(h=0)
text(x=b3, y=ifelse(p$rotation[,4] > 0, -0.01, p$rotation[,4]-0.01), labels=names(p$rotation[,4]), adj=1, srt=90, xpd=NA)

# Close the pdf file
dev.off() 
#save this to file

# Plot variable lines and data on chosen PC axes biplots ----------------------------------
fviz_pca_var(p, axes = c(3,2), col.var = "contrib",
             gradient.cols = c("brown", "green"), #"Blues",
             ggtheme = theme_minimal(),
             title = "PCA - Variable Plot")

# Final graph -------------------------------------------------------------

# Biplot of individuals and variables
# Keep only the labels for variables
# Change the color by groups, add ellipses
pca1 <- fviz_pca_biplot(p, axes = c(3,2), label = "var", col.var = "blue", 
                        #gradient.cols = c("brown", "blue"), #"Blues",
                        habillage=kta$Eruption_ID,
                        addEllipses=TRUE, ellipse.level=0.68, 
                        palette = erup_col,
                        #c("blue", "grey50", "black", "salmon"),
                        title = "PC2/PC3 and Variables (68% CI)",
                        ggtheme = theme_classic()) + 
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

pca1
# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS11C.pdf", height = c(16), width = c(16), dpi = 600, units = "cm")

# Plot Scree and PCA  ------------------------------------

# Figure 11B - not used 
ggarrange(pca_scree, pca1,
          common.legend = FALSE, 
          align = c("h"),
          legend = "right", 
          labels = c("A", "B"),
          font.label = list(size = 16, color = "black"),
          ncol = 1, nrow = 2) 
# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS11B.pdf", height = c(29.7), width = c(21), dpi = 600, units = "cm")


# Cluster analysis --------------------------------------------------------

#clear plot window
dev.off()
#reset parameters to default
.pardefault <- par()
par(.pardefault)

# k-means clustering - Dim 1 vs Dim 2 ------------------------------------------------------

# Visualize kmeans clustering using repel = TRUE to avoid overplotting
kta.sqrt.scaled <- kta.sqrt.Z[ME_n]
head(kta.sqrt.scaled)
km.res <- kmeans(kta.sqrt.scaled, 5, nstart = 5)
fviz_cluster(km.res, geom = c("point"), kta.sqrt.scaled, ellipse.type = "euclid",
             palette = "RdRlBu", ggtheme = theme_minimal())

# Hierachical clustering - Vegan - Dim plot ------------------------------------------

# Use hcut() which computes hclust and cut the tree
hc.cut <- hcut(kta.sqrt.scaled, k = 5, hc_method = "wards") #can be complete, single, average or wards
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = TRUE, labelsize = 2, rect = TRUE) + 
  theme(text = element_text(size = 7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
fviz_cluster(hc.cut, ellipse.type = "convex")


# Hierarchical clustering in dendextend and gplot  -----------------------------

#clear plot window
dev.off()
#reset parameters to default
.pardefault <- par()
par(.pardefault)

# determine the number of clusters to use 
# optimal number of clusters usng average silhouette method (https://uc-r.github.io/kmeans_clustering#silo) 
fviz_nbclust(kta.sqrt.Z[, ME_n], FUN = hcut, method = "silhouette")

# Figure 12A ---------------------------------------------------------------

# USE THIS ONE - optimal number of clusters using Gap Statistic Method - bootstrapping MCMC approach (Tibshirani et al., 2001)
gap_stat <- clusGap(kta.sqrt.Z[, ME_n], FUN = kmeans, nstart = 25,
                    K.max = 15, B = 50)
p9 <- fviz_gap_stat(gap_stat) 
p9 +theme(text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)
)

# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS12A.pdf", height = c(12), width = c(12), dpi = 600, units = "cm")

# Figure 7B ---------------------------------------------------------------

#clear plot window
dev.off()

# Color as a function of the cluster - using hc as defined above - need to define colours to match number of groups
par(mar=c(4,1,1,7))
# Run clustering with just data matrix
kta.matrix <- as.matrix(kta.sqrt.Z[, ME_n])
head(kta.matrix)
# Link row names from original matrix
rownames(kta.matrix) <- kta.sqrt.Z$LongID
dend <- as.dendrogram(hclust(d = dist(x = kta.matrix), "ward.D2"))

#par(mfrow = c(1,2), mar = c(2,2,1,0))
dend <- dend %>%
  color_branches(k = 11) %>%
  set("branches_lwd", 0.75) %>%
  set("branches_lty", 1) %>%
  set("labels_cex", 0.05) %>%
  set("labels_col", value = BrBG2, k=11) %>%
  set("branches_k_color", value = BrBG2, k=10)
# use this if havent set the colours: dend <- color_labels(dend, k = 10)
# code is the same as:labels_colors(dend)  <- get_leaves_branches_col(dend)
p10 <- plot(dend, horiz=TRUE, axes=TRUE)

#get a list of cluster group members
dend_list <- get_subdendrograms(dend, 11)
lapply(dend_list, labels)

#10 colour scheme 
#c("skyblue", "steelblue3", "blue3", "burlywood3", "chocolate4", "salmon2", "red", "darkred", "darkorange",  "darkmagenta")

# save plot in working directory - A4 dimensions
ggsave("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_et_al_2020_Pato/Figures/FigureS12B.pdf", height = c(29.7), width = c(21), dpi = 600, units = "cm")

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

