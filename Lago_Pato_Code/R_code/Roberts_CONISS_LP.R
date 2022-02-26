
#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Load libaries -----------------------------------------------------------

library(ggplot2)
library(ggsci)
library(ggpubr)
library(ellipse)  # for PCA and cluster
library(dplyr)  # for transforming data
library(factoextra) # for PCA and cluster
library(reshape2)
library(cowplot) # for plotting
library(vegan) # for paleo analysis
library(rioja) # for paleo analysis
library(mgcv)
library(tidyr)
library(tidypaleo)

# Set working directory ---------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Roberts_2021")
#check working directory
getwd()

# CONISS - CONSTRAINED CLUSTER ANALYSIS -----------------------------------

# Import pre-transformed data file from this folder
# To run quickly / smoothly - smooth data to 1 or 2 mm from 200 um
# https://www.dropbox.com/sh/ubzleyty1kt63zx/AACyCkY3KRkhKjHwInPBjWZRa?dl=0 

x <- read.csv("Data/CONISS/LP08_2013_2mm_TSN_sqrt_CONISS.csv", 
              header=TRUE, row.names=1, sep=",", check.names=FALSE)

# OR

x <- read.csv("Data/CONISS/LP08_2018_2mm_TSN_sqrt_CONISS.csv", 
              header=TRUE, row.names=1, sep=",", check.names=FALSE)

# OR

x <- read.csv("Data/CONISS/LP16_2020_1mm_TSN_sqrt_CONISS.csv", 
              header=TRUE, row.names=1, sep=",", check.names=FALSE)
head(x)

##create % sum dataset if not as % data already or sums to 100% if not already 
totals<-apply(x, 1, sum)
x.pc<-x/totals*100
mx<-apply(x.pc, 2, max)
x.pc.5<-x.pc[, mx >0.5]

#plots % stratigraphic data using first column as either depth or age
depth <-as.numeric(rownames(x))
strat.plot(x.pc.5, y.rev=TRUE, cex.lab=0.5)
strat.plot(x.pc.5, y.rev=TRUE, scale.percent=TRUE, cex.lab=0.5)
strat.plot(x.pc.5, y.rev=TRUE, scale.percent=TRUE, yvar=depth, cex.lab=0.5)

# Run stratigraphically constrained cluster analysis (after Grimm 1987) 
# by the method of incremental sum of squares with Hellingers transformation
# based on Borcard D., Gillet F., Legendre P. (2018) Cluster Analysis. 
# In: Numerical Ecology with R. Use R!. Springer for details 
# doi.org/10.1007/978-3-319-71404-2_4

# Hellinger transformation is used to balance the effects of rare taxa. 
# Constrained Incremental Sum of Squares (CONISS) and broken stick modelling identifies the maximum number 
# of sediment units that significantly differ from a random model
x.hel<-decostand(x, method="hellinger")
diss<-dist(x.hel)
clust<-chclust(diss, method="coniss")

# broken stick analysis on dataset and plots results - ng=number of groups on x axis
# where the lines first meet or cross is the statistically signifcant 
# number of groups - in this example = 4 groups
# not always that clear with a large number of groups in, for example, ITRAX data - smooth data first
# and/or where sum of squares first major dip occurs?

bstick(clust, ng=20) #plots the broken stick graph - with 10 groups shown - can change this if needed

#plots strat plot with CONISS cluster zones
x <-strat.plot(x.pc.5, y.rev=TRUE, scale.percent=TRUE, yvar=depth,clust=clust)

# defines and adds the number of cluster zones
# need to enter this data based on broken stick graph after clust,
# provides depth or age of the zone boundary - can be imported into C2

# 5 groups for LP08 2013 dataset
# 7 groups for LP08 2018 dataset
# 11 groups for LP16 2020 dataset

addClustZone(x, clust, 11, col="blue")
which(diff(cutree(clust,k=11))>0)


# TO DO -------------------------------------------------------------------

# https://easystats.github.io/parameters/reference/n_clusters.html

#Number of clusters to extract
#Source: R/n_clusters.R
#This function runs many existing procedures for determining how many clusters are present in data. 
#It returns the number of clusters based on the maximum consensus. In case of ties, it will select the solution with fewer clusters.


n_clusters(
  x,
  standardize = TRUE,
  force = FALSE,
  package = c("NbClust", "mclust", "cluster", "M3C"),
  fast = TRUE,
  ...
)
