#### Timeseries Analysis

## Loading Packages

library(dplyr)
library(tidyr)
library(readxl)
library(ggpubr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(RColorBrewer)
library(gridExtra)
library(ggbiplot)
library(stats)

## setting directory

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N")

## Sourcing Data

source("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/1.0 Calculations - TP and SigmaV.R")
csiaadataclean_transpose <- data.frame(read_excel("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/AllCoralBulkandAA.xlsx")) %>%
  dplyr::select(Coral, bulk_15_n, Time) %>%
  dplyr::rename(bulk.15n = bulk_15_n, age = Time) %>%
  dplyr::mutate(Coral = ifelse(Coral == "64344", "P", Coral),
                Coral = ifelse(Coral == "35104", "T", Coral),
                Coral = ifelse(Coral == "15131", "F", Coral),
                Coral = ifelse(Coral == "47996", "M", Coral))
## Plotting Functions

# Making Individual Plots

timeseries_plot <- function(dataset_1, dataset_2, coral, AA){
  timeseries_1 <- dataset_1 %>%
    #dplyr::select(age, Coral,	AA) %>%
    filter(Coral == coral)
  timeseries_2 <- dataset_2 %>%
    dplyr::select(age, bulk.15n, Coral) %>%
    filter(Coral == coral)
  
  time_series_bulk <- xyplot(bulk.15n ~ age, timeseries_2,
                             type = "l", ylab = "Bulk 15n")
  timeseries_AA <- xyplot(get(AA) ~ age, timeseries_1,
                                type = "p", ylab = {{AA}})
  plot1 <- doubleYScale(time_series_bulk, timeseries_AA, add.ylab2 = TRUE)
}

# Grouping All Plots Together

allAA_plot <- function(dataset_1, dataset_2, coral) {
  ala_plot <- timeseries_plot(dataset_1, dataset_2, coral, "ala")
  asp_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Asp")
  glu_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Glu")
  gly_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Gly")
  ile_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Ile")
  leu_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Leu")
  lys_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Lys")
  phe_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Phe")
  pro_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Pro")
  thr_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Thr")
  val_plot <- timeseries_plot(dataset_1, dataset_2, coral, "Val")
  tp_plot <- timeseries_plot(dataset_1, dataset_2, coral, "TP")
  sigmaV_plot <- timeseries_plot(dataset_1, dataset_2, coral, "sigmaV")
  all_plots <- grid.arrange(ala_plot, asp_plot, glu_plot, gly_plot,
                            ile_plot, leu_plot, lys_plot, phe_plot,
                            pro_plot,  val_plot, tp_plot, sigmaV_plot,
                            ncol=4)
}

### Looking at all plots by coral

plots_64344 <- allAA_plot(all_corals_Niso, csiaadataclean_transpose, "P")

plots_47996 <- allAA_plot(all_corals_Niso, csiaadataclean_transpose, "M")

plots_35104 <- allAA_plot(all_corals_Niso, csiaadataclean_transpose, "T")

plots_15131 <- allAA_plot(all_corals_Niso, csiaadataclean_transpose, "F")

