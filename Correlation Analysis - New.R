## New Correlation Analysis of CSIAA and bulk N data

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

## setting directory

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N")

## Functions

csiaa_corr <- function(coral_data) {
  all_aa <- c("Ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Lys", "Phe", "Pro", "Thr", "Val")
  bulk <- coral_data[,1]
  coral_data <- coral_data[,2:11]
  mat <- matrix(nrow = 10, ncol = 2)
  for (i in 1:10) {
    res <- cor.test(bulk, coral_data[,i], 
                    method = "pearson")
    mat[i,1] <- res$estimate
    mat[i,2] <- res$p.value
    rm(res)
  }
  mat <- as.data.frame(mat)
  mat <- cbind(all_aa, mat) %>%
    dplyr::mutate(significance = case_when(V2 < 0.1 ~ "significant"))
  return(mat)
}

## Reading Data
all_corals_bulkAA <- data.frame(read_excel("AllCoralBulkandAA.xlsx", sheet = "Sheet1")) %>%
  dplyr::filter(!is.na(Phe)) %>%
  dplyr::select(-Depth, -Time, -Uncertainty)

## Subsetting M

coral_m <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "47996") %>%
  dplyr::select(-Coral)

## Subsetting P

coral_p <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "64344") %>%
  dplyr::select(-Coral)

## Subsetting T

coral_t <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "35104") %>%
  dplyr::select(-Coral)

## Subsetting All Corals

all_corals_corr <- all_corals %>%
  dplyr::select(-sample.id, -depth, -Coral)

## Correlations
corr_all <- csiaa_corr(all_corals_corr)
corr_m <- csiaa_corr(coral_m)
corr_p <- csiaa_corr(coral_p)
corr_t <- csiaa_corr(coral_t)
