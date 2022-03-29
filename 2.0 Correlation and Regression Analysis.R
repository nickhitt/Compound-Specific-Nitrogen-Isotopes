## Final Correlation Analysis of CSIAA and bulk N data

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

## Function for Correlation

csiaa_corr <- function(coral_data) {
  all_aa <- c("Ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Lys", "Phe", "Pro", "Thr", "Val")
  bulk <- coral_data[,1] ## pass only bulk through AA data only
  coral_data <- coral_data[,2:12]
  mat <- matrix(nrow = 11, ncol = 2)
  for (i in 1:11) {
    res <- cor.test(bulk, as.numeric(coral_data[,i]), 
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

source("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/1.0 Calculations - TP and SigmaV.R")
all_corals_bulkAA <- all_corals_Niso %>%
  dplyr::select(-sample.id, -depth, -age)

## Subsetting M

coral_m <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "M") %>%
  dplyr::select(-Coral)

## Subsetting P

coral_p <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "P") %>%
  dplyr::select(-Coral)

## Subsetting T

coral_t <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "T") %>%
  dplyr::select(-Coral)

## Subsetting F

coral_f <- all_corals_bulkAA %>%
  dplyr::filter(Coral == "F") %>%
  dplyr::select(-Coral)

## Correlations
corr_all <- csiaa_corr(all_corals_bulkAA)
corr_m <- csiaa_corr(coral_m)
corr_p <- csiaa_corr(coral_p)
corr_t <- csiaa_corr(coral_t)
corr_f <- csiaa_corr(coral_f)


#################### Regression Figures by AA

## Plotting function
regress_plot <- function(data, amino){
  regress_plot <- data %>%
    ggplot(mapping = aes(bulk.15n, {{amino}}, group = Coral)) +
    facet_wrap(~ Coral, scales = "free") + 
    geom_point() + geom_smooth(method = "lm")
}

# making figures
all_aa <- c("ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Lys", "Phe", "Pro", "Thr", "Val")
ala <- regress_plot(all_corals_bulkAA, ala)

asp <- regress_plot(all_corals_bulkAA, Asp)

glu <- regress_plot(all_corals_bulkAA, Glu)

gly <- regress_plot(all_corals_bulkAA, Gly)

ile <- regress_plot(all_corals_bulkAA, Ile)

leu <- regress_plot(all_corals_bulkAA, Leu)

lys <- regress_plot(all_corals_bulkAA, Lys)

phe <- regress_plot(all_corals_bulkAA, Phe)

pro <- regress_plot(all_corals_bulkAA, Pro)

thr <- regress_plot(all_corals_bulkAA, Thr)

val <- regress_plot(all_corals_bulkAA, Val)

#################### Regression Figures by Coral

## Formatting Data

data_reshaped <- reshape(all_corals_bulkAA, 
        varying = c(all_aa), 
        v.names = "Value",
        timevar = c("AA"), 
        times = c(all_aa), 
        direction = "long") %>%
  dplyr::select(-id) %>%
  dplyr::mutate(Value = round(as.numeric(Value), 2))

row.names(data_reshaped) <- 1:nrow(data_reshaped)

## Y Axis Rounding Function
scaleFUN <- function(x) {sprintf("%.2f", x)}

## Plotting function
regress_plot_coral <- function(data, coral){
  regress_plot <- data %>%
    dplyr::filter(Coral == coral) %>%
    ggplot(mapping = aes(bulk.15n, Value, group = AA)) +
    facet_wrap(~ AA, scales = "free") + 
    geom_point() + geom_smooth(method = "lm") +
    theme(text = element_text(size=20))
}

p <- regress_plot_coral(data_reshaped, "P")
p

ggsave("Bulk vs. AAs - 64344", device = png, dpi = 320, 
       path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Updated Regression Figures")

t <- regress_plot_coral(data_reshaped, "T")
t

ggsave("Bulk vs. AAs - 35104", device = "png", dpi = 320, 
       path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Updated Regression Figures")

m <- regress_plot_coral(data_reshaped, "M")
m

ggsave("Bulk vs. AAs - 47996", device = "png", dpi = 320, 
       path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Updated Regression Figures")


f <- regress_plot_coral(data_reshaped, "F")
f

ggsave("Bulk vs. AAs - 15131", device = "png", dpi = 320, 
       path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Updated Regression Figures")


all_corals <- data_reshaped %>%
  ggplot(mapping = aes(bulk.15n, Value, group = AA)) +
  facet_wrap(~ AA, scales = "free") + 
  geom_point() + geom_smooth(method = "lm") +
  theme(text = element_text(size=20))
all_corals

ggsave("Bulk vs. AAs - All Corals", device = "png", dpi = 320, 
       path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Updated Regression Figures")




