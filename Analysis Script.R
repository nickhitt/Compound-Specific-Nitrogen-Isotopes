## Correlation Analysis of CSIAA and bulk N data

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

## reading data
all_corals <- data.frame(read_excel("csiaadataclean.xlsx", sheet = "All Data"))
all_corals_transpose <- data.frame(read_excel("csiaadataclean_transpose.xlsx", sheet = "All Data"))
north <- data.frame(read_excel("csiaadataclean.xlsx", sheet = "North"))
east <- data.frame(read_excel("csiaadataclean.xlsx", sheet = "East"))
t <- read.csv("csiaadataclean.csv")
m <- data.frame(read_excel("csiaadataclean.xlsx", sheet = "47996"))
p <- data.frame(read_excel("csiaadataclean.xlsx", sheet = "64344"))
all_corals_bulk <- data.frame(read_excel("AllCoralBulk.xlsx", sheet = "Sheet1"))
all_corals_bulkAA <- data.frame(read_excel("AllCoralBulkandAA.xlsx", sheet = "Sheet1"))


## plotting boxplot of all bulk and amino acid data by coral
trophic <- c("ala", "Asp", "Glu", "Ile", "Leu", "Pro", "Val")
source <- c("Gly", "Phe")
bulk <- c("bulk 15n")

all_corals_transpose <- all_corals_transpose %>% 
  dplyr::select(sample.id, depth, Value, Isotope, Coral.ID, Time) %>%
  mutate(Coral_name = case_when(Coral.ID == "M" ~ "47996",
                          Coral.ID == "T" ~ "35104",
                          Coral.ID == "P" ~ "64344",
                          Coral.ID == "F" ~ "15131")) %>%
  mutate(type = case_when(Isotope %in% source ~ "Source",
                          Isotope %in% trophic ~ "Trophic",
                          Isotope %in% bulk ~ "Bulk")) %>%
  mutate(location = case_when(Coral.ID == "M" ~ "STF",
                              Coral.ID == "T" ~ "BOP",
                              Coral.ID == "P" ~ "BOP",
                              Coral.ID == "F" ~ "STF"))

## Plotting coral isotope boxplots by coral and location

boxplots_corals <- all_corals_transpose %>%
  filter(Isotope != "Thr") %>%
  ggplot(aes(Isotope, Value, color = type) )+ geom_boxplot() +
  facet_wrap(Coral_name~location)

boxplots_corals
# Run to save figure
 #ggsave("AAll Amino Acid Boxplots By Coral.png", boxplots_corals,
       #path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures")

## Plotting boxplots of coral isotopes by isotope

boxplots_isotope <- all_corals_transpose %>%
  filter(Isotope != "Thr") %>%
  ggplot(aes(Coral_name, Value, color = Coral_name) )+ geom_boxplot() +
  facet_grid(~Isotope) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

boxplots_isotope
# Run to save figure
 #ggsave("Individual Amino Acid Boxplots By Coral.png", boxplots_corals,
        #path = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures")


#############3 Plotting individual amino acids against bulk by coral

####### Ala

## 35104
timeseries_ala1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	ala) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_ala1,
                          type = "l", ylab = "Bulk 15n")
timeseries_ala_BOP2 <- xyplot(ala ~ Time, timeseries_ala1,
                              type = "p", ylab = "Ala 15N")
ala_35104 <- doubleYScale(time_series_bulk, timeseries_ala_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_ala2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	ala) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_ala2,
                           type = "l", ylab = "Bulk 15n")
timeseries_ala_STF2 <- xyplot(ala ~ Time, timeseries_ala2,
                              type = "p", ylab = "Ala 15N")
ala_15131 <- doubleYScale(time_series_bulk, timeseries_ala_STF2, add.ylab2 = TRUE)


## 47996
timeseries_ala3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	ala) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_ala3,
                           type = "l", ylab = "Bulk 15n")
timeseries_ala_STF1 <- xyplot(ala ~ Time, timeseries_ala3,
                              type = "p", ylab = "Ala 15N")
ala_47996 <- doubleYScale(time_series_bulk, timeseries_ala_STF1, add.ylab2 = TRUE)

## 64344
timeseries_ala4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	ala) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_ala4,
                           type = "l", ylab = "Bulk 15n")
timeseries_ala_EAC1 <- xyplot(ala ~ Time, timeseries_ala4,
                              type = "p", ylab = "Ala 15N")
ala_64344 <- doubleYScale(time_series_bulk, timeseries_ala_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

ala_plots <- grid.arrange(ala_64344, ala_47996, 
             ala_15131, ala_35104, ncol=2)

######## Asp
## 35104
timeseries_asp1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Asp) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_asp1,
                           type = "l", ylab = "Bulk 15n")
timeseries_asp_BOP2 <- xyplot(Asp ~ Time, timeseries_asp1,
                              type = "p", ylab = "Asp 15N")
asp_35104 <- doubleYScale(time_series_bulk, timeseries_asp_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_asp2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Asp) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_asp2,
                           type = "l", ylab = "Bulk 15n")
timeseries_asp_STF2 <- xyplot(Asp ~ Time, timeseries_asp2,
                              type = "p", ylab = "Asp 15N")
asp_15131 <- doubleYScale(time_series_bulk, timeseries_asp_STF2, add.ylab2 = TRUE)


## 47996
timeseries_asp3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Asp) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_asp3,
                           type = "l", ylab = "Bulk 15n")
timeseries_asp_STF1 <- xyplot(Asp ~ Time, timeseries_asp3,
                              type = "p", ylab = "Asp 15N")
asp_47996 <- doubleYScale(time_series_bulk, timeseries_asp_STF1, add.ylab2 = TRUE)

## 64344
timeseries_asp4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Asp) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_asp4,
                           type = "l", ylab = "Bulk 15n")
timeseries_asp_EAC1 <- xyplot(Asp ~ Time, timeseries_asp4,
                              type = "p", ylab = "Asp 15N")
asp_64344 <- doubleYScale(time_series_bulk, timeseries_asp_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

asp_plots <- grid.arrange(asp_64344, asp_47996, 
             asp_15131, asp_35104, ncol=2)

########### Glu
## 35104
timeseries_Glu1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Glu) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Glu1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Glu_BOP2 <- xyplot(Glu ~ Time, timeseries_Glu1,
                              type = "p", ylab = "Glu 15N")
Glu_35104 <- doubleYScale(time_series_bulk, timeseries_Glu_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Glu2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Glu) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Glu2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Glu_STF2 <- xyplot(Glu ~ Time, timeseries_Glu2,
                              type = "p", ylab = "Glu 15N")
Glu_15131 <- doubleYScale(time_series_bulk, timeseries_Glu_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Glu3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Glu) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Glu3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Glu_STF1 <- xyplot(Glu ~ Time, timeseries_Glu3,
                              type = "p", ylab = "Glu 15N")
Glu_47996 <- doubleYScale(time_series_bulk, timeseries_Glu_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Glu4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Glu) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Glu4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Glu_EAC1 <- xyplot(Glu ~ Time, timeseries_Glu4,
                              type = "p", ylab = "Glu 15N")
Glu_64344 <- doubleYScale(time_series_bulk, timeseries_Glu_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

glu_plots <- grid.arrange(Glu_64344, Glu_47996, 
             Glu_15131, Glu_35104, ncol=2)
#Gly
## 35104
timeseries_Gly1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Gly) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Gly1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Gly_BOP2 <- xyplot(Gly ~ Time, timeseries_Gly1,
                              type = "p", ylab = "Gly 15N")
Gly_35104 <- doubleYScale(time_series_bulk, timeseries_Gly_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Gly2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Gly) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Gly2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Gly_STF2 <- xyplot(Gly ~ Time, timeseries_Gly2,
                              type = "p", ylab = "Gly 15N")
Gly_15131 <- doubleYScale(time_series_bulk, timeseries_Gly_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Gly3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Gly) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Gly3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Gly_STF1 <- xyplot(Gly ~ Time, timeseries_Gly3,
                              type = "p", ylab = "Gly 15N")
Gly_47996 <- doubleYScale(time_series_bulk, timeseries_Gly_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Gly4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Gly) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Gly4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Gly_EAC1 <- xyplot(Gly ~ Time, timeseries_Gly4,
                              type = "p", ylab = "Gly 15N")
Gly_64344 <- doubleYScale(time_series_bulk, timeseries_Gly_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

gly_plots <- grid.arrange(Gly_64344, Gly_47996, 
             Gly_15131, Gly_35104, ncol=2)

#Ile
## 35104
timeseries_Ile1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Ile) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Ile1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Ile_BOP2 <- xyplot(Ile ~ Time, timeseries_Ile1,
                              type = "p", ylab = "Ile 15N")
Ile_35104 <- doubleYScale(time_series_bulk, timeseries_Ile_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Ile2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Ile) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Ile2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Ile_STF2 <- xyplot(Ile ~ Time, timeseries_Ile2,
                              type = "p", ylab = "Ile 15N")
Ile_15131 <- doubleYScale(time_series_bulk, timeseries_Ile_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Ile3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Ile) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Ile3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Ile_STF1 <- xyplot(Ile ~ Time, timeseries_Ile3,
                              type = "p", ylab = "Ile 15N")
Ile_47996 <- doubleYScale(time_series_bulk, timeseries_Ile_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Ile4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Ile) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Ile4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Ile_EAC1 <- xyplot(Ile ~ Time, timeseries_Ile4,
                              type = "p", ylab = "Ile 15N")
Ile_64344 <- doubleYScale(time_series_bulk, timeseries_Ile_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

ile_plots <- grid.arrange(Ile_64344, Ile_47996, 
             Ile_15131, Ile_35104, ncol=2)
#Leu
## 35104
timeseries_Leu1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Leu) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Leu1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Leu_BOP2 <- xyplot(Leu ~ Time, timeseries_Leu1,
                              type = "p", ylab = "Leu 15N")
Leu_35104 <- doubleYScale(time_series_bulk, timeseries_Leu_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Leu2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Leu) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Leu2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Leu_STF2 <- xyplot(Leu ~ Time, timeseries_Leu2,
                              type = "p", ylab = "Leu 15N")
Leu_15131 <- doubleYScale(time_series_bulk, timeseries_Leu_STF2, add.ylab2 = TRUE)

## 47996
timeseries_Leu3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Leu) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Leu3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Leu_STF1 <- xyplot(Leu ~ Time, timeseries_Leu3,
                              type = "p", ylab = "Leu 15N")
Leu_47996 <- doubleYScale(time_series_bulk, timeseries_Leu_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Leu4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Leu) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Leu4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Leu_EAC1 <- xyplot(Leu ~ Time, timeseries_Leu4,
                              type = "p", ylab = "Leu 15N")
Leu_64344 <- doubleYScale(time_series_bulk, timeseries_Leu_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

leu_plots <- grid.arrange(Leu_64344, Leu_47996, 
             Leu_15131, Leu_35104, ncol=2)

#Phe
## 35104
timeseries_Phe1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Phe) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Phe1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Phe_BOP2 <- xyplot(Phe ~ Time, timeseries_Phe1,
                              type = "p", ylab = "Phe 15N")
Phe_35104 <- doubleYScale(time_series_bulk, timeseries_Phe_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Phe2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Phe) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Phe2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Phe_STF2 <- xyplot(Phe ~ Time, timeseries_Phe2,
                              type = "p", ylab = "Phe 15N")
Phe_15131 <- doubleYScale(time_series_bulk, timeseries_Phe_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Phe3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Phe) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Phe3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Phe_STF1 <- xyplot(Phe ~ Time, timeseries_Phe3,
                              type = "p", ylab = "Phe 15N")
Phe_47996 <- doubleYScale(time_series_bulk, timeseries_Phe_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Phe4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Phe) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Phe4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Phe_EAC1 <- xyplot(Phe ~ Time, timeseries_Phe4,
                              type = "p", ylab = "Phe 15N")
Phe_64344 <- doubleYScale(time_series_bulk, timeseries_Phe_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

phe_plots <- grid.arrange(Phe_64344, Phe_47996, 
             Phe_15131, Phe_35104, ncol=2)

### Bay of Plenty
timeseries_Phe5 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Phe) %>%
  filter(Coral == "64344" | Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Phe5,
                           type = "l", ylab = "Bulk 15n")
timeseries_Phe_EAC3 <- xyplot(Phe ~ Time, timeseries_Phe5,
                              type = "p", ylab = "Phe 15N")
Phe_BOP <- doubleYScale(time_series_bulk, timeseries_Phe_EAC3, add.ylab2 = TRUE)


#Pro
## 35104
timeseries_Pro1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Pro) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Pro1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Pro_BOP2 <- xyplot(Pro ~ Time, timeseries_Pro1,
                              type = "p", ylab = "Pro 15N")
Pro_35104 <- doubleYScale(time_series_bulk, timeseries_Pro_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Pro2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Pro) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Pro2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Pro_STF2 <- xyplot(Pro ~ Time, timeseries_Pro2,
                              type = "p", ylab = "Pro 15N")
Pro_15131 <- doubleYScale(time_series_bulk, timeseries_Pro_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Pro3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Pro) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Pro3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Pro_STF1 <- xyplot(Pro ~ Time, timeseries_Pro3,
                              type = "p", ylab = "Pro 15N")
Pro_47996 <- doubleYScale(time_series_bulk, timeseries_Pro_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Pro4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Pro) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Pro4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Pro_EAC1 <- xyplot(Pro ~ Time, timeseries_Pro4,
                              type = "p", ylab = "Pro 15N")
Pro_64344 <- doubleYScale(time_series_bulk, timeseries_Pro_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure

pro_plots <- grid.arrange(Pro_64344, Pro_47996, 
             Pro_15131, Pro_35104, ncol=2)

#Val
## 35104
timeseries_Val1 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Val) %>%
  filter(Coral == "35104")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Val1,
                           type = "l", ylab = "Bulk 15n")
timeseries_Val_BOP2 <- xyplot(Val ~ Time, timeseries_Val1,
                              type = "p", ylab = "Val 15N")
Val_35104 <- doubleYScale(time_series_bulk, timeseries_Val_BOP2, add.ylab2 = TRUE)

## 15131

timeseries_Val2 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Val) %>%
  filter(Coral == "15131")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Val2,
                           type = "l", ylab = "Bulk 15n")
timeseries_Val_STF2 <- xyplot(Val ~ Time, timeseries_Val2,
                              type = "p", ylab = "Val 15N")
Val_15131 <- doubleYScale(time_series_bulk, timeseries_Val_STF2, add.ylab2 = TRUE)


## 47996
timeseries_Val3 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Val) %>%
  filter(Coral == "47996")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Val3,
                           type = "l", ylab = "Bulk 15n")
timeseries_Val_STF1 <- xyplot(Val ~ Time, timeseries_Val3,
                              type = "p", ylab = "Val 15N")
Val_47996 <- doubleYScale(time_series_bulk, timeseries_Val_STF1, add.ylab2 = TRUE)

## 64344
timeseries_Val4 <- all_corals_bulkAA %>%
  dplyr::select(Time, bulk_15_n, Coral,	Val) %>%
  filter(Coral == "64344")

time_series_bulk <- xyplot(bulk_15_n ~ Time, timeseries_Val4,
                           type = "l", ylab = "Bulk 15n")
timeseries_Val_EAC1 <- xyplot(Val ~ Time, timeseries_Val4,
                              type = "p", ylab = "Val 15N")
Val_64344 <- doubleYScale(time_series_bulk, timeseries_Val_EAC1, add.ylab2 = TRUE)

## Combining Figures into 1 Figure
png(file = "~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N/Figures/Val_Figures.png")

val_plots <- grid.arrange(Val_64344, Val_47996, 
             Val_15131, Val_35104, ncol=2)
dev.off()
    

#### Principal Components Analysis

## 35104
pca_35104 <- all_corals_bulkAA %>%
  filter(Coral == "35104") %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

pca.35104 <- prcomp(pca_35104, center = TRUE,scale. = TRUE)

summary(pca.35104)

ggbiplot(pca.35104)

## 64344
pca_64344 <- all_corals_bulkAA %>%
  filter(Coral == "64344") %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

pca.64344 <- prcomp(pca_64344, center = TRUE,scale. = TRUE)

summary(pca.64344)

ggbiplot(pca.64344)

## 47996
pca_47996 <- all_corals_bulkAA %>%
  filter(Coral == "47996") %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

pca.47996 <- prcomp(pca_47996, center = TRUE,scale. = TRUE)

summary(pca.47996)

ggbiplot(pca.47996)

## 15131
pca_15131 <- all_corals_bulkAA %>%
  filter(Coral == "15131") %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

pca.15131 <- prcomp(pca_15131, center = TRUE,scale. = TRUE)

summary(pca.15131)

ggbiplot(pca.15131)

## All Corals
pca_all <- all_corals_bulkAA %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

pca.all <- prcomp(pca_all, center = TRUE,scale. = TRUE)

summary(pca.all)

ggbiplot(pca.all)

###### Linear Model Fits

### 35104
# All AAs
lm_35104 <- lm(bulk_15_n ~ ala + Asp + Glu + Gly + Ile + Leu + Phe + Pro + Val,
   pca_35104)
summary(lm_35104)

# Only Bulk
lm_b_35104 <- all_corals_bulkAA %>%
  filter(Coral == "35104") %>%
  select(-Depth, -Uncertainty, -Coral) 

lm_35104_bulk <- lm(Time ~ bulk_15_n, lm_b_35104)
summary(lm_35104_bulk)

lm_source_35104 <- all_corals_bulkAA %>%
  filter(Coral == "35104") %>%
  select(-Depth, -Uncertainty, -Coral, Time, Phe, Gly) 

lm_35104_phe <- lm(Time ~ Phe + Gly,
               lm_source_35104)
summary(lm_35104_phe)


### 47996
lm_47996 <- lm(bulk_15_n ~ ala + Asp + Glu + Gly + Ile + Leu + Phe + Pro + Val,
               pca_47996)
summary(lm_47996)

### 64344
lm_64344 <- lm(bulk_15_n ~ ala + Asp + Glu + Gly + Ile + Leu + Phe + Pro + Val,
               pca_64344)
summary(lm_64344)


#### Bay of Plenty
lm_bop_dat <- all_corals_bulkAA %>%
  filter(Coral == "64344" | Coral == "35104") %>%
  select(-Depth, -Uncertainty, -Coral, -Time) %>%
  filter(!is.na(Phe))

lm_bop <- lm(bulk_15_n ~ ala + Asp + Glu + Gly + Ile + Leu + Phe + Pro + Val,
               lm_bop_dat)
summary(lm_bop)





## function for correlation matrix

coral_correl <- function(data){
  data <- data[,4:14]
  cors <- c()
  unc_cor <- c()
  slope <- c()
  p <- c()
  r2 <- c()
  for (i in 1:10) {
    data_1 <- c()
    data_2 <- c()
    data_1 <- data[,1]
    data_2 <- data[,i+1]
    dat <- data.frame(data_1, data_2)
    d <- cor.test(data_1, data_2, method = "pearson")
    l <- lm(data_1 ~ data_2, dat)
    l <- summary(l)
    slope[i] <- l$coefficients[2]
    p[i] <- l$coefficients[2,4]
    r2[i] <- l$r.squared
    cors[i] <- d$estimate
    unc_cor[i] <- d$p.value
  }
  name <- c("ala", "asp", "Glu",	"Gly",	"Ile",	"Leu",	"Phe",	"Pro",	"Thr",	"Val")
  cbind(name,cors,unc_cor, slope, p, r2)
}
  
p_correls <- coral_correl(p)
m_correls <- coral_correl(m)
t_correls <- coral_correl(t)
east_correls <- coral_correl(east)
north_correls <- coral_correl(north)
all_correls <- coral_correl(all_corals)


##Estimating Changes in N_fix

nfix <- function(data){
  bulkdata <- data[,4]
  glydata <- data[,8]
  n_fix <- c()
  for (i in 1:length(bulkdata)){
    n_fix[i] <- 1 - ((glydata[i] - -1)/(6.5 - glydata[i]))
  }
  n_fix <- (n_fix/n_fix[1])-1
}

nfix_t <- nfix(t)
nfix_m <- nfix(m)
nfix_p <- nfix(p)
nfix_east <- nfix(east)
nfix_north <- nfix(north)
#nfix_all <- nfix(all)
all <- data.frame(nfix_t,nfix_m,nfix_p,nfix_north,nfix_east)


## Estimating Changes in Trophic Level

ntroph <- function(data){
  gludata <- data[,7]
  phedata <- data[,11]
  n_troph <- c()
  for (i in 1:length(gludata)){
    n_troph[i] <- (((gludata[i] - phedata[i])-3.4)/7.6)+1
  }
}

ntroph_t <- ntroph(t)
ntroph_m <- ntroph(m)
ntroph_p <- ntroph(p)
ntroph_east <- ntroph(east)
ntroph_north <- ntroph(north)
