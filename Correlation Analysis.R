## Correlation Analysis of CSIAA and bulk C data

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

## Reading Data
all_corals_bulkAA <- data.frame(read_excel("AllCoralBulkandAA.xlsx", sheet = "Sheet1"))
all_corals <- data.frame(read_excel("csiaadataclean_corr.xlsx", sheet = "All Data"))
all_corals_transpose <- data.frame(read_excel("csiaadataclean_transpose.xlsx", sheet = "All Data"))

## Bulk, Essential and Non-Essential IDs 
trophic <- c("Ala", "Asp", "Glu", "Ile", "Leu", "Pro", "Val")
source <- c("Gly", "Phe")
bulk <- c("bulk15n")
all_aa <- c("Ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Phe", "Pro", "Val")

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

# all_corals <- all_corals %>% 
#   dplyr::select(sample.id, depth, Value, Coral) %>%
#   mutate(Coral_name = case_when(Coral == "M" ~ "47996",
#                                 Coral == "T" ~ "35104",
#                                 Coral == "P" ~ "64344",
#                                 Coral == "F" ~ "15131")) %>%
#   mutate(location = case_when(Coral == "M" ~ "STF",
#                               Coral == "T" ~ "BOP",
#                               Coral == "P" ~ "BOP",
#                               Coral == "F" ~ "STF"))

## All Corals
data_all <- all_corals %>%
  subset(Coral == "M" | Coral == "P" | Coral == "T") 

lm_ala_all <- lm(bulk15n ~ Ala,
                 data_all)

ala_all <- ggplot(data_all, aes(bulk15n, Ala)) +
  geom_point() + geom_smooth(method = "lm") 

lm_asp_all <- lm(bulk15n ~ Asp ,
                 data_all)
summary(lm_asp_all)

asp_all <- ggplot(data_all, aes(bulk15n, Asp)) +
  geom_point() + geom_smooth(method = "lm")

lm_glu_all <- lm(bulk15n ~ Glu,
                 data_all)
summary(lm_glu_all)

glu_all <- ggplot(data_all, aes(bulk15n, Glu)) +
  geom_point() + geom_smooth(method = "lm")

lm_gly_all <- lm(bulk15n ~ Gly,
                 data_all)
summary(lm_gly_all)

gly_all <- ggplot(data_all, aes(bulk15n, Gly)) +
  geom_point() + geom_smooth(method = "lm")

lm_ile_all <- lm(bulk15n ~ Ile,
                 data_all)
summary(lm_ile_all)

ile_all <- ggplot(data_all, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

lm_leu_all <- lm(bulk15n ~ Leu,
                 data_all)
summary(lm_leu_all)

leu_all <- ggplot(data_all, aes(bulk15n, Leu)) +
  geom_point() + geom_smooth(method = "lm")

lm_phe_all <- lm(bulk15n ~ Phe,
                 data_all)
summary(lm_phe_all)

phe_all <- ggplot(data_all, aes(bulk15n, Phe)) +
  geom_point() + geom_smooth(method = "lm")

lm_pro_all <- lm(bulk15n ~ Pro,
                 data_all)
summary(lm_pro_all)

pro_all <- ggplot(data_all, aes(bulk15n, Pro)) +
  geom_point() + geom_smooth(method = "lm")

lm_val_all <- lm(bulk15n ~ Val,
                 data_all)
summary(lm_val_all)

val_all <- ggplot(data_all, aes(bulk15n, Val)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_all, asp_all, glu_all, gly_all,
             ile_all, leu_all, phe_all, pro_all,
             val_all, nrow = 3, ncol = 3)

r2_all <- c(summary(lm_ala_all)$r.squared, summary(lm_asp_all)$r.squared,
            summary(lm_glu_all)$r.squared, summary(lm_gly_all)$r.squared,
            summary(lm_ile_all)$r.squared, summary(lm_leu_all)$r.squared,
            summary(lm_phe_all)$r.squared, summary(lm_pro_all)$r.squared,
            summary(lm_val_all)$r.squared) %>%
  as.data.frame()

## BOP Corals

data_BOP <- all_corals %>%
  subset(Coral == "T" | Coral == "P") 

ala_BOP <- ggplot(data_BOP, aes(bulk15n, Ala)) +
  geom_point() + geom_smooth(method = "lm")

asp_BOP <- ggplot(data_BOP, aes(bulk15n, Asp)) +
  geom_point() + geom_smooth(method = "lm")

glu_BOP <- ggplot(data_BOP, aes(bulk15n, Glu)) +
  geom_point() + geom_smooth(method = "lm")

gly_BOP <- ggplot(data_BOP, aes(bulk15n, Gly)) +
  geom_point() + geom_smooth(method = "lm")

ile_BOP <- ggplot(data_BOP, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

leu_BOP <- ggplot(data_BOP, aes(bulk15n, Leu)) +
  geom_point() + geom_smooth(method = "lm")

ile_BOP <- ggplot(data_BOP, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

phe_BOP <- ggplot(data_BOP, aes(bulk15n, Phe)) +
  geom_point() + geom_smooth(method = "lm")

pro_BOP <- ggplot(data_BOP, aes(bulk15n, Pro)) +
  geom_point() + geom_smooth(method = "lm")

val_BOP <- ggplot(data_BOP, aes(bulk15n, Val)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_BOP, asp_BOP, glu_BOP, gly_BOP,
             ile_BOP, leu_BOP, phe_BOP, pro_BOP,
             val_BOP, nrow = 3, ncol = 3)

## 35104 and 47996

data_eac_stf <- all_corals %>%
  subset(Coral == "M" | Coral == "T") 

lm_ala_eac_stf <- lm(bulk15n ~ Ala,
                     data_eac_stf)

ala_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Ala, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm") 

lm_asp_eac_stf <- lm(bulk15n ~ Asp ,
                     data_eac_stf)
summary(lm_asp_eac_stf)

asp_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Asp, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_glu_eac_stf <- lm(bulk15n ~ Glu,
                     data_eac_stf)
summary(lm_glu_eac_stf)

glu_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Glu, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_gly_eac_stf <- lm(bulk15n ~ Gly,
                     data_eac_stf)
summary(lm_gly_eac_stf)

gly_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Gly, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_ile_eac_stf <- lm(bulk15n ~ Ile,
                     data_eac_stf)
summary(lm_ile_eac_stf)

ile_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Ile, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_leu_eac_stf <- lm(bulk15n ~ Leu,
                     data_eac_stf)
summary(lm_leu_eac_stf)

leu_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Leu, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_phe_eac_stf <- lm(bulk15n ~ Phe,
                     data_eac_stf)
summary(lm_phe_eac_stf)

phe_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Phe, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_pro_eac_stf <- lm(bulk15n ~ Pro,
                     data_eac_stf)
summary(lm_pro_eac_stf)

pro_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Pro, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

lm_val_eac_stf <- lm(bulk15n ~ Val,
                     data_eac_stf)
summary(lm_val_eac_stf)

val_eac_stf <- ggplot(data_eac_stf, aes(bulk15n, Val, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_eac_stf, asp_eac_stf, glu_eac_stf, gly_eac_stf,
             ile_eac_stf, leu_eac_stf, phe_eac_stf, pro_eac_stf,
             val_eac_stf, nrow = 3, ncol = 3)

r2_eac_stf <- c(summary(lm_ala_eac_stf)$r.squared, summary(lm_asp_eac_stf)$r.squared,
                summary(lm_glu_eac_stf)$r.squared, summary(lm_gly_eac_stf)$r.squared,
                summary(lm_ile_eac_stf)$r.squared, summary(lm_leu_eac_stf)$r.squared,
                summary(lm_phe_eac_stf)$r.squared, summary(lm_pro_eac_stf)$r.squared,
                summary(lm_val_eac_stf)$r.squared) %>%
  as.data.frame()

## 35104
data_35104 <- all_corals %>%
  subset(Coral == "T") 

lm_ala_35104 <- lm(bulk15n ~ Ala,
                   data_35104)

ala_35104 <- ggplot(data_35104, aes(bulk15n, Ala)) +
  geom_point() + geom_smooth(method = "lm") 

lm_asp_35104 <- lm(bulk15n ~ Asp ,
                   data_35104)
summary(lm_asp_35104)

asp_35104 <- ggplot(data_35104, aes(bulk15n, Asp)) +
  geom_point() + geom_smooth(method = "lm")

lm_glu_35104 <- lm(bulk15n ~ Glu,
                   data_35104)
summary(lm_glu_35104)

glu_35104 <- ggplot(data_35104, aes(bulk15n, Glu)) +
  geom_point() + geom_smooth(method = "lm")

lm_gly_35104 <- lm(bulk15n ~ Gly,
                   data_35104)
summary(lm_gly_35104)

gly_35104 <- ggplot(data_35104, aes(bulk15n, Gly)) +
  geom_point() + geom_smooth(method = "lm")

lm_ile_35104 <- lm(bulk15n ~ Ile,
                   data_35104)
summary(lm_ile_35104)

ile_35104 <- ggplot(data_35104, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

lm_leu_35104 <- lm(bulk15n ~ Leu,
                   data_35104)
summary(lm_leu_35104)

leu_35104 <- ggplot(data_35104, aes(bulk15n, Leu)) +
  geom_point() + geom_smooth(method = "lm")

lm_phe_35104 <- lm(bulk15n ~ Phe,
                   data_35104)
summary(lm_phe_35104)

phe_35104 <- ggplot(data_35104, aes(bulk15n, Phe)) +
  geom_point() + geom_smooth(method = "lm")

lm_pro_35104 <- lm(bulk15n ~ Pro,
                   data_35104)
summary(lm_pro_35104)

pro_35104 <- ggplot(data_35104, aes(bulk15n, Pro)) +
  geom_point() + geom_smooth(method = "lm")

lm_val_35104 <- lm(bulk15n ~ Val,
                   data_35104)
summary(lm_val_35104)

val_35104 <- ggplot(data_35104, aes(bulk15n, Val)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_35104, asp_35104, glu_35104, gly_35104,
              ile_35104, leu_35104, phe_35104, pro_35104,
              val_35104, nrow = 3, ncol = 3)

r2_35104 <- c(summary(lm_ala_35104)$r.squared, summary(lm_asp_35104)$r.squared,
              summary(lm_glu_35104)$r.squared, summary(lm_gly_35104)$r.squared,
              summary(lm_ile_35104)$r.squared, summary(lm_leu_35104)$r.squared,
              summary(lm_phe_35104)$r.squared, summary(lm_pro_35104)$r.squared,
              summary(lm_val_35104)$r.squared) %>%
  as.data.frame()

## 47996

data_47996 <- all_corals %>%
  subset(Coral == "M") 

lm_ala_47996 <- lm(bulk15n ~ Ala,
                   data_47996)

ala_47996 <- ggplot(data_47996, aes(bulk15n, Ala)) +
  geom_point() + geom_smooth(method = "lm") 

lm_asp_47996 <- lm(bulk15n ~ Asp ,
                   data_47996)
summary(lm_asp_47996)

asp_47996 <- ggplot(data_47996, aes(bulk15n, Asp)) +
  geom_point() + geom_smooth(method = "lm")

lm_glu_47996 <- lm(bulk15n ~ Glu,
                   data_47996)
summary(lm_glu_47996)

glu_47996 <- ggplot(data_47996, aes(bulk15n, Glu)) +
  geom_point() + geom_smooth(method = "lm")

lm_gly_47996 <- lm(bulk15n ~ Gly,
                   data_47996)
summary(lm_gly_47996)

gly_47996 <- ggplot(data_47996, aes(bulk15n, Gly)) +
  geom_point() + geom_smooth(method = "lm")

lm_ile_47996 <- lm(bulk15n ~ Ile,
                   data_47996)
summary(lm_ile_47996)

ile_47996 <- ggplot(data_47996, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

lm_leu_47996 <- lm(bulk15n ~ Leu,
                   data_47996)
summary(lm_leu_47996)

leu_47996 <- ggplot(data_47996, aes(bulk15n, Leu)) +
  geom_point() + geom_smooth(method = "lm")

lm_phe_47996 <- lm(bulk15n ~ Phe,
                   data_47996)
summary(lm_phe_47996)

phe_47996 <- ggplot(data_47996, aes(bulk15n, Phe)) +
  geom_point() + geom_smooth(method = "lm")

lm_pro_47996 <- lm(bulk15n ~ Pro,
                   data_47996)
summary(lm_pro_47996)

pro_47996 <- ggplot(data_47996, aes(bulk15n, Pro)) +
  geom_point() + geom_smooth(method = "lm")

lm_val_47996 <- lm(bulk15n ~ Val,
                   data_47996)
summary(lm_val_47996)

val_47996 <- ggplot(data_47996, aes(bulk15n, Val)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_47996, asp_47996, glu_47996, gly_47996,
             ile_47996, leu_47996, phe_47996, pro_47996,
             val_47996, nrow = 3, ncol = 3)

r2_47996 <- c(summary(lm_ala_47996)$r.squared, summary(lm_asp_47996)$r.squared,
              summary(lm_glu_47996)$r.squared, summary(lm_gly_47996)$r.squared,
              summary(lm_ile_47996)$r.squared, summary(lm_leu_47996)$r.squared,
              summary(lm_phe_47996)$r.squared, summary(lm_pro_47996)$r.squared,
              summary(lm_val_47996)$r.squared) %>%
  as.data.frame()

## 64344
data_64344 <- all_corals %>%
  subset(Coral == "P") 

lm_ala_64344 <- lm(bulk15n ~ Ala,
                   data_64344)

ala_64344 <- ggplot(data_64344, aes(bulk15n, Ala)) +
  geom_point() + geom_smooth(method = "lm") 

lm_asp_64344 <- lm(bulk15n ~ Asp ,
                   data_64344)
summary(lm_asp_64344)

asp_64344 <- ggplot(data_64344, aes(bulk15n, Asp)) +
  geom_point() + geom_smooth(method = "lm")

lm_glu_64344 <- lm(bulk15n ~ Glu,
                   data_64344)
summary(lm_glu_64344)

glu_64344 <- ggplot(data_64344, aes(bulk15n, Glu)) +
  geom_point() + geom_smooth(method = "lm")

lm_gly_64344 <- lm(bulk15n ~ Gly,
                   data_64344)
summary(lm_gly_64344)

gly_64344 <- ggplot(data_64344, aes(bulk15n, Gly)) +
  geom_point() + geom_smooth(method = "lm")

lm_ile_64344 <- lm(bulk15n ~ Ile,
                   data_64344)
summary(lm_ile_64344)

ile_64344 <- ggplot(data_64344, aes(bulk15n, Ile)) +
  geom_point() + geom_smooth(method = "lm")

lm_leu_64344 <- lm(bulk15n ~ Leu,
                   data_64344)
summary(lm_leu_64344)

leu_64344 <- ggplot(data_64344, aes(bulk15n, Leu)) +
  geom_point() + geom_smooth(method = "lm")

lm_phe_64344 <- lm(bulk15n ~ Phe,
                   data_64344)
summary(lm_phe_64344)

phe_64344 <- ggplot(data_64344, aes(bulk15n, Phe)) +
  geom_point() + geom_smooth(method = "lm")

lm_pro_64344 <- lm(bulk15n ~ Pro,
                   data_64344)
summary(lm_pro_64344)

pro_64344 <- ggplot(data_64344, aes(bulk15n, Pro)) +
  geom_point() + geom_smooth(method = "lm")

lm_val_64344 <- lm(bulk15n ~ Val,
                   data_64344)
summary(lm_val_64344)

val_64344 <- ggplot(data_64344, aes(bulk15n, Val)) +
  geom_point() + geom_smooth(method = "lm")

grid.arrange(ala_64344, asp_64344, glu_64344, gly_64344,
             ile_64344, leu_64344, phe_64344, pro_64344,
             val_64344, nrow = 3, ncol = 3)

r2_64344 <- c(summary(lm_ala_64344)$r.squared, summary(lm_asp_64344)$r.squared,
              summary(lm_glu_64344)$r.squared, summary(lm_gly_64344)$r.squared,
              summary(lm_ile_64344)$r.squared, summary(lm_leu_64344)$r.squared,
              summary(lm_phe_64344)$r.squared, summary(lm_pro_64344)$r.squared,
              summary(lm_val_64344)$r.squared) %>%
  as.data.frame()