### Lys Phe Comparison

library(dplyr)
library(ggplot2)
library(readxl)
library(stats)
library(RColorBrewer)

## loading data

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N")

all_corals_Niso <- data.frame(read_excel("csiaadataclean_with_lys.xlsx", sheet = "All Data"))

corrs_m <- all_corals_Niso %>%
  dplyr::filter(Coral == "M") 

corrs_m <- cor.test(corrs_m$Phe, corrs_m$Lys)

corrs_p <- all_corals_Niso %>%
  dplyr::filter(Coral == "P") 

corrs_p <- cor.test(corrs_p$Phe, corrs_p$Lys)

corrs_t <- all_corals_Niso %>%
  dplyr::filter(Coral == "T") 

corrs_t <- cor.test(corrs_t$Phe, corrs_t$Lys)

corrs_f <- all_corals_Niso %>%
  dplyr::filter(Coral == "F") 

corrs_f <- cor.test(corrs_f$Phe, corrs_f$Lys)

p1 <- all_corals_Niso %>%
  #dplyr::filter(Coral != "F") %>%
  ggplot(mapping = aes(Phe, Lys, group = Coral)) +
  facet_wrap(~ Coral, scales = "free") + 
  geom_point() + geom_smooth(method = "lm")
p1


p2 <- all_corals_Niso %>%
  #dplyr::filter(Coral != "F") %>%
  ggplot(mapping = aes(Phe, Lys)) +
  #facet_wrap(~ Coral, scales = "free") + 
  geom_point() + geom_smooth(method = "lm")
p2

model <- lm(Lys ~ Phe, data = all_corals_Niso, method = "qr")

corrs_all <- all_corals_Niso 

corrs_all <- cor.test(corrs_all$Phe, corrs_all$Lys)



names <- c("35104", "47996", "64344", "15131")
cors <- c(corrs_t$estimate, corrs_m$estimate, corrs_p$estimate, corrs_f$estimate)
pvals <- c(corrs_t$p.value, corrs_m$p.value, corrs_p$p.value, corrs_f$p.value)

corr_matrix <- cbind(names, cors, pvals)



boxplots_corals <- all_corals_Niso %>%
  ggplot(aes(Value, color = type) )+ geom_boxplot() +
  facet_wrap(~Coral)

boxplots_corals





