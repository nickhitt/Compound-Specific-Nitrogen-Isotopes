### Trophic Position Comparison

library(dplyr)
library(ggplot2)
library(readxl)
library(stats)
library(RColorBrewer)

## loading data

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N")

## Calculating TP

all_corals_Niso <- data.frame(read_excel("csiaadataclean_with_lys.xlsx", sheet = "All Data")) %>%
  dplyr::mutate(TP = 1 + ((Glu + 3.4) - Phe - 3.4)/7.6) 

## Calculating Sigma V -- NEED TO FINISH

sigmaV <- function(dataset, coral){
  mat <- c()
  data <- dataset %>%
    dplyr::filter(Coral == coral)  %>%
    dplyr::mutate(ala_1 = abs(ala - mean(ala)), 
                  asp_1 = abs(Asp - mean(Asp)),
                  glu_1 = abs(Glu - mean(Glu)),
                  ile_1 = abs(Ile - mean(Ile)),
                  leu_1 = abs(Leu - mean(Leu)),
                  pro_1 = abs(Pro - mean(Pro)),
                  sigmaV = 1/(6)*(ala_1 + asp_1 + glu_1 +
                                    ile_1 + leu_1 + pro_1)) %>%
    dplyr::select(-ala_1, -asp_1, -glu_1, -ile_1, -leu_1, -pro_1)
}

sigmaV_m <- sigmaV(all_corals_Niso, "M")
sigmaV_p <- sigmaV(all_corals_Niso, "P")
sigmaV_t <- sigmaV(all_corals_Niso, "T")
sigmaV_f <- sigmaV(all_corals_Niso, "F")

## Plotting TP

p1 <- all_corals_Niso %>%
  ggplot(mapping = aes(sample.id, TP, group = Coral)) +
  facet_wrap(~ Coral) +
  geom_point() + geom_smooth(method = "lm")
p1

corrs_m <- all_corals_Niso %>%
  dplyr::filter(Coral == "M") 

corrs_m <- cor.test(corrs_m$Phe, corrs_m$TP)
corrs_m_bulk <- cor.test(corrs_m$bulk.15n, corrs_m$TP)

corrs_p <- all_corals_Niso %>%
  dplyr::filter(Coral == "P") 

corrs_p <- cor.test(corrs_p$Phe, corrs_p$TP)
corrs_p_bulk <- cor.test(corrs_p$bulk.15n, corrs_p$TP)

corrs_t <- all_corals_Niso %>%
  dplyr::filter(Coral == "T") 

corrs_t <- cor.test(corrs_t$Phe, corrs_t$TP)
corrs_t_bulk <- cor.test(corrs_t$bulk.15n, corrs_t$TP)


names <- c("35104", "47996", "64344")
cors <- c(corrs_t$estimate, corrs_m$estimate, corrs_p$estimate)
pvals <- c(corrs_t$p.value, corrs_m$p.value, corrs_p$p.value)

corr_matrix <- cbind(names, cors, pvals)

cors <- c(corrs_t_bulk$estimate, corrs_m_bulk$estimate, corrs_p_bulk$estimate)
pvals <- c(corrs_t_bulk$p.value, corrs_m_bulk$p.value, corrs_p_bulk$p.value)

corr_matrix_bulk_tp <- cbind(names, cors, pvals)

