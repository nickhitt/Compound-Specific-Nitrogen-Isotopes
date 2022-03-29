### Trophic Position Comparison

library(dplyr)
library(ggplot2)
library(readxl)
library(stats)
library(RColorBrewer)

## setting directory

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/CSIAA_N")

## Calculating TP

all_corals_Niso <- data.frame(read_excel("csiaadataclean_with_lys.xlsx", sheet = "All Data")) %>%
  dplyr::mutate(TP = 1 + ((Glu + 3.4) - Phe - 3.4)/7.6) 

## Calculating Sigma V 

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

all_corals_Niso <- rbind(sigmaV_f, sigmaV_m, sigmaV_p, sigmaV_t)



