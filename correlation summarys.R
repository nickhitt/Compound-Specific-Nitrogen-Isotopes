## Correlation Analysis of CSIAA and bulk N data

library(dplyr)
library(tidyr)
library(readxl)
library(ggpubr)
library(ggplot2)

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

## plotting boxplot of all bulk and amino acid data by coral
trophic <- c("ala", "Asp", "Glu", "Ile", "Leu", "Pro", "Val")
source <- c("Gly", "Phe")

all_corals_transpose <- all_corals_transpose %>% 
  dplyr::select(sample.id, depth, Value, Isotope, Coral.ID) %>%
  mutate(Coral_name = case_when(Coral.ID == "M" ~ "47996",
                          Coral.ID == "T" ~ "35104",
                          Coral.ID == "P" ~ "64344",
                          Coral.ID == "F" ~ "15131")) %>%
  mutate(type = case_when(Isotope %in% source ~ "Source",
                           Isotope %in% trophic ~ "Trophic")) %>%
  mutate(location = case_when(Coral.ID == "M" ~ "STF",
                              Coral.ID == "T" ~ "BOP",
                              Coral.ID == "P" ~ "BOP",
                              Coral.ID == "F" ~ "STF"))

boxplots <- all_corals_transpose %>%
  filter(Isotope != "Thr") %>%
  ggplot(aes(Isotope, Value, color = type) )+ geom_boxplot() +
  facet_wrap(Coral_name~location)

boxplots

boxplots_35104 <- all_corals_transpose[,which(all_corals_transpose$Coral.ID == "T")] %>%
  ggplot(aes(Isotope, Value)) + geom_boxplot() 




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
