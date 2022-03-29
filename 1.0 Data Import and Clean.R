## Loading Packages

library(dplyr)
library(readxl)
library(data.table)   

setwd("~/Dropbox/Marsden Black Coral Project/R Codes/Bacon/New Age Models")

raw <- data.frame(read_excel("~/Dropbox/Marsden Black Coral Project/R Codes/Bacon/New Age Models/Raw Dates.xlsx")) %>%
  dplyr::mutate(Depth = Depth * 0.05) %>%
  dplyr::rename(lab_ID = Sample) %>%
  dplyr::rename(Age = AgeCorr.1950.BP) %>%
  dplyr::rename(error = X2s) %>%
  dplyr::rename(depth = Depth) %>%
  dplyr::arrange(lab_ID, Age, error, depth)

f_coral <- raw[raw$Sample %like% "F", ]
m_coral <- raw[raw$Sample %like% "M", ]
p_coral <- raw[raw$Sample %like% "P", ]
t_coral <- raw[raw$Sample %like% "T", ]

write.csv(f_coral, file = "F.csv")
write.csv(m_coral, file = "M.csv")
write.csv(p_coral, file = "P.csv")
write.csv(t_coral, file = "T.csv")

