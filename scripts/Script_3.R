#' Supplementary code to the article:
#' Klinkovska K, Sperandii MG, Travnicek B, Chytry M (2023). 
#' Significant decline in habitat specialists in semi-dry grasslands over four decades. 
#' Biodiversity and Conservation.
#' 
#' author Klara Klinkovska, klinkovska.klara@gmail.com
#' R version 4.3.0

library(adespatial) # version 0.3-23
library (tidyverse) # version 2.0.0
#' dplyr v. 1.1.2
#' readr v. 2.1.4
#' tibble v. 3.2.1

# import data -------------------------------------------------------------

# header data 
head <- read_csv("data/Klinkovska_et_al_semi_dry_grasslands_S_Moravia_head.csv") 

head80 <- head %>% filter(Rs_observ == "1980")
head22 <- head %>% filter(Rs_observ == "2022")

# species data + square-root transformation
spe <- read_csv("data/Klinkovska_et_al_semi_dry_grasslands_S_Moravia_species.csv") 

spe80 <- spe %>% semi_join(head80) %>% sqrt() %>% select(-Releve_number) %>% as.data.frame()
spe22 <- spe %>% semi_join(head22) %>% sqrt() %>% select(-Releve_number) %>% as.data.frame()

# species exchange ratio --------------------------------------------------
tbi.jac <- TBI(spe80, spe22, method = "jaccard", test.t.perm = T, nperm = 999)
tbi.ruz <- TBI(spe80, spe22, method = "ruzicka", test.t.perm = T, nperm = 999)

# proportion of plots where dissimilarity due to species losses prevailed 
# over the dissimilarity due to species gains
table(tbi.jac$BCD.mat$Change)
(58/90)*100

table(tbi.ruz$BCD.mat$Change)
(47/90)*100
