#' Supplementary code to the article:
#' Klinkovska K, Sperandii MG, Travnicek B, Chytry M (2023). 
#' Significant decline in habitat specialists in semi-dry grasslands over four decades. 
#' Biodiversity and Conservation.
#' 
#' author Klara Klinkovska, klinkovska.klara@gmail.com
#' R version 4.3.0

library(patchwork) # version 1.1.2
library(vegan) # version 2.6-4
library(ggrepel) # version 0.9.3
library(corrplot) # version 0.92
library(tidyverse) # version 2.0.0
#' dplyr v. 1.1.2
#' ggplot2 v. 3.4.2
#' readr v. 2.1.4
#' stringr v. 1.5.0
#' tibble v. 3.2.1

# species data + square-root transformation
spe <- read_csv("data/Klinkovska_et_al_semi_dry_grasslands_S_Moravia_species.csv") %>% 
  select(-Releve_number) %>% 
  sqrt()

# species characteristics
spe.traits <- read_csv("data/Klinkovska_et_al_semi_dry_grasslands_S_Moravia_species_data.csv") %>% 
  mutate(across(c(starts_with("Ellenberg"), starts_with("dist")), ~as.numeric(.)))

# header data 
head <- read_csv("data/Klinkovska_et_al_semi_dry_grasslands_S_Moravia_head.csv")  

# unconstrained ordination ------------------------------------------------

pcoa <- capscale(spe ~ 1, distance = "bray", sqrt.dist = T)

head.ord <- bind_cols(head, as_tibble(scores(pcoa)$sites))

# ordination plot -----------------------------------------------------------

plot1 <- head.ord %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = as.factor(Rs_observ),
                 shape = as.factor(Rs_observ)),
             size = 3) +
  scale_shape_manual(values = 21:22, 
                     labels=c("1980s", "2022")) +
  scale_fill_manual(values = c(4, 7), labels=c("1980s", "2022"))+
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        plot.title = element_text(size = 11))+
  geom_path(aes(MDS1, MDS2, group=Rs_plot), arrow = arrow(length=unit(0.3, "cm"))) +
  ggtitle("(a) 1980s-2022, PCoA, sites")+
  labs(x="Axis 1", y="Axis 2")

# ordinations with species unconstrained ----------------------------------

# scores for species
spe.fit.pcoa <- envfit(pcoa, spe[colSums(spe!=0)>5]) 

spe.pcoa <- scores(spe.fit.pcoa, "vectors") %>%
  as_tibble(rownames = "species") %>% 
  mutate(MDS1=MDS1*3, 
         MDS2=MDS2*3,
         p = spe.fit.pcoa$vectors$pvals, 
         r = spe.fit.pcoa$vectors$r, 
         short_name = str_c(str_split_i(species, "\\s", 1) %>% str_sub(., 1, 3), 
                            str_split_i(species, "\\s", 2) %>% str_sub(., 1, 3), sep = "."))

# ordination plot with species
plot1.spe <- spe.pcoa %>% 
  filter(p < 0.05) %>% 
  slice_max(r, n = 80) %>% 
  ggplot(aes(MDS1, MDS2)) + 
  geom_text_repel(aes(label = short_name), max.overlaps = Inf,
                  alpha = 0.6, cex = 3.5)+
  theme_bw()+
  theme(plot.title = element_text(size=11))+
  ggtitle("(b) 1980s-2022, PCoA, species")+
  labs(x="Axis1", y="Axis2")

# bind together ordination plots for sites and species
plot1 + plot1.spe
ggsave("plots/Fig_3_ordination.png", width = 10, height = 5)

# constrained ordination --------------------------------------------------

dord <- capscale(spe ~ head$Rs_observ+Condition(head$Rs_plot), distance = "bray", sqrt.dist = T)

# significance test
anova(dord, permutations=how(blocks=as.factor(head$Rs_plot), nperm=999))

# significantly increasing and decreasing species
dord.spe <- envfit(dord, spe, choices = 1, 
                    permutations = how(blocks=as.factor(head$Rs_plot), nperm=999), display="lc")

sig.spe <- as_tibble(dord.spe$vectors$arrows, rownames = "species") %>% 
  bind_cols(p = dord.spe$vectors$pvals, r = dord.spe$vectors$r) %>% 
  filter(p < 0.05) %>% 
  arrange(CAP1, p)

# correlation of species scores with their Ellenberg-type indicator values
spe.score <- as.data.frame(scores(dord, display = "species"))

# difference through time
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Light), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Temperature), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Moisture), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Soil_Reaction), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Nutrients), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$dist_freq), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$dist_sev), method = "spearman")

# testing differences between protected and non-protected -----------------

# constrained ordination
dord.prot <- capscale(spe ~ head$protected*head$Rs_observ + Condition(head$protected+head$Rs_observ), 
                  distance = "bray", sqrt.dist = T)

# significance test
anova(dord.prot, permutations=how(plots=Plots(as.factor(head$Rs_plot), type="free"), 
                              within=Within("none"), nperm=999))

# correlations of species scores with their indicator values
spe.score <- as.data.frame(scores(dord.prot, display = "species"))

cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Light), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Temperature), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Moisture), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Soil_Reaction), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$Ellenberg_Nutrients), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$dist_freq), method = "spearman")
cor.test(spe.score$CAP1, as.numeric(spe.traits$dist_sev), method = "spearman")

# correlations between EIVs -------------------------------------------------
corr <- cor(spe.traits %>% select(starts_with("Ellenberg"), starts_with("dist")), 
              method = "spearman", use = "pairwise.complete.obs")

rownames(corr) <- c("Light", "Temperature", "Moisture", "Soil Reaction", 
                    "Nutrients", "Disturbance frequency", "Disturbance severity")

colnames(corr) <- c("Light", "Temperature", "Moisture", "Soil Reaction", 
                    "Nutrients", "Disturbance frequency", "Disturbance severity")

corrplot(corr, method = "color", type = "upper", addCoef.col = "black", diag = F, 
         tl.col = "black", tl.srt = 45)

png(filename="plots/Appendix_S5_corrplot.png", width=8.5, height=6, units="cm", res=1024)
corrplot(corr, method = "color", type = "upper", addCoef.col = "black", diag = F, 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5, cl.cex = 0.5, number.cex = 0.5)

dev.off()

