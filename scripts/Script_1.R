#' Supplementary code to the article:
#' Klinkovska K, Sperandii MG, Trávníček B, Chytrý M (2023). 
#' Significant decline in habitat specialists in semi-dry grasslands over four decades. 
#' Biodiversity and Conservation.
#' 
#' author Klara Klinkovska, klinkovska.klara@gmail.com
#' R version 4.3.0

library(geepack) # version 1.3.9
library(nlme) # version 3.1-162
library(broom) # version 1.0.5
library(patchwork) # version 1.1.2
library(ggalluvial) # version 0.12.5
library(readxl) # version 1.4.2
library(tidyverse) # version 2.0.0
#' dplyr v. 1.1.2
#' ggplot2 v. 3.4.2
#' readr v. 2.1.4
#' stringr v. 1.5.0
#' tibble v. 3.2.1
#' tidyr v. 1.3.0

# import data ------------------------------------------------------------

# species data + species characteristics
spe.traits <- read_csv("data/Klinkovska_et_al_basiphilous_grasslands_S_Moravia_species.csv") %>%
  pivot_longer(cols = -Releve_number, values_to = "cover_perc", names_to = "SpecName_Layer") %>% 
  filter(cover_perc > 0) %>% 
  left_join(read_csv("data/Klinkovska_et_al_basiphilous_grasslands_S_Moravia_species_data.csv"))

# header data + calculate species richness, number of alien species etc.
head <- read_csv("data/Klinkovska_et_al_basiphilous_grasslands_S_Moravia_head.csv") %>% 
  left_join(spe.traits %>% group_by(Releve_number) %>% count(name = "spe.nr")) %>% 
  left_join(spe.traits %>% filter(!is.na(alien)) %>% group_by(Releve_number) %>% count(name = "alien.nr")) %>% 
  left_join(spe.traits %>% filter(str_detect(conserv_stat, "CR|VU|EN")) %>% group_by(Releve_number) %>% count(name = "end.nr")) %>% 
  left_join(spe.traits %>% filter(THE_THF == "TRUE") %>% group_by(Releve_number) %>% count(name = "diag.nr")) %>% 
  mutate(across(c(spe.nr, end.nr, alien.nr, diag.nr), ~ifelse(is.na(.), 0, .)))
  
# changes in vegetation types ---------------------------------------------
head %>% 
  filter(Rs_observ == "2022") %>% 
  mutate(across(c(ESY0, ESY), ~ifelse(str_detect(., "THE"), ., str_sub(., 1, 3)))) %>% 
  mutate(across(c(ESY0, ESY), ~str_replace(.x, "\\?|\\+|LCA", "Unclassified"))) %>% 
  mutate(ESY = factor(ESY, levels = c("TD", "TDA", "THD", "THE", "THE01", "THE03", 
                           "THF", "THH", "THG", "XCB", "XCC", "Unclassified")), 
         ESY0 = factor(ESY0, levels = c("TDA", "THD", "THE", "THE01", "THE03", 
                                      "THF", "THH", "Unclassified"))) %>% 
  select(Releve_number, ESY, ESY0) %>% 
  ggplot(aes(axis1 = ESY0, axis2 = as.factor(ESY)))+
  geom_alluvium(aes(fill = ESY0))+
  stat_stratum(width = 1/3)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()+
  scale_x_discrete(labels = c("1980s", "2022"))+
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.x = element_blank())

# alluvial plot was further modified in Inscape 1.2.2
ggsave("plots/Fig_2_alluvial.svg", width = 12, height = 8)

# species richness --------------------------------------------------------

# model 
m2 <- gls(spe.nr~Rs_observ*protected, cor=corCompSymm(form = ~1|Rs_plot), data = head) # with interaction
m2 <- gls(spe.nr~Rs_observ+protected, cor=corCompSymm(form = ~1|Rs_plot), data = head) # interaction not significant -> model without interaction
summary(m2)
plot(m2)
qqnorm(m2)

# confidence intervals
both <- paste(head$Rs_observ, head$protected)
c.int <- intervals(gls(spe.nr~factor(both)-1, cor=corCompSymm(form=~1|Rs_plot), data = head))[["coef"]] %>% 
  as_tibble(rownames = "both") %>% 
  mutate(Rs_observ = str_sub(both, 13, 16), protected = str_sub(both, 18, 18))

# plot
plot.a <- c.int %>% 
  ggplot(aes(protected, est., fill = Rs_observ))+
  scale_fill_manual(values = c(4, 7), labels = c("1980s", "2022"))+
  geom_pointrange(aes(ymin = lower, ymax = upper), pch = 21, size = 0.8, 
                  position = position_dodge(width = 0.4), fatten = 5)+
  theme_bw()+
  theme(legend.title = element_blank(), legend.text = element_text(size = 11), 
        axis.text.x = element_text(size = 11))+
  xlab("")+
  ggtitle("(a) Number of species per plot",
          subtitle = ("time: p < 0.001; protection status: p = 0.175; \ninteraction: p= 0.053")) + 
  ylab(expression("Number of species"))+
  scale_x_discrete(labels = c("Non-protected", "Protected"))+
  scale_y_continuous(breaks = c(32, 34, 36, 38, 40, 42))

# threatened species ------------------------------------------------------

# model
m1 <- geeglm(cbind(end.nr, spe.nr-end.nr)~Rs_observ*protected, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
m1 <- geeglm(cbind(end.nr, spe.nr-end.nr)~Rs_observ+protected, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
summary(m1)

# confidence intervals
fact<-as.factor(paste(head$Rs_observ, head$protected))
m2 <- geeglm(cbind(end.nr, spe.nr-end.nr)~fact-1, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
summary(m2)

c.int <- tidy(m2, conf.int = T, exponentiate = T) %>% 
  mutate(across(c(conf.low, conf.high, estimate), ~.*100)) %>% 
  mutate(Rs_observ = str_sub(term, 5, 8), protected = str_sub(term, 10, 10))

# plot
plot.b <- c.int %>% 
  ggplot(aes(protected, y = estimate, fill = Rs_observ)) +
  scale_fill_manual(values = c(4, 7), labels = c("1980s", "2022"))+
  geom_pointrange(aes(ymin= conf.low, 
                      ymax = conf.high), pch = 21, position = position_dodge(0.4), size = 0.8, 
                  fatten = 5)+
  theme_bw()+
  xlab("")+
  ggtitle("(b) % of threatened species per plot",
          subtitle = ("time: p < 0.001; protection status: p = 0.083; \ninteraction: p = 0.856"))+ 
  ylab(expression("% of threatened species"))+
  theme(legend.title = element_blank(), legend.text = element_text(size = 11), 
        axis.text.x = element_text(size = 11))+
  scale_x_discrete(labels = c("Non-protected", "Protected"))

# proportion of specialists -----------------------------------------------

# filter only plots classified in 1980s into the Cirsio-Brachypodion pinnati or Bromion erecti alliances
head.spec <- head %>% 
  filter(str_detect(ESY0, "THE|THF"))

# model
m1 <- geeglm(cbind(diag.nr, spe.nr-diag.nr)~Rs_observ*protected, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head.spec)
summary(m1)

# confidence interval
both<-as.factor(paste(head.spec$Rs_observ, head.spec$protected))
m2 <- geeglm(cbind(diag.nr, spe.nr-diag.nr)~both-1, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head.spec)

summary(m2)

c.int <- tidy(m2, conf.int = T, exponentiate = T) %>% 
  mutate(across(c(conf.low, conf.high, estimate), ~.*100)) %>% 
  mutate(Rs_observ = str_sub(term, 5, 8), protected = str_sub(term, 10, 10))

# plot
plot.c <- c.int %>% 
  ggplot(aes(protected, y = estimate, fill = Rs_observ)) +
  scale_fill_manual(values = c(4, 7), labels = c("1980s", "2022"))+
  geom_pointrange(aes(ymin= conf.low, 
                      ymax = conf.high), pch = 21, position = position_dodge(0.4), size = 0.8, 
                  fatten = 5)+
  theme_bw()+
  xlab("")+
  ggtitle(" (c) % of semi-dry grassland specialists per plot",
          subtitle = ("time: p < 0.001; protection status: p = 0.005; \ninteraction: p = 0.005")) + 
  ylab(expression("% of semi-dry grassland specialists"))+
  theme(legend.title = element_blank(), legend.text = element_text(size = 11), 
        axis.text.x = element_text(size = 11))+
  scale_x_discrete(labels = c("Non-protected", "Protected"))

# proportion of alien species ---------------------------------------------

# model
m1 <- geeglm(cbind(alien.nr, spe.nr-alien.nr)~Rs_observ*protected, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
m1 <- geeglm(cbind(alien.nr, spe.nr-alien.nr)~Rs_observ+protected, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
summary(m1)

# cofidence intervals
fact<-as.factor(paste(head$Rs_observ, head$protected))
m2 <- geeglm(cbind(alien.nr, spe.nr-alien.nr)~fact-1, id = Rs_plot, 
             family = binomial(link = "logit"), corstr = "exchangeable", data = head)
summary(m2)

c.int <- tidy(m2, conf.int = T, exponentiate = T) %>% 
  mutate(across(c(conf.low, conf.high, estimate), ~.*100)) %>% 
  mutate(Rs_observ = str_sub(term, 5, 8), protected = str_sub(term, 10, 10))

# plot
plot.d <- c.int %>% 
  ggplot(aes(protected, y = estimate, fill = Rs_observ)) +
  scale_fill_manual(values = c(4, 7), labels = c("1980s", "2022"))+
  geom_pointrange(aes(ymin= conf.low, 
                      ymax = conf.high), pch = 21, position = position_dodge(0.4), size = 0.8, 
                  fatten = 5)+
  theme_bw()+
  xlab("")+
  ggtitle("(d) % of alien species per plot",
          subtitle = ("time: p < 0.001; protection status: p < 0.001; \ninteraction: p = 0.12")) + 
  ylab(expression("% of alien species"))+
  theme(legend.title = element_blank(), legend.text = element_text(size = 11), 
        axis.text.x = element_text(size = 11))+
  scale_x_discrete(labels = c("Non-protected", "Protected"))

#### numbers of invasive species in old and new plots
spe.traits %>% 
  filter(inv_status == "inv") %>% 
  left_join(head %>% select(Releve_number, Rs_observ)) %>% 
  group_by(Taxon, Rs_observ) %>% 
  count()

# combine all plots together ----------------------------------------------
plot.a + plot.b + plot.c + plot.d + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("plots/Fig_4_pointrange.png", width = 9.5, height = 9.5)
