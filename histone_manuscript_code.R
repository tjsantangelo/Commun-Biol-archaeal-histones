###The code here was authored by Breanna R. Wenck to analyze the data in 
###the manuscript:
###Archaeal histone-based chromatin structures regulate 
###transcription elongation rates. Commun Biol. (2024)

library(tidyverse)
library(viridis)
library(viridisLite)

###Code used for Figures 2c and 4d###

am <- read_csv("excel/avg_RNAP_movement.csv")

for (i in 1:nrow(am)){
  if(am$TFS[i] == 'no') {
    am$norm_am = (am$avg_mvmt/1.8433360)
  } else if(am$TFS[i] == 'yes') {
    am$norm_am = am$avg_mvmt/1.8441235
  }
}

am_noTFS <- am %>% 
  filter(TFS == 'no') %>% 
  select(type:norm_am)

am_ss <- am %>% 
  group_by(type, TFS) %>%
  summarise(mean_am = mean(avg_mvmt),
            min_am = min(avg_mvmt),
            max_am = max(avg_mvmt),
            n_am = n(),
            sd_am = sd(avg_mvmt),
            se_am = (sd_am/sqrt(n_am)),
            UL = (mean_am + se_am),
            LL = (mean_am - se_am))

for (i in 1:nrow(am_ss)){
  if(am_ss$TFS[i] == 'no') {
    am_ss$norm_mean_am = (am_ss$mean_am/1.8433360)
    am_ss$norm_sd_am = (am_ss$sd_am/1.8433360)
    am_ss$norm_se_am = (am_ss$norm_sd_am/sqrt(am_ss$n_am))
    am_ss$norm_UL = (am_ss$norm_mean_am + am_ss$norm_se_am)
    am_ss$norm_LL = (am_ss$norm_mean_am - am_ss$norm_se_am)
  } else if(am_ss$TFS[i] == 'yes') {
    am_ss$norm_mean_am = am_ss$mean_am/1.8441235
    am_ss$norm_sd_am = am_ss$sd_am/1.8441235
    am_ss$norm_se_am = (am_ss$norm_sd_am/sqrt(am_ss$n_am))
    am_ss$norm_UL = (am_ss$norm_mean_am + am_ss$norm_se_am)
    am_ss$norm_LL = (am_ss$norm_mean_am - am_ss$norm_se_am)
  }
}

am_ss_noTFS <- am_ss %>% 
  filter(TFS == 'no') %>% 
  select(type:norm_LL)

###Combined plot +/-TFS for publication
ggplot() +
  geom_point(data = am, 
             mapping = aes(x = norm_am,
                           y = type,
                           color = type,
                           shape = TFS),
             alpha = 0.3) +
  geom_point(data = am_ss, 
             mapping = aes(x = norm_mean_am,
                           y = type,
                           color = type,
                           shape = TFS,
                           size = TFS),
             alpha = 0.8,
             stroke = 2.0) +
  scale_size_manual(values = c(5,4)) +
  scale_shape_manual(values = c(16,1), 
                     labels = c("- TFS", "+ TFS"),
                     guide = guide_legend(override.aes = list(size = c(5, 4)))) +
  guides(color = "none", size = "none") +
  geom_errorbarh(inherit.aes = FALSE, 
                 data = am_ss, 
                 mapping = aes(y = type,
                               color = type,
                               group = TFS,
                               linetype = TFS,
                               xmin = norm_LL, 
                               xmax = norm_UL, 
                               height = 0.2),
                 alpha = 0.3,
                 show.legend = FALSE) +
  labs(x = "Average RNAP movement relative \nto HTkA-free conditions", 
       y = "HTkA condition of template") +
  scale_y_discrete(limits = c("E19K/G52K", "E19K", "G17D", "E34A", "G52K","R11A", 
                              "WT", "T55L", "E3A","R20S","HTkA-free"),
                   labels = c("E19K/G52K"="E19K/G52K", "E19K"="E19K", "G17D"="G17D", 
                              "E34A"="E34A", "G52K"="G52K","R11A"="R11A",
                              "WT"="WT", "T55L"="T55L", "E3A"="E3A",
                              "R20S"="R20S","HTkA-free"="HTkA-free")) +
  scale_colour_manual(breaks = c("HTkA-free", "R20S","E3A","T55L", "WT","R11A",
                                 "G52K", "E34A", "G17D", "E19K", "E19K/G52K"), 
                      values = c("#004134", "greenyellow", "steelblue3", 
                                 "springgreen4", "grey68", "navy",
                                 "lightslateblue", "lightsteelblue", 
                                 "#FF6633", "purple4","darkviolet")) +
  scale_x_continuous(position = "bottom",
                     guide = guide_axis(position = "top"),
                     breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", face = 'italic',
                                 size = 10, color = "#263238"),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.85),
        legend.text = element_text(family = '', face = 'bold', size = 11))

###No TFS plot for publication
ggplot() +
  geom_point(data = am_noTFS, 
             mapping = aes(x = norm_am,
                           y = type,
                           color = type),
             alpha = 0.3) +
  geom_point(data = am_ss_noTFS, 
             mapping = aes(x = norm_mean_am,
                           y = type,
                           color = type,
                           size = 5),
             alpha = 0.8) +
  guides(color = "none", size = "none") +
  geom_errorbarh(inherit.aes = FALSE, 
                 data = am_ss_noTFS, 
                 mapping = aes(y = type,
                               color = type,
                               xmin = norm_LL, 
                               xmax = norm_UL, 
                               height = 0.2),
                 alpha = 0.3,
                 show.legend = FALSE) +
  labs(x = "Average RNAP movement relative \nto HTkA-free conditions", 
       y = "HTkA condition of template") +
  scale_y_discrete(limits = c("E19K/G52K", "E19K", "G17D", "E34A", "G52K","R11A", 
                              "WT", "T55L", "E3A","R20S","HTkA-free"),
                   labels = c("E19K/G52K"="E19K/G52K", "E19K"="E19K", "G17D"="G17D", 
                              "E34A"="E34A", "G52K"="G52K","R11A"="R11A",
                              "WT"="WT", "T55L"="T55L", "E3A"="E3A",
                              "R20S"="R20S","HTkA-free"="HTkA-free")) +
  scale_colour_manual(breaks = c("HTkA-free", "R20S","E3A","T55L", "WT","R11A",
                                 "G52K", "E34A", "G17D", "E19K", "E19K/G52K"), 
                      values = c("#004134", "greenyellow", "steelblue3", 
                                 "springgreen4", "grey68", "navy",
                                 "lightslateblue", "lightsteelblue", 
                                 "#FF6633", "purple4","darkviolet")) +
  scale_x_continuous(position = "bottom",
                     guide = guide_axis(position = "top"),
                     breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", face = 'italic',
                                 size = 10, color = "#263238"),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.85),
        legend.text = element_text(family = '', face = 'bold', size = 11))

###Code used for Figures 2a and 4b ###

trx_all <- read_csv("excel/hist_meas_all.csv")

trx_all <- trx_all %>%
  mutate(time = as.factor(time),
         type = as.factor(type),
         type = fct_inorder(type),
         bin = as.factor(bin),
         bin = fct_rev(bin))

trx_all_ss <- trx_all %>% 
  group_by(type, time, bin) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

trx_all_ss %>% ggplot(aes(x = time,
                          y = mean_per,
                          fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (sec)", y = "Percentage of transcripts") +
  facet_wrap(~ type, scales = "free_x") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))

TFS_trx_all <- read_csv("TFS/TFS_meas_all.csv")

TFS_trx_all <- TFS_trx_all %>%
  mutate(time = as.factor(time),
         type = as.factor(type),
         type = fct_inorder(type),
         bin = as.factor(bin),
         bin = fct_rev(bin))

TFS_trx_all_ss <- TFS_trx_all %>% 
  group_by(type, time, bin) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per),
            .groups = "drop")

TFS_trx_all_ss %>% ggplot(aes(x = time,
                              y = mean_per,
                              fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (sec)", y = "Percentage of transcripts (+TFS)") +
  facet_wrap(~ type, scales = "free_x") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))


###Code used for Figures 2b and 4c###

bins <- read_csv("excel/bins1_7.csv")

bins <- bins %>%
  mutate(time = as.factor(time),
         type = as.factor(type),
         type = fct_inorder(type),
         bin = as.factor(bin),
         TFS = as.factor(TFS))

bins_ss <- bins %>% 
  group_by(type, bin, time, TFS) %>%
  summarise(prct = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (prct + se_per),
            LL = (prct - se_per),
            .groups = "drop")

bin7_noTFS <- bins %>% 
  filter(TFS == "no") %>% 
  filter(bin == "7") %>% 
  select(type:TFS) %>% 
  mutate(prct = percent)

bin7_ss_noTFS <- bins_ss %>% 
  filter(TFS == "no") %>% 
  filter(bin == "7") %>% 
  select(type:LL)

bin7_TFS <- bins %>% 
  filter(TFS == "yes") %>% 
  filter(bin == "7") %>% 
  select(type:TFS) %>% 
  mutate(prct = percent)

bin7_ss_TFS <- bins_ss %>% 
  filter(TFS == "yes") %>%
  filter(bin == "7") %>% 
  select(type:LL)

###No TFS plot for publication
ggplot() +
  geom_point(data = bin7_noTFS,  
             mapping = aes(x = time,
                           y = prct,
                           color = type),
             alpha = 0.15) +
  geom_point(data = bin7_ss_noTFS,
             mapping = aes(x = time,
                           y = prct,
                           color = type,
                           size = 4),
             alpha = 0.6) +
  guides(size = "none") +
  geom_line(data = bin7_ss_noTFS,
            mapping = aes(x = time,
                          y = prct,
                          color = type,
                          group = type),
            alpha = 0.3) +
  geom_errorbar(inherit.aes = FALSE, 
                data = bin7_ss_noTFS,
                mapping = aes(x = time,
                              ymin = LL, 
                              ymax = UL, 
                              width = 0.1,
                              color = type), 
                alpha = 0.4) +
  xlab("Time (sec)") +
  ylab(bquote(bold('% TECs'['+209-231']))) +
  expand_limits(y = 0) +
  ylim(0, 65) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue"),
                      name = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 10, color = "#263238"),
        legend.title = element_text(family = '', face = 'bold', size = 12),
        legend.text = element_text(family = '', size = 10, 
                                   face = 'italic', color = "#263238"),
        legend.title.align = 0.5)


###TFS plot for publication 

ggplot() +
  geom_point(data = bin7_TFS,  
             mapping = aes(x = time,
                           y = prct,
                           color = type),
             shape = 1, 
             stroke = 0.9,
             alpha = 0.2) +
  geom_point(data = bin7_ss_TFS,
             mapping = aes(x = time,
                           y = prct,
                           color = type),
             size = 3,
             shape = 1, 
             stroke = 1,
             alpha = 0.8) +
  guides(size = "none") +
  geom_line(data = bin7_ss_TFS,
            mapping = aes(x = time,
                          y = prct,
                          color = type,
                          group = type),
            alpha = 0.3) +
  geom_errorbar(inherit.aes = FALSE, 
                data = bin7_ss_TFS,
                mapping = aes(x = time,
                              ymin = LL, 
                              ymax = UL, 
                              width = 0.1,
                              color = type), 
                alpha = 0.2) +
  xlab("Time (sec)") +
  ylab(bquote(bold('% TECs'['+209-231']))) +
  expand_limits(y = 0) +
  ylim(0, 75) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue"),
                      name = "Condition (+TFS)") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 10, color = "#263238"),
        legend.title = element_text(family = '', face = 'bold', size = 12),
        legend.text = element_text(family = '', size = 10, 
                                   face = 'italic', color = "#263238"),
        legend.title.align = 0.5)

###Code used for Figures 3d and 4e###

ppp <- read_csv("excel/all_reps_ppp.csv")

ppp <- ppp %>% 
  mutate(Time = as.factor(Time),
         Type = as.factor(Type))

ppp_ss <- ppp %>% 
  group_by(Time, Type) %>% 
  summarise(mean_per = mean(Percent),
            sd_per = sd(Percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per),
            .groups = "drop")

ppp_TFS <- read_csv("TFS/all_TFS_ppp.csv")

ppp_TFS <- ppp_TFS %>% 
  mutate(time = as.factor(time),
         type = as.factor(type))

ppp_TFS_ss <- ppp_TFS %>% 
  group_by(time, type) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

###No TFS plot for publication
ggplot() +
  geom_point(data = ppp,  
             mapping = aes(x = Time,
                           y = Percent,
                           color = Type),
             alpha = 0.2) +
  geom_point(data = ppp_ss,
             mapping = aes(x = Time,
                           y = mean_per,
                           color = Type,
                           size = 3),
             alpha = 0.6) +
  geom_line(data = ppp_ss,
            mapping = aes(x = Time,
                          y = mean_per,
                          color = Type,
                          group = Type),
            alpha = 0.3) +
  guides(size = "none") +
  geom_errorbar(inherit.aes = FALSE, 
                data = ppp_ss,
                mapping = aes(x = Time,
                              ymin = LL, 
                              ymax = UL, 
                              width = 0.1,
                              color = Type),
                alpha = 0.2) +
  labs(x = "Time (sec)", y = bquote(bold('% TECs'['+70']))) +
  scale_y_continuous(trans = "log10") +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue"),
                      name = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 10, color = "#263238"),
        legend.title = element_text(family = '', face = 'bold', size = 12),
        legend.text = element_text(family = '', size = 10, 
                                   face = 'italic', color = "#263238"),
        legend.title.align = 0.5)

###TFS plot for publication 
ggplot() +
  geom_point(data = ppp_TFS,  
             mapping = aes(x = time,
                           y = percent,
                           color = type),
             shape = 1, 
             stroke = 0.9,
             alpha = 0.2) +
  geom_point(data = ppp_TFS_ss,
             mapping = aes(x = time,
                           y = mean_per,
                           color = type),
             size = 3,
             shape = 1, 
             stroke = 1,
             alpha = 0.8) +
  geom_line(data = ppp_TFS_ss,
            mapping = aes(x = time,
                          y = mean_per,
                          color = type,
                          group = type),
            alpha = 0.3) +
  guides(size = "none") +
  geom_errorbar(inherit.aes = FALSE, 
                data = ppp_TFS_ss,
                mapping = aes(x = time,
                              ymin = LL, 
                              ymax = UL, 
                              width = 0.1,
                              color = type),
                alpha = 0.2) +
  labs(x = "Time (sec)", y = bquote(bold('% TECs'['+70']))) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 10, 100),
                     expand = expansion(mult = 0.1)) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue"),
                      name = "Condition (+TFS)") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 10, color = "#263238"),
        legend.title = element_text(family = '', face = 'bold', size = 12),
        legend.text = element_text(family = '', size = 10, 
                                   face = 'italic', color = "#263238"),
        legend.title.align = 0.5)

###Code used for Figure 3f###

exchange <- read_csv("excel/exchange_msmts_long.csv")

exchange <- exchange %>% mutate(time = as.factor(time),
                                mins = as.factor(mins),
                                type = as.factor(type),
                                type = fct_inorder(type),
                                bin = as.factor(bin),
                                bin = fct_rev(bin))

exchange %>% ggplot(aes(x = mins,
                        y = percent,
                        fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.9) +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (min)", y = "Percentage of transcripts") +
  facet_wrap(~ type, scales = "free_x", ncol = 4, shrink = TRUE) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 30),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))

###Code used for Figure 3c###

ec <- read_csv("excel/ext_chase_msmts_long.csv")

ec <- ec %>% mutate(time = as.factor(time),
                    mins = as.factor(mins),
                    type = as.factor(type),
                    type = fct_inorder(type),
                    bin = as.factor(bin),
                    bin = fct_rev(bin))

ec %>% ggplot(aes(x = mins,
                  y = percent,
                  fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (min)", y = "Percentage of transcripts") +
  facet_wrap(~ type) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))

###Code used for Supplementary Figure 4###

decay_1 <- read_csv("excel/no_TFS_decay_1.csv")

decay_1 <- decay_1 %>% 
  mutate(TFS = as.factor(TFS),
         type = as.factor(type),
         type = fct_inorder(type),
         time = as.factor(time))

decay_noTFS <- decay_1 %>% 
  filter(TFS == 'no') %>% 
  select(type:TFS)

decay_TFS <- decay_1 %>% 
  filter(TFS == 'yes') %>% 
  select(type:TFS)

decay_noTFS_ss <- decay_noTFS %>% 
  group_by(type, time) %>% 
  summarize(mean_per = mean(percent),
            min_per = min(percent),
            max_per = max(percent),
            n_per = n(),
            sd_per = sd(percent),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per),
            .groups = "drop")

decay_TFS_ss <- decay_TFS %>% 
  group_by(type, time) %>% 
  summarize(mean_per = mean(percent),
            min_per = min(percent),
            max_per = max(percent),
            n_per = n(),
            sd_per = sd(percent),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per),
            .groups = "drop")

###No TFS plot for publication

decay_noTFS_ss %>% ggplot(aes(x = time,
                              y = mean_per,
                              fill = type)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 0.2),
                alpha = 0.4) +
  scale_fill_manual(breaks = c("HTkA-free", "WT","G17D","E19K", "G52K",
                               "E19K/G52K", "R20S", "T55L", "E3A", "R11A", 
                               "E34A"),
                    values = c("#004134", "grey68", "#FF6633", "purple4",
                               "lightslateblue", "darkviolet", "greenyellow", 
                               "springgreen4", "steelblue3", "navy",
                               "lightsteelblue")) +
  xlab("Time (sec)") +
  ylab(bquote(bold('% TECs'['+58']))) +
  facet_wrap(~ type, scales = "free_x") +
  geom_point(data = decay_noTFS,
             mapping = aes(x = time, 
                           y = percent, 
                           color = type),
             alpha = 0.2) +
  guides(color = FALSE) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  scale_y_continuous(breaks = c(10, 30, 50),
                     limits = c(0,65)) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "none",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))

###TFS plot for publication 

decay_TFS_ss %>% ggplot(aes(x = time,
                            y = mean_per,
                            fill = type)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 0.2),
                alpha = 0.4) +
  scale_fill_manual(breaks = c("HTkA-free", "WT","G17D","E19K", "G52K",
                               "E19K/G52K", "R20S", "T55L", "E3A", "R11A", 
                               "E34A"),
                    values = c("#004134", "grey68", "#FF6633", "purple4",
                               "lightslateblue", "darkviolet", "greenyellow", 
                               "springgreen4", "steelblue3", "navy",
                               "lightsteelblue")) +
  xlab("Time (sec)") +
  ylab(bquote(bold('% TECs'['+58']~'(+TFS)'))) +
  facet_wrap(~ type, scales = "free_x") +
  geom_point(data = decay_TFS,
             mapping = aes(x = time, 
                           y = percent, 
                           color = type),
             alpha = 0.2) +
  guides(color = FALSE) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  scale_y_continuous(breaks = c(10, 30, 50),
                     limits = c(0, 65)) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        axis.text = element_text(family = "", size = 8, color = "black"),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(family = '', face = 'bold', size = 10),
        legend.text = element_text(size = 8),
        legend.position = "none",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 10))
