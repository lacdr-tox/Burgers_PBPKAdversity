rm(list = ls())
library(tidyverse)
library(ggnewscale)
library(patchwork)

SaveFigures = 0

data1 = read.csv("ALLdata_KDs.csv") %>%
  mutate(concentration = round(concentration))

#check low cell counts
plot = ggplot(data = data1) + theme_bw() +
  geom_point(aes(x = TimePoint, y = obj_cell_count, color = replicate, shape = as.factor(well_tech))) +
  xlab("Time (hours)") + ylab("Cell Count") + labs(shape = "Well position", color = "Replicate") +
  facet_grid(cols = vars(concentration), rows = vars(siRNA))
# print(plot)
check = filter(data1, obj_cell_count<100)
low_reps = unique(check$REP_si_conc)
data = filter(data1, REP_si_conc != low_reps)
check2 = filter(data, obj_cell_count<100)

# compute mean of well positions 
data_meanT = data %>%
  group_by(concentration, replicate, TimePoint, siRNA) %>%
  summarize(count_mt = mean(obj_cell_count, na.rm = TRUE), 
            count_sd = sd(obj_cell_count, na.rm = TRUE),
            AnVPos_mt = mean(AnVPos_0.1, na.rm = TRUE), 
            AnVPos_sdt = sd(AnVPos_0.1, na.rm = TRUE),
            PIPos_mt = mean(PIPos_0.1, na.rm = TRUE), 
            PIPos_sdt = sd(PIPos_0.1, na.rm = TRUE)
            )

# plot = ggplot(data = data_meanT) + theme_bw() +
#   geom_point(aes(x = TimePoint, y = count_mt, color = replicate)) +
#   xlab("Time (hours)") + ylab("Cell Count") + labs(color = "Replicate") +
#   facet_grid(cols = vars(concentration), rows = vars(siRNA))
# print(plot)

# compute mean of biological replicates
data_mean = data_meanT %>%
  group_by(concentration, TimePoint, siRNA) %>%
  summarize(count_M = mean(count_mt, na.rm = TRUE), 
            count_SD = sd(count_mt, na.rm = TRUE),
            AnVPos_M = mean(AnVPos_mt, na.rm = TRUE), 
            AnVPos_SD = sd(AnVPos_mt, na.rm = TRUE),
            PIPos_M = mean(PIPos_mt, na.rm = TRUE), 
            PIPos_SD = sd(PIPos_mt, na.rm = TRUE)
            )

# Normalise cell count
data_T0 = filter(data_mean, TimePoint == 1) %>%
  rename(count_T0 = count_M) %>%
  ungroup() %>%
  select(concentration, siRNA, count_T0)

data_mean_norm = data_mean %>%
  left_join(data_T0) %>%
  mutate(count_norm_M = count_M/count_T0,
         count_norm_SD = count_SD/count_T0) %>%
  mutate(Apop_cells = AnVPos_M*count_norm_M,
         Healthy_cells = (1-AnVPos_M)*count_norm_M)

# ### PI ###
# plot = ggplot(data = data_mean_norm) +
#   #facet with concentration and siRNA
#   geom_point(data = data_meanT, aes(x = TimePoint, y = PIPos_mt, color = replicate)) +
#   geom_line(aes(x = TimePoint, y = PIPos_M), linewidth = 1.5) +
#   facet_grid(cols = vars(concentration), rows = vars(siRNA)) +
#   #labels
#   xlab("Time (hours)") + ylab("Fraction PI positive cells")
# print(plot)
# 
# #facet with concentration
# plot = ggplot(data = data_mean_norm) +
#   geom_ribbon(aes(x = TimePoint, ymin = PIPos_M-PIPos_SD, ymax = PIPos_M+PIPos_SD, fill = siRNA), alpha = 0.2) +
#   geom_point(aes(x = TimePoint, y = PIPos_M, color = siRNA)) +
#   facet_wrap(vars(concentration)) +
#   xlab("Time (hours)") + ylab("Fraction PI positive cells")
# print(plot)
# 
# # facet_wrap(vars(siRNA))
# plot = ggplot(data = data_mean_norm) +
#   #facet with siRNA
#   geom_ribbon(aes(x = TimePoint, ymin = PIPos_M-PIPos_SD, ymax = PIPos_M+PIPos_SD, fill = as.factor(concentration)), alpha = 0.2) +
#   geom_point(aes(x = TimePoint, y = PIPos_M, color = as.factor(concentration))) +
#   facet_wrap(vars(siRNA)) +
#   labs(fill = "concentration", color = "concentration") +
#   xlab("Time (hours)") + ylab("Fraction PI positive cells")
# print(plot)

# label for facets
dose_names = as_labeller(c(
  `0`   = "0~mu*M~NIT",
  `8`   = "8~mu*M~NIT",
  `15`  = "15~mu*M~NIT",
  `20`  = "20~mu*M~NIT",
  `30`  = "30~mu*M~NIT",
  `60`  = "60~mu*M~NIT",
  `120` = "120~mu*M~NIT"),
  label_parsed)


### Normalised cell count ###
count_plot = ggplot(data = data_mean_norm) + theme_bw() +
  geom_ribbon(aes(x = TimePoint, ymin = count_norm_M-count_norm_SD, ymax = count_norm_M+count_norm_SD, fill = as.factor(siRNA)), alpha = 0.2) +
  geom_point(aes(x = TimePoint, y = count_norm_M, color = as.factor(siRNA))) +
  labs(fill = "siRNA", color = "siRNA") + 
  facet_wrap(vars(concentration), labeller = labeller(concentration = dose_names)) +
  xlab("Time (hours)") + ylab("Normalised cell count")
print(count_plot)
 
data_mean_norm_siPOOL = data_mean_norm %>%
  filter(siRNA == "siKinasePool") %>%
  select(-siRNA)

### AnV ###
plot_anv = ggplot() + theme_bw() +
  # siKinasePool & mock
  geom_ribbon(data = filter(data_mean_norm_siPOOL), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = "siKinasePool"), alpha = 0.2) +
  geom_point(data = filter(data_mean_norm_siPOOL), 
             aes(x = TimePoint, y = AnVPos_M, color = "siKinasePool")) +
  # others
  geom_ribbon(data = data_mean_norm %>% filter(siRNA != "siKinasePool", siRNA != "mock" ), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = siRNA), alpha = 0.2) +
  geom_point(data = data_mean_norm %>% filter(siRNA != "siKinasePool", siRNA != "mock" ), 
             aes(x = TimePoint, y = AnVPos_M, color = siRNA)) +
  facet_grid(cols = vars(siRNA), rows = vars(concentration), labeller = labeller(concentration = dose_names)) + 
  labs(fill = "", color = "") + theme(legend.position = "bottom") +
  xlab("Time (hours)") + ylab("Fraction AnV positive cells")
print(plot_anv)

### one dose selection ###
dose = 120
plot_120 = ggplot() + theme_bw() +
  # siKinasePool & mock
  geom_ribbon(data = filter(data_mean_norm_siPOOL, concentration == dose), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = "siKinasePool"), alpha = 0.2) +
  geom_point(data = filter(data_mean_norm_siPOOL, concentration == dose), 
             aes(x = TimePoint, y = AnVPos_M, color = "siKinasePool")) +
  # others
  geom_ribbon(data = data_mean_norm %>% filter(concentration == dose, siRNA != "siKinasePool", siRNA != "mock" ), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = siRNA), alpha = 0.2) +
  geom_point(data = data_mean_norm %>% filter(concentration == dose, siRNA != "siKinasePool", siRNA != "mock" ), 
               aes(x = TimePoint, y = AnVPos_M, color = siRNA)) +
  facet_wrap(vars(siRNA)) + labs(fill = "", color = "") + theme(legend.position = "bottom") +
  xlab("Time (hours)") + ylab("Fraction AnV positive cells")

print(plot_120)

plot_120_zoom = ggplot() + theme_bw() +
  # siKinasePool & mock
  geom_ribbon(data = filter(data_mean_norm_siPOOL, concentration == dose), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = "siKinasePool"), alpha = 0.2) +
  geom_point(data = filter(data_mean_norm_siPOOL, concentration == dose), 
             aes(x = TimePoint, y = AnVPos_M, color = "siKinasePool")) +
  # others
  geom_ribbon(data = data_mean_norm %>% filter(concentration == dose, siRNA != "siKinasePool", siRNA != "mock" ), 
              aes(x = TimePoint, ymin = AnVPos_M-AnVPos_SD, ymax = AnVPos_M+AnVPos_SD, fill = "siRNA"), alpha = 0.2) +
  geom_point(data = data_mean_norm %>% filter(concentration == dose, siRNA != "siKinasePool", siRNA != "mock" ), 
             aes(x = TimePoint, y = AnVPos_M, color = "siRNA")) +
  facet_wrap(vars(siRNA)) + labs(fill = "", color = "") + theme(legend.position = "bottom") +
  scale_color_viridis_d(begin = 0.25, end = 0.75, aesthetics = c("color", "fill")) +
  xlab("Time (hours)") + ylab("Fraction AnV positive cells") + xlim(0,20) + ylim(0,0.6)

print(plot_120_zoom)


## Save Figures ##
if(SaveFigures == 1){
  ggsave("Figure_KD_anv.pdf", plot = plot_anv, path = "Figures", width = 20,
         height = 25, units = "cm", dpi = 300)
  ggsave("Figure_KD120_anv.svg", plot = plot_120_zoom, path = "Figures", width = 13,
         height = 13, units = "cm", dpi = 300)
}

