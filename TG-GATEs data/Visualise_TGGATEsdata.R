library(tidyverse)
rm(list=ls())

SaveFigures = 1

treatment = as_labeller(c(
  `nitrofurantoin` = "NIT",  
  `diclofenac`     = "DIC",
  `ketoconazole`   = "KET"),
  label_parsed)

# Load data
pathology_data = rbind(read.csv("open_tggates_pathology_NIT.csv"),
                       read.csv("open_tggates_pathology_DIC.csv"),
                       read.csv("open_tggates_pathology_KET.csv")) %>%
  filter(FINDING_TYPE %in% c("Necrosis", "Cellular infiltration", "Hypertrophy", "Vacuolization, cytoplasmic")) %>%
  filter(!SACRIFICE_PERIOD %in% c("6 hr", "9 hr")) %>%
  select(COMPOUND_NAME, INDIVIDUAL_ID, DOSE_LEVEL, SACRIFICE_PERIOD, FINDING_TYPE, GRADE_TYPE) %>%
  mutate(SACRIFICE_PERIOD = as.numeric(gsub(" day", "", SACRIFICE_PERIOD))) %>%
  mutate(FINDING_TYPE = ifelse(FINDING_TYPE == "Vacuolization, cytoplasmic", "Cytoplasmic vacuolisation", FINDING_TYPE))


full_data = pathology_data %>%
  expand(COMPOUND_NAME, INDIVIDUAL_ID, DOSE_LEVEL, SACRIFICE_PERIOD, FINDING_TYPE) %>%
  left_join(pathology_data) %>%
  mutate(GRADE_TYPE = ifelse(is.na(GRADE_TYPE), "absent", GRADE_TYPE))

# Plot data
plot = ggplot(full_data) + theme_bw() +
  geom_jitter(aes(x = SACRIFICE_PERIOD, y = GRADE_TYPE, colour = DOSE_LEVEL), height = 0.3, width = 0.6) +
  facet_grid(rows = vars(COMPOUND_NAME), cols = vars(FINDING_TYPE)) + expand_limits(y=0) +
  ylab("Severity") + xlab("Time (days)") +
  scale_y_discrete(limits = c("absent", "minimal", "slight", "moderate", "severe")) +
  scale_color_viridis_d(name = "Dose", limits = c("Control", "Low", "Middle", "High"))

# print(plot)                       


plot_necrosis = ggplot(filter(full_data, FINDING_TYPE == "Necrosis")) + theme_bw() +
  geom_jitter(aes(x = SACRIFICE_PERIOD, y = GRADE_TYPE, colour = DOSE_LEVEL), height = 0.3, width = 0.5) +
  facet_grid(cols = vars(COMPOUND_NAME), labeller = labeller(COMPOUND_NAME = treatment)) + 
  expand_limits(y=0) + ylab("Severity necrosis") + xlab("Time (days)") +
  scale_y_discrete(limits = c("absent", "minimal", "slight", "moderate", "severe")) +
  scale_color_viridis_d(name = "Dose", limits = c("Control", "Low", "Middle", "High")) +
  theme(legend.position = "bottom")
print(plot_necrosis) 

plot_infiltration = ggplot(filter(full_data, FINDING_TYPE == "Cellular infiltration")) + theme_bw() +
  geom_jitter(aes(x = SACRIFICE_PERIOD, y = GRADE_TYPE, colour = DOSE_LEVEL), height = 0.3, width = 0.5) +
  facet_grid(cols = vars(COMPOUND_NAME), labeller = labeller(COMPOUND_NAME = treatment)) + 
  expand_limits(y=0) + ylab("Severity cellular infiltration") + xlab("Time (days)") +
  scale_y_discrete(limits = c("absent", "minimal", "slight", "moderate", "severe")) +
  scale_color_viridis_d(name = "Dose", limits = c("Control", "Low", "Middle", "High"))
print(plot_infiltration) 

plot_vacuolisation = ggplot(filter(full_data, FINDING_TYPE == "Cytoplasmic vacuolisation")) + theme_bw() +
  geom_jitter(aes(x = SACRIFICE_PERIOD, y = GRADE_TYPE, colour = DOSE_LEVEL), height = 0.3, width = 0.5) +
  facet_grid(cols = vars(COMPOUND_NAME), labeller = labeller(COMPOUND_NAME = treatment)) + 
  expand_limits(y=0) + ylab("Severity cytoplasmic vacuolisation") + xlab("Time (days)") +
  scale_y_discrete(limits = c("absent", "minimal", "slight", "moderate", "severe")) +
  scale_color_viridis_d(name = "Dose", limits = c("Control", "Low", "Middle", "High"))
print(plot_vacuolisation) 

plot_hypertrophy = ggplot(filter(full_data, FINDING_TYPE == "Hypertrophy")) + theme_bw() +
  geom_jitter(aes(x = SACRIFICE_PERIOD, y = GRADE_TYPE, colour = DOSE_LEVEL), height = 0.3, width = 0.5) +
  facet_grid(cols = vars(COMPOUND_NAME), labeller = labeller(COMPOUND_NAME = treatment)) + 
  expand_limits(y=0) + ylab("Severity hypertrophy") + xlab("Time (days)") +
  scale_y_discrete(limits = c("absent", "minimal", "slight", "moderate", "severe")) +
  scale_color_viridis_d(name = "Dose", limits = c("Control", "Low", "Middle", "High"))
print(plot_hypertrophy) 

# Save Figures
if(SaveFigures){
  ggsave("Figure_pathology_necrosis.svg", plot = plot_necrosis, path = "Figures", 
         width = 18, height = 8, units = "cm")
  ggsave("Figure_pathology_infiltration.svg", plot = plot_infiltration, path = "Figures", 
         width = 20, height = 7, units = "cm")
  ggsave("Figure_pathology_vacuolisation.svg", plot = plot_vacuolisation, path = "Figures", 
         width = 20, height = 7, units = "cm")
  ggsave("Figure_pathology_hypertrophy.svg", plot = plot_hypertrophy, path = "Figures", 
         width = 20, height = 7, units = "cm")
}
