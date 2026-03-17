# First load the deSolve package. If this gives an error, you probably need 
# to install it first with install.packages("deSolve")
library(deSolve)
library(tidyverse)

rm(list=ls())

dose_names = as_labeller(c(
  `0.2`   = "0~mu*M~NIT",
  `7.5`   = "7.5~mu*M~NIT",
  `30`  = "30~mu*M~NIT",
  `60`  = "60~mu*M~NIT",
  `120`  = "120~mu*M~NIT"),
  label_parsed)

var_names = as_labeller(c(
  `Apop_cells`   = "Apoptotic~cells",
  `Live_cells`   = "Live~cells",
  `Total_cells`  = "Total~cells"),
  label_parsed)

SaveFigures = 0

source("model_functions.R")

# Load  measurement data
measurements_POOL = read.csv(
  file = "POOLdata_max20_NIT_20240529.csv", header=TRUE) %>%
  rename(dose_uMadj = concentration) %>%
  filter(dose_uMadj>0 & dose_uMadj!=15 & dose_uMadj!=20)
 
mean_POOL = measurements_POOL %>%
  group_by(StateVar, dose_uMadj, TimePoint)%>%
  summarise(mean = mean(data4modelReal, na.rm=TRUE), 
            sd = sd(data4modelReal, na.rm=TRUE))

measurements_siNRF2 = read.csv(
  file = "siNRF2data_max20_NIT_20240529.csv", header=TRUE) %>%
  rename(dose_uMadj = concentration) %>%
  filter(dose_uMadj>0 & dose_uMadj!=15 & dose_uMadj!=20)

mean_siNRF2 = measurements_siNRF2 %>%
  group_by(StateVar, dose_uMadj, TimePoint) %>%
  summarise(mean = mean(data4modelReal, na.rm=TRUE), 
            sd = sd(data4modelReal, na.rm=TRUE))

doses = unique(measurements_POOL$dose_uMadj)

# Function that defines the model including the stress input function
model <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 =  - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2* NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    dCompound = S - degrade*Compound
    dGSH = buildGBase + Vmax_G * ((NRF2**nG)/(K_G**nG + NRF2**nG)) - degradG * (degGA4 * ATF4 + 1) * GSH
    dcells = (growth_rate) * cells - (death_byC_rate/(1+q*GSH)) * (CHOP) *cells
    dapop = (death_byC_rate/(1+q*GSH)) * (CHOP) * cells - float_rate * apop
    list(c(dS, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1, dCompound, dGSH, dcells, dapop))
  })
}

tspan <-  seq(0, 20, by = 1)

baseparms_NIT = read.csv(file = "NITparms_EB.csv", header=TRUE)

#first we fitted the model for siNRF2
simu_pars = baseparms_NIT %>%
  rbind(read.csv("..//ParameterCalibration//Adversity_calibration//20240715_091310_EB_KDModel_parameterEstimates_parmset9_cost_881.56.csv") %>% select(-init_value)) %>%
  add_row(X = "q", est_value = 0)
  
inistate = c(ATF4 = 0.30702814, 
             CHOP = 0.378908704, 
             KEAP1 = 1, 
             modified_K1 = 0, 
             NRF2 = 0.557889857, 
             SRXN1 = 0.467017171, 
             Compound = 0,
             GSH = 1,
             cells = simu_pars$est_value[simu_pars$X=="cells_init"],
             apop = simu_pars$est_value[simu_pars$X=="apop_init"])

ECs = c(0, baseparms_NIT$est_value[baseparms_NIT$X=="EC_6"],
        baseparms_NIT$est_value[baseparms_NIT$X=="EC_30"],
        baseparms_NIT$est_value[baseparms_NIT$X=="EC_60"],
        baseparms_NIT$est_value[baseparms_NIT$X=="EC_120"]
)

# simulation for siNRF2
simulation_siNRF2 = run_model(model, tspan, doses, ECs, inistate, simu_pars$X, simu_pars$est_value) %>%
  mutate(StateVar = ifelse(StateVar == "cells", "Live_cells", ifelse(StateVar == "apop", "Apop_cells", StateVar))) %>%
  pivot_wider(values_from = value, names_from = StateVar) %>%
  mutate(Total_cells = Live_cells + Apop_cells) %>%
  pivot_longer(c(Total_cells, Live_cells, Apop_cells), values_to = "value", names_to = "StateVar")

if (0) {
  max_i = 100
  costs = data.frame(simulation_number = 1:max_i, cost = 0)

  for (i in 1:max_i){
    simu_pars$est_value[simu_pars$X == "q"] = i/max_i #simu_pars_POOL$est_value[simu_pars_POOL$X == "q"]
  
    simulation_POOL = run_model(model, tspan, doses, ECs, inistate, simu_pars$X, simu_pars$est_value) %>%
      mutate(StateVar = ifelse(StateVar == "cells", "Live_cells", ifelse(StateVar == "apop", "Apop_cells", StateVar))) %>%
      pivot_wider(values_from = value, names_from = StateVar) %>%
      mutate(Total_cells = Live_cells + Apop_cells) %>%
      pivot_longer(c(Total_cells, Live_cells, Apop_cells), values_to = "value", names_to = "StateVar")
    
    cost_data = simulation_POOL %>%
      select(dose_uMadj, time, StateVar, value) %>%
      rename(TimePoint = time, simulation = value) %>%
      full_join(measurements_POOL) %>%
      mutate(difference = abs(data4modelReal - simulation))
  
    COST = sum(cost_data$difference, na.rm = TRUE)
  
    costs$cost[i] = COST
}
  plot_cost = ggplot(data = costs) + theme_bw() +
    geom_point(aes(x = simulation_number/max_i, y = cost)) +
    xlab("q") + ylab("Cost")
  
  print(plot_cost)
  if(SaveFigures == 1){
    ggsave("Figure_cost_KDmodels.svg", plot = plot_cost, path = "Figures", width = 12,
           height = 10, units = "cm", dpi = 300)
  }
}

simu_pars2 = simu_pars
simu_pars2$est_value[simu_pars2$X == "q"] = 0.33

simulation_POOL = run_model(model, tspan, doses, ECs, inistate, simu_pars2$X, simu_pars2$est_value) %>%
  mutate(StateVar = ifelse(StateVar == "cells", "Live_cells", ifelse(StateVar == "apop", "Apop_cells", StateVar))) %>%
  pivot_wider(values_from = value, names_from = StateVar) %>%
  mutate(Total_cells = Live_cells + Apop_cells) %>%
  pivot_longer(c(Total_cells, Live_cells, Apop_cells), values_to = "value", names_to = "StateVar")


# combined - all conc
plot_siNRF2 = ggplot(NULL) + theme_bw() +
  # data mean & sd
  geom_ribbon(data = filter(mean_siNRF2), aes(x = TimePoint, ymax = mean+sd, ymin = mean-sd, fill = "siNFE2L2 data"), alpha = 0.2) +
  geom_point(data = filter(mean_siNRF2), aes(x = TimePoint, y = mean, colour = "siNFE2L2 data")) +
  #simulation
  geom_line(data = simulation_siNRF2, aes(x = time, y = value, col = "siNFE2L2 simulation"), linewidth = 1.5) +
  # data mean & sd
  # geom_ribbon(data = filter(mean_POOL), aes(x = TimePoint, ymax = mean+sd, ymin = mean-sd, fill = "siPool data"), alpha = 0.2) +
  # geom_point(data = filter(mean_POOL), aes(x = TimePoint, y = mean, col = "siPool data")) +
  # #simulation
  # geom_line(data = simulation_POOL, aes(x = time, y = value, col = "siPool simulation"), linewidth = 1.5) +
  scale_colour_viridis_d(name = "", aesthetics = c("colour", "fill"), begin = 0.2, end = 0.4) +
  facet_grid(rows = vars(StateVar), cols = vars(dose_uMadj), scales = "free",
             labeller = labeller(dose_uMadj = dose_names, StateVar = var_names)) +
  theme(legend.position="bottom") + xlab("Time (hours)") + ylab("Relative number of cells")
print(plot_siNRF2)

plot_siPOOL = ggplot(NULL) + theme_bw() +
  # data mean & sd
  # geom_ribbon(data = filter(mean_siNRF2), aes(x = TimePoint, ymax = mean+sd, ymin = mean-sd, fill = "siNRF2 data"), alpha = 0.2) +
  # geom_point(data = filter(mean_siNRF2), aes(x = TimePoint, y = mean, colour = "siNRF2 data")) +
  # #simulation
  # geom_line(data = simulation_siNRF2, aes(x = time, y = value, col = "siNRF2 simulation"), linewidth = 1.5) +
  # data mean & sd
  geom_ribbon(data = filter(mean_POOL), aes(x = TimePoint, ymax = mean+sd, ymin = mean-sd, fill = "siPool data"), alpha = 0.2) +
  geom_point(data = filter(mean_POOL), aes(x = TimePoint, y = mean, col = "siPool data")) +
  #simulation
  geom_line(data = simulation_POOL, aes(x = time, y = value, col = "siPool simulation"), linewidth = 1.5) +
  scale_colour_viridis_d(name = "", aesthetics = c("colour", "fill"), begin = 0.6, end = 0.8) +
  facet_grid(rows = vars(StateVar), cols = vars(dose_uMadj), scales = "free",
             labeller = labeller(dose_uMadj = dose_names, StateVar = var_names)) +   
               theme(legend.position="bottom") + xlab("Time (hours)") + ylab("Relative number of cells")
print(plot_siPOOL)

if(SaveFigures == 1){
  ggsave("Figure_KDmodels_siNRF2.svg", plot = plot_siNRF2, path = "Figures", width = 20,
         height = 14, units = "cm", dpi = 300)
  ggsave("Figure_KDmodels_siPOOL.svg", plot = plot_siPOOL, path = "Figures", width = 20,
         height = 14, units = "cm", dpi = 300)
}
