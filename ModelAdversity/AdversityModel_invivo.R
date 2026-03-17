# First load the deSolve package. If this gives an error, you probably need 
# to install it first with install.packages("deSolve")
library(deSolve)
library(tidyverse)

rm(list=ls())

SaveFigures = 0

var_names = as_labeller(c(
  `apoptotic_cells` = "Apoptotic~cells",
  `live_cells`      = "Live~cells",
  `total_cells`     = "Total~cells"),
  label_parsed)

source("model_functions.R")

# Function that defines the model including the stress input function
model <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = 0
    dcells = (growth_rate/(1 + (Ki*(ATF4-0.30702814)))) * (1-((cells)/CM))*cells - (death_byC_rate/(1 + (q*GSH))) * (CHOP-0.378908704) *cells - death_rate*cells
    dapop = (death_byC_rate/(1 + (q*GSH))) * (CHOP-0.378908704) * cells + death_rate*cells - x*float_rate * apop
    list(c(dS, dcells, dapop))
  })
}

tspan <-  seq(0, 15*24, by = 1)

adv_parms = read.csv("..//ParameterCalibration//Adversity_calibration//20250605_135045_EB_NITDMSOModel_parameterEstimates_parmset2_cost_56.46.csv") %>%
  select(X, est_value) %>%
  add_row(X = "q", est_value = 0.33)
adv_parms$est_value[adv_parms$X == "CM"] = 1 
adv_parms$est_value[adv_parms$X == "death_rate"] = 0


EC = c(0)
doses = c(0)
  
ATF4 = 0.30702814 
CHOP = 0.378908704
GSH = 1

# x = 1/3

simulation = data.frame()
for (X in c(1, 2)){
  for(y in c(0.1, 0.25, 0.67)){
    x = c(1/3,1)[X]
    simulation1 = run_model(model, tspan, doses, EC, c(cells = 1-y, apop = y), 
                        adv_parms$X, adv_parms$est_value) %>%
    mutate(StateVar = ifelse(StateVar == "cells", "live_cells", ifelse(StateVar == "apop", "apoptotic_cells", StateVar))) %>%
    pivot_wider(values_from = value, names_from = StateVar) %>%
    pivot_longer(c(live_cells, apoptotic_cells), values_to = "value", names_to = "StateVar")

    simulation1$damage = paste(y*100, "% damage", sep = "")
    simulation1$label = c("Model with reduced value for rd", "Model with in vitro value for rd")[X]
    
    simulation = rbind(simulation, simulation1)
  }}


# plot results
plot = ggplot(NULL) + theme_bw() +
  geom_line(data = simulation, aes(x = time/24, y = value, col = label, 
            linetype = label), linewidth = 1.5) +
  facet_grid(cols = vars(damage), rows = vars(StateVar),
             labeller = labeller(StateVar = var_names)) + 
  scale_color_viridis_d(name = "", aesthetics = c("colour", "fill")) +
  theme(legend.position="bottom") + expand_limits(y = 0) +
  xlab("Time (days)") + ylab("Relative number of cells") + labs(linetype = "")

print(plot)

if(SaveFigures == 1){
  ggsave("Figure_invivo.svg", plot = plot, path = "Figures", width = 20,
         height = 12, units = "cm", dpi = 300)
}
