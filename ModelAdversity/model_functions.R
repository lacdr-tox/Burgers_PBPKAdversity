modeldose <- function(modelname, stress, tspan, initial_states, param_names, param_values, conc){
  pars.nostim = stats::setNames(as.list(param_values), param_names)
  inistate = initial_states
  out <-  deSolve::ode(
    y = c(S = stress, inistate),
    times = tspan,
    func = match.fun(modelname),
    parms = c(pars.nostim, dose = conc)
  )
  return(out)
}

# This function simulates the ODE model with the selected parameters and returns the simulation.
run_model <- function(modelname, tspan, doses, ECs, initial_states, param_names, param_values){
  for (i in 1:length(doses)){
    xx <- modeldose(modelname, stress = ECs[i], tspan, initial_states, param_names, param_values, doses[i])
    dose_uMadj = doses[i]
    
    xx = cbind(dose_uMadj, xx)
    if (i == 1) {
      modeldata = data.frame(xx)
    }
    else{
      modeldata = data.frame(rbind(modeldata, xx))
    }
  }
  modeldata = pivot_longer(modeldata, cols = 3:ncol(modeldata), names_to = "StateVar")
  return(modeldata)
}

include_p53 <- function(simulation) {
  simulation_P53 = filter(simulation, StateVar == "P53") %>% 
    rename(value_p53 = value) %>%
    mutate(StateVar = "p53")
  
  simulation_P53P = filter(simulation, StateVar == "P53P") %>% 
    rename(value_p53p = value) %>%
    mutate(StateVar = "p53")
  
  simulation_p53 = full_join(simulation_P53P, simulation_P53) %>%
    mutate(value = value_p53 + value_p53p) %>%
    select(-c("value_p53p", "value_p53"))
  
  simulation = rbind(simulation, simulation_p53)
  
  return(simulation)
}

sum_G1earlyG0 <- function(simulation) {
  simulation_G1early = filter(simulation, StateVar == "G1early") %>% 
    rename(value_G1early = value) %>%
    mutate(StateVar = "G1early_G0")
  
  simulation_G0 = filter(simulation, StateVar == "G0") %>% 
    rename(value_G0 = value) %>%
    mutate(StateVar = "G1early_G0")
  
  simulation_G1earlyG0 = full_join(simulation_G1early, simulation_G0) %>%
    mutate(value = value_G1early + value_G0) %>%
    select(-c("value_G1early", "value_G0"))
  
  simulation = rbind(simulation, simulation_G1earlyG0)
  
  return(simulation)
}

