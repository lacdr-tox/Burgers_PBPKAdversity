#PBPK modeling 
library(tidyverse)
library(deSolve)
rm(list=ls(all=TRUE))

### Choose what to model
Tstep = 0.1       # Time step for simulation
Species = c("Human", "Rat")[2]  # select which species to simulate
Compound = c("NIT", "DIC", "KET")[3] # select which compound to simulate
multiple = TRUE   # Currently multiple dosing is only option for oral intake
dose_times = c(0) #, 8, 16) #c(0,6,12,18)  # #Time of doses taken in hours after t0 = 0, always 0 when multiple == FALSE
doses_per_day = length(dose_times) # Doses per day
gender = c("male", "female")[1]
elderly = c(FALSE, TRUE)[1] # only works for humans
severe_GFR = c(FALSE, TRUE)[1] # only works for humans
  
### Data
file_name = "DATA_oral_plasma.csv"
Expdata <- read.csv(paste0(getwd(), "/", file_name)) %>%
  filter(species == Species & treatment == Compound) %>%
  # mutate(sd = 0)    # use when you want to plot without errorbar
  mutate(sd = sd) # use when you want to plot with errorbar (if sd available)


choose_Tmax = TRUE
if(length(unique(Expdata$dose)) == 0 | choose_Tmax == TRUE) {
  Tmax = 1*24# Indicate timespan you want to simulate
  } else{
  Tmax = max(Expdata$time)  # Extract max time from experimental data
}

choose_dose = TRUE
if(length(unique(Expdata$dose)) == 0 | choose_dose == TRUE) {
  doses = c(150) # Indicate which doses you want to simulate
  Expdata = Expdata %>% filter(dose %in% doses)
} else{
  doses = unique(Expdata$dose)  # Extract doses from experimental data
}

Times = seq(0, Tmax, Tstep)

### Load parameters
source(paste(Species,"_parms.R", sep = ""))
source(paste(Compound,"_parms.R", sep = ""))

if(Species == "Human"){
  Final_Rdose1 = doses
} else {
  Final_Rdose1 = doses*BW
}

### Physiological parameters scaling 
#Blood flow to organs Scaling
QCblood = QCC*BW; 		 			        #Initial cardiac output for blood L/h  
QCplasma = QCblood;   					    #Adjust initial cardiac output for plasma flow  
Qliver= FQliver* QCplasma;   	 			#Plasma flow to liver 									  	
Qgut= FQgut* QCplasma;    	 			  #Plasma flow to gut 
Qkidney= FQkidney*QCplasma; 				#Plasma flow to kidney
Qfat = FQfat*QCplasma;              #plasma flow to fat  
Qrestbody = QCplasma - (Qliver + Qkidney + Qfat + Qgut);  #plasma flow to rest of the body 
vgut = Fgut*BW; 								    #volume of gut  
vliver = Fliver*BW;  								#volume of liver  
vkidney = Fkidney*BW;								#volume of kidney
vfilterate =Ffilterate*BW;          #volume of tubules  
vfat = Ffat*BW;                   	#volume of fat 
vplasma = Fplasma*BW;               #volume of plasma 
vrestbody = 0.84*BW - vliver- vkidney- vfat- vplasma -vgut -vfilterate;
Qurine = QurineC*BW;
VmaxM = VmaxMC*vliver;

Tr = Trc*BW
Tm = Tmc*BW
Vehr = Vehrc*BW
kbile = kbile*BW

### Compile parameters
para  <-unlist(c(data.frame(
  QCblood ,
  QCplasma ,
  Qliver ,
  Qgut,
  Qkidney ,
  Qfilterate,
  Qfat,
  Qrestbody,
  vliver,
  vkidney,
  vfilterate,
  vfat,
  vplasma,
  vgut,
  vrestbody,
  fu,
  fugut,
  Kgut_plasma,            
  Kliver_plasma,            
  Kkidney_plasma,       
  Kfat_plasma,        
  Krestbody_plasma,
  VmaxM,
  KmM,
  kbp,
  kgutabs,
  Qurine,
  Vehr,
  Kehr,
  kbile,
  Trc,
  Tmc,
  Kt,
  multiple, 
  doses_per_day
)))

### Initial conditions
yini <- unlist(c(data.frame(
  Agutlumen = 0,
  Agut = 0,
  Aliver = 0,
  Abile = 0,
  Akidney = 0,
  Afilterate = 0,
  Adelay = 0, 
  Afat = 0,
  Arestbody = 0,
  Aurine = 0,  
  Afeces = 0,
  Aplasma = 0
)))

### Dosing function  
anaus <- function(t,t0,t1){
  y <- (tanh(100*(t-t0)) - tanh(100*(t-t1)))/2
  return(y)
  }

### PBPK model
PBTKmod <- function(para,dose) {
  derivs <- function(t, y, para) {
    with (as.list(c(y, para)),
          {
            period <- 3/60                                              # (h)       |uptake period
            koa <- dose/period                                          # (mg/h)    |uptake rate  assumed
            if(multiple){
              t0 = dose_times + rep((0:floor(Tmax/24))*24, each=doses_per_day)}
            else{t0 = 0}
            t1 <- t0 + period
            OD <- data.frame(t0, t1) 
            
            ANAUS.O   <- sum(apply(OD, 1, function(x) anaus(t, x[1], x[2])))
            
            Input <-  koa*ANAUS.O 
            
            cgut = Agut/vgut
            
            cliver = Aliver/vliver
            
            cfat = Afat/vfat
            
            ckidney = Akidney/vkidney
            
            cfilterate = Afilterate/vfilterate
            
            crestbody = Arestbody/vrestbody
            
            cplasma = Aplasma/vplasma
            
            RAM = VmaxM*cliver*(fu)/(cliver*(fu)+KmM)
            
            Filteration = Qfilterate*ckidney*((fu/kbp)/Kkidney_plasma)

            dAgutlumen = -kgutabs*Agutlumen + Input - kfeces*Agutlumen + kbile*Abile;

            dAgut = kgutabs*Agutlumen + Qgut*cplasma*(fu/kbp)- Qgut*cgut*((fugut/kbp)/Kgut_plasma);                             

            dAliver = Qliver*cplasma*(fu/kbp) + Qgut*cgut*((fugut/kbp)/Kgut_plasma) - (Qliver +Qgut)*cliver*((fu/kbp)/Kliver_plasma) - RAM - (Vehr*cliver*fu)/(Kehr + cliver*fu) 
            
            dAbile = (Vehr*cliver*fu)/(Kehr + cliver*fu) - kbile*Abile
            
            dAkidney = Qkidney*cplasma*(fu/kbp) - Qkidney*ckidney*((fu/kbp)/Kkidney_plasma) - Filteration + Tr*cfilterate - (Tm*(ckidney*((fu/kbp)/Kkidney_plasma)))/(Kt+(ckidney*((fu/kbp)/Kkidney_plasma)))
            
            dAfilterate = Filteration - Qfilterate*cfilterate - Tr*cfilterate + (Tm*(ckidney*((fu/kbp)/Kkidney_plasma)))/(Kt+(ckidney*((fu/kbp)/Kkidney_plasma)))
            
            dAdelay  = Qfilterate*cfilterate - Qurine *Adelay;
            
            dAfat = Qfat*cplasma*(fu/kbp) - Qfat*cfat*((fu/kbp)/Kfat_plasma)    
            
            dArestbody = Qrestbody*cplasma*(fu/kbp) - Qrestbody*crestbody*((fu/kbp)/Krestbody_plasma)
            
            dAurine = Qurine *Adelay;
            
            dAfeces = kfeces*Agutlumen
            
            dAplasma = (Qliver + Qgut)*cliver*((fu/kbp)/Kliver_plasma)+ Qkidney*ckidney*((fu/kbp)/Kkidney_plasma)+ 
              Qfat*cfat*((fu/kbp)/Kfat_plasma)+ Qrestbody*crestbody*((fu/kbp)/Krestbody_plasma)- QCplasma*cplasma*(fu/kbp)
            
            dydt = c(dAgutlumen, dAgut, dAliver, dAbile, dAkidney, dAfilterate, dAdelay, dAfat,dArestbody,dAurine, dAfeces, dAplasma)
            
            conc <- c(cgut = cgut, cliver = cliver, ckidney = ckidney, cfilterate = cfilterate, cfat = cfat, cplasma=cplasma)
            
            res  <- list(dydt, conc)
            return(res)
          })}
  
  times = seq(0, Tmax, Tstep)
  
  return(ode(y=yini, func=derivs, times=times, parms=para, method="lsoda")) 
}


### Simulate model
z = list()
for(i in 1:length(Final_Rdose1)){

  z[[i]] <- PBTKmod(para, dose = Final_Rdose1[[i]])
  
  }

v1<- data.frame(do.call(rbind,z))

experiment = length(Final_Rdose1)
df <- data.frame(doses)
df1 = df
colnames(df1)= c("dose")
v1$Simulation = (rep(paste0(rep("Simu",experiment),1:experiment),each = length(Times)))  # as we have 3 variables
v1$Chemical = rep(paste0(Compound), each = nrow(v1))
v1$dose = rep(df1$dose,each = length(Times))
Simulationfile = v1 %>%  gather(variable, value, Agutlumen:cplasma)

plot_plasmaliver = function (save) {
  Simulation = Simulationfile %>% filter(.,variable %in% c("cplasma", "cliver"))
  exp_data = Expdata %>% filter(variable %in% c("cplasma"))
  
  plot = ggplot() + theme_bw() +
    geom_line(data =Simulation,aes(x = time, y = value, color="Simulation"), linewidth =1) +
    geom_errorbar(data =exp_data,aes(x = time, ymin = conc-sd, ymax = conc+sd,
                                     color=as.factor(Study)), size = 1) +
    geom_point(data =exp_data,aes(x = time, y = conc, color=as.factor(Study)), size = 2) +
    labs(x = "Time (hours)", y = "Concentration (mg/L)",
         title = ifelse(Species == "Rat", 
                        paste(Species, " - ", doses, " mg/kg ", Compound, " (", round(doses*BW, digits = 2), " mg)", sep = ""),
                        paste(Species, " - ", doses, " mg ", Compound, " (", round(doses/BW, digits = 2), " mg/kg)", sep = ""))) +
    facet_wrap(~factor(variable, levels = c("cplasma", "cliver"), labels = c("Plasma", "Liver")), scales = "free") +
    # scale_y_continuous(trans = 'log10')+
    scale_color_viridis_d(breaks = c("Simulation", sort(unique(Expdata$Study))), direction = -1, end = 0.5) +
    # scale_color_manual(breaks = c("Simulation", sort(unique(Expdata$Study))), values = c("#440154", "#21918c", "#90d743")) +
    theme(legend.title = element_blank(), legend.position="bottom") + 
    xlim(0,Tmax)
  
  print(plot)
  
  if(save == 1){
    ggsave(paste("Figure", Species, doses, Compound,"oral.svg", sep = "_"), plot = plot, path = "Figures", width = 15,
           height = 8, units = "cm", dpi = 300)
  }
}

plot_compartments = function(save){
  Simulation = Simulationfile %>% filter(variable %in% c("cgut", "cliver", "ckidney", "cfilterate", "cfat", "cplasma"))

  plot = ggplot() + theme_bw() +
    geom_line(data =Simulation,aes(x = time, y = value, color="Simulation"), linewidth =1) +
    labs(x = "Time (hours)", y = "Concentration (mg/L)",
         title = paste(Species, "-", doses, ifelse(Species == "Rat", "mg/kg", "mg"),Compound)) +
    facet_wrap(~variable, scales = "free") +
    # scale_y_continuous(trans = 'log10')+
    scale_color_viridis_d(breaks = c("Simulation", sort(unique(Expdata$Study))), direction = 1, end = 0.5)+
    theme(legend.position="none") + xlim(0,Tmax)
  
  print(plot)
  
  if(save == 1){
    ggsave(paste("Figure_compartments", Species, doses, Compound,"oral.svg", sep = "_"), plot = plot, path = "Figures", width = 20,
           height = 11, units = "cm", dpi = 300)
  }
}

plot_plasmaliver(save = 0)

# plot_compartments(save = 0)

if (0) {
  if(elderly){
    name = paste("simu_cliver//", paste(Compound, doses, "mg", Tmax/24, "days", Species, gender,"elderly.csv", sep = "_"), sep = "")
  }else{
    name = paste("simu_cliver//", paste(Compound, doses, "mg", Tmax/24, "days", Species, gender,"regular.csv", sep = "_"), sep = "")
  }
  write.csv(filter(Simulationfile, variable == "cliver"), file = name)
  
}

# name = paste("simu_cliver//", paste(Compound, doses, "mg", Tmax/24, "days", Species, gender,"severeGFR.csv", sep = "_"), sep = "")
# write.csv(filter(Simulationfile, variable == "cliver"), file = name)

