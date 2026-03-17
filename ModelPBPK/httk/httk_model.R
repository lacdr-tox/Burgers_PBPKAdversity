library(tidyverse)
library(httk)
rm(list=ls(all=TRUE))

save = 0

data_plasma = read.csv("../DATA_oral_plasma.csv")

parameters_NIT <- parameterize_pbtk(chem.name = "nitrofurantoin", species = "human")
# parameters_NIT_rat <- parameterize_pbtk(chem.name = "nitrofurantoin", species = "rat")

parameters_DIC <- parameterize_pbtk(chem.name = "diclofenac", species = "human")
parameters_DIC_rat <- parameterize_pbtk(chem.name = "diclofenac", species = "rat")

parameters_KET <- parameterize_pbtk(chem.name = "ketoconazole", species = "human")
# parameters_KET_rat <- parameterize_pbtk(chem.name = "ketoconazole", species = "rat")

out_NIT <- solve_pbtk(parameters = parameters_NIT, 
                      species        = "Human",
                      daily.dose     = 50/70,
                      doses.per.day  = 1,
                      days           = 1
)

out_DIC <- solve_pbtk(parameters = parameters_DIC, 
                      species        = "Human",
                      daily.dose     = 50/70,
                      doses.per.day  = 1,
                      days           = 1
)

out_KET <- solve_pbtk(parameters = parameters_KET, 
                      species        = "Human",
                      daily.dose     = 200/70,
                      doses.per.day  = 1,
                      days           = 1
)


Data_NIT = filter(data_plasma, treatment == "NIT" & dose == 50 & species == "Human")
plot_NIT = ggplot() + theme_bw() +
  geom_line(data = out_NIT, aes(x = time*24, y = Cplasma*parameters_NIT$MW*0.001, 
                                colour = "Simulation httk"), linewidth = 1) +
  geom_errorbar(data = Data_NIT, aes(x = time, ymin = conc-sd, ymax = conc+sd, 
                                     color=as.factor(Study)), linewidth = 1) +
  geom_point(data = Data_NIT, aes(x = time, y = conc, color=as.factor(Study)), 
             size = 2) +
  scale_color_viridis_d(breaks = c("Simulation httk", sort(unique(Data_NIT$Study))), 
                        direction = -1, end = 0.5) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  xlab("Time (hours)") + ylab("plasma concentration (mg/L)") + 
  labs(title = "Nitrofurantoin (50 mg)")

Data_DIC = filter(data_plasma, treatment == "DIC" & dose == 50 & species == "Human")
plot_DIC = ggplot() + theme_bw() +
  geom_line(data = out_DIC, aes(x = time*24, y = Cplasma*parameters_DIC$MW*0.001, 
                                colour = "Simulation httk"), linewidth = 1) +
  geom_errorbar(data = Data_DIC, aes(x = time, ymin = conc-sd, ymax = conc+sd, 
                                     color=as.factor(Study)), linewidth = 1) +
  geom_point(data = Data_DIC, aes(x = time, y = conc, color=as.factor(Study)), 
             size = 2) +
  scale_color_manual(breaks = c("Simulation httk", sort(unique(Data_DIC$Study))), values = c("#440154", "#21918c", "#90d743")) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  xlab("Time (hours)") + ylab("plasma concentration (mg/L)") + labs(title = "Dicolfenac (50 mg)")


Data_KET = filter(data_plasma, treatment == "KET" & dose == 200 & species == "Human" & time<25)
plot_KET = ggplot() + theme_bw() +
  geom_line(data = out_KET, aes(x = time*24, y = Cplasma*parameters_KET$MW*0.001, 
                                colour = "Simulation httk"), linewidth = 1) +
  geom_errorbar(data = Data_KET, aes(x = time, ymin = conc-sd, ymax = conc+sd, 
                                     color=as.factor(Study)), linewidth = 1) +
  geom_point(data = Data_KET, aes(x = time, y = conc, color=as.factor(Study)), 
             size = 2) +
  scale_color_manual(breaks = c("Simulation httk", sort(unique(data_plasma$Study))), values = c("#440154", "#21918c", "#90d743")) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  xlab("Time (hours)") + ylab("plasma concentration (mg/L)") + labs(title = "Ketoconazole (200 mg)")

print(plot_NIT)
print(plot_DIC)
print(plot_KET)

if(save == 1){
  ggsave("httk_figure_NIT.svg", plot = plot_NIT, path = "../Figures", width = 9,
         height = 9, units = "cm", dpi = 300)
  ggsave("httk_figure_DIC.svg", plot = plot_DIC, path = "../Figures", width = 9,
         height = 9, units = "cm", dpi = 300)
  ggsave("httk_figure_KET.svg", plot = plot_KET, path = "../Figures", width = 9,
         height = 9, units = "cm", dpi = 300)
}
