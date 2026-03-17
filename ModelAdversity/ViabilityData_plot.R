library(tidyverse)

rm(list=ls())

SaveFigures = 1

for (n in 3){
  compound = c("NIT", "DIC", "KET")[n]

  data= (read.csv(file = paste("ViabilityData_", compound, ".csv", sep = ""), header=TRUE))

  Figure_count = ggplot(data = data) + theme_bw() +
    geom_ribbon(aes(x = mean_time, ymin = mean_count_norm-sd_count_norm, 
                    ymax = mean_count_norm+sd_count_norm, 
                    fill = as.factor(dose_uM)), alpha = 0.2) +
    geom_line(aes(x = mean_time, y = mean_count_norm, color = as.factor(dose_uM)), 
              linewidth=1) +
    xlab("Time (hours)") + ylab("Normalised cell count") +
    labs(color = expression(paste("Concentration KET (", mu, "M)")), 
         fill = expression(paste("Concentration KET (", mu, "M)")))

  print(Figure_count)

  Figure_AnV = ggplot(data = data) + theme_bw() +
    geom_ribbon(aes(x = mean_time, ymin = mean_AnV-sd_AnV, ymax = mean_AnV+sd_AnV, 
                    fill = as.factor(dose_uM)), alpha = 0.2) +
    geom_line(aes(x = mean_time, y = mean_AnV, color = as.factor(dose_uM)), 
              linewidth=1) + 
    xlab("Time (hours)") + ylab("Fraction AnV positive cells") +
    labs(color = expression(paste("Concentration KET (", mu, "M)")), 
         fill = expression(paste("Concentration KET (", mu, "M)")))

  print(Figure_AnV)


  Figure_PI = ggplot(data = data) + theme_bw() +
    geom_ribbon(aes(x = mean_time, ymin = mean_PI-sd_PI, ymax = mean_PI+sd_PI, 
                    fill = as.factor(dose_uM)), alpha = 0.2) +
    geom_line(aes(x = mean_time, y = mean_PI, color = as.factor(dose_uM)), 
              linewidth=1) + 
    xlab("Time (hours)") + ylab("Fraction PI positive cells") +
    labs(color = expression(paste("Concentration KET (", mu, "M)")), 
         fill = expression(paste("Concentration KET (", mu, "M)")))

  print(Figure_PI)


  if(SaveFigures == 1){
    ggsave(paste("Figure_count_", compound, ".svg", sep= ""), plot = Figure_count, path = "Figures", width = 15.5,
           height = 10, units = "cm", dpi = 300)
    ggsave(paste("Figure_AnV_", compound, ".svg", sep= ""), plot = Figure_AnV, path = "Figures", width = 15.5,
           height = 10, units = "cm", dpi = 300)
    ggsave(paste("Figure_PI_", compound, ".svg", sep= ""), plot = Figure_PI, path = "Figures", width = 15.5,
           height = 10, units = "cm", dpi = 300)
  }
}
