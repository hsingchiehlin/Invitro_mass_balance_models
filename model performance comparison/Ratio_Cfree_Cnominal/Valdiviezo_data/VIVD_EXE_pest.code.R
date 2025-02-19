library(ggplot2)
library(scales)
library(patchwork)

source("MCSim/function.R")
set_PATH()
makemod()

VIVD.Che.list <- read.csv("VIVD.Che.para_pest.csv")

model <- "INERIS_VIVD_dynamic.model"
makemcsim(model)

for(i in 1:20){
  
  file.name <- paste("VIVD_dynamic_", i, ".in", sep="")
  Q_total <- paste("Q_total = ", 300*0.000001*0.000001*10, ";\n", sep="")
  
  cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
  sink(file.name, append=TRUE)
  
  cat('OutputFile("sim.out");\n\n')
  
  cat('# Parameters

temp         = 37;    # culture temperature (degrees celsius)
Vmedium      = 300;    # (microL)
Sysdiameter  = 15;   # (mm)
Vwell        = 4000;   # (microL), 1mm3 = 1microL
Celldiameter = 30;    # (microm)
ncell        = 1e-30;\n\n')
  
  cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NL    = 0.0348;  # fraction (v/v) neutral lipid 
f_NLrbc = 0.0017;  # ~        in red blood cells

f_NP    = 0.0252;  # fraction (v/v) neutral phospholipid
f_NPrbc = 0.0029;  # ~        in red blood cells

f_IW    = 0.586;   # fration (v/v) intracellular water
f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n\n')
  
  cat('# Culture medium composition

pH_medium     = 7.4;      # culture medium pH
f_serum       = 0.1;      # fraction (v/v) of serum in total media
pH_serum      = 7.0;      # FBS pH
f_FBS_albumin = 0.0176;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.00046;  # fraction (v/v) of neutral lipid in FBS 


# Physicochemical properties of test substance (at 25 degrees Celsius)

# cmp_types: 1:neutral, 2:monobase, 3:dibase, 4:monoacid, 5:diacid, 6:ampholyte\n')
  
  
  cmp_type <- paste("cmp_type = ", VIVD.Che.list$Compound.type[i], ";\n", sep="")
  pKa1 <- paste("pKa1 = ", VIVD.Che.list$pKa[i], ";\n", sep="")
  cat(cmp_type)
  cat(pKa1)
  
  cat('pKa2        = 6;
BP          = 0.1;
fu          = -1; # unknown, to be computed\n')
  
  MW <- paste("MW = ", VIVD.Che.list$MW[i], ";\n", sep="")
  logPow_Tref <- paste("logPow_Tref = ", VIVD.Che.list$LogKOW[i], ";\n", sep="")
  logKAW_Tref <- paste("logKAW_Tref = ", VIVD.Che.list$LogKAW[i], ";\n", sep="")
  cat(MW)
  cat(logPow_Tref)
  cat(logKAW_Tref)
  
  
  cat('K_met = 0; # metabolic clearance per cell, L/min


Simulation {\n\n')
  
  cat(Q_total)
  
  cat('PrintStep(Q_met, 
            C_cells,
            C_air,
            C_plastic,
            C_free_diss,
            C_IW,
            C_lyso,
            C_mito,
            C_parent, 0, 1440, 60);

  PrintStep(Q_parent,
            Q_media_free,
            Q_cells,
            Q_media,
            Q_air,
            Q_plastic,
            Q_IW,
            Q_lyso,
            Q_mito, 0, 1440, 60);

}

End.\n')      
  
  sink()
  
  out <- mcsim(model, file.name)
  
  VIVD.Che.list[i, c(9:27)] <- out[25, ] # put %abs into the "mc_sim_results"
  
  write.csv(VIVD.Che.list, "VIVD.Che.sim.results_pest_10uM.csv")
  file.remove(paste("MCSim/VIVD_dynamic_", i, ".in", sep=""))
  print(i)
}


MassBalan.results <- read.csv("VIVD.Che.sim.results.csv")
MassBalan.results$C_free_VIVD <- MassBalan.results$C_free_VIVD*1000000
MassBalan.results$C_free_Fischer <- MassBalan.results$C_free_Fischer*1000000
MassBalan.results$C_free_IVBMB <- MassBalan.results$C_free_IVBMB*1000000
MassBalan.results$Class <- factor(MassBalan.results$Class)


F.1 <- 
ggplot(MassBalan.results, aes(x = C_free_VIVD, y = C_free_Fischer, group = Class))+
  geom_point(aes(color = Class), alpha=0.6, size=3) +
  xlab("C_free calculated by VIVD model (uM)") +
  ylab("C_free calculated by Fischer model (uM)")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        legend.position = c(0.2, 0.8),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))

F.2 <- 
ggplot(MassBalan.results, aes(x = C_free_VIVD, y = C_free_IVBMB, group = Class))+
  geom_point(aes(color = Class), alpha=0.6, size=3) +
  xlab("C_free calculated by VIVD model (uM)") +
  ylab("C_free calculated by IVBMB model (uM)")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))
F.3 <- 
ggplot(MassBalan.results, aes(x = C_free_IVBMB, y = C_free_Fischer, group = Class))+
  geom_point(aes(color = Class), alpha=0.6, size=3) +
  xlab("C_free calculated by IVBMB model (uM)") +
  ylab("C_free calculated by Fischer model (uM)")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))

tiff(width = 13.5, height = 4, "Figure_Scatter plot_MassbalaModel.tiff", units = "in", res=300)
F.1 + F.3 + F.2  +
  plot_layout(nrow = 1, widths = c(1/3, 1/3, 1/3))
dev.off()


# Bar chart
MassBalan.plot.dat <- data.frame(Chemical = rep(MassBalan.results$Compound, 3),
                                 Class = rep(MassBalan.results$Class, 3),
                                 C_free = c(MassBalan.results$C_free_Fischer, MassBalan.results$C_free_VIVD, MassBalan.results$C_free_IVBMB),
                                 Model = rep(c("Fischer", "VIVD", "IVBMB"), each = 40))
MassBalan.plot.dat$Model <- factor(MassBalan.plot.dat$Model, levels = c("Fischer", "VIVD", "IVBMB"))
ggplot(data=MassBalan.plot.dat, aes(x=Chemical, y=C_free, fill = Model)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + 
  facet_wrap(~ Class)+
  scale_y_log10(lim = c(1E-7, 1E+1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("Chemcials") +
  ylab("C_free (uM)")+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        legend.position = c(0.2, 0.8),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))