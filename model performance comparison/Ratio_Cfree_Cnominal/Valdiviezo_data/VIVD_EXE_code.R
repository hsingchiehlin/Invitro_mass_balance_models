library(ggplot2)
library(scales)
library(patchwork)

source("MCSim/function.R")
set_PATH()
makemod()

VIVD.Che.list <- read.csv("VIVD.Che.para.csv")

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

# cmp_types: 1:neutral, 2:monobase, 3:dibase, 4:monoacid, 5:diacid, 6:ampholyte

cmp_type = 1; # neutral

pKa1        = 5;
pKa2        = 6;
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
  
  write.csv(VIVD.Che.list, "VIVD.Che.sim.results_10uM.csv")
  file.remove(paste("MCSim/VIVD_dynamic_", i, ".in", sep=""))
  print(i)
}


for(i in 1:20){
  
  file.name <- paste("VIVD_dynamic_", i, ".in", sep="")
  Q_total <- paste("Q_total = ", 300*0.000001*0.000001*1, ";\n", sep="")
  
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

# cmp_types: 1:neutral, 2:monobase, 3:dibase, 4:monoacid, 5:diacid, 6:ampholyte

cmp_type = 1; # neutral

pKa1        = 5;
pKa2        = 6;
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
  
  write.csv(VIVD.Che.list, "VIVD.Che.sim.results_1uM.csv")
  file.remove(paste("MCSim/VIVD_dynamic_", i, ".in", sep=""))
  print(i)
}