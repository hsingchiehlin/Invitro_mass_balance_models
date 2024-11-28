source("MCSim/function.R")
set_PATH()
makemod()

VIVD.Che.list <- read.csv("VIVD.Che.para.csv")
VIVD.Cell.list <- read.csv("VIVD.Cell.para.csv")


model <- "INERIS_VIVD_dynamic.model"
makemcsim(model)

VIVD.results <- c()
VIVD.results.1um <- VIVD.results
#2%FBS_96WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 2,
                        Well = 96,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
                        )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 150*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_2FBS_96WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 150;    # (microL)
Sysdiameter  = 6.4;    # (mm)
Vwell        = 360;    # (microL), 1mm3 = 1microL
ncell        = 20000;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.02;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]

    print(i)
    
  }
  
VIVD.results <- rbind(VIVD.results, results)
  
}


#20%FBS_96WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 20,
                        Well = 96,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 150*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_20FBS_96WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 150;    # (microL)
Sysdiameter  = 6.4;    # (mm)
Vwell        = 360;    # (microL), 1mm3 = 1microL
ncell        = 20000;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.2;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_2FBS_96WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#2%FBS_384WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 2,
                        Well = 384,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 40*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_2FBS_384WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 40;    # (microL)
Sysdiameter  = 2.7;    # (mm)
Vwell        = 112;    # (microL), 1mm3 = 1microL
ncell        = 5600;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.02;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_2FBS_384WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#20%FBS_384WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 20,
                        Well = 384,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 40*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_20FBS_384WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 40;    # (microL)
Sysdiameter  = 2.7;    # (mm)
Vwell        = 112;    # (microL), 1mm3 = 1microL
ncell        = 5600;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.2;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_20FBS_384WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#VIVD.results$Conc <- rep(c("1 uM", "0.01 uM"), each = 4316)
VIVD.results$Cell.type[which(VIVD.results$Cell.type == "Ratcerebellargranulecell")] <- "Rat cerebellar granule cell"
VIVD.results$Q_cells <- VIVD.results$Q_cells * 1e6
VIVD.results$Q_media_free <- VIVD.results$Q_media_free * 1e6
VIVD.results$Q_media <- VIVD.results$Q_media * 1e6

write.csv(VIVD.results, "VIVD_results.csv")



# 0.01 um
VIVD.Che.list$CNOM_uM <- VIVD.Che.list$CNOM_uM*0.01
#2%FBS_96WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 2,
                        Well = 96,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 150*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_2FBS_96WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 150;    # (microL)
Sysdiameter  = 6.4;    # (mm)
Vwell        = 360;    # (microL), 1mm3 = 1microL
ncell        = 20000;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.02;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#20%FBS_96WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 20,
                        Well = 96,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 150*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_20FBS_96WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 150;    # (microL)
Sysdiameter  = 6.4;    # (mm)
Vwell        = 360;    # (microL), 1mm3 = 1microL
ncell        = 20000;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.2;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_2FBS_96WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#2%FBS_384WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 2,
                        Well = 384,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 40*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_2FBS_384WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 40;    # (microL)
Sysdiameter  = 2.7;    # (mm)
Vwell        = 112;    # (microL), 1mm3 = 1microL
ncell        = 5600;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.02;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_2FBS_384WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}


#20%FBS_384WELL
for (j in 1:length(VIVD.Cell.list$Cell)){
  
  Cell <- VIVD.Cell.list$Cell[j]
  Cell.Diameter <- VIVD.Cell.list$Diameter[j]
  Cell.Storage_lipids <- VIVD.Cell.list$Storage_lipids[j]
  Cell.Membrane_lipids <- VIVD.Cell.list$Membrane_lipids[j]
  Cell.Water_content <- VIVD.Cell.list$Water_content[j]
  
  results <- data.frame(Chemical = VIVD.Che.list$Chemical,
                        CAS = VIVD.Che.list$CAS,
                        Cell.type = Cell,
                        FBS = 20,
                        Well = 384,
                        Q_cells = NA,
                        Q_media_free = NA,
                        Q_media = NA,
                        F_free_media = NA
  )
  
  
  for(i in 1:length(VIVD.Che.list$Chemical)){
    
    Q_total <- paste("Q_total = ", 40*0.000001*0.000001*VIVD.Che.list$CNOM_uM[i], ";\n", sep="")
    
    
    file.name <- paste("VIVD_dynamic_20FBS_384WELL_", Cell, "_", i, ".in", sep="")
    
    cat("# Units:
# quantities:     mol
# volumes:        L
# concentrations: mol/L
# time:           minutes\n", file = file.name)
    sink(file.name, append=F)
    
    cat('OutputFile("sim.out");\n\n')
    
    cat('# Parameters

temp         = 37;     # culture temperature (degrees celsius)
Vmedium      = 40;    # (microL)
Sysdiameter  = 2.7;    # (mm)
Vwell        = 112;    # (microL), 1mm3 = 1microL
ncell        = 5600;\n')
    
    Celldiameter <- paste("Celldiameter = ", Cell.Diameter, ";\n\n", sep="")
    cat(Celldiameter)
    
    cat('# Cell composition data

pH_IW   = 7.0;     # intracellular water pH

f_NLrbc = 0.0017;  # ~        in red blood cells

f_NPrbc = 0.0029;  # ~        in red blood cells

f_IWrbc = 0.666;   # ~       for red blood cells

AP      = 5.09;    # acidic phospholipid concentration (mg/g)
AP_rbc  = 0.44;    # ~      in red blood cells

pH_mito = 8.0;     # pH mitochondria
pH_lyso = 4.0;     # pH lysosome

f_lyso  = 0.01;    # fraction (v/v) lysosome
f_mito  = 0.1;     # fraction (v/v) mitochondria\n')
    
    
    f_NL <- paste("f_NL    = ", Cell.Storage_lipids, ";\n", sep="")
    f_NP <- paste("f_NP    = ", Cell.Membrane_lipids, ";\n", sep="")
    f_IW <- paste("f_IW    = ", Cell.Water_content, ";\n\n", sep="")
    cat(f_NL)
    cat(f_NP)
    cat(f_IW)
    
    cat('# Culture medium composition

pH_medium     = 7.4;       # culture medium pH
f_serum       = 0.2;      # fraction (v/v) of serum in total media
pH_serum      = 7.4;       # FBS pH
f_FBS_albumin = 0.024;   # fraction (v/v) of albumin in FBS
f_FBS_NL      = 0.0019;   # fraction (v/v) of neutral lipid in FBS 


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
            Q_lyso, 0, 1440, 60);

}

End.\n')      
    
    sink()
    
    out <- mcsim(model, file.name)
    
    results$Q_cells[i] <- out$Q_cells[25]
    results$Q_media_free[i] <- out$Q_media_free[25]
    results$Q_media[i] <- out$Q_media[25]
    results$F_free_media[i] <- out$C_free_diss[25]/out$C_parent[25]
    
    print(i)
    file.remove(paste("VIVD_dynamic_20FBS_384WELL_", Cell, "_", i, ".in", sep=""))
    
  }

VIVD.results <- rbind(VIVD.results, results)

}
