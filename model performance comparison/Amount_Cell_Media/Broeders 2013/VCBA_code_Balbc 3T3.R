vcba.che.list <- read.csv("VCBA.che.list_Balbc 3T3.csv")

vcba.che.list$x1 <- NA
vcba.che.list$x2 <- NA
vcba.che.list$x3 <- NA
vcba.che.list$Cdcomp <- NA
vcba.che.list$ncells <- NA
vcba.che.list$Lt <- NA
vcba.che.list$St <- NA

cellType <- "3T3"

for(i in 1){
  
  ######################DESCRIPTION OF xdot function=DIFFERENTIAL EQUATIONS #######################
  xdot<-function(t,state,parameters) {
    with(as.list(c(state,parameters)), {
      #x1 total concentration in the medium
      #x2 concentration in the headspace
      #x3 concentration inside the cells
      rmax<-0
      Ksat<-0
      Cdcomp<-xx1/(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V)) #Concentration of dissolved compound in mol.l-1. Eq. 4. 
      #CONDITION 1) there are via cells; 
      # 1) i) there is no compound entering then: xdot3<- -kmet*xx3-weight_change*xx3
      # 1) ii) there is compound entering and so xdot3 is as represented in equation
      if(ncells>1) {
        W<-mcells/ncells #Wet Weight (g)
        Vcell<-Vcells*1E6/ncells #cm3
        DeltaC<-Cdcomp-(xx3/(MWcomp*BCF))
        if(DeltaC==0.0) {
          rexchange<-0.0
        }else {
          rexchange<-MWcomp*((Vcell ^(2/3))*rda*DeltaC+rmax*DeltaC/(DeltaC+Ksat))/W
        }
        xdot3<-rexchange-kmet*xx3-weight_change*xx3 #Function for concentration inside the cell (g.gww-1) ###
        #compound uptake or sendit back to the medium by cells in mol.l-1.s-1
        cells_up<--chemdead+(rexchange-kmet*xx3)*mcells/MWcomp/((1e3)*V) 
      } else {
        #CONDITION 2) cells are all dead
        xdot3<-0
        cells_up<--chemdead
      }
      # kgcomp mass transfer coefficient on the air (m.s-1)
      # klcomp mass transfer coefficient on the water film (m.s-1)
      # KGLcomp dimensionless gas-liquid distribution coefficient.
      Fawcomp<-kaw*(-Cdcomp+xx2/KGLcomp) # diffusive air-water exchange (mol.m-2.s-1) 
      Fdecomp<-kdeccomp*Cdcomp #medium decomposition in mol.l-1.s-1
      Fdecompa<-kdecacomp*xx2 #headspace decomposition in mol.l-1.s-1
      Flosses<-Fexch*xx2 #headspace losses. For Fexch is set to 0, this is 0, but can be used when modeling cross contamination
      ########## Equation for total concentration the medium (mol.l-1.s1) #########
      xdot1<-(P*Fawcomp-Fdecomp-cells_up)
      ########### Equation for concentration in the headspace (mol.l1.s-1) #########
      xdot2<-(-Ph*Fawcomp-Fdecompa-Flosses)
      list(c(xdot1,xdot2,xdot3))
    })
  }
  ########################################################################## 
  #END OF DIFFERENTIAL EQUATIONS 
  ######################################################################################
  library(deSolve) #Program to solve differential equations
  ##################################### DESCRIPTION OF CONSTANTS ##########################
  ######## WELL GEOMETRY ############
  # To return AssayVol, P, Ph and SP
  PI<-pi
  #Well Description in meters
  diaTop<-0.02228 #Diameter of the top of the well
  diaBot<-0.0222 #Dimater of the bottom of the well
  Depth<-0.01753#Depth of the well
  V<-2.0E-6 #Volume of medium in the well 0.0003 (m3)
  ##Equations
  H<-diaTop*Depth/(diaTop-diaBot)
  h<-H-Depth
  #Total Volume well:
  vVessel<-0.0000069 
  vt<-PI*h*((0.5*diaBot)^2)/3
  hh<-((V+vt)*12*h*h/(PI*diaBot*diaBot))^(1/3)
  xx<-diaBot*hh*0.5/h
  llt<-(hh^2+xx^2)^(1/2)
  lli<-(h^2+(diaBot*0.5)^2)^(1/2)
  SPt<-PI*xx*llt
  SPi<-PI*diaBot*0.5*lli
  #surface of the plastic in contact with the medium
  SP<-SPt-SPi+PI*((diaBot*0.5)^2) 
  #cell assay surface:
  AS<-PI*(xx^2) 
  P<-AS/V #element to simplify differential equation xdot1
  #head space volume (m3)
  Vh=vVessel-V
  Ph<-AS/Vh #element to simplify differential equation xdot2
  ####### AIR-WATER EXCHANGE #########
  # To return kdeccomp, kdecacomp, kgcomp,KGLcomp, Klcomp and kaw
  # Constant values introduced knime
  kdeccomp<-0 #Compound's Rate of Degradation in water in s-1
  kdecacomp<-0 #Compound's Rate of Degradation in water in s-1
  Te<- 37+273.15 #Air temperature from Celsius to K. Experiment temperature constant at 37 C.
  MWair<-28.8 #Air molecular weight
  MWcomp<-vcba.che.list$MW[i] #Compound Molecular Weight
  Press<-1 #Air Pressure in atm.
  Svair<-20.1 #atomic diffusion air
  Svcomp<-90.96 #atomic diffusion
  RR<-8.3144 #Universal gas constant kJ(mol.K)-1
  H37<-vcba.che.list$H37[i] #Henry??s constant at 37 C
  vb<-vcba.che.list$vb[i] #molar volume at its normal boiling point (cm^3/g mol)
  kgH2O<-3.0e-3 #mass transfer coefficient for water m*s-1
  kLCO2<-4.1e-2 #mass transfer coefficient of CO2 in the Water side
  fi<-2.6 #association factor or organic solutes diffusing into water
  #Parameter for the exchange of compound between wells. Null for now.
  Fexch<-0 
  ##Evaporation and Degradation equations
  #Organic compound gas phase diffusion coefficient (m2/s). From Fuller 1966 
  valgas<-1E-7*((MWair+MWcomp)/(MWair*MWcomp))^0.5 
  valgas<-valgas/(Press*((Svair^0.33)+(Svcomp^0.33))^2) 
  DGcomp<-valgas*Te^(1.75) #Diffusion coefficient in the air.
  DGw<-1.23655e-9*Te^(1.75) #Diffusion coefficients for water in air 
  kgcomp<-kgH2O*(DGcomp/DGw)^0.67 #Gas phase mass transfer coefficient.
  KGLcomp<-H37/(RR*Te) #Gas-Liquid partition Coefficient (dimensionless)
  muw<-0.6913 #simplified water viscosity in cP at 37C. 
  denw<-0.933 #simplified water density in gr/cm3 at 37C. 
  Fi<-7.4E-12*(fi*MWcomp)^0.5/(vb)^0.6 
  DLcomp<-Fi*Te/muw #Liquid phase diffusion coefficients (m2/s). 
  SCAcomp<-muw*1e-3/(1e3*denw*DLcomp) # Schimdt number. 
  klcomp<-kLCO2*(SCAcomp/600)^(-0.5) #mass transfer coefficient on the water film .
  kaw<-(kgcomp*KGLcomp*klcomp)/(klcomp+kgcomp*KGLcomp) #mass transfer coefficient. 
  ######### BCF AND PARTITION #########
  # To return mass fractions and densities fo each component Kp, Ks, Kl, KL, KP, BCF and rad
  logkow<-vcba.che.list$logkow[i] #logarithm of partition octanol-water
  #mass fractions in % weight
  faq<-0.614 #3T3 cell line aqueous fraction
  fL<-0.142 #3T3 cell line lipidic fraction
  fP<-0.244 #3T3 cell line proteic fraction
  #densities in g/L
  rhoaq<-1000 #aqueous phase density
  rhoL<-900 #lipid phase density
  rhoP<-1350 #proteic phase density
  #Lipid and protein in medium serum. Mind the % of serum supplementing the medium
  S0<-0.0234*1 #protein content for 5 % serum -> 10% *2
  L0<-0.08*1 #lipid content for 5% serum
  ##Equations for partition
  #plastic well partition
  Kp<-10^(0.97*logkow-6.94)
  #serum well partition
  if (logkow <1.09) {
    vals<-1.31
  } else if (logkow>=1.09&&logkow<=4.6) {
    vals<-0.57*logkow+0.69
  } else if (logkow>4.6) {
    vals<-logkow-1.3
  }
  Ks<-10^(vals-1.178)
  #lipid partition
  Kl<-10^(1.25*logkow-3.70)
  ##cell partition##
  LC<-170.4
  SC<-4.4 
  KL<-LC*Kl
  KP<-SC*Ks
  #Bioconcentration factor 
  #Does not contain Lipid and protein concentration as these parameters are already contained in KL and KP
  BCF<-(faq/rhoaq)+(fL*KL/rhoL)+(fP*KP/rhoP)
  alldiss<-function(Lt=Lt,St=St) {
    disss<-(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V))
    disscs<-(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V))/(Ks*St)
    disscl<-(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V))/(Kl*Lt)
    disscp<-(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V))/(Kp*SP/V)
    allinone<-c(disss,disscs,disscl,disscp)
    return(allinone)
  }
  ###Cell Permeability Equations###
  logp<--1.1711+0.98*logkow-0.0011*MWcomp # QSAR for log cell permeability(cm*h-1)
  Per<-10^logp # Cell permeability still in cm*h-1
  rda<-Per/36000 # Rate of uptake=Cell permeability from cm*h-1 to L*cm-2*s-1 
  #Herein is considered rate of uptake is the same of elimination. If this assumption change parameters have to be revised
  ####Metabolism####
  kmet<-0 #3T3 are considered metabolic imcompetent cells
  ################################## END CONSTANTS DESCRIPTION##############################
  coreModel<-function(ci,cell_i,numCells,nec,kt) {
    #Parameters for leslie matrix. Already structured in arrays of 1 row and 4 columns
    di<-c(2880,2880,2880,2880) # duration at each stage in min (9.63,3.65,3.45,2.26)
    zi<-c(1e-15,1e-15,1e-15,1e-15) #mortality at each stage (0.005,0.005,0.04,0.04)
    di<-array(di,dim=c(1,4)) 
    zi<-array(zi,dim=c(1,4))
    ##Initial values for the time loop. Once the program runs once the loop, these values stop being used and actualised for the time run(kl)
    numCells<-15000 ##Initial cell number
    ncells<-numCells
    N<-c(0,100,0,0)*ncells/100 #fraction of cells in each cell cycle phase (50.7,19.2,18.1798,11.9202)
    N<-array(N,dim=c(4))
    mcells_old<-2.08E-9*N[1]+2.4E-9*N[2]+2.4E-9*N[3]+2.4E-9*N[4] #Cell Mass per cell cycle phase In g
    mcells<-mcells_old
    Vcells<-1.73E-15*N[1]+2.4E-15*N[2]+2.4E-15*N[3]+2.4E-15*N[4] #Cell volume per cell cycle phase in m^3
    Lt<-L0
    St<-S0
    weight_change<-0
    chemdead<-0
    x1<-ci
    x2<-0
    x3<-cell_i
    #Definition of the Duration
    totalTime<-48 #Time in h
    kl_fin<-totalTime #in h
    val<-array(kl_fin) #array for cell growth
    dead_cells<-array(kl_fin) #array for cell death
    ###################################### TIME CYCLE ##############################
    for (kl in 1:kl_fin) {
      t<-seq((kl-1)*3600,kl*3600,by=1) #time in seconds but making steps of 60 s
      ##Differential equations solving
      parameters<-
        c(Ks=Ks,St=St,Kl=Kl,Kp=Kp,Lt=Lt,SP=SP,V=V,kgcomp=kgcomp,KGLcomp=KGLcomp,klcomp=klcomp,ncells=ncells,mcells=mcells,Vcells=Vcells,faq=faq,rhoaq=rhoaq,fL=fL,
          KL=KL,rhoL=rhoL,fP=fP,KP=KP,rhoP=rhoP,MWcomp=MWcomp,rda=rda,kmet=kmet,weight_change=weight_change,chemdead=chemdead,kdeccomp=kdeccomp,P=P,
          kdecacomp=kdecacomp,Fexch=Fexch,Ph=Ph)
      state<-c(xx1=x1,xx2=x2,xx3=x3)
      out<-ode(y=state,times=t,func=xdot,parms=parameters,method="radau",atol=1e-4,rtol=1e-4,hmax=1)
      out_row<-length(out[,1])
      out_col<-length(out[1,])
      out_<-matrix(data=out,nrow=out_row,ncol=out_col)
      ##Mortality NEC and Kt
      cq<-mean(out_[,4])
      val2 <- cq - nec
      val1<-max(0,val2)
      za<-zi+kt*val1
      pii<-exp(-za)
      gamma<-(1-pii)*(pii^(di-1))/(1-pii^di)
      PS<-pii*(1-gamma)
      GS<-pii*gamma
      ##density dependence HepG2. Needs to be revised for other cell lines may also have density dependence. To run take the commen-out
      if(cellType =="HepG2") {
        Fii<-3.8062
        dvol<-35.0676e-5
        F<-Fii*exp(-ncells/dvol)
      } else {
        F<-1e-15 #optimised fecundity 
      }
      ##Leslie Matrix##
      numSteps<-4 #how many cell phases the cell line has; Ex:3T3 has 4 and HepaRG 1.
      if(numSteps==1) {
        dead_cells[kl]<-N[1]*(1-(PS[1]+GS[1]))
        LE<-c(PS[1])
        LE<-array(LE,dim=c(1,1))
        N[1]=LE*N[1]
        #The PS and GS for other cell phases are set to 0 so in valdead it will not count with these phases
        PS[2:4]<-0 
        GS[1:4]<-0
      } 
      #Because the cell we are using here is 3T3 and it has 4 phases, then it will run this option
      if(numSteps==4) { 
        dead_cells[kl]<-N[1]*(1-(PS[1]+GS[1]))+N[2]*(1-(PS[2]+GS[2]))+N[3]*(1-(PS[3]+GS[3]))+N[4]*(1-(PS[4]+GS[4]))
        LE<-c(PS[1],GS[1],0,0,0,PS[2],GS[2],0,0,0,PS[3],GS[3],F,0,0,PS[4])
        LE<-array(LE,dim=c(4,4))
        N<-LE%*%N
      }
      ##Number of cells at each moment##
      val[kl]<-sum(N)
      ncells<-val[kl]
      ##When all cells die
      if(ncells<1) {
        ncells<-0
        out_[,4]<-0
        Wcells<-0
        mcells<-0
        Vcells<-0
      } else {
        mcells<-2.08E-9*N[1]+2.4E-9*N[2]+2.4E-9*N[3]+2.4E-9*N[4] # recalculate the total mass of cells
        Vcells<-1.73E-15*N[1]+2.4E-15*N[2]+2.4E-15*N[3]+2.4E-15*N[4] # recalculate the total volume of cells 
        Wcell<-mcells/ncells
        weight_change<-(mcells-mcells_old)/mcells/ncells/3600 
        mcells_old<-mcells #So for the new cycle the herein calculated mcells will be the mcells(kl-1)=mcells_old
        ##Return of chemical, lipid and protein to medium from cells death##
        chemdead<-dead_cells[kl]*out_[out_row,4]*Wcell*1E-3/MWcomp/V/3600 
      }
      valdead<-N[1]*(1-(PS[1]+GS[1]))*2.08E-9+N[2]*(1-(PS[2]+GS[2]))*2.40E-9+N[3]*(1-(PS[3]+GS[3]))*2.4E-9+N[4]*(1-
                                                                                                                   (PS[4]+GS[4]))*2.4E-9
      Lt<-Lt+(fL*valdead*1E-3)/V
      #albumin MW= 66400 g/mol
      St<-St+(fP*valdead)/66400/V
      #output of the differential equations and for cycle (concentration in air, medium and intracellularly)
      x0<-out_[out_row,1]
      x1<-out_[out_row,2]
      x2<-out_[out_row,3]
      x3<-out_[out_row,4]
      Cdcomp<-x1/(1+(Ks*St)+(Kl*Lt)+(Kp*SP/V))
    }
    # Start new cycle 
    return(c(x1,x2,x3,Cdcomp,ncells,Lt,St))
  }
  ################################# END LOOP #################################
  ###Description of the procedures for different Modes, can only run one mode at a time. 
  
  vcba.che.list[i, 8:14] <- 
    coreModel(ci = vcba.che.list$CNOM_M[i], cell_i = 0, numCells = 15000, nec = 1, kt = 0) #ci unit: M
  
  print(i)
}

V_cell <- 2.4e-12 #L (2.4e-15 m3)

vcba.che.list$Cell_nmole <- vcba.che.list$x3*V_cell*vcba.che.list$ncells*(1E+6)

write.csv(vcba.che.list, "vcba.che.list_results_Balbc 3T3.csv")



