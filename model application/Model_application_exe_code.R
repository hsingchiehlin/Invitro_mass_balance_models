library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(weights)
library(wCorr)
library(httk) # VER. 2.3.0
library(ggh4x)
library(patchwork)
library(spatstat)

# ----------- load data -----------
invivo_fup <- read.csv("invivo_fup.csv")

#original in vivo POD
origpod <- read.csv(file.path("POD database_15che","15 chemicals orig POD x Css.csv"))

#original in vivo POD + DAF
poddaf <- read.csv(file.path("POD database_15che","original POD x DAF x Css for 15 chemicals.csv"))

#BBMD
BMD <- read.csv(file.path("POD database_15che","BMD x Css for 15 chemicals.csv"))
BMD$Chemical <- origpod$Chemical
BMD$Classification <- origpod$Classification

#HED
HED <- read.csv(file.path("POD database_15che","HED x Css for 15 chemicals.csv"))
HED$Chemical <- origpod$Chemical
HED$Classification <- origpod$Classification

for(i in 1:15){
  
  origpod$Css_h[i] <- calc_analytic_css(chem.cas = origpod$CAS[i], model = "3compartmentss", parameterize.args = list(adjusted.Funbound.plasma=FALSE,adjusted.Clint = FALSE))
  poddaf$Css_h[i] <- calc_analytic_css(chem.cas = poddaf$CAS[i], model = "3compartmentss", parameterize.args = list(adjusted.Funbound.plasma=FALSE,adjusted.Clint = FALSE))
  BMD$Css_h[i] <- calc_analytic_css(chem.cas = BMD$CAS[i], model = "3compartmentss", parameterize.args = list(adjusted.Funbound.plasma=FALSE,adjusted.Clint = FALSE))
  HED$Css[i] <- calc_analytic_css(chem.cas = HED$CAS[i], model = "3compartmentss", parameterize.args = list(adjusted.Funbound.plasma=FALSE,adjusted.Clint = FALSE))
}

origpod$PODCss <- origpod$value
poddaf$poddaf.Css <- poddaf$VALUE

BMD$BMD50.Css <- BMD$X50.
BMD$BMD5.Css <- BMD$X5.
BMD$BMD95.Css <- BMD$X95.

HED$HED50.Css <- HED$X50.
HED$HED5.Css <- HED$X5.
HED$HED95.Css <- HED$X95.

invitro_fub <- read.csv("invitro_fub_cfree_nonminal.csv")
#ToxCast active
toxcast <- read.csv(file.path("POD database_15che","toxcast ac50.csv"))
tox.5perc <- data.frame(sapply(toxcast, quantile, prob=0.05, na.rm = T))
tox.5perc$CAS <- origpod$CAS
colnames(tox.5perc)[1] <- "ToxCast"

#iPSC+LCL
cell.t <- read.csv(file.path("POD database_15che","cellline POD distribution for 15 chemicals.csv"))
cell.min <- data.frame(sapply(cell.t[2:16], min, na.rm=T))
cell.min$CAS <- origpod$CAS
colnames(cell.min)[1] <- "iPSCLCL"

constant.POD <- data.frame(CAS = invitro_fub$CAS[1:15],
                           iPSCLCL = 1)

constant.POD.2 <- data.frame(CAS = invitro_fub$CAS[1:15],
                             iPSCLCL = 0.01)

#Cardio
#cardiomyocytes (single donor)
cardio.df <- read.csv(file.path("POD database_15che","Cardio POD for 15 chemicals.csv"))
cardio <- full_join(cardio.df, origpod[,c(1,2)], by="CAS")
cardio <- cardio[-which(is.na(cardio$Chemical.y)),]


for(i in 1:15){
  
  tox.5perc$ToxCast[i] <- tox.5perc$ToxCast[i]/origpod$Css_h[which(tox.5perc$CAS[i] == origpod$CAS)]
  
  cell.min$iPSCLCL[i] <- cell.min$iPSCLCL[i]/origpod$Css_h[which(cell.min$CAS[i] == origpod$CAS)]
  
  constant.POD$iPSCLCL[i] <- constant.POD$iPSCLCL[i]/origpod$Css_h[which(constant.POD$CAS[i] == origpod$CAS)]
  
  constant.POD.2$iPSCLCL[i] <- constant.POD.2$iPSCLCL[i]/origpod$Css_h[which(constant.POD.2$CAS[i] == origpod$CAS)]
  
  cardio$Min[i] <- cardio$Min[i]/origpod$Css_h[which(cardio$CAS[i] == origpod$CAS)]
  
}


# -------------- without any adjustment --------------
##traditional POD
pod.tox.df <- full_join(origpod, tox.5perc, by="CAS") # x=ToxCast, y=PODCss, xlab("ToxCast POD (ÂµM)")+ylab(bquote(Reg~POD[A]~x~Css~(ÂµM)))+
pod.tox.s <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Spearman"))), digits=2)
pod.tox.p <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Pearson"))), digits=2)
pod.tox.mae <- mean(abs(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast)))
pod.tox.me <- mean(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast))

##traditional POD x DAF
poddaf.tox.df <- full_join(poddaf, tox.5perc, by="CAS") # x=ToxCast, y=poddaf.Css, xlab("ToxCast POD (ÂµM)")+ylab(bquote(Reg~POD[H]~x~Css~(ÂµM)))+
poddaf.tox.s <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Spearman"))), digits=2)
poddaf.tox.p <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Pearson"))), digits=2)
poddaf.tox.mae <- mean(abs(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast)))
poddaf.tox.me <- mean(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast))

##BBMD
bbmd.tox.df <- full_join(BMD, tox.5perc, by="CAS") 
wt.bmd.tox <- 1/(log10(bbmd.tox.df$BMD95.Css)-log10(bbmd.tox.df$BMD5.Css))^2
bbmd.tox.w.s <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Spearman"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.w.p <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Pearson"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.mae <- weighted.mean(abs(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast)), wt.bmd.tox)
bbmd.tox.me <- weighted.mean(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast), wt.bmd.tox)

##HED
hed.tox.df <- full_join(HED, tox.5perc, by="CAS") 
wt.hed.tox <- 1/(log10(hed.tox.df$HED95.Css)-log10(hed.tox.df$HED5.Css))^2
hed.tox.w.s <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Spearman"), weights=wt.hed.tox)), digits=2)
hed.tox.w.p <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Pearson"), weights=wt.hed.tox)), digits=2)
hed.tox.mae <- weighted.mean(abs(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast)), wt.hed.tox)
hed.tox.me <- weighted.mean(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast), wt.hed.tox)


##in vivo POD vs iPSC+LCL POD scatter plot
##traditional POD
pod.cell.df <- full_join(origpod, cell.min, by="CAS") # x=iPSCLCL, y=PODCss, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(Reg~POD[A]~x~Css~(ÂµM)))+
pod.cell.s <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.cell.p <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.cell.mae <- mean(abs(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL)))
pod.cell.me <- mean(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL))

##traditional POD x DAF
poddaf.cell.df <- full_join(poddaf, cell.min, by="CAS") # x=iPSCLCL, y=poddaf.Css, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(Reg~POD[H]~x~Css~(ÂµM)))+
poddaf.cell.s <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.cell.p <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.cell.mae <- mean(abs(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL)))
poddaf.cell.me <- mean(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL))

##BBMD
bbmd.cell.df <- full_join(BMD, cell.min, by="CAS") 
wt.bmd.cell <- 1/(log10(bbmd.cell.df$BMD95.Css)-log10(bbmd.cell.df$BMD5.Css))^2
bbmd.cell.w.s <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.w.p <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.mae <- weighted.mean(abs(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL)),wt.bmd.cell)
bbmd.cell.me <- weighted.mean(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL),wt.bmd.cell)

##HED
hed.cell.df <- full_join(HED, cell.min, by="CAS") 
wt.hed.cell <- 1/(log10(hed.cell.df$HED95.Css)-log10(hed.cell.df$HED5.Css))^2
hed.cell.w.s <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.cell)), digits=2)
hed.cell.w.p <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.cell)), digits=2)
hed.cell.mae <- weighted.mean(abs(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL)),wt.hed.cell)
hed.cell.me <- weighted.mean(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL),wt.hed.cell)

##in vivo POD vs Cardio POD scatter plot
##traditional POD
pod.cardio.df <- full_join(origpod, cardio, by="CAS") 
pod.car.s <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Spearman")), digits=2)
pod.car.p <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Pearson")), digits=2)
pod.car.mae <- mean(abs(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min)))
pod.car.me <- mean(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min))

##traditional POD x DAF
poddaf.cardio.df <- full_join(poddaf, cardio, by="CAS") 
poddaf.car.s <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Spearman")), digits=2)
poddaf.car.p <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Pearson")), digits=2)
poddaf.car.mae <- mean(abs(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min)))
poddaf.car.me <- mean(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min))

##BBMD
bbmd.cardio.df <- full_join(BMD, cardio, by="CAS") 
wt.bmd.cardio <- 1/(log10(bbmd.cardio.df$BMD95.Css)-log10(bbmd.cardio.df$BMD5.Css))^2
bbmd.car.w.s <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Spearman"), weights=wt.bmd.cardio), digits=2)
bbmd.car.w.p <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Pearson"), weights=wt.bmd.cardio), digits=2)
bbmd.car.mae <- weighted.mean(abs(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min)),wt.bmd.cardio)
bbmd.car.me <- weighted.mean(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min),wt.bmd.cardio)

##HED
hed.cardio.df <- full_join(HED, cardio, by="CAS")
wt.hed.cardio <- 1/(log10(hed.cardio.df$HED95.Css)-log10(hed.cardio.df$HED5.Css))^2
hed.car.w.s <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Spearman"), weights=wt.hed.cardio), digits=2)
hed.car.w.p <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Pearson"), weights=wt.hed.cardio), digits=2)
hed.car.mae <- weighted.mean(abs(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min)),wt.hed.cardio)
hed.car.me <- weighted.mean(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min),wt.hed.cardio)


##in vivo POD vs constant POD scatter plot
##traditional POD
pod.constant.df <- full_join(origpod, constant.POD, by="CAS") 
pod.constant.s <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae <- mean(abs(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL)))
pod.constant.me <- mean(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df <- full_join(poddaf, constant.POD, by="CAS") 
poddaf.constant.s <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae <- mean(abs(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL)))
poddaf.constant.me <- mean(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL))

##BBMD
bbmd.constant.df <- full_join(BMD, constant.POD, by="CAS") 
wt.bmd.constant <- 1/(log10(bbmd.constant.df$BMD95.Css)-log10(bbmd.constant.df$BMD5.Css))^2
bbmd.constant.w.s <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.w.p <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.mae <- weighted.mean(abs(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL)),wt.bmd.constant)
bbmd.constant.me <- weighted.mean(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL),wt.bmd.constant)

##HED
hed.constant.df <- full_join(HED, constant.POD, by="CAS") # 
wt.hed.constant <- 1/(log10(hed.constant.df$HED95.Css)-log10(hed.constant.df$HED5.Css))^2
hed.constant.w.s <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant)), digits=2)
hed.constant.w.p <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant)), digits=2)
hed.constant.mae <- weighted.mean(abs(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL)),wt.hed.constant)
hed.constant.me <- weighted.mean(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL),wt.hed.constant)


##in vivo POD vs constant POD scatter plot (conc = 0.01)
##traditional POD
pod.constant.df.2 <- full_join(origpod, constant.POD.2, by="CAS") 
pod.constant.s.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae.2 <- mean(abs(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL)))
pod.constant.me.2 <- mean(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df.2 <- full_join(poddaf, constant.POD.2, by="CAS") 
poddaf.constant.s.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae.2 <- mean(abs(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL)))
poddaf.constant.me.2 <- mean(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL))

##BBMD
bbmd.constant.df.2 <- full_join(BMD, constant.POD.2, by="CAS") 
wt.bmd.constant.2 <- 1/(log10(bbmd.constant.df.2$BMD95.Css)-log10(bbmd.constant.df.2$BMD5.Css))^2
bbmd.constant.w.s.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.w.p.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.mae.2 <- weighted.mean(abs(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL)),wt.bmd.constant.2)
bbmd.constant.me.2 <- weighted.mean(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL),wt.bmd.constant.2)

##HED
hed.constant.df.2 <- full_join(HED, constant.POD.2, by="CAS") 
wt.hed.constant.2 <- 1/(log10(hed.constant.df.2$HED95.Css)-log10(hed.constant.df.2$HED5.Css))^2
hed.constant.w.s.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant.2)), digits=2)
hed.constant.w.p.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant.2)), digits=2)
hed.constant.mae.2 <- weighted.mean(abs(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL)),wt.hed.constant.2)
hed.constant.me.2 <- weighted.mean(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL),wt.hed.constant.2)

POD.data.all_withoutadj <- data.frame(Group = "Without adjustment",
                                      CAS = c(pod.tox.df$CAS, poddaf.tox.df$CAS, bbmd.tox.df$CAS, hed.tox.df$CAS,
                                              pod.cell.df$CAS, poddaf.cell.df$CAS, bbmd.cell.df$CAS, hed.cell.df$CAS,
                                              pod.cardio.df$CAS, poddaf.cardio.df$CAS, bbmd.cardio.df$CAS, hed.cardio.df$CAS,
                                              pod.constant.df$CAS, poddaf.constant.df$CAS, bbmd.constant.df$CAS, hed.constant.df$CAS,
                                              pod.constant.df.2$CAS, poddaf.constant.df.2$CAS, bbmd.constant.df.2$CAS, hed.constant.df.2$CAS),
                                      x.val = c(pod.tox.df$ToxCast, poddaf.tox.df$ToxCast, bbmd.tox.df$ToxCast, hed.tox.df$ToxCast,
                                                pod.cell.df$iPSCLCL, poddaf.cell.df$iPSCLCL, bbmd.cell.df$iPSCLCL, hed.cell.df$iPSCLCL,
                                                pod.cardio.df$Min, poddaf.cardio.df$Min, bbmd.cardio.df$Min, hed.cardio.df$Min,
                                                pod.constant.df$iPSCLCL, poddaf.constant.df$iPSCLCL, bbmd.constant.df$iPSCLCL, hed.constant.df$iPSCLCL,
                                                pod.constant.df.2$iPSCLCL, poddaf.constant.df.2$iPSCLCL, bbmd.constant.df.2$iPSCLCL, hed.constant.df.2$iPSCLCL),
                                      x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 60),
                                      y.val = c(pod.tox.df$PODCss, poddaf.tox.df$poddaf.Css, bbmd.tox.df$BMD50.Css, hed.tox.df$HED50.Css,
                                                pod.cell.df$PODCss, poddaf.cell.df$poddaf.Css, bbmd.cell.df$BMD50.Css, hed.cell.df$HED50.Css,
                                                pod.cardio.df$PODCss, poddaf.cardio.df$poddaf.Css, bbmd.cardio.df$BMD50.Css, hed.cardio.df$HED50.Css,
                                                pod.constant.df$PODCss, poddaf.constant.df$poddaf.Css, bbmd.constant.df$BMD50.Css, hed.constant.df$HED50.Css,
                                                pod.constant.df.2$PODCss, poddaf.constant.df.2$poddaf.Css, bbmd.constant.df.2$BMD50.Css, hed.constant.df.2$HED50.Css),
                                      y.val.min = c(rep(NA, 30), bbmd.tox.df$BMD5.Css, hed.tox.df$HED5.Css,
                                                    rep(NA, 30), bbmd.cell.df$BMD5.Css, hed.cell.df$HED5.Css,
                                                    rep(NA, 30), bbmd.cardio.df$BMD5.Css, hed.cardio.df$HED5.Css,
                                                    rep(NA, 30), bbmd.constant.df$BMD5.Css, hed.constant.df$HED5.Css,
                                                    rep(NA, 30), bbmd.constant.df.2$BMD5.Css, hed.constant.df.2$HED5.Css),
                                      y.val.max = c(rep(NA, 30), bbmd.tox.df$BMD95.Css, hed.tox.df$HED95.Css,
                                                    rep(NA, 30), bbmd.cell.df$BMD95.Css, hed.cell.df$HED95.Css,
                                                    rep(NA, 30), bbmd.cardio.df$BMD95.Css, hed.cardio.df$HED95.Css,
                                                    rep(NA, 30), bbmd.constant.df$BMD95.Css, hed.constant.df$HED95.Css,
                                                    rep(NA, 30), bbmd.constant.df.2$BMD95.Css, hed.constant.df.2$HED95.Css),
                                      y.type = rep(rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"), each =15), 5))
RHO.data.all_withoutadj <- data.frame(Group = "Without adjustment",
                                      x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                      y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                      rho.value = c(pod.tox.s, poddaf.tox.s, bbmd.tox.w.s, hed.tox.w.s,
                                                    pod.cell.s, poddaf.cell.s, bbmd.cell.w.s, hed.cell.w.s,
                                                    pod.car.s, poddaf.car.s, bbmd.car.w.s, hed.car.w.s,
                                                    pod.constant.s, poddaf.constant.s, bbmd.constant.w.s, hed.constant.w.s,
                                                    pod.constant.s.2, poddaf.constant.s.2, bbmd.constant.w.s.2, hed.constant.w.s.2))
R.data.all_withoutadj <- data.frame(Group = "Without adjustment",
                                    x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                    y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                    rho.value = c(pod.tox.p, poddaf.tox.p, bbmd.tox.w.p, hed.tox.w.p,
                                                  pod.cell.p, poddaf.cell.p, bbmd.cell.w.p, hed.cell.w.p,
                                                  pod.car.p, poddaf.car.p, bbmd.car.w.p, hed.car.w.p,
                                                  pod.constant.p, poddaf.constant.p, bbmd.constant.w.p, hed.constant.w.p,
                                                  pod.constant.p.2, poddaf.constant.p.2, bbmd.constant.w.p.2, hed.constant.w.p.2))

mae.data.all_withoutadj <- data.frame(Group = "Without adjustment",
                                      x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                      y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                      rho.value = c(pod.tox.mae, poddaf.tox.mae, bbmd.tox.mae, hed.tox.mae,
                                                    pod.cell.mae, poddaf.cell.mae, bbmd.cell.mae, hed.cell.mae,
                                                    pod.car.mae, poddaf.car.mae, bbmd.car.mae, hed.car.mae,
                                                    pod.constant.mae, poddaf.constant.mae, bbmd.constant.mae, hed.constant.mae,
                                                    pod.constant.mae.2, poddaf.constant.mae.2, bbmd.constant.mae.2, hed.constant.mae.2))
me.data.all_withoutadj <- data.frame(Group = "Without adjustment",
                                     x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                     y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                     rho.value = c(pod.tox.me, poddaf.tox.me, bbmd.tox.me, hed.tox.me,
                                                   pod.cell.me, poddaf.cell.me, bbmd.cell.me, hed.cell.me,
                                                   pod.car.me, poddaf.car.me, bbmd.car.me, hed.car.me,
                                                   pod.constant.me, poddaf.constant.me, bbmd.constant.me, hed.constant.me,
                                                   pod.constant.me.2, poddaf.constant.me.2, bbmd.constant.me.2, hed.constant.me.2))



# -------------- with adjustment: fup and fub from Fischer --------------
invitro_fub <- read.csv("invitro_fub_cfree_nonminal.csv")
#ToxCast active
toxcast <- read.csv(file.path("POD database_15che","toxcast ac50.csv"))
tox.5perc <- data.frame(sapply(toxcast, quantile, prob=0.05, na.rm = T))
tox.5perc$CAS <- origpod$CAS
colnames(tox.5perc)[1] <- "ToxCast"

#iPSC+LCL
cell.t <- read.csv(file.path("POD database_15che","cellline POD distribution for 15 chemicals.csv"))
cell.min <- data.frame(sapply(cell.t[2:16], min, na.rm=T))
cell.min$CAS <- origpod$CAS
colnames(cell.min)[1] <- "iPSCLCL"

constant.POD <- data.frame(CAS = invitro_fub$CAS[1:15],
                           iPSCLCL = 1)

constant.POD.2 <- data.frame(CAS = invitro_fub$CAS[1:15],
                             iPSCLCL = 0.01)

#Cardio
#cardiomyocytes (single donor)
cardio.df <- read.csv(file.path("POD database_15che","Cardio POD for 15 chemicals.csv"))
cardio <- full_join(cardio.df, origpod[,c(1,2)], by="CAS")
cardio <- cardio[-which(is.na(cardio$Chemical.y)),]


for(i in 1:15){
  
  tox.5perc$ToxCast[i] <- (tox.5perc$ToxCast[i]*invitro_fub$Fischer[which(tox.5perc$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(tox.5perc$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(tox.5perc$CAS[i] == invivo_fup$CAS)])
  
  cell.min$iPSCLCL[i] <- (cell.min$iPSCLCL[i]*invitro_fub$Fischer[which(cell.min$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cell.min$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cell.min$CAS[i] == invivo_fup$CAS)])
  
  constant.POD$iPSCLCL[i] <- (constant.POD$iPSCLCL[i]*invitro_fub$Fischer[which(constant.POD$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD$CAS[i] == invivo_fup$CAS)])
  
  constant.POD.2$iPSCLCL[i] <- (constant.POD.2$iPSCLCL[i]*invitro_fub$Fischer[which(constant.POD.2$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD.2$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD.2$CAS[i] == invivo_fup$CAS)])
  
  cardio$Min[i] <- (cardio$Min[i]*invitro_fub$Fischer[which(cardio$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cardio$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cardio$CAS[i] == invivo_fup$CAS)])
  
}



##traditional POD
pod.tox.df <- full_join(origpod, tox.5perc, by="CAS") # x=ToxCast, y=PODCss, xlab("ToxCast POD (ÂµM)")+ylab(bquote(Reg~POD[A]~x~Css~(ÂµM)))+
pod.tox.df <- pod.tox.df[-(which(is.na(pod.tox.df$PODCss*pod.tox.df$ToxCast))),]
pod.tox.s <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Spearman"))), digits=2)
pod.tox.p <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Pearson"))), digits=2)
pod.tox.mae <- mean(abs(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast)))
pod.tox.me <- mean(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast))

##traditional POD x DAF
poddaf.tox.df <- full_join(poddaf, tox.5perc, by="CAS") # x=ToxCast, y=poddaf.Css, xlab("ToxCast POD (ÂµM)")+ylab(bquote(Reg~POD[H]~x~Css~(ÂµM)))+
poddaf.tox.df <- poddaf.tox.df[-(which(is.na(poddaf.tox.df$poddaf.Css*poddaf.tox.df$ToxCast))),]
poddaf.tox.s <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Spearman"))), digits=2)
poddaf.tox.p <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Pearson"))), digits=2)
poddaf.tox.mae <- mean(abs(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast)))
poddaf.tox.me <- mean(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast))

##BBMD
bbmd.tox.df <- full_join(BMD, tox.5perc, by="CAS") # x=ToxCast, y=BMD50.Css, xlab("ToxCast POD (ÂµM)")+ylab(bquote(BMA~BMD[A]~x~Css~(ÂµM)))+
bbmd.tox.df <- bbmd.tox.df[-(which(is.na(bbmd.tox.df$BMD50.Css*bbmd.tox.df$ToxCast))),]
wt.bmd.tox <- 1/(log10(bbmd.tox.df$BMD95.Css)-log10(bbmd.tox.df$BMD5.Css))^2
bbmd.tox.w.s <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Spearman"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.w.p <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Pearson"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.mae <- weighted.mean(abs(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast)), wt.bmd.tox)
bbmd.tox.me <- weighted.mean(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast), wt.bmd.tox)

##HED
hed.tox.df <- full_join(HED, tox.5perc, by="CAS") # x=ToxCast, y=HED50.Css, xlab("ToxCast POD (ÂµM)")+ylab(bquote(BMA~BMD[H]~x~Css~(ÂµM)))+
hed.tox.df <- hed.tox.df[-(which(is.na(hed.tox.df$HED50.Css*hed.tox.df$ToxCast))),]
wt.hed.tox <- 1/(log10(hed.tox.df$HED95.Css)-log10(hed.tox.df$HED5.Css))^2
hed.tox.w.s <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Spearman"), weights=wt.hed.tox)), digits=2)
hed.tox.w.p <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Pearson"), weights=wt.hed.tox)), digits=2)
hed.tox.mae <- weighted.mean(abs(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast)), wt.hed.tox)
hed.tox.me <- weighted.mean(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast), wt.hed.tox)


##in vivo POD vs iPSC+LCL POD scatter plot
##traditional POD
pod.cell.df <- full_join(origpod, cell.min, by="CAS") # x=iPSCLCL, y=PODCss, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(Reg~POD[A]~x~Css~(ÂµM)))+
pod.cell.df <- pod.cell.df[-(which(is.na(pod.cell.df$PODCss*pod.cell.df$iPSCLCL))),]
pod.cell.s <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.cell.p <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.cell.mae <- mean(abs(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL)))
pod.cell.me <- mean(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL))

##traditional POD x DAF
poddaf.cell.df <- full_join(poddaf, cell.min, by="CAS") # x=iPSCLCL, y=poddaf.Css, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(Reg~POD[H]~x~Css~(ÂµM)))+
poddaf.cell.df <- poddaf.cell.df[-(which(is.na(poddaf.cell.df$poddaf.Css*poddaf.cell.df$iPSCLCL))),]
poddaf.cell.s <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.cell.p <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.cell.mae <- mean(abs(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL)))
poddaf.cell.me <- mean(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL))

##BBMD
bbmd.cell.df <- full_join(BMD, cell.min, by="CAS") # x=iPSCLCL, y=BMD50.Css, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(BMA~BMD[A]~x~Css~(ÂµM)))+
bbmd.cell.df <- bbmd.cell.df[-(which(is.na(bbmd.cell.df$BMD50.Css*bbmd.cell.df$iPSCLCL))),]
wt.bmd.cell <- 1/(log10(bbmd.cell.df$BMD95.Css)-log10(bbmd.cell.df$BMD5.Css))^2
bbmd.cell.w.s <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.w.p <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.mae <- weighted.mean(abs(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL)),wt.bmd.cell)
bbmd.cell.me <- weighted.mean(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL),wt.bmd.cell)

##HED
hed.cell.df <- full_join(HED, cell.min, by="CAS") # x=iPSCLCL, y=HED50.Css, xlab("six-cell-type POD (ÂµM)")+ylab(bquote(BMA~BMD[H]~x~Css~(ÂµM)))+
hed.cell.df <- hed.cell.df[-(which(is.na(hed.cell.df$HED50.Css*hed.cell.df$iPSCLCL))),]
wt.hed.cell <- 1/(log10(hed.cell.df$HED95.Css)-log10(hed.cell.df$HED5.Css))^2
hed.cell.w.s <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.cell)), digits=2)
hed.cell.w.p <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.cell)), digits=2)
hed.cell.mae <- weighted.mean(abs(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL)),wt.hed.cell)
hed.cell.me <- weighted.mean(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL),wt.hed.cell)

##in vivo POD vs Cardio POD scatter plot
##traditional POD
pod.cardio.df <- full_join(origpod, cardio, by="CAS") # x=Min, y=PODCss, xlab("hiPSC-CM POD (ÂµM)")+ylab(bquote(Reg~POD[A]~x~Css~(ÂµM)))
pod.cardio.df <- pod.cardio.df[-(which(is.na(pod.cardio.df$Min * pod.cardio.df$PODCss))),]
pod.car.s <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Spearman")), digits=2)
pod.car.p <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Pearson")), digits=2)
pod.car.mae <- mean(abs(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min)))
pod.car.me <- mean(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min))

##traditional POD x DAF
poddaf.cardio.df <- full_join(poddaf, cardio, by="CAS") # x=Min, y=poddaf.Css, xlab("hiPSC-CM POD (ÂµM)")+ylab(bquote(Reg~POD[H]~x~Css~(ÂµM)))
poddaf.cardio.df <- poddaf.cardio.df[-(which(is.na(poddaf.cardio.df$Min * poddaf.cardio.df$poddaf.Css))),]
poddaf.car.s <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Spearman")), digits=2)
poddaf.car.p <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Pearson")), digits=2)
poddaf.car.mae <- mean(abs(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min)))
poddaf.car.me <- mean(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min))

##BBMD
bbmd.cardio.df <- full_join(BMD, cardio, by="CAS") 
bbmd.cardio.df <- bbmd.cardio.df[-(which(is.na(bbmd.cardio.df$Min * bbmd.cardio.df$BMD50.Css))),]
wt.bmd.cardio <- 1/(log10(bbmd.cardio.df$BMD95.Css)-log10(bbmd.cardio.df$BMD5.Css))^2
bbmd.car.w.s <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Spearman"), weights=wt.bmd.cardio), digits=2)
bbmd.car.w.p <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Pearson"), weights=wt.bmd.cardio), digits=2)
bbmd.car.mae <- weighted.mean(abs(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min)),wt.bmd.cardio)
bbmd.car.me <- weighted.mean(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min),wt.bmd.cardio)

##HED
hed.cardio.df <- full_join(HED, cardio, by="CAS")
hed.cardio.df <- hed.cardio.df[-(which(is.na(hed.cardio.df$Min * hed.cardio.df$HED50.Css))),]
wt.hed.cardio <- 1/(log10(hed.cardio.df$HED95.Css)-log10(hed.cardio.df$HED5.Css))^2
hed.car.w.s <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Spearman"), weights=wt.hed.cardio), digits=2)
hed.car.w.p <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Pearson"), weights=wt.hed.cardio), digits=2)
hed.car.mae <- weighted.mean(abs(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min)),wt.hed.cardio)
hed.car.me <- weighted.mean(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min),wt.hed.cardio)


##in vivo POD vs constant POD scatter plot
##traditional POD
pod.constant.df <- full_join(origpod, constant.POD, by="CAS")
pod.constant.df <- pod.constant.df[-(which(is.na(pod.constant.df$PODCss*pod.constant.df$iPSCLCL))),]
pod.constant.s <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae <- mean(abs(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL)))
pod.constant.me <- mean(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df <- full_join(poddaf, constant.POD, by="CAS")
poddaf.constant.df <- poddaf.constant.df[-(which(is.na(poddaf.constant.df$poddaf.Css*poddaf.constant.df$iPSCLCL))),]
poddaf.constant.s <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae <- mean(abs(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL)))
poddaf.constant.me <- mean(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL))

##BBMD
bbmd.constant.df <- full_join(BMD, constant.POD, by="CAS") 
bbmd.constant.df <- bbmd.constant.df[-(which(is.na(bbmd.constant.df$BMD50.Css*bbmd.constant.df$iPSCLCL))),]
wt.bmd.constant <- 1/(log10(bbmd.constant.df$BMD95.Css)-log10(bbmd.constant.df$BMD5.Css))^2
bbmd.constant.w.s <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.w.p <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.mae <- weighted.mean(abs(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL)),wt.bmd.constant)
bbmd.constant.me <- weighted.mean(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL),wt.bmd.constant)

##HED
hed.constant.df <- full_join(HED, constant.POD, by="CAS") 
hed.constant.df <- hed.constant.df[-(which(is.na(hed.constant.df$HED50.Css*hed.constant.df$iPSCLCL))),]
wt.hed.constant <- 1/(log10(hed.constant.df$HED95.Css)-log10(hed.constant.df$HED5.Css))^2
hed.constant.w.s <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant)), digits=2)
hed.constant.w.p <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant)), digits=2)
hed.constant.mae <- weighted.mean(abs(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL)),wt.hed.constant)
hed.constant.me <- weighted.mean(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL),wt.hed.constant)


##in vivo POD vs constant POD scatter plot (conc = 0.01)
##traditional POD
pod.constant.df.2 <- full_join(origpod, constant.POD.2, by="CAS") 
pod.constant.df.2 <- pod.constant.df.2[-(which(is.na(pod.constant.df.2$PODCss*pod.constant.df.2$iPSCLCL))),]
pod.constant.s.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae.2 <- mean(abs(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL)))
pod.constant.me.2 <- mean(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df.2 <- full_join(poddaf, constant.POD.2, by="CAS") 
poddaf.constant.df.2 <- poddaf.constant.df.2[-(which(is.na(poddaf.constant.df.2$poddaf.Css*poddaf.constant.df.2$iPSCLCL))),]
poddaf.constant.s.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae.2 <- mean(abs(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL)))
poddaf.constant.me.2 <- mean(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL))

##BBMD
bbmd.constant.df.2 <- full_join(BMD, constant.POD.2, by="CAS") 
bbmd.constant.df.2 <- bbmd.constant.df.2[-(which(is.na(bbmd.constant.df.2$BMD50.Css*bbmd.constant.df.2$iPSCLCL))),]
wt.bmd.constant.2 <- 1/(log10(bbmd.constant.df.2$BMD95.Css)-log10(bbmd.constant.df.2$BMD5.Css))^2
bbmd.constant.w.s.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.w.p.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.mae.2 <- weighted.mean(abs(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL)),wt.bmd.constant.2)
bbmd.constant.me.2 <- weighted.mean(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL),wt.bmd.constant.2)

##HED
hed.constant.df.2 <- full_join(HED, constant.POD.2, by="CAS") 
hed.constant.df.2 <- hed.constant.df.2[-(which(is.na(hed.constant.df.2$HED50.Css*hed.constant.df.2$iPSCLCL))),]
wt.hed.constant.2 <- 1/(log10(hed.constant.df.2$HED95.Css)-log10(hed.constant.df.2$HED5.Css))^2
hed.constant.w.s.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant.2)), digits=2)
hed.constant.w.p.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant.2)), digits=2)
hed.constant.mae.2 <- weighted.mean(abs(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL)),wt.hed.constant.2)
hed.constant.me.2 <- weighted.mean(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL),wt.hed.constant.2)

POD.data.all_Fischer <- data.frame(Group = "Using Fischer model",
                                   CAS = c(pod.tox.df$CAS, poddaf.tox.df$CAS, bbmd.tox.df$CAS, hed.tox.df$CAS,
                                           pod.cell.df$CAS, poddaf.cell.df$CAS, bbmd.cell.df$CAS, hed.cell.df$CAS,
                                           pod.cardio.df$CAS, poddaf.cardio.df$CAS, bbmd.cardio.df$CAS, hed.cardio.df$CAS,
                                           pod.constant.df$CAS, poddaf.constant.df$CAS, bbmd.constant.df$CAS, hed.constant.df$CAS,
                                           pod.constant.df.2$CAS, poddaf.constant.df.2$CAS, bbmd.constant.df.2$CAS, hed.constant.df.2$CAS),
                                   x.val = c(pod.tox.df$ToxCast, poddaf.tox.df$ToxCast, bbmd.tox.df$ToxCast, hed.tox.df$ToxCast,
                                             pod.cell.df$iPSCLCL, poddaf.cell.df$iPSCLCL, bbmd.cell.df$iPSCLCL, hed.cell.df$iPSCLCL,
                                             pod.cardio.df$Min, poddaf.cardio.df$Min, bbmd.cardio.df$Min, hed.cardio.df$Min,
                                             pod.constant.df$iPSCLCL, poddaf.constant.df$iPSCLCL, bbmd.constant.df$iPSCLCL, hed.constant.df$iPSCLCL,
                                             pod.constant.df.2$iPSCLCL, poddaf.constant.df.2$iPSCLCL, bbmd.constant.df.2$iPSCLCL, hed.constant.df.2$iPSCLCL),
                                   x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 48),
                                   y.val = c(pod.tox.df$PODCss, poddaf.tox.df$poddaf.Css, bbmd.tox.df$BMD50.Css, hed.tox.df$HED50.Css,
                                             pod.cell.df$PODCss, poddaf.cell.df$poddaf.Css, bbmd.cell.df$BMD50.Css, hed.cell.df$HED50.Css,
                                             pod.cardio.df$PODCss, poddaf.cardio.df$poddaf.Css, bbmd.cardio.df$BMD50.Css, hed.cardio.df$HED50.Css,
                                             pod.constant.df$PODCss, poddaf.constant.df$poddaf.Css, bbmd.constant.df$BMD50.Css, hed.constant.df$HED50.Css,
                                             pod.constant.df.2$PODCss, poddaf.constant.df.2$poddaf.Css, bbmd.constant.df.2$BMD50.Css, hed.constant.df.2$HED50.Css),
                                   y.val.min = c(rep(NA, 24), bbmd.tox.df$BMD5.Css, hed.tox.df$HED5.Css,
                                                 rep(NA, 24), bbmd.cell.df$BMD5.Css, hed.cell.df$HED5.Css,
                                                 rep(NA, 24), bbmd.cardio.df$BMD5.Css, hed.cardio.df$HED5.Css,
                                                 rep(NA, 24), bbmd.constant.df$BMD5.Css, hed.constant.df$HED5.Css,
                                                 rep(NA, 24), bbmd.constant.df.2$BMD5.Css, hed.constant.df.2$HED5.Css),
                                   y.val.max = c(rep(NA, 24), bbmd.tox.df$BMD95.Css, hed.tox.df$HED95.Css,
                                                 rep(NA, 24), bbmd.cell.df$BMD95.Css, hed.cell.df$HED95.Css,
                                                 rep(NA, 24), bbmd.cardio.df$BMD95.Css, hed.cardio.df$HED95.Css,
                                                 rep(NA, 24), bbmd.constant.df$BMD95.Css, hed.constant.df$HED95.Css,
                                                 rep(NA, 24), bbmd.constant.df.2$BMD95.Css, hed.constant.df.2$HED95.Css),
                                   y.type = rep(rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"), each =12), 5))
RHO.data.all_Fischer <- data.frame(Group = "Using Fischer model",
                                   x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                   y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                   rho.value = c(pod.tox.s, poddaf.tox.s, bbmd.tox.w.s, hed.tox.w.s,
                                                 pod.cell.s, poddaf.cell.s, bbmd.cell.w.s, hed.cell.w.s,
                                                 pod.car.s, poddaf.car.s, bbmd.car.w.s, hed.car.w.s,
                                                 pod.constant.s, poddaf.constant.s, bbmd.constant.w.s, hed.constant.w.s,
                                                 pod.constant.s.2, poddaf.constant.s.2, bbmd.constant.w.s.2, hed.constant.w.s.2))
R.data.all_Fischer <- data.frame(Group = "Using Fischer model",
                                 x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                 y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                 rho.value = c(pod.tox.p, poddaf.tox.p, bbmd.tox.w.p, hed.tox.w.p,
                                               pod.cell.p, poddaf.cell.p, bbmd.cell.w.p, hed.cell.w.p,
                                               pod.car.p, poddaf.car.p, bbmd.car.w.p, hed.car.w.p,
                                               pod.constant.p, poddaf.constant.p, bbmd.constant.w.p, hed.constant.w.p,
                                               pod.constant.p.2, poddaf.constant.p.2, bbmd.constant.w.p.2, hed.constant.w.p.2))

mae.data.all_Fischer <- data.frame(Group = "Using Fischer model",
                                   x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                   y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                   rho.value = c(pod.tox.mae, poddaf.tox.mae, bbmd.tox.mae, hed.tox.mae,
                                                 pod.cell.mae, poddaf.cell.mae, bbmd.cell.mae, hed.cell.mae,
                                                 pod.car.mae, poddaf.car.mae, bbmd.car.mae, hed.car.mae,
                                                 pod.constant.mae, poddaf.constant.mae, bbmd.constant.mae, hed.constant.mae,
                                                 pod.constant.mae.2, poddaf.constant.mae.2, bbmd.constant.mae.2, hed.constant.mae.2))
me.data.all_Fischer <- data.frame(Group = "Using Fischer model",
                                  x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                  y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                  rho.value = c(pod.tox.me, poddaf.tox.me, bbmd.tox.me, hed.tox.me,
                                                pod.cell.me, poddaf.cell.me, bbmd.cell.me, hed.cell.me,
                                                pod.car.me, poddaf.car.me, bbmd.car.me, hed.car.me,
                                                pod.constant.me, poddaf.constant.me, bbmd.constant.me, hed.constant.me,
                                                pod.constant.me.2, poddaf.constant.me.2, bbmd.constant.me.2, hed.constant.me.2))



# -------------- with adjustment: fup and fub from IVMBM --------------
invitro_fub <- read.csv("invitro_fub_cfree_nonminal.csv")
#ToxCast active
toxcast <- read.csv(file.path("POD database_15che","toxcast ac50.csv"))
tox.5perc <- data.frame(sapply(toxcast, quantile, prob=0.05, na.rm = T))
tox.5perc$CAS <- origpod$CAS
colnames(tox.5perc)[1] <- "ToxCast"

#iPSC+LCL
cell.t <- read.csv(file.path("POD database_15che","cellline POD distribution for 15 chemicals.csv"))
cell.min <- data.frame(sapply(cell.t[2:16], min, na.rm=T))
cell.min$CAS <- origpod$CAS
colnames(cell.min)[1] <- "iPSCLCL"

constant.POD <- data.frame(CAS = invitro_fub$CAS[1:15],
                           iPSCLCL = 1)

constant.POD.2 <- data.frame(CAS = invitro_fub$CAS[1:15],
                             iPSCLCL = 0.01)

#Cardio
#cardiomyocytes (single donor)
cardio.df <- read.csv(file.path("POD database_15che","Cardio POD for 15 chemicals.csv"))
cardio <- full_join(cardio.df, origpod[,c(1,2)], by="CAS")
cardio <- cardio[-which(is.na(cardio$Chemical.y)),]


for(i in 1:15){
  
  tox.5perc$ToxCast[i] <- (tox.5perc$ToxCast[i]*invitro_fub$IVBMB[which(tox.5perc$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(tox.5perc$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(tox.5perc$CAS[i] == invivo_fup$CAS)])
  
  cell.min$iPSCLCL[i] <- (cell.min$iPSCLCL[i]*invitro_fub$IVBMB[which(cell.min$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cell.min$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cell.min$CAS[i] == invivo_fup$CAS)])
  
  constant.POD$iPSCLCL[i] <- (constant.POD$iPSCLCL[i]*invitro_fub$IVBMB[which(constant.POD$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD$CAS[i] == invivo_fup$CAS)])
  
  constant.POD.2$iPSCLCL[i] <- (constant.POD.2$iPSCLCL[i]*invitro_fub$IVBMB[which(constant.POD.2$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD.2$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD.2$CAS[i] == invivo_fup$CAS)])
  
  cardio$Min[i] <- (cardio$Min[i]*invitro_fub$IVBMB[which(cardio$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cardio$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cardio$CAS[i] == invivo_fup$CAS)])
  
}

##traditional POD
pod.tox.df <- full_join(origpod, tox.5perc, by="CAS") 
pod.tox.s <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Spearman"))), digits=2)
pod.tox.p <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Pearson"))), digits=2)
pod.tox.mae <- mean(abs(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast)))
pod.tox.me <- mean(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast))

##traditional POD x DAF
poddaf.tox.df <- full_join(poddaf, tox.5perc, by="CAS") 
poddaf.tox.s <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Spearman"))), digits=2)
poddaf.tox.p <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Pearson"))), digits=2)
poddaf.tox.mae <- mean(abs(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast)))
poddaf.tox.me <- mean(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast))

##BBMD
bbmd.tox.df <- full_join(BMD, tox.5perc, by="CAS") 
wt.bmd.tox <- 1/(log10(bbmd.tox.df$BMD95.Css)-log10(bbmd.tox.df$BMD5.Css))^2
bbmd.tox.w.s <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Spearman"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.w.p <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Pearson"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.mae <- weighted.mean(abs(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast)), wt.bmd.tox)
bbmd.tox.me <- weighted.mean(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast), wt.bmd.tox)

##HED
hed.tox.df <- full_join(HED, tox.5perc, by="CAS") 
wt.hed.tox <- 1/(log10(hed.tox.df$HED95.Css)-log10(hed.tox.df$HED5.Css))^2
hed.tox.w.s <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Spearman"), weights=wt.hed.tox)), digits=2)
hed.tox.w.p <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Pearson"), weights=wt.hed.tox)), digits=2)
hed.tox.mae <- weighted.mean(abs(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast)), wt.hed.tox)
hed.tox.me <- weighted.mean(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast), wt.hed.tox)


##in vivo POD vs iPSC+LCL POD scatter plot
##traditional POD
pod.cell.df <- full_join(origpod, cell.min, by="CAS") 
pod.cell.s <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.cell.p <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.cell.mae <- mean(abs(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL)))
pod.cell.me <- mean(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL))

##traditional POD x DAF
poddaf.cell.df <- full_join(poddaf, cell.min, by="CAS") 
poddaf.cell.s <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.cell.p <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.cell.mae <- mean(abs(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL)))
poddaf.cell.me <- mean(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL))

##BBMD
bbmd.cell.df <- full_join(BMD, cell.min, by="CAS") 
wt.bmd.cell <- 1/(log10(bbmd.cell.df$BMD95.Css)-log10(bbmd.cell.df$BMD5.Css))^2
bbmd.cell.w.s <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.w.p <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.mae <- weighted.mean(abs(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL)),wt.bmd.cell)
bbmd.cell.me <- weighted.mean(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL),wt.bmd.cell)

##HED
hed.cell.df <- full_join(HED, cell.min, by="CAS") 
wt.hed.cell <- 1/(log10(hed.cell.df$HED95.Css)-log10(hed.cell.df$HED5.Css))^2
hed.cell.w.s <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.cell)), digits=2)
hed.cell.w.p <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.cell)), digits=2)
hed.cell.mae <- weighted.mean(abs(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL)),wt.hed.cell)
hed.cell.me <- weighted.mean(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL),wt.hed.cell)

##in vivo POD vs Cardio POD scatter plot
##traditional POD
pod.cardio.df <- full_join(origpod, cardio, by="CAS") 
pod.car.s <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Spearman")), digits=2)
pod.car.p <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Pearson")), digits=2)
pod.car.mae <- mean(abs(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min)))
pod.car.me <- mean(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min))

##traditional POD x DAF
poddaf.cardio.df <- full_join(poddaf, cardio, by="CAS") 
poddaf.car.s <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Spearman")), digits=2)
poddaf.car.p <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Pearson")), digits=2)
poddaf.car.mae <- mean(abs(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min)))
poddaf.car.me <- mean(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min))

##BBMD
bbmd.cardio.df <- full_join(BMD, cardio, by="CAS") 
wt.bmd.cardio <- 1/(log10(bbmd.cardio.df$BMD95.Css)-log10(bbmd.cardio.df$BMD5.Css))^2
bbmd.car.w.s <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Spearman"), weights=wt.bmd.cardio), digits=2)
bbmd.car.w.p <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Pearson"), weights=wt.bmd.cardio), digits=2)
bbmd.car.mae <- weighted.mean(abs(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min)),wt.bmd.cardio)
bbmd.car.me <- weighted.mean(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min),wt.bmd.cardio)

##HED
hed.cardio.df <- full_join(HED, cardio, by="CAS")
wt.hed.cardio <- 1/(log10(hed.cardio.df$HED95.Css)-log10(hed.cardio.df$HED5.Css))^2
hed.car.w.s <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Spearman"), weights=wt.hed.cardio), digits=2)
hed.car.w.p <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Pearson"), weights=wt.hed.cardio), digits=2)
hed.car.mae <- weighted.mean(abs(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min)),wt.hed.cardio)
hed.car.me <- weighted.mean(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min),wt.hed.cardio)


##in vivo POD vs constant POD scatter plot
##traditional POD
pod.constant.df <- full_join(origpod, constant.POD, by="CAS") 
pod.constant.s <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae <- mean(abs(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL)))
pod.constant.me <- mean(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df <- full_join(poddaf, constant.POD, by="CAS") 
poddaf.constant.s <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae <- mean(abs(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL)))
poddaf.constant.me <- mean(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL))

##BBMD
bbmd.constant.df <- full_join(BMD, constant.POD, by="CAS") 
wt.bmd.constant <- 1/(log10(bbmd.constant.df$BMD95.Css)-log10(bbmd.constant.df$BMD5.Css))^2
bbmd.constant.w.s <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.w.p <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.mae <- weighted.mean(abs(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL)),wt.bmd.constant)
bbmd.constant.me <- weighted.mean(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL),wt.bmd.constant)

##HED
hed.constant.df <- full_join(HED, constant.POD, by="CAS") 
wt.hed.constant <- 1/(log10(hed.constant.df$HED95.Css)-log10(hed.constant.df$HED5.Css))^2
hed.constant.w.s <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant)), digits=2)
hed.constant.w.p <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant)), digits=2)
hed.constant.mae <- weighted.mean(abs(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL)),wt.hed.constant)
hed.constant.me <- weighted.mean(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL),wt.hed.constant)


##in vivo POD vs constant POD scatter plot (conc = 0.01)
##traditional POD
pod.constant.df.2 <- full_join(origpod, constant.POD.2, by="CAS") 
pod.constant.s.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae.2 <- mean(abs(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL)))
pod.constant.me.2 <- mean(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df.2 <- full_join(poddaf, constant.POD.2, by="CAS") 
poddaf.constant.s.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae.2 <- mean(abs(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL)))
poddaf.constant.me.2 <- mean(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL))

##BBMD
bbmd.constant.df.2 <- full_join(BMD, constant.POD.2, by="CAS") 
wt.bmd.constant.2 <- 1/(log10(bbmd.constant.df.2$BMD95.Css)-log10(bbmd.constant.df.2$BMD5.Css))^2
bbmd.constant.w.s.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.w.p.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.mae.2 <- weighted.mean(abs(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL)),wt.bmd.constant.2)
bbmd.constant.me.2 <- weighted.mean(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL),wt.bmd.constant.2)

##HED
hed.constant.df.2 <- full_join(HED, constant.POD.2, by="CAS") 
wt.hed.constant.2 <- 1/(log10(hed.constant.df.2$HED95.Css)-log10(hed.constant.df.2$HED5.Css))^2
hed.constant.w.s.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant.2)), digits=2)
hed.constant.w.p.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant.2)), digits=2)
hed.constant.mae.2 <- weighted.mean(abs(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL)),wt.hed.constant.2)
hed.constant.me.2 <- weighted.mean(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL),wt.hed.constant.2)

POD.data.all_IVMBM <- data.frame(Group = "Using IVMBM model",
                                 CAS = c(pod.tox.df$CAS, poddaf.tox.df$CAS, bbmd.tox.df$CAS, hed.tox.df$CAS,
                                         pod.cell.df$CAS, poddaf.cell.df$CAS, bbmd.cell.df$CAS, hed.cell.df$CAS,
                                         pod.cardio.df$CAS, poddaf.cardio.df$CAS, bbmd.cardio.df$CAS, hed.cardio.df$CAS,
                                         pod.constant.df$CAS, poddaf.constant.df$CAS, bbmd.constant.df$CAS, hed.constant.df$CAS,
                                         pod.constant.df.2$CAS, poddaf.constant.df.2$CAS, bbmd.constant.df.2$CAS, hed.constant.df.2$CAS),
                                 x.val = c(pod.tox.df$ToxCast, poddaf.tox.df$ToxCast, bbmd.tox.df$ToxCast, hed.tox.df$ToxCast,
                                           pod.cell.df$iPSCLCL, poddaf.cell.df$iPSCLCL, bbmd.cell.df$iPSCLCL, hed.cell.df$iPSCLCL,
                                           pod.cardio.df$Min, poddaf.cardio.df$Min, bbmd.cardio.df$Min, hed.cardio.df$Min,
                                           pod.constant.df$iPSCLCL, poddaf.constant.df$iPSCLCL, bbmd.constant.df$iPSCLCL, hed.constant.df$iPSCLCL,
                                           pod.constant.df.2$iPSCLCL, poddaf.constant.df.2$iPSCLCL, bbmd.constant.df.2$iPSCLCL, hed.constant.df.2$iPSCLCL),
                                 x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 60),
                                 y.val = c(pod.tox.df$PODCss, poddaf.tox.df$poddaf.Css, bbmd.tox.df$BMD50.Css, hed.tox.df$HED50.Css,
                                           pod.cell.df$PODCss, poddaf.cell.df$poddaf.Css, bbmd.cell.df$BMD50.Css, hed.cell.df$HED50.Css,
                                           pod.cardio.df$PODCss, poddaf.cardio.df$poddaf.Css, bbmd.cardio.df$BMD50.Css, hed.cardio.df$HED50.Css,
                                           pod.constant.df$PODCss, poddaf.constant.df$poddaf.Css, bbmd.constant.df$BMD50.Css, hed.constant.df$HED50.Css,
                                           pod.constant.df.2$PODCss, poddaf.constant.df.2$poddaf.Css, bbmd.constant.df.2$BMD50.Css, hed.constant.df.2$HED50.Css),
                                 y.val.min = c(rep(NA, 30), bbmd.tox.df$BMD5.Css, hed.tox.df$HED5.Css,
                                               rep(NA, 30), bbmd.cell.df$BMD5.Css, hed.cell.df$HED5.Css,
                                               rep(NA, 30), bbmd.cardio.df$BMD5.Css, hed.cardio.df$HED5.Css,
                                               rep(NA, 30), bbmd.constant.df$BMD5.Css, hed.constant.df$HED5.Css,
                                               rep(NA, 30), bbmd.constant.df.2$BMD5.Css, hed.constant.df.2$HED5.Css),
                                 y.val.max = c(rep(NA, 30), bbmd.tox.df$BMD95.Css, hed.tox.df$HED95.Css,
                                               rep(NA, 30), bbmd.cell.df$BMD95.Css, hed.cell.df$HED95.Css,
                                               rep(NA, 30), bbmd.cardio.df$BMD95.Css, hed.cardio.df$HED95.Css,
                                               rep(NA, 30), bbmd.constant.df$BMD95.Css, hed.constant.df$HED95.Css,
                                               rep(NA, 30), bbmd.constant.df.2$BMD95.Css, hed.constant.df.2$HED95.Css),
                                 y.type = rep(rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"), each =15), 5))
RHO.data.all_IVMBM <- data.frame(Group = "Using IVMBM model",
                                 x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                 y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                 rho.value = c(pod.tox.s, poddaf.tox.s, bbmd.tox.w.s, hed.tox.w.s,
                                               pod.cell.s, poddaf.cell.s, bbmd.cell.w.s, hed.cell.w.s,
                                               pod.car.s, poddaf.car.s, bbmd.car.w.s, hed.car.w.s,
                                               pod.constant.s, poddaf.constant.s, bbmd.constant.w.s, hed.constant.w.s,
                                               pod.constant.s.2, poddaf.constant.s.2, bbmd.constant.w.s.2, hed.constant.w.s.2))
R.data.all_IVMBM <- data.frame(Group = "Using IVMBM model",
                               x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                               y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                               rho.value = c(pod.tox.p, poddaf.tox.p, bbmd.tox.w.p, hed.tox.w.p,
                                             pod.cell.p, poddaf.cell.p, bbmd.cell.w.p, hed.cell.w.p,
                                             pod.car.p, poddaf.car.p, bbmd.car.w.p, hed.car.w.p,
                                             pod.constant.p, poddaf.constant.p, bbmd.constant.w.p, hed.constant.w.p,
                                             pod.constant.p.2, poddaf.constant.p.2, bbmd.constant.w.p.2, hed.constant.w.p.2))

mae.data.all_IVMBM <- data.frame(Group = "Using IVMBM model",
                                 x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                 y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                 rho.value = c(pod.tox.mae, poddaf.tox.mae, bbmd.tox.mae, hed.tox.mae,
                                               pod.cell.mae, poddaf.cell.mae, bbmd.cell.mae, hed.cell.mae,
                                               pod.car.mae, poddaf.car.mae, bbmd.car.mae, hed.car.mae,
                                               pod.constant.mae, poddaf.constant.mae, bbmd.constant.mae, hed.constant.mae,
                                               pod.constant.mae.2, poddaf.constant.mae.2, bbmd.constant.mae.2, hed.constant.mae.2))
me.data.all_IVMBM <- data.frame(Group = "Using IVMBM model",
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                rho.value = c(pod.tox.me, poddaf.tox.me, bbmd.tox.me, hed.tox.me,
                                              pod.cell.me, poddaf.cell.me, bbmd.cell.me, hed.cell.me,
                                              pod.car.me, poddaf.car.me, bbmd.car.me, hed.car.me,
                                              pod.constant.me, poddaf.constant.me, bbmd.constant.me, hed.constant.me,
                                              pod.constant.me.2, poddaf.constant.me.2, bbmd.constant.me.2, hed.constant.me.2))





# -------------- with adjustment: fup and fub from VIVD --------------
invitro_fub <- read.csv("invitro_fub_cfree_nonminal.csv")
#ToxCast active
toxcast <- read.csv(file.path("POD database_15che","toxcast ac50.csv"))
tox.5perc <- data.frame(sapply(toxcast, quantile, prob=0.05, na.rm = T))
tox.5perc$CAS <- origpod$CAS
colnames(tox.5perc)[1] <- "ToxCast"

#iPSC+LCL
cell.t <- read.csv(file.path("POD database_15che","cellline POD distribution for 15 chemicals.csv"))
cell.min <- data.frame(sapply(cell.t[2:16], min, na.rm=T))
cell.min$CAS <- origpod$CAS
colnames(cell.min)[1] <- "iPSCLCL"

constant.POD <- data.frame(CAS = invitro_fub$CAS[1:15],
                           iPSCLCL = 1)

constant.POD.2 <- data.frame(CAS = invitro_fub$CAS[1:15],
                             iPSCLCL = 0.01)

#Cardio
#cardiomyocytes (single donor)
cardio.df <- read.csv(file.path("POD database_15che","Cardio POD for 15 chemicals.csv"))
cardio <- full_join(cardio.df, origpod[,c(1,2)], by="CAS")
cardio <- cardio[-which(is.na(cardio$Chemical.y)),]


for(i in 1:15){
  
  tox.5perc$ToxCast[i] <- (tox.5perc$ToxCast[i]*invitro_fub$VIVD[which(tox.5perc$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(tox.5perc$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(tox.5perc$CAS[i] == invivo_fup$CAS)])
  
  cell.min$iPSCLCL[i] <- (cell.min$iPSCLCL[i]*invitro_fub$VIVD[which(cell.min$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cell.min$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cell.min$CAS[i] == invivo_fup$CAS)])
  
  constant.POD$iPSCLCL[i] <- (constant.POD$iPSCLCL[i]*invitro_fub$VIVD[which(constant.POD$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD$CAS[i] == invivo_fup$CAS)])
  
  constant.POD.2$iPSCLCL[i] <- (constant.POD.2$iPSCLCL[i]*invitro_fub$VIVD[which(constant.POD.2$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD.2$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD.2$CAS[i] == invivo_fup$CAS)])
  
  cardio$Min[i] <- (cardio$Min[i]*invitro_fub$VIVD[which(cardio$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cardio$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cardio$CAS[i] == invivo_fup$CAS)])
  
}

##traditional POD
pod.tox.df <- full_join(origpod, tox.5perc, by="CAS") 
pod.tox.s <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Spearman"))), digits=2)
pod.tox.p <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Pearson"))), digits=2)
pod.tox.mae <- mean(abs(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast)))
pod.tox.me <- mean(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast))

##traditional POD x DAF
poddaf.tox.df <- full_join(poddaf, tox.5perc, by="CAS") 
poddaf.tox.s <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Spearman"))), digits=2)
poddaf.tox.p <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Pearson"))), digits=2)
poddaf.tox.mae <- mean(abs(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast)))
poddaf.tox.me <- mean(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast))

##BBMD
bbmd.tox.df <- full_join(BMD, tox.5perc, by="CAS") 
wt.bmd.tox <- 1/(log10(bbmd.tox.df$BMD95.Css)-log10(bbmd.tox.df$BMD5.Css))^2
bbmd.tox.w.s <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Spearman"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.w.p <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Pearson"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.mae <- weighted.mean(abs(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast)), wt.bmd.tox)
bbmd.tox.me <- weighted.mean(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast), wt.bmd.tox)

##HED
hed.tox.df <- full_join(HED, tox.5perc, by="CAS") 
wt.hed.tox <- 1/(log10(hed.tox.df$HED95.Css)-log10(hed.tox.df$HED5.Css))^2
hed.tox.w.s <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Spearman"), weights=wt.hed.tox)), digits=2)
hed.tox.w.p <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Pearson"), weights=wt.hed.tox)), digits=2)
hed.tox.mae <- weighted.mean(abs(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast)), wt.hed.tox)
hed.tox.me <- weighted.mean(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast), wt.hed.tox)


##in vivo POD vs iPSC+LCL POD scatter plot
##traditional POD
pod.cell.df <- full_join(origpod, cell.min, by="CAS") 
pod.cell.s <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.cell.p <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.cell.mae <- mean(abs(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL)))
pod.cell.me <- mean(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL))

##traditional POD x DAF
poddaf.cell.df <- full_join(poddaf, cell.min, by="CAS") 
poddaf.cell.s <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.cell.p <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.cell.mae <- mean(abs(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL)))
poddaf.cell.me <- mean(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL))

##BBMD
bbmd.cell.df <- full_join(BMD, cell.min, by="CAS") 
wt.bmd.cell <- 1/(log10(bbmd.cell.df$BMD95.Css)-log10(bbmd.cell.df$BMD5.Css))^2
bbmd.cell.w.s <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.w.p <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.mae <- weighted.mean(abs(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL)),wt.bmd.cell)
bbmd.cell.me <- weighted.mean(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL),wt.bmd.cell)

##HED
hed.cell.df <- full_join(HED, cell.min, by="CAS") 
wt.hed.cell <- 1/(log10(hed.cell.df$HED95.Css)-log10(hed.cell.df$HED5.Css))^2
hed.cell.w.s <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.cell)), digits=2)
hed.cell.w.p <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.cell)), digits=2)
hed.cell.mae <- weighted.mean(abs(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL)),wt.hed.cell)
hed.cell.me <- weighted.mean(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL),wt.hed.cell)

##in vivo POD vs Cardio POD scatter plot
##traditional POD
pod.cardio.df <- full_join(origpod, cardio, by="CAS") 
pod.car.s <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Spearman")), digits=2)
pod.car.p <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Pearson")), digits=2)
pod.car.mae <- mean(abs(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min)))
pod.car.me <- mean(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min))

##traditional POD x DAF
poddaf.cardio.df <- full_join(poddaf, cardio, by="CAS") 
poddaf.car.s <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Spearman")), digits=2)
poddaf.car.p <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Pearson")), digits=2)
poddaf.car.mae <- mean(abs(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min)))
poddaf.car.me <- mean(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min))

##BBMD
bbmd.cardio.df <- full_join(BMD, cardio, by="CAS") 
wt.bmd.cardio <- 1/(log10(bbmd.cardio.df$BMD95.Css)-log10(bbmd.cardio.df$BMD5.Css))^2
bbmd.car.w.s <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Spearman"), weights=wt.bmd.cardio), digits=2)
bbmd.car.w.p <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Pearson"), weights=wt.bmd.cardio), digits=2)
bbmd.car.mae <- weighted.mean(abs(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min)),wt.bmd.cardio)
bbmd.car.me <- weighted.mean(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min),wt.bmd.cardio)

##HED
hed.cardio.df <- full_join(HED, cardio, by="CAS")
wt.hed.cardio <- 1/(log10(hed.cardio.df$HED95.Css)-log10(hed.cardio.df$HED5.Css))^2
hed.car.w.s <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Spearman"), weights=wt.hed.cardio), digits=2)
hed.car.w.p <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Pearson"), weights=wt.hed.cardio), digits=2)
hed.car.mae <- weighted.mean(abs(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min)),wt.hed.cardio)
hed.car.me <- weighted.mean(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min),wt.hed.cardio)


##in vivo POD vs constant POD scatter plot
##traditional POD
pod.constant.df <- full_join(origpod, constant.POD, by="CAS") 
pod.constant.s <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae <- mean(abs(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL)))
pod.constant.me <- mean(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df <- full_join(poddaf, constant.POD, by="CAS") 
poddaf.constant.s <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae <- mean(abs(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL)))
poddaf.constant.me <- mean(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL))

##BBMD
bbmd.constant.df <- full_join(BMD, constant.POD, by="CAS") 
wt.bmd.constant <- 1/(log10(bbmd.constant.df$BMD95.Css)-log10(bbmd.constant.df$BMD5.Css))^2
bbmd.constant.w.s <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.w.p <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.mae <- weighted.mean(abs(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL)),wt.bmd.constant)
bbmd.constant.me <- weighted.mean(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL),wt.bmd.constant)

##HED
hed.constant.df <- full_join(HED, constant.POD, by="CAS") 
wt.hed.constant <- 1/(log10(hed.constant.df$HED95.Css)-log10(hed.constant.df$HED5.Css))^2
hed.constant.w.s <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant)), digits=2)
hed.constant.w.p <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant)), digits=2)
hed.constant.mae <- weighted.mean(abs(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL)),wt.hed.constant)
hed.constant.me <- weighted.mean(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL),wt.hed.constant)


##in vivo POD vs constant POD scatter plot (conc = 0.01)
##traditional POD
pod.constant.df.2 <- full_join(origpod, constant.POD.2, by="CAS") 
pod.constant.s.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae.2 <- mean(abs(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL)))
pod.constant.me.2 <- mean(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df.2 <- full_join(poddaf, constant.POD.2, by="CAS") 
poddaf.constant.s.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae.2 <- mean(abs(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL)))
poddaf.constant.me.2 <- mean(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL))

##BBMD
bbmd.constant.df.2 <- full_join(BMD, constant.POD.2, by="CAS") 
wt.bmd.constant.2 <- 1/(log10(bbmd.constant.df.2$BMD95.Css)-log10(bbmd.constant.df.2$BMD5.Css))^2
bbmd.constant.w.s.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.w.p.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.mae.2 <- weighted.mean(abs(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL)),wt.bmd.constant.2)
bbmd.constant.me.2 <- weighted.mean(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL),wt.bmd.constant.2)

##HED
hed.constant.df.2 <- full_join(HED, constant.POD.2, by="CAS") 
wt.hed.constant.2 <- 1/(log10(hed.constant.df.2$HED95.Css)-log10(hed.constant.df.2$HED5.Css))^2
hed.constant.w.s.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant.2)), digits=2)
hed.constant.w.p.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant.2)), digits=2)
hed.constant.mae.2 <- weighted.mean(abs(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL)),wt.hed.constant.2)
hed.constant.me.2 <- weighted.mean(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL),wt.hed.constant.2)

POD.data.all_VIVD <- data.frame(Group = "Using VIVD model",
                                CAS = c(pod.tox.df$CAS, poddaf.tox.df$CAS, bbmd.tox.df$CAS, hed.tox.df$CAS,
                                        pod.cell.df$CAS, poddaf.cell.df$CAS, bbmd.cell.df$CAS, hed.cell.df$CAS,
                                        pod.cardio.df$CAS, poddaf.cardio.df$CAS, bbmd.cardio.df$CAS, hed.cardio.df$CAS,
                                        pod.constant.df$CAS, poddaf.constant.df$CAS, bbmd.constant.df$CAS, hed.constant.df$CAS,
                                        pod.constant.df.2$CAS, poddaf.constant.df.2$CAS, bbmd.constant.df.2$CAS, hed.constant.df.2$CAS),
                                x.val = c(pod.tox.df$ToxCast, poddaf.tox.df$ToxCast, bbmd.tox.df$ToxCast, hed.tox.df$ToxCast,
                                          pod.cell.df$iPSCLCL, poddaf.cell.df$iPSCLCL, bbmd.cell.df$iPSCLCL, hed.cell.df$iPSCLCL,
                                          pod.cardio.df$Min, poddaf.cardio.df$Min, bbmd.cardio.df$Min, hed.cardio.df$Min,
                                          pod.constant.df$iPSCLCL, poddaf.constant.df$iPSCLCL, bbmd.constant.df$iPSCLCL, hed.constant.df$iPSCLCL,
                                          pod.constant.df.2$iPSCLCL, poddaf.constant.df.2$iPSCLCL, bbmd.constant.df.2$iPSCLCL, hed.constant.df.2$iPSCLCL),
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 60),
                                y.val = c(pod.tox.df$PODCss, poddaf.tox.df$poddaf.Css, bbmd.tox.df$BMD50.Css, hed.tox.df$HED50.Css,
                                          pod.cell.df$PODCss, poddaf.cell.df$poddaf.Css, bbmd.cell.df$BMD50.Css, hed.cell.df$HED50.Css,
                                          pod.cardio.df$PODCss, poddaf.cardio.df$poddaf.Css, bbmd.cardio.df$BMD50.Css, hed.cardio.df$HED50.Css,
                                          pod.constant.df$PODCss, poddaf.constant.df$poddaf.Css, bbmd.constant.df$BMD50.Css, hed.constant.df$HED50.Css,
                                          pod.constant.df.2$PODCss, poddaf.constant.df.2$poddaf.Css, bbmd.constant.df.2$BMD50.Css, hed.constant.df.2$HED50.Css),
                                y.val.min = c(rep(NA, 30), bbmd.tox.df$BMD5.Css, hed.tox.df$HED5.Css,
                                              rep(NA, 30), bbmd.cell.df$BMD5.Css, hed.cell.df$HED5.Css,
                                              rep(NA, 30), bbmd.cardio.df$BMD5.Css, hed.cardio.df$HED5.Css,
                                              rep(NA, 30), bbmd.constant.df$BMD5.Css, hed.constant.df$HED5.Css,
                                              rep(NA, 30), bbmd.constant.df.2$BMD5.Css, hed.constant.df.2$HED5.Css),
                                y.val.max = c(rep(NA, 30), bbmd.tox.df$BMD95.Css, hed.tox.df$HED95.Css,
                                              rep(NA, 30), bbmd.cell.df$BMD95.Css, hed.cell.df$HED95.Css,
                                              rep(NA, 30), bbmd.cardio.df$BMD95.Css, hed.cardio.df$HED95.Css,
                                              rep(NA, 30), bbmd.constant.df$BMD95.Css, hed.constant.df$HED95.Css,
                                              rep(NA, 30), bbmd.constant.df.2$BMD95.Css, hed.constant.df.2$HED95.Css),
                                y.type = rep(rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"), each =15), 5))
RHO.data.all_VIVD <- data.frame(Group = "Using VIVD model",
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                rho.value = c(pod.tox.s, poddaf.tox.s, bbmd.tox.w.s, hed.tox.w.s,
                                              pod.cell.s, poddaf.cell.s, bbmd.cell.w.s, hed.cell.w.s,
                                              pod.car.s, poddaf.car.s, bbmd.car.w.s, hed.car.w.s,
                                              pod.constant.s, poddaf.constant.s, bbmd.constant.w.s, hed.constant.w.s,
                                              pod.constant.s.2, poddaf.constant.s.2, bbmd.constant.w.s.2, hed.constant.w.s.2))
R.data.all_VIVD <- data.frame(Group = "Using VIVD model",
                              x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                              y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                              rho.value = c(pod.tox.p, poddaf.tox.p, bbmd.tox.w.p, hed.tox.w.p,
                                            pod.cell.p, poddaf.cell.p, bbmd.cell.w.p, hed.cell.w.p,
                                            pod.car.p, poddaf.car.p, bbmd.car.w.p, hed.car.w.p,
                                            pod.constant.p, poddaf.constant.p, bbmd.constant.w.p, hed.constant.w.p,
                                            pod.constant.p.2, poddaf.constant.p.2, bbmd.constant.w.p.2, hed.constant.w.p.2))

mae.data.all_VIVD <- data.frame(Group = "Using VIVD model",
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                rho.value = c(pod.tox.mae, poddaf.tox.mae, bbmd.tox.mae, hed.tox.mae,
                                              pod.cell.mae, poddaf.cell.mae, bbmd.cell.mae, hed.cell.mae,
                                              pod.car.mae, poddaf.car.mae, bbmd.car.mae, hed.car.mae,
                                              pod.constant.mae, poddaf.constant.mae, bbmd.constant.mae, hed.constant.mae,
                                              pod.constant.mae.2, poddaf.constant.mae.2, bbmd.constant.mae.2, hed.constant.mae.2))
me.data.all_VIVD <- data.frame(Group = "Using VIVD model",
                               x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                               y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                               rho.value = c(pod.tox.me, poddaf.tox.me, bbmd.tox.me, hed.tox.me,
                                             pod.cell.me, poddaf.cell.me, bbmd.cell.me, hed.cell.me,
                                             pod.car.me, poddaf.car.me, bbmd.car.me, hed.car.me,
                                             pod.constant.me, poddaf.constant.me, bbmd.constant.me, hed.constant.me,
                                             pod.constant.me.2, poddaf.constant.me.2, bbmd.constant.me.2, hed.constant.me.2))




# -------------- with adjustment: fup and fub from VCBA --------------
invitro_fub <- read.csv("invitro_fub_cfree_nonminal.csv")
#ToxCast active
toxcast <- read.csv(file.path("POD database_15che","toxcast ac50.csv"))
tox.5perc <- data.frame(sapply(toxcast, quantile, prob=0.05, na.rm = T))
tox.5perc$CAS <- origpod$CAS
colnames(tox.5perc)[1] <- "ToxCast"

#iPSC+LCL
cell.t <- read.csv(file.path("POD database_15che","cellline POD distribution for 15 chemicals.csv"))
cell.min <- data.frame(sapply(cell.t[2:16], min, na.rm=T))
cell.min$CAS <- origpod$CAS
colnames(cell.min)[1] <- "iPSCLCL"

constant.POD <- data.frame(CAS = invitro_fub$CAS[1:15],
                           iPSCLCL = 1)

constant.POD.2 <- data.frame(CAS = invitro_fub$CAS[1:15],
                             iPSCLCL = 0.01)

#Cardio
#cardiomyocytes (single donor)
cardio.df <- read.csv(file.path("POD database_15che","Cardio POD for 15 chemicals.csv"))
cardio <- full_join(cardio.df, origpod[,c(1,2)], by="CAS")
cardio <- cardio[-which(is.na(cardio$Chemical.y)),]


for(i in 1:15){
  
  tox.5perc$ToxCast[i] <- (tox.5perc$ToxCast[i]*invitro_fub$VCBA[which(tox.5perc$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(tox.5perc$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(tox.5perc$CAS[i] == invivo_fup$CAS)])
  
  cell.min$iPSCLCL[i] <- (cell.min$iPSCLCL[i]*invitro_fub$VCBA[which(cell.min$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cell.min$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cell.min$CAS[i] == invivo_fup$CAS)])
  
  constant.POD$iPSCLCL[i] <- (constant.POD$iPSCLCL[i]*invitro_fub$VCBA[which(constant.POD$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD$CAS[i] == invivo_fup$CAS)])
  
  constant.POD.2$iPSCLCL[i] <- (constant.POD.2$iPSCLCL[i]*invitro_fub$VCBA[which(constant.POD.2$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(constant.POD.2$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(constant.POD.2$CAS[i] == invivo_fup$CAS)])
  
  cardio$Min[i] <- (cardio$Min[i]*invitro_fub$VCBA[which(cardio$CAS[i] == invitro_fub$CAS)])/(origpod$Css_h[which(cardio$CAS[i] == origpod$CAS)]*invivo_fup$fup[which(cardio$CAS[i] == invivo_fup$CAS)])
  
}

##traditional POD
pod.tox.df <- full_join(origpod, tox.5perc, by="CAS") 
pod.tox.s <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Spearman"))), digits=2)
pod.tox.p <- round((weightedCorr(log10(pod.tox.df$PODCss),log10(pod.tox.df$ToxCast), method=c("Pearson"))), digits=2)
pod.tox.mae <- mean(abs(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast)))
pod.tox.me <- mean(log10(pod.tox.df$PODCss)-log10(pod.tox.df$ToxCast))

##traditional POD x DAF
poddaf.tox.df <- full_join(poddaf, tox.5perc, by="CAS") 
poddaf.tox.s <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Spearman"))), digits=2)
poddaf.tox.p <- round((weightedCorr(log10(poddaf.tox.df$poddaf.Css),log10(poddaf.tox.df$ToxCast), method=c("Pearson"))), digits=2)
poddaf.tox.mae <- mean(abs(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast)))
poddaf.tox.me <- mean(log10(poddaf.tox.df$poddaf.Css)-log10(poddaf.tox.df$ToxCast))

##BBMD
bbmd.tox.df <- full_join(BMD, tox.5perc, by="CAS") 
wt.bmd.tox <- 1/(log10(bbmd.tox.df$BMD95.Css)-log10(bbmd.tox.df$BMD5.Css))^2
bbmd.tox.w.s <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Spearman"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.w.p <- round((weightedCorr(log10(bbmd.tox.df$BMD50.Css),log10(bbmd.tox.df$ToxCast), method=c("Pearson"), weights=wt.bmd.tox)), digits=2)
bbmd.tox.mae <- weighted.mean(abs(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast)), wt.bmd.tox)
bbmd.tox.me <- weighted.mean(log10(bbmd.tox.df$BMD50.Css)-log10(bbmd.tox.df$ToxCast), wt.bmd.tox)

##HED
hed.tox.df <- full_join(HED, tox.5perc, by="CAS") 
wt.hed.tox <- 1/(log10(hed.tox.df$HED95.Css)-log10(hed.tox.df$HED5.Css))^2
hed.tox.w.s <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Spearman"), weights=wt.hed.tox)), digits=2)
hed.tox.w.p <- round((weightedCorr(log10(hed.tox.df$HED50.Css),log10(hed.tox.df$ToxCast), method=c("Pearson"), weights=wt.hed.tox)), digits=2)
hed.tox.mae <- weighted.mean(abs(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast)), wt.hed.tox)
hed.tox.me <- weighted.mean(log10(hed.tox.df$HED50.Css)-log10(hed.tox.df$ToxCast), wt.hed.tox)


##in vivo POD vs iPSC+LCL POD scatter plot
##traditional POD
pod.cell.df <- full_join(origpod, cell.min, by="CAS") 
pod.cell.s <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.cell.p <- round((weightedCorr(log10(pod.cell.df$PODCss),log10(pod.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.cell.mae <- mean(abs(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL)))
pod.cell.me <- mean(log10(pod.cell.df$PODCss)-log10(pod.cell.df$iPSCLCL))

##traditional POD x DAF
poddaf.cell.df <- full_join(poddaf, cell.min, by="CAS") 
poddaf.cell.s <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.cell.p <- round((weightedCorr(log10(poddaf.cell.df$poddaf.Css),log10(poddaf.cell.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.cell.mae <- mean(abs(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL)))
poddaf.cell.me <- mean(log10(poddaf.cell.df$poddaf.Css)-log10(poddaf.cell.df$iPSCLCL))

##BBMD
bbmd.cell.df <- full_join(BMD, cell.min, by="CAS") 
wt.bmd.cell <- 1/(log10(bbmd.cell.df$BMD95.Css)-log10(bbmd.cell.df$BMD5.Css))^2
bbmd.cell.w.s <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.w.p <- round((weightedCorr(log10(bbmd.cell.df$BMD50.Css),log10(bbmd.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.cell)), digits=2)
bbmd.cell.mae <- weighted.mean(abs(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL)),wt.bmd.cell)
bbmd.cell.me <- weighted.mean(log10(bbmd.cell.df$BMD50.Css)-log10(bbmd.cell.df$iPSCLCL),wt.bmd.cell)

##HED
hed.cell.df <- full_join(HED, cell.min, by="CAS") 
wt.hed.cell <- 1/(log10(hed.cell.df$HED95.Css)-log10(hed.cell.df$HED5.Css))^2
hed.cell.w.s <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.cell)), digits=2)
hed.cell.w.p <- round((weightedCorr(log10(hed.cell.df$HED50.Css),log10(hed.cell.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.cell)), digits=2)
hed.cell.mae <- weighted.mean(abs(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL)),wt.hed.cell)
hed.cell.me <- weighted.mean(log10(hed.cell.df$HED50.Css)-log10(hed.cell.df$iPSCLCL),wt.hed.cell)

##in vivo POD vs Cardio POD scatter plot
##traditional POD
pod.cardio.df <- full_join(origpod, cardio, by="CAS") 
pod.car.s <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Spearman")), digits=2)
pod.car.p <- round(weightedCorr(log10(pod.cardio.df$PODCss),log10(pod.cardio.df$Min), method=c("Pearson")), digits=2)
pod.car.mae <- mean(abs(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min)))
pod.car.me <- mean(log10(pod.cardio.df$PODCss)-log10(pod.cardio.df$Min))

##traditional POD x DAF
poddaf.cardio.df <- full_join(poddaf, cardio, by="CAS") 
poddaf.car.s <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Spearman")), digits=2)
poddaf.car.p <- round(weightedCorr(log10(poddaf.cardio.df$poddaf.Css),log10(poddaf.cardio.df$Min), method=c("Pearson")), digits=2)
poddaf.car.mae <- mean(abs(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min)))
poddaf.car.me <- mean(log10(poddaf.cardio.df$poddaf.Css)-log10(poddaf.cardio.df$Min))

##BBMD
bbmd.cardio.df <- full_join(BMD, cardio, by="CAS") 
wt.bmd.cardio <- 1/(log10(bbmd.cardio.df$BMD95.Css)-log10(bbmd.cardio.df$BMD5.Css))^2
bbmd.car.w.s <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Spearman"), weights=wt.bmd.cardio), digits=2)
bbmd.car.w.p <- round(weightedCorr(log10(bbmd.cardio.df$BMD50.Css),log10(bbmd.cardio.df$Min), method=c("Pearson"), weights=wt.bmd.cardio), digits=2)
bbmd.car.mae <- weighted.mean(abs(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min)),wt.bmd.cardio)
bbmd.car.me <- weighted.mean(log10(bbmd.cardio.df$BMD50.Css)-log10(bbmd.cardio.df$Min),wt.bmd.cardio)

##HED
hed.cardio.df <- full_join(HED, cardio, by="CAS")
wt.hed.cardio <- 1/(log10(hed.cardio.df$HED95.Css)-log10(hed.cardio.df$HED5.Css))^2
hed.car.w.s <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Spearman"), weights=wt.hed.cardio), digits=2)
hed.car.w.p <- round(weightedCorr(log10(hed.cardio.df$HED50.Css),log10(hed.cardio.df$Min), method=c("Pearson"), weights=wt.hed.cardio), digits=2)
hed.car.mae <- weighted.mean(abs(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min)),wt.hed.cardio)
hed.car.me <- weighted.mean(log10(hed.cardio.df$HED50.Css)-log10(hed.cardio.df$Min),wt.hed.cardio)


##in vivo POD vs constant POD scatter plot
##traditional POD
pod.constant.df <- full_join(origpod, constant.POD, by="CAS") 
pod.constant.s <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p <- round((weightedCorr(log10(pod.constant.df$PODCss),log10(pod.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae <- mean(abs(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL)))
pod.constant.me <- mean(log10(pod.constant.df$PODCss)-log10(pod.constant.df$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df <- full_join(poddaf, constant.POD, by="CAS") 
poddaf.constant.s <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p <- round((weightedCorr(log10(poddaf.constant.df$poddaf.Css),log10(poddaf.constant.df$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae <- mean(abs(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL)))
poddaf.constant.me <- mean(log10(poddaf.constant.df$poddaf.Css)-log10(poddaf.constant.df$iPSCLCL))

##BBMD
bbmd.constant.df <- full_join(BMD, constant.POD, by="CAS") 
wt.bmd.constant <- 1/(log10(bbmd.constant.df$BMD95.Css)-log10(bbmd.constant.df$BMD5.Css))^2
bbmd.constant.w.s <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.w.p <- round((weightedCorr(log10(bbmd.constant.df$BMD50.Css),log10(bbmd.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant)), digits=2)
bbmd.constant.mae <- weighted.mean(abs(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL)),wt.bmd.constant)
bbmd.constant.me <- weighted.mean(log10(bbmd.constant.df$BMD50.Css)-log10(bbmd.constant.df$iPSCLCL),wt.bmd.constant)

##HED
hed.constant.df <- full_join(HED, constant.POD, by="CAS") 
wt.hed.constant <- 1/(log10(hed.constant.df$HED95.Css)-log10(hed.constant.df$HED5.Css))^2
hed.constant.w.s <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant)), digits=2)
hed.constant.w.p <- round((weightedCorr(log10(hed.constant.df$HED50.Css),log10(hed.constant.df$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant)), digits=2)
hed.constant.mae <- weighted.mean(abs(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL)),wt.hed.constant)
hed.constant.me <- weighted.mean(log10(hed.constant.df$HED50.Css)-log10(hed.constant.df$iPSCLCL),wt.hed.constant)


##in vivo POD vs constant POD scatter plot (conc = 0.01)
##traditional POD
pod.constant.df.2 <- full_join(origpod, constant.POD.2, by="CAS") 
pod.constant.s.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
pod.constant.p.2 <- round((weightedCorr(log10(pod.constant.df.2$PODCss),log10(pod.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
pod.constant.mae.2 <- mean(abs(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL)))
pod.constant.me.2 <- mean(log10(pod.constant.df.2$PODCss)-log10(pod.constant.df.2$iPSCLCL))

##traditional POD x DAF
poddaf.constant.df.2 <- full_join(poddaf, constant.POD.2, by="CAS") 
poddaf.constant.s.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Spearman"))), digits=2)
poddaf.constant.p.2 <- round((weightedCorr(log10(poddaf.constant.df.2$poddaf.Css),log10(poddaf.constant.df.2$iPSCLCL), method=c("Pearson"))), digits=2)
poddaf.constant.mae.2 <- mean(abs(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL)))
poddaf.constant.me.2 <- mean(log10(poddaf.constant.df.2$poddaf.Css)-log10(poddaf.constant.df.2$iPSCLCL))

##BBMD
bbmd.constant.df.2 <- full_join(BMD, constant.POD.2, by="CAS") 
wt.bmd.constant.2 <- 1/(log10(bbmd.constant.df.2$BMD95.Css)-log10(bbmd.constant.df.2$BMD5.Css))^2
bbmd.constant.w.s.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.w.p.2 <- round((weightedCorr(log10(bbmd.constant.df.2$BMD50.Css),log10(bbmd.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.bmd.constant.2)), digits=2)
bbmd.constant.mae.2 <- weighted.mean(abs(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL)),wt.bmd.constant.2)
bbmd.constant.me.2 <- weighted.mean(log10(bbmd.constant.df.2$BMD50.Css)-log10(bbmd.constant.df.2$iPSCLCL),wt.bmd.constant.2)

##HED
hed.constant.df.2 <- full_join(HED, constant.POD.2, by="CAS") 
wt.hed.constant.2 <- 1/(log10(hed.constant.df.2$HED95.Css)-log10(hed.constant.df.2$HED5.Css))^2
hed.constant.w.s.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Spearman"), weights=wt.hed.constant.2)), digits=2)
hed.constant.w.p.2 <- round((weightedCorr(log10(hed.constant.df.2$HED50.Css),log10(hed.constant.df.2$iPSCLCL), method=c("Pearson"), weights=wt.hed.constant.2)), digits=2)
hed.constant.mae.2 <- weighted.mean(abs(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL)),wt.hed.constant.2)
hed.constant.me.2 <- weighted.mean(log10(hed.constant.df.2$HED50.Css)-log10(hed.constant.df.2$iPSCLCL),wt.hed.constant.2)

POD.data.all_VCBA <- data.frame(Group = "Using VCBA model",
                                CAS = c(pod.tox.df$CAS, poddaf.tox.df$CAS, bbmd.tox.df$CAS, hed.tox.df$CAS,
                                        pod.cell.df$CAS, poddaf.cell.df$CAS, bbmd.cell.df$CAS, hed.cell.df$CAS,
                                        pod.cardio.df$CAS, poddaf.cardio.df$CAS, bbmd.cardio.df$CAS, hed.cardio.df$CAS,
                                        pod.constant.df$CAS, poddaf.constant.df$CAS, bbmd.constant.df$CAS, hed.constant.df$CAS,
                                        pod.constant.df.2$CAS, poddaf.constant.df.2$CAS, bbmd.constant.df.2$CAS, hed.constant.df.2$CAS),
                                x.val = c(pod.tox.df$ToxCast, poddaf.tox.df$ToxCast, bbmd.tox.df$ToxCast, hed.tox.df$ToxCast,
                                          pod.cell.df$iPSCLCL, poddaf.cell.df$iPSCLCL, bbmd.cell.df$iPSCLCL, hed.cell.df$iPSCLCL,
                                          pod.cardio.df$Min, poddaf.cardio.df$Min, bbmd.cardio.df$Min, hed.cardio.df$Min,
                                          pod.constant.df$iPSCLCL, poddaf.constant.df$iPSCLCL, bbmd.constant.df$iPSCLCL, hed.constant.df$iPSCLCL,
                                          pod.constant.df.2$iPSCLCL, poddaf.constant.df.2$iPSCLCL, bbmd.constant.df.2$iPSCLCL, hed.constant.df.2$iPSCLCL),
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 60),
                                y.val = c(pod.tox.df$PODCss, poddaf.tox.df$poddaf.Css, bbmd.tox.df$BMD50.Css, hed.tox.df$HED50.Css,
                                          pod.cell.df$PODCss, poddaf.cell.df$poddaf.Css, bbmd.cell.df$BMD50.Css, hed.cell.df$HED50.Css,
                                          pod.cardio.df$PODCss, poddaf.cardio.df$poddaf.Css, bbmd.cardio.df$BMD50.Css, hed.cardio.df$HED50.Css,
                                          pod.constant.df$PODCss, poddaf.constant.df$poddaf.Css, bbmd.constant.df$BMD50.Css, hed.constant.df$HED50.Css,
                                          pod.constant.df.2$PODCss, poddaf.constant.df.2$poddaf.Css, bbmd.constant.df.2$BMD50.Css, hed.constant.df.2$HED50.Css),
                                y.val.min = c(rep(NA, 30), bbmd.tox.df$BMD5.Css, hed.tox.df$HED5.Css,
                                              rep(NA, 30), bbmd.cell.df$BMD5.Css, hed.cell.df$HED5.Css,
                                              rep(NA, 30), bbmd.cardio.df$BMD5.Css, hed.cardio.df$HED5.Css,
                                              rep(NA, 30), bbmd.constant.df$BMD5.Css, hed.constant.df$HED5.Css,
                                              rep(NA, 30), bbmd.constant.df.2$BMD5.Css, hed.constant.df.2$HED5.Css),
                                y.val.max = c(rep(NA, 30), bbmd.tox.df$BMD95.Css, hed.tox.df$HED95.Css,
                                              rep(NA, 30), bbmd.cell.df$BMD95.Css, hed.cell.df$HED95.Css,
                                              rep(NA, 30), bbmd.cardio.df$BMD95.Css, hed.cardio.df$HED95.Css,
                                              rep(NA, 30), bbmd.constant.df$BMD95.Css, hed.constant.df$HED95.Css,
                                              rep(NA, 30), bbmd.constant.df.2$BMD95.Css, hed.constant.df.2$HED95.Css),
                                y.type = rep(rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"), each =15), 5))
RHO.data.all_VCBA <- data.frame(Group = "Using VCBA model",
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                rho.value = c(pod.tox.s, poddaf.tox.s, bbmd.tox.w.s, hed.tox.w.s,
                                              pod.cell.s, poddaf.cell.s, bbmd.cell.w.s, hed.cell.w.s,
                                              pod.car.s, poddaf.car.s, bbmd.car.w.s, hed.car.w.s,
                                              pod.constant.s, poddaf.constant.s, bbmd.constant.w.s, hed.constant.w.s,
                                              pod.constant.s.2, poddaf.constant.s.2, bbmd.constant.w.s.2, hed.constant.w.s.2))
R.data.all_VCBA <- data.frame(Group = "Using VCBA model",
                              x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                              y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                              rho.value = c(pod.tox.p, poddaf.tox.p, bbmd.tox.w.p, hed.tox.w.p,
                                            pod.cell.p, poddaf.cell.p, bbmd.cell.w.p, hed.cell.w.p,
                                            pod.car.p, poddaf.car.p, bbmd.car.w.p, hed.car.w.p,
                                            pod.constant.p, poddaf.constant.p, bbmd.constant.w.p, hed.constant.w.p,
                                            pod.constant.p.2, poddaf.constant.p.2, bbmd.constant.w.p.2, hed.constant.w.p.2))

mae.data.all_VCBA <- data.frame(Group = "Using VCBA model",
                                x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                                y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                                rho.value = c(pod.tox.mae, poddaf.tox.mae, bbmd.tox.mae, hed.tox.mae,
                                              pod.cell.mae, poddaf.cell.mae, bbmd.cell.mae, hed.cell.mae,
                                              pod.car.mae, poddaf.car.mae, bbmd.car.mae, hed.car.mae,
                                              pod.constant.mae, poddaf.constant.mae, bbmd.constant.mae, hed.constant.mae,
                                              pod.constant.mae.2, poddaf.constant.mae.2, bbmd.constant.mae.2, hed.constant.mae.2))
me.data.all_VCBA <- data.frame(Group = "Using VCBA model",
                               x.type = rep(c("ToxCast POD", "Six-cell-type POD", "hiPSC-CM", "Constant POD (1 uM)", "Constant POD  (0.01 uM)"), each = 4),
                               y.type = rep(c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"),5),
                               rho.value = c(pod.tox.me, poddaf.tox.me, bbmd.tox.me, hed.tox.me,
                                             pod.cell.me, poddaf.cell.me, bbmd.cell.me, hed.cell.me,
                                             pod.car.me, poddaf.car.me, bbmd.car.me, hed.car.me,
                                             pod.constant.me, poddaf.constant.me, bbmd.constant.me, hed.constant.me,
                                             pod.constant.me.2, poddaf.constant.me.2, bbmd.constant.me.2, hed.constant.me.2))



# --------- Plot -----------
POD.data.all <- rbind(POD.data.all_withoutadj, POD.data.all_Fischer, POD.data.all_IVMBM, POD.data.all_VIVD, POD.data.all_VCBA)
POD.data.all$Group <- factor(POD.data.all$Group, levels = c("Without adjustment", "Using Fischer model", "Using IVMBM model", "Using VIVD model", "Using VCBA model"),
                             labels = c("No Adjustmemt","Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
POD.data.all$x.type <- factor(POD.data.all$x.type, levels = c("Constant POD  (0.01 uM)", "Constant POD (1 uM)", "ToxCast POD", "hiPSC-CM", "Six-cell-type POD"))
POD.data.all$y.type <- factor(POD.data.all$y.type, levels = c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"))
POD.data.all <- POD.data.all[-which(POD.data.all$x.type == "hiPSC-CM"),]


RHO.data.all <- rbind(RHO.data.all_withoutadj, RHO.data.all_Fischer, RHO.data.all_IVMBM, RHO.data.all_VIVD, RHO.data.all_VCBA)
RHO.data.all$Group <- factor(RHO.data.all$Group, levels = c("Without adjustment", "Using Fischer model", "Using IVMBM model", "Using VIVD model", "Using VCBA model"),
                             labels = c("No Adjustmemt","Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
RHO.data.all$x.type <- factor(RHO.data.all$x.type, levels = c("Constant POD  (0.01 uM)", "Constant POD (1 uM)", "ToxCast POD", "hiPSC-CM", "Six-cell-type POD"))
RHO.data.all$y.type <- factor(RHO.data.all$y.type, levels = c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"))
RHO.data.all <- RHO.data.all[-which(RHO.data.all$x.type == "hiPSC-CM"),]

R.data.all <- rbind(R.data.all_withoutadj, R.data.all_Fischer, R.data.all_IVMBM, R.data.all_VIVD, R.data.all_VCBA)
R.data.all$Group <- factor(R.data.all$Group, levels = c("Without adjustment", "Using Fischer model", "Using IVMBM model", "Using VIVD model", "Using VCBA model"),
                           labels = c("No Adjustmemt","Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
R.data.all$x.type <- factor(R.data.all$x.type, levels = c("Constant POD  (0.01 uM)", "Constant POD (1 uM)", "ToxCast POD", "hiPSC-CM", "Six-cell-type POD"))
R.data.all$y.type <- factor(R.data.all$y.type, levels = c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"))
R.data.all <- R.data.all[-which(R.data.all$x.type == "hiPSC-CM"),]

mae.data.all <- rbind(mae.data.all_withoutadj, mae.data.all_Fischer, mae.data.all_IVMBM, mae.data.all_VIVD, mae.data.all_VCBA)
mae.data.all$Group <- factor(mae.data.all$Group, levels = c("Without adjustment", "Using Fischer model", "Using IVMBM model", "Using VIVD model", "Using VCBA model"),
                             labels = c("No Adjustmemt","Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
mae.data.all$x.type <- factor(mae.data.all$x.type, levels = c("Constant POD  (0.01 uM)", "Constant POD (1 uM)", "ToxCast POD", "hiPSC-CM", "Six-cell-type POD"))
mae.data.all$y.type <- factor(mae.data.all$y.type, levels = c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"))
mae.data.all <- mae.data.all[-which(mae.data.all$x.type == "hiPSC-CM"),]

me.data.all <- rbind(me.data.all_withoutadj, me.data.all_Fischer, me.data.all_IVMBM, me.data.all_VIVD, me.data.all_VCBA)
me.data.all$Group <- factor(me.data.all$Group, levels = c("Without adjustment", "Using Fischer model", "Using IVMBM model", "Using VIVD model", "Using VCBA model"),
                            labels = c("No Adjustmemt","Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
me.data.all$x.type <- factor(me.data.all$x.type, levels = c("Constant POD  (0.01 uM)", "Constant POD (1 uM)", "ToxCast POD", "hiPSC-CM", "Six-cell-type POD"))
me.data.all$y.type <- factor(me.data.all$y.type, levels = c("Reg POD_A", "Reg POD_H", "BMA BMD_A", "BMA BMD_H"))
me.data.all <- me.data.all[-which(me.data.all$x.type == "hiPSC-CM"),]

scatter <- 
  ggplot(POD.data.all[which(POD.data.all$y.type == "BMA BMD_H"),], aes(x=x.val, y=y.val, color=Group))+
  geom_point()+
  geom_errorbar(aes(ymin=y.val.min, ymax=y.val.max))+
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 3, color = "grey")+
  geom_smooth(method=lm, formula = y~offset(x), se=FALSE, size = 0.5, linetype = 4, color = "#EE5044")+
  geom_smooth(method=lm, formula = y~x, se=FALSE, size = 0.5, linetype = 2, color = "#3B7BE3")+
  geom_text(
    data    = mae.data.all[which(mae.data.all$y.type == "BMA BMD_H"),],
    mapping = aes(x=30, y=1E-4, label = paste0("MAE ==", round(rho.value, digits = 2))), parse = TRUE,
    color = "#EE5044"
  )+
  
  geom_text(
    data    = me.data.all[which(me.data.all$y.type == "BMA BMD_H"),],
    mapping = aes(x=30, y=1E-6, label = paste0("ME ==", round(rho.value, digits = 2))), parse = TRUE,
    color = "#EE5044"
  )+
  
  geom_text(
    data    = RHO.data.all[which(RHO.data.all$y.type == "BMA BMD_H"),],
    mapping = aes(x=1E-5, y=1E+1, label = paste0("rho ==", rho.value)), parse = TRUE,
    color = "#3B7BE3"
  )+
  
  geom_text(
    data    = R.data.all[which(R.data.all$y.type == "BMA BMD_H"),],
    mapping = aes(x=1E-5, y=1E+3, label = paste0("r ==", rho.value)), parse = TRUE,
    color = "#3B7BE3"
  )+
  
  scale_color_manual(values = c("black", "black", "black", "black", "black"))+
  facet_nested(Group~x.type, nest_line = element_line(linetype = 1))+
  scale_x_log10(limits=c(1e-7, 1e+4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(limits=c(1e-7, 1e+4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("OED from in vitro (mg/kg/day)")+ylab("Reg OED (mg/kg/day)")+
  theme_bw()+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        ggh4x.facet.nestline = element_line(size = 0.5),
        legend.position = "none", 
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", size=12),
        axis.text = element_text(color = "black", size=9))
scatter
ggsave("Figure 5.jpeg", width=11, height=10) #Saving 8.57 x 8.96 in image
ggsave("Figure 5.pdf", width=11, height=10)


library(plyr)
library(ggnewscale)
bar.plot.data <- rbind(R.data.all, RHO.data.all, mae.data.all, me.data.all)
bar.plot.data$Group_2 <- as.character(bar.plot.data$Group)
bar.plot.data$Group_2[which(bar.plot.data$Group_2 != "No Adjustmemt")] <- "Fub Adjustmemt"
bar.plot.data$Group_2 <- factor(bar.plot.data$Group_2, levels = c("No Adjustmemt", "Fub Adjustmemt"))
bar.plot.data$index <- rep(c("Pearson's r", "Spearman's rho", "MAE", "ME"), each = 80)
bar.plot.data$index <- factor(bar.plot.data$index, levels = c("Pearson's r", "Spearman's rho", "MAE", "ME"))
bar.plot.data_2 <- ddply(bar.plot.data, .(x.type, y.type, Group_2, index), summarize, median=median(rho.value))

ggplot() + 
  geom_col(data =bar.plot.data_2, aes(x=x.type, y=median, fill = Group_2, color = Group_2), position ="dodge2") + 
  scale_fill_manual(values = c("white", "#B8BEF5"))+
  scale_color_manual(values = c("black", "#B8BEF5"))+
  new_scale_colour()+
  geom_point(data =bar.plot.data, aes(x=x.type, y=rho.value, shape = Group, group = Group_2), position = position_dodge(width = 0.9))+
  #scale_color_manual(values = c("#00000000", "black", "black", "black", "black", "black"))+
  geom_hline(yintercept = 0)+
  xlab("In vitro POD")+ylab("Value")+
  facet_grid(index~y.type, scale = "free_y") + #facets to add gaps
  
  facetted_pos_scales(
    y = list(
      index == "Pearson's r" ~ scale_y_continuous(limits = c(0, 1)),
      index == "Spearman's rho" ~ scale_y_continuous(limits = c(0, 1)),
      index == "MAE" ~ scale_y_continuous(limits = c(0, 3.5)),
      index == "ME" ~ scale_y_continuous(limits = c(-3.5, 3.5))
    )
  )+
  theme_bw()+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        ggh4x.facet.nestline = element_line(size = 0.5),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=10),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black", size=12),
        axis.text.y = element_text(color = "black", size=10),
        axis.text.x = element_text(angle = 45, color = "black", size=10, hjust=1))
ggsave("Figure S11.jpeg", width = 12.8, height = 6.29) #Saving 12.8 x 6.29 in image
ggsave("Figure S11.pdf", width = 12.8, height = 6.29) 


