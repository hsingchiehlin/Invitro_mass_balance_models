library(patchwork)
library(car)
library(plyr)
library(compositions)
library(ggh4x)
library(tidyverse)
library(scatterplot3d)
library(plotly)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(xlsx)

# --------------- che space -------------
chemical_data1 <- read.csv("chemical_dataset1.csv")
chemical_data2 <- read.csv("chemical_dataset2.csv")
chemical_data3 <- read.csv("chemical_dataset3.csv")

chemical_data1$Dataset <- "116 chemicals"
chemical_data2$Dataset <- "CERAPP chemicals"
chemical_data3$Dataset <- "15 chemicals"
colnames(chemical_data1)[1] <- "List"
colnames(chemical_data2)[1] <- "List"

ALL_Che_PCAdat <- rbind(chemical_data1, chemical_data2, chemical_data3)
ALL_Che_PCAdat_clean <- na.omit(ALL_Che_PCAdat)

Che_pc <- prcomp(ALL_Che_PCAdat_clean[,5:13], scale  = TRUE, center = TRUE)
summary(Che_pc)
pca_data <- as.data.frame(Che_pc$x[, 1:3])
pca_data$Dataset <- ALL_Che_PCAdat_clean$Dataset

jpeg(width = 9.08, height = 7.68, "Figure 2.jpeg", units = "in", res=300)
pdf(width = 9.08, height = 7.68, "Figure 2.pdf")
layout(matrix(c(1,3,2,4), 2))
colors <- c("#D53A0B","#5281E000","#99999903")
colors <- colors[as.factor(ALL_Che_PCAdat$Dataset)]
scatterplot3d(ALL_Che_PCAdat$OPERAOctanolWaterPartitionCoefficient,
              ALL_Che_PCAdat$OPERAWaterSolubility,
              ALL_Che_PCAdat$OPERAVaporPressure,
              xlab="LogKow",
              ylab="",
              zlab="Vapor Pressure",
              color = colors, main="A", pch = 16, angle=30)
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,expression("Solubility"),srt=30)

colors <- c("#D53A0B","#5281E000","#99999903")
colors <- colors[as.factor(pca_data$Dataset)]
scatterplot3d(pca_data[,1:3], pch = 16, color=colors, angle=30, main="B", ylab="")
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,expression(PC2),srt=30)


colors <- c("#D53A0B00","#5281E0","#99999903")
colors <- colors[as.factor(ALL_Che_PCAdat$Dataset)]
scatterplot3d(ALL_Che_PCAdat$OPERAOctanolWaterPartitionCoefficient,
              ALL_Che_PCAdat$OPERAWaterSolubility,
              ALL_Che_PCAdat$OPERAVaporPressure,
              xlab="LogKow",
              ylab="",
              zlab="Vapor Pressure",
              color = colors, main="C", pch = 16, angle=30)
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,expression("Solubility"),srt=30)

colors <- c("#D53A0B00","#5281E0","#99999903")
colors <- colors[as.factor(pca_data$Dataset)]
scatterplot3d(pca_data[,1:3], pch = 16, color=colors, angle=30, main="D", ylab="")
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,expression(PC2),srt=30)
dev.off()


# 2D
p.1 <- 
  ggplot(ALL_Che_PCAdat, aes(OPERAOctanolWaterPartitionCoefficient,
                             OPERAWaterSolubility,color=Dataset))+
  geom_point()+
  xlab("LogKow")+ylab("Solubility")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())
p.2 <- 
  ggplot(ALL_Che_PCAdat, aes(OPERAOctanolWaterPartitionCoefficient,
                             OPERAVaporPressure,color=Dataset))+
  geom_point()+
  xlab("LogKow")+ylab("Vapor Pressure")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())
p.3 <- 
  ggplot(ALL_Che_PCAdat, aes(OPERAWaterSolubility,
                             OPERAVaporPressure,color=Dataset))+
  geom_point()+
  xlab("Solubility")+ylab("Vapor Pressure")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())


p.4 <- 
  ggplot(pca_data, aes(PC1,PC2,color=Dataset))+
  geom_point()+
  #xlab("LogKow")+ylab("Solubility")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())
p.5 <- 
  ggplot(pca_data, aes(PC1,PC3,color=Dataset))+
  geom_point()+
  #xlab("Solubility")+ylab("Vapor Pressure")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())
p.6 <- 
  ggplot(pca_data, aes(PC2,PC3,color=Dataset))+
  geom_point()+
  #xlab("LogKow")+ylab("Vapor Pressure")+
  scale_color_manual(values = c("#D53A0B","#5281E0","#99999908"))+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())

jpeg(width = 10.6, height = 6.5, "Figure.S1.jpeg", units = "in", res=300)
p.1 + p.2 + p.3 + p.4 + p.5 + p.6 + plot_layout(widths = c(1/3,1/3,1/3),
                                                heights = c(1/2,1/2))
dev.off()

pdf(width = 10.6, height = 6.5, "Figure.S1.pdf")
p.1 + p.2 + p.3 + p.4 + p.5 + p.6 + plot_layout(widths = c(1/3,1/3,1/3),
                                                heights = c(1/2,1/2))
dev.off()

write.csv(ALL_Che_PCAdat, "ALL_Che_PCAdat.csv")
write.xlsx(ALL_Che_PCAdat, "ALL_Che_PCAdat.xlsx")


#  ------------- Model comparesion: Figure S2, S4  -------------
Val.results.pest_single <- read.csv("Val.results.pest.forplot.csv")
Val.results.pest_mixture <- read.csv("Val.results.pest.forplot_mixture.csv")
Val.results.hus <- read.csv("Val.results.hus.forplot.csv")
Val.results.sch <- read.csv("Val.results.sch.forplot.csv")
Val.results.tan <- read.csv("Val.results.tan.forplot.csv")
Val.results.pfas <- read.csv("Val.results.pfas.forplot.csv")
Val.results.unilever <- read.csv("Val.results.unilever.forplot.csv")
Val.results.Blanchette <- read.csv("Val.results.Blanchette.forplot.csv")

diff.data_single <- aggregate(.~Chemical, data = Val.results.pest_single[,c(1:11)], mean)
Val.results.pest_single <- diff.data_single[, c(1:6)]
Val.results.pest_single[c(21:40),] <- diff.data_single[, c(1,7:11)]

diff.data_mixture <- aggregate(.~Chemical, data = Val.results.pest_mixture[,c(1:11)], mean)
Val.results.pest_mixture <- diff.data_mixture[, c(1:6)]
Val.results.pest_mixture[c(21:40),] <- diff.data_mixture[, c(1,7:11)]
Val.results.pest <- rbind(Val.results.pest_single, Val.results.pest_mixture)
colnames(Val.results.pest)[2] <- "Obs"
Val.results.pest_CAS <- read.csv("Val.results.pest.forplot.csv")$CAS[1:20]
Val.results.pest$ref <- rep(c("Valdiviezo et al.(10 uM - Single)", "Valdiviezo et al.(1 uM - Single)",
                              "Valdiviezo et al.(10 uM- Mixture)", "Valdiviezo et al.(1 uM- Mixture)"), each = 20)
Val.results.pest$IOC <- "N"
Val.results.pest$IOC[which(Val.results.pest$Chemical == "2,4-Dinitrophenol")] <- "A"
Val.results.pest$CAS <- rep(Val.results.pest_CAS, 4)


diff.data_pfas <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,2)], mean)
diff.data_pfas[, c(3:6)] <- Val.results.pfas[c(1:14), c(3:6)]
diff.data_pfas$Obs.10uM <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,7)], mean)$Obs.10uM
diff.data_pfas[, c(8:11)] <- Val.results.pfas[c(1:14), c(8:11)]
Val.results.pfas <- diff.data_pfas[, c(1:6)]
Val.results.pfas[c(15:28),] <- diff.data_pfas[, c(1,7:11)]
colnames(Val.results.pfas)[2] <- "Obs"
Val.results.pfas$ref <- rep(c("PFAS (100 uM)", "PFAS (10 uM)"), each = 14)
Val.results.pfas$IOC <- "A"
Val.results.pfas$CAS <- read.csv("Val.results.pfas.forplot.csv")$CAS[1:28]


Val.results.tan$ref <- "Tanneberger et al.(Fish RTgill)"
Val.results.sch$ref <- rep(c("Schug et al.(Fish RTgut - High dose)", "Schug et al.(Fish RTgut - Low dose)"), each = 9)

comb.data <- rbind(Val.results.pest, Val.results.hus, Val.results.tan, Val.results.sch, Val.results.pfas, Val.results.Blanchette)
comb.data$MassBalanceIssue <- 0
comb.data <- rbind(comb.data, Val.results.unilever)

comb.data$'Cell type' <- c(rep("Cell-free", 80), rep("Human", 16), rep("Fish", 45),rep("Cell-free", 28), rep("Cell-free", 90), rep("Cell-free", 10))
comb.data$'Data source' <- c(rep("Our data", 80), rep("Literature", 16), rep("Literature", 45),rep("Our data", 28), rep("Literature", 90), rep("Our data", 10))
comb.data$conc <- c(rep(10,20),rep(1,20),rep(10,20),rep(1,20), 
                    c(3680,863,89.4,316,411,263,1100,215,2440,439,2.12,49.2,17.7,316,184,242,
                      0.835586139,0.049152126,14037.07518,98787.4465,39.34624697,0.075103267,4.1281912,712.8514056,2.883625129,310.5310531,93.3125972,23.00025556,4749.731472,9.110787172,1913.591359,10.12658228,14359.72461,33.3562586, 1509.328358,47.83950617,2.399314482,219.1887676,78.2369146,182.9931973,17656.06596,195.4161641,1000,
                      6.090133983,165.3439153,11.34429949,5.200208008,0.50942435,5.042864347,0.449640288,1.057529611,0.03869969,1827.040195,5291.005291,680.6579694,104.0041602,61.13092206,75.6429652	,26.97841727,21.15059222,7.73993808), rep(100,14),rep(10,14),rep(10,10),rep(10,90))
che.prop.list <- read.csv("che.prop.list.csv")
plot.data <- data.frame(Chemical = rep(comb.data$Chemical, 4),
                        CAS = rep(comb.data$CAS, 4),
                        Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 269),
                        Obs.val = rep(comb.data$Obs, 4),
                        pre.val = c(comb.data$X.free_Fischer,
                                    comb.data$X.free_VIVD,
                                    comb.data$X.free_IVMBM,
                                    comb.data$X.free_VCBA),
                        ref = rep(comb.data$ref, 4),
                        `Cell type` = rep(comb.data$`Cell type`, 4),
                        `Data source` = rep(comb.data$`Data source`, 4),
                        MassBalanceIssue = rep(comb.data$MassBalanceIssue, 4),
                        IOC = rep(comb.data$IOC, 4),
                        conc = rep(comb.data$conc, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}

#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Cell.type <- factor(plot.data$Cell.type, levels = c("Human", "Fish", "Cell-free"))
plot.data$Data.source <- factor(plot.data$Data.source, levels = c("Our data", "Literature"))
plot.data$MassBalanceIssue <- factor(plot.data$MassBalanceIssue)
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))

plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)" & plot.data$IOC != "N"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_3A <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(color = MassBalanceIssue,shape = Data.source), alpha=0.8) +
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/8, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/4, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 1-8/5, y = -7+8/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 1-8/5, y = -7+8/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 1-8/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental free percentage)") +
  ylab("log10(Predicted free percentage)")+
  ggtitle("A. Fraction of free chemical concentration in media")+# (all datapoints)
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 1)+
  ylim(-7, 1)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("black", "black"))+ #"#FB684B", "#058ED9"
  scale_shape_manual(values = c(21, 16))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))

Fig_S2A <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(aes(fill = IOC.lab), position = "stack", color = "black") + 
  #geom_text(
  #  data    = MAD,
  #  mapping = aes(x = -4, y = 50, label = paste("MAE = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  #geom_text(
  #  data    = ME,
  #  mapping = aes(x = -4, y = 40, label = paste("ME = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  xlim(-5.5,3)+
  facet_grid(~Model)+ #, scales = "free_x"
  geom_vline(xintercept = 0,
             linetype="dashed", color = "dark red")+
  theme_bw() + 
  xlab("log10(Predicted F_free) - log10(Observed F_free)") +
  ylab("Count")+
  ggtitle("A. Fraction of free chemical concentration in media")+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))

plot.data_freeratio <- plot.data
write.csv(plot.data_freeratio, "freeratio_modelcomp_data.csv")



Fig_S0A <- # mass balance
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(color = MassBalanceIssue,shape = Data.source), alpha=0.8) +
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/8, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/4, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 1-8/5, y = -7+8/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 1-8/5, y = -7+8/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 1-8/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental free percentage)") +
  ylab("log10(Predicted free percentage)")+
  ggtitle("A. All datapoints")+# (all datapoints)
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 1)+
  ylim(-7, 1)+
  scale_color_manual(values = c("black", "#FB684B"))+ #"#FB684B", "#058ED9"
  scale_shape_manual(values = c(21, 16))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))




plot.data <- plot.data[which(plot.data$MassBalanceIssue == 0),]


corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S0B <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(color = MassBalanceIssue,shape = Data.source), alpha=0.8) +
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/8, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/4, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 1-8/5, y = -7+8/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 1-8/5, y = -7+8/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 1-8/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  xlab("log10(Experimental free percentage)") +
  ylab("log10(Predicted free percentage)")+
  ggtitle("B. Removing datapoints with mass balance issue")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 1)+
  ylim(-7, 1)+
  scale_color_manual(values = c("black", "#FB684B"))+ #"#FB684B", "#058ED9"
  scale_shape_manual(values = c(21, 16))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))

jpeg(width = 10, height = 5.8, "Figure S_massbalance.jpeg", units = "in", res=300)
Fig_S0A + Fig_S0B + plot_layout(height = c(1/2,1/2))
dev.off()
pdf(width = 10, height = 5.8, "Figure S_massbalance.pdf")
Fig_S0A + Fig_S0B + plot_layout(height = c(1/2,1/2))
dev.off()




#prepare data (SOT poster data)
Val.results <- read.csv("Val.results.forplot_medconc.csv")

plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 27),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)" & plot.data$IOC != "N"),]


corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))


MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_3B <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species),color = "black", alpha=0.8) +#, shape=Cell.type
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/8, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/4, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0.5-5.5/5, y = -5, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in media)") +
  ylab("log10(Predicted chemical amount in media)")+
  ggtitle("B. Amount of chemicals in media")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-5, 0.5)+
  ylim(-5, 0.5)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))


Fig_S2B <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(aes(fill = IOC.lab), position = "stack", color = "black") + 
  #geom_text(
  #  data    = MAD,
  #  mapping = aes(x = -4, y = 5, label = paste("MAE = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  #geom_text(
  #  data    = ME,
  #  mapping = aes(x = -4, y = 4.2, label = paste("ME = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  facet_grid(~Model)+ #, scales = "free_x"
  geom_vline(xintercept = 0,
             linetype="dashed", color = "dark red")+
  theme_bw() + 
  xlab("log10(Predicted chemical amount in media) - log10(Observed chemical amount in media)") +
  ylab("Count")+
  xlim(-5.5,3)+
  ggtitle("B. Amount of chemicals in media")+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
plot.data_media <- plot.data
write.csv(plot.data_media, "media_modelcomp_data.csv")


Val.results <- read.csv("Val.results.forplot_cellconc.csv")
plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 43),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)" & plot.data$IOC != "N"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_3C <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species), color = "black", alpha=0.8) + #aes(color = Species, shape=Cell.type)
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/8, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/4, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0-6/5, y = -7+6/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0-6/5, y = -7+6/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0-6/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in cells)") +
  ylab("log10(Predicted chemical amount in cells)")+
  ggtitle("C. Amount of chemicals in cells")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 0)+
  ylim(-7, 0)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))


Fig_S2C <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(aes(fill = IOC.lab), position = "stack", color = "black") + 
  #geom_text(
  #  data    = MAD,
  #  mapping = aes(x = -4, y = 8, label = paste("MAE = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  #geom_text(
  #  data    = ME,
  #  mapping = aes(x = -4, y = 6.8, label = paste("ME = ", round(mean, 2), sep = "")),
  #  color = "black"
  #)+
  xlim(-5.5,3)+
  facet_grid(~Model)+ #, scales = "free_x"
  geom_vline(xintercept = 0,
             linetype="dashed", color = "dark red")+
  theme_bw() + 
  xlab("log10(Predicted chemical amount in cells) - log10(Observed chemical amount in cells)") +
  ylab("Count")+
  ggtitle("C. Amount of chemicals in cells")+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
plot.data_cell <- plot.data
write.csv(plot.data_cell, "cell_modelcomp_data.csv")

jpeg(width = 10, height = 8, "Figure.S2.jpeg", units = "in", res=300)
Fig_S2A + Fig_S2B + Fig_S2C + plot_layout(height = c(1/3,1/3,1/3))
dev.off()

pdf(width = 10, height = 8, "Figure.S2.pdf")
Fig_S2A + Fig_S2B + Fig_S2C + plot_layout(height = c(1/3,1/3,1/3))
dev.off()

#Figure S4-S5:error trend
plot.data_media$obs.frac <- plot.data_media$Obs.val/plot.data_media$amount
plot.data_media$pre.frac <- plot.data_media$pre.val/plot.data_media$amount
plot.data_cell$obs.frac <- plot.data_cell$Obs.val/plot.data_cell$amount
plot.data_cell$pre.frac <- plot.data_cell$pre.val/plot.data_cell$amount

plot.data_media$diff <- log10(plot.data_media$pre.frac)-log10(plot.data_media$obs.frac)
plot.data_cell$diff <- log10(plot.data_cell$pre.frac)-log10(plot.data_cell$obs.frac)

Corr.plot.data <- rbind(plot.data_freeratio[,c(1:3,21,10:18)],
                        plot.data_media[,c(1:3,20,10:17,7)])



Corr.plot.data <- data.frame(Chemical = rep(Corr.plot.data$Chemical, 8),
                             Model = rep(Corr.plot.data$Model, 8),
                             Diff = rep(Corr.plot.data$diff, 8),
                             #Type = rep(Corr.plot.data$Type, 7),
                             Che.prop.value = c(log10(Corr.plot.data$H37),
                                                Corr.plot.data$MP,
                                                Corr.plot.data$log.KOW,
                                                Corr.plot.data$log.KAW,
                                                log10(Corr.plot.data$Solubility),
                                                Corr.plot.data$MW,
                                                log10(Corr.plot.data$conc),
                                                Corr.plot.data$pKa),
                             Che.prop.name = rep(c("log10_H37", "MP", "logkow", "logkaw", "log10_Solubility", "MW", "log10_Conc", "pKa"), each = 1047))

Corr.plot.data <- Corr.plot.data[-which(Corr.plot.data$Chemical == "cyclosporine A"),]

jpeg(width = 15, height = 9.1, "Figure S4.jpeg", units = "in", res=300)
pdf(width = 15, height = 9.1, "Figure S4.pdf")
ggplot(Corr.plot.data)+
  geom_smooth(aes(x = Che.prop.value, y = Diff), method = "lm", colour = "#7BD389", fill = "#7BD38920") + 
  geom_point(aes(x = Che.prop.value, y = Diff), color = "#004777") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3)+
  xlab("Chemical property value") +
  ylab("Log10(Predicted % of chemical in media / Measured % of chemical in media")+
  facet_grid(Model~Che.prop.name, scales = "free")+
  scale_shape_manual(values = c(16, 21))+
  #ylim(-6, 3)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        title = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        legend.position = "none",
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        #legend.position = c(0.055, 0.8), 
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
dev.off()

Corr.plot.data <- Corr.plot.data[-which(is.na(Corr.plot.data$Che.prop.value)),]
lm_eqn <- 
  function(Corr.plot.data){
    m <- lm(Diff ~ Che.prop.value, Corr.plot.data);
    
    Intercept <- summary(m)[["coefficients"]][1,1];
    Intercept.se <- summary(m)[["coefficients"]][1,2];
    Intercept.pval <- summary(m)[["coefficients"]][1,4];
    
    slope <- summary(m)[["coefficients"]][2,1];
    slope.se <- summary(m)[["coefficients"]][2,2];
    slope.pval <- summary(m)[["coefficients"]][2,4];
    
    r.squared <- summary(m)[["r.squared"]];
    adj.r.squared <- summary(m)[["r.squared"]];
    r.squared.pval <- pf(summary(m)$fstatistic[1], summary(m)$fstatistic[2], summary(m)$fstatistic[3], lower.tail=F);
    
    eq <- c(Intercept, Intercept.se, Intercept.pval, slope, slope.se, slope.pval, r.squared, adj.r.squared, r.squared.pval);                 
  }


corr.results.diff_media <- ddply(Corr.plot.data,.(Che.prop.name, Model),lm_eqn)
colnames(corr.results.diff_media)[3:11] <- c("Intercept", "Intercept.se", "Intercept.pval", 
                                             "slope", "slope.se", "slope.pval", 
                                             "r.squared", "adj.r.squared", "r.squared.pval")
write.csv(corr.results.diff_media, "corr.results.diff_media.csv")



Error_cell_prec.results <- plot.data_cell


Corr.plot.data <- data.frame(Chemical = rep(Error_cell_prec.results$Chemical, 8),
                             Model = rep(Error_cell_prec.results$Model, 8),
                             Diff = rep(Error_cell_prec.results$diff, 8),
                             Che.prop.value = c(log10(Error_cell_prec.results$H37),
                                                Error_cell_prec.results$MP,
                                                Error_cell_prec.results$log.KOW,
                                                Error_cell_prec.results$log.KAW,
                                                log10(Error_cell_prec.results$Solubility),
                                                Error_cell_prec.results$MW,
                                                log10(Error_cell_prec.results$conc),
                                                Error_cell_prec.results$pKa),
                             Che.prop.name = rep(c("log10_H37", "MP", "logkow", "logkaw", "log10_Solubility", "MW", "log10_Conc", "pKa"), each = 151))

Corr.plot.data <- Corr.plot.data[-which(Corr.plot.data$Chemical == "cyclosporine A"),]

jpeg(width = 15, height = 9.1, "Figure S5.jpeg", units = "in", res=300)
pdf(width = 15, height = 9.1, "Figure S5.pdf")
ggplot(Corr.plot.data)+
  geom_smooth(aes(x = Che.prop.value, y = Diff), method = "lm", colour = "#7BD389", fill = "#7BD38920") + 
  geom_point(aes(x = Che.prop.value, y = Diff), color = "#AB2C21") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3)+
  xlab("Chemical property value") +
  ylab("Log10(Predicted % of chemical in cell / Measured% of chemical in cell")+
  facet_grid(Model~Che.prop.name, scales = "free")+
  ylim(-4, 4)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        title = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        #legend.position = c(0.055, 0.8), 
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
dev.off()


lm_eqn <- 
  function(Corr.plot.data){
    m <- lm(Diff ~ Che.prop.value, Corr.plot.data);
    
    Intercept <- summary(m)[["coefficients"]][1,1];
    Intercept.se <- summary(m)[["coefficients"]][1,2];
    Intercept.pval <- summary(m)[["coefficients"]][1,4];
    
    slope <- summary(m)[["coefficients"]][2,1];
    slope.se <- summary(m)[["coefficients"]][2,2];
    slope.pval <- summary(m)[["coefficients"]][2,4];
    
    r.squared <- summary(m)[["r.squared"]];
    adj.r.squared <- summary(m)[["r.squared"]];
    r.squared.pval <- pf(summary(m)$fstatistic[1], summary(m)$fstatistic[2], summary(m)$fstatistic[3], lower.tail=F);
    
    eq <- c(Intercept, Intercept.se, Intercept.pval, slope, slope.se, slope.pval, r.squared, adj.r.squared, r.squared.pval);     
  }
Corr.plot.data <- Corr.plot.data[-which(is.na(Corr.plot.data$Che.prop.value)),]
corr.results.diff_cell <- ddply(Corr.plot.data,.(Che.prop.name, Model),lm_eqn)
colnames(corr.results.diff_cell)[3:11] <- c("Intercept", "Intercept.se", "Intercept.pval", 
                                            "slope", "slope.se", "slope.pval", 
                                            "r.squared", "adj.r.squared", "r.squared.pval")
write.csv(corr.results.diff_cell, "corr.results.diff_cell.csv")

#  ------------- Model comparision: Figure S3A  -------------
Val.results.pest_single <- read.csv("Val.results.pest.forplot.csv")
Val.results.pest_mixture <- read.csv("Val.results.pest.forplot_mixture.csv")
Val.results.hus <- read.csv("Val.results.hus.forplot.csv")
Val.results.sch <- read.csv("Val.results.sch.forplot.csv")
Val.results.tan <- read.csv("Val.results.tan.forplot.csv")
Val.results.pfas <- read.csv("Val.results.pfas.forplot.csv")
Val.results.unilever <- read.csv("Val.results.unilever.forplot.csv")
Val.results.Blanchette <- read.csv("Val.results.Blanchette.forplot.csv")

diff.data_single <- aggregate(.~Chemical, data = Val.results.pest_single[,c(1:11)], mean)
Val.results.pest_single <- diff.data_single[, c(1:6)]
Val.results.pest_single[c(21:40),] <- diff.data_single[, c(1,7:11)]

diff.data_mixture <- aggregate(.~Chemical, data = Val.results.pest_mixture[,c(1:11)], mean)
Val.results.pest_mixture <- diff.data_mixture[, c(1:6)]
Val.results.pest_mixture[c(21:40),] <- diff.data_mixture[, c(1,7:11)]
Val.results.pest <- rbind(Val.results.pest_single, Val.results.pest_mixture)
colnames(Val.results.pest)[2] <- "Obs"
Val.results.pest_CAS <- read.csv("Val.results.pest.forplot.csv")$CAS[1:20]
Val.results.pest$ref <- rep(c("Valdiviezo et al.(10 uM - Single)", "Valdiviezo et al.(1 uM - Single)",
                              "Valdiviezo et al.(10 uM- Mixture)", "Valdiviezo et al.(1 uM- Mixture)"), each = 20)
Val.results.pest$IOC <- "N"
Val.results.pest$IOC[which(Val.results.pest$Chemical == "2,4-Dinitrophenol")] <- "A"
Val.results.pest$CAS <- rep(Val.results.pest_CAS, 4)


diff.data_pfas <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,2)], mean)
diff.data_pfas[, c(3:6)] <- Val.results.pfas[c(1:14), c(3:6)]
diff.data_pfas$Obs.10uM <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,7)], mean)$Obs.10uM
diff.data_pfas[, c(8:11)] <- Val.results.pfas[c(1:14), c(8:11)]
Val.results.pfas <- diff.data_pfas[, c(1:6)]
Val.results.pfas[c(15:28),] <- diff.data_pfas[, c(1,7:11)]
colnames(Val.results.pfas)[2] <- "Obs"
Val.results.pfas$ref <- rep(c("PFAS (100 uM)", "PFAS (10 uM)"), each = 14)
Val.results.pfas$IOC <- "A"
Val.results.pfas$CAS <- read.csv("Val.results.pfas.forplot.csv")$CAS[1:28]


Val.results.tan$ref <- "Tanneberger et al.(Fish RTgill)"
Val.results.sch$ref <- rep(c("Schug et al.(Fish RTgut - High dose)", "Schug et al.(Fish RTgut - Low dose)"), each = 9)

comb.data <- rbind(Val.results.pest, Val.results.hus, Val.results.tan, Val.results.sch, Val.results.pfas, Val.results.Blanchette)
comb.data$MassBalanceIssue <- 0
comb.data <- rbind(comb.data, Val.results.unilever)

comb.data$'Cell type' <- c(rep("Cell-free", 80), rep("Human", 16), rep("Fish", 45),rep("Cell-free", 28), rep("Cell-free", 90), rep("Cell-free", 10))
comb.data$'Data source' <- c(rep("Our data", 80), rep("Literature", 16), rep("Literature", 45),rep("Our data", 28), rep("Literature", 90), rep("Our data", 10))
comb.data$conc <- c(rep(10,20),rep(1,20),rep(10,20),rep(1,20), 
                    c(3680,863,89.4,316,411,263,1100,215,2440,439,2.12,49.2,17.7,316,184,242,
                      0.835586139,0.049152126,14037.07518,98787.4465,39.34624697,0.075103267,4.1281912,712.8514056,2.883625129,310.5310531,93.3125972,23.00025556,4749.731472,9.110787172,1913.591359,10.12658228,14359.72461,33.3562586, 1509.328358,47.83950617,2.399314482,219.1887676,78.2369146,182.9931973,17656.06596,195.4161641,1000,
                      6.090133983,165.3439153,11.34429949,5.200208008,0.50942435,5.042864347,0.449640288,1.057529611,0.03869969,1827.040195,5291.005291,680.6579694,104.0041602,61.13092206,75.6429652	,26.97841727,21.15059222,7.73993808), rep(100,14),rep(10,14),rep(10,10),rep(10,90))
che.prop.list <- read.csv("che.prop.list.csv")
plot.data <- data.frame(Chemical = rep(comb.data$Chemical, 4),
                        CAS = rep(comb.data$CAS, 4),
                        Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 269),
                        Obs.val = rep(comb.data$Obs, 4),
                        pre.val = c(comb.data$X.free_Fischer,
                                    comb.data$X.free_VIVD,
                                    comb.data$X.free_IVMBM,
                                    comb.data$X.free_VCBA),
                        ref = rep(comb.data$ref, 4),
                        `Cell type` = rep(comb.data$`Cell type`, 4),
                        `Data source` = rep(comb.data$`Data source`, 4),
                        MassBalanceIssue = rep(comb.data$MassBalanceIssue, 4),
                        IOC = rep(comb.data$IOC, 4),
                        conc = rep(comb.data$conc, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}

#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Cell.type <- factor(plot.data$Cell.type, levels = c("Human", "Fish", "Cell-free"))
plot.data$Data.source <- factor(plot.data$Data.source, levels = c("Our data", "Literature"))
plot.data$MassBalanceIssue <- factor(plot.data$MassBalanceIssue)
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))

plot.data <- plot.data[-which(plot.data$IOC != "N"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3A1 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(color = MassBalanceIssue,shape = Data.source), alpha=0.8) +
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 1-8/5, y = -7+8/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 1-8/5, y = -7+8/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 1-8/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental free percentage)") +
  ylab("log10(Predicted free percentage)")+
  ggtitle("A1. Fraction of free chemical concentration in media")+# (all datapoints)
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 1)+
  ylim(-7, 1)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("black", "black"))+ #"#FB684B", "#058ED9"
  scale_shape_manual(values = c(21, 16))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))

#prepare data (SOT poster data)
Val.results <- read.csv("Val.results.forplot_medconc.csv")

plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 27),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$IOC != "N"),]


corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))


MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3A2 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species),color = "black", alpha=0.8) +#, shape=Cell.type
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0.5-5.5/5, y = -5, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in media)") +
  ylab("log10(Predicted chemical amount in media)")+
  ggtitle("A2. Amount of chemicals in media")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-5, 0.5)+
  ylim(-5, 0.5)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))





Val.results <- read.csv("Val.results.forplot_cellconc.csv")
plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 43),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$IOC != "N"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3A3 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species), color = "black", alpha=0.8) + #aes(color = Species, shape=Cell.type)
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0-6/5, y = -7+6/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0-6/5, y = -7+6/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0-6/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in cells)") +
  ylab("log10(Predicted chemical amount in cells)")+
  ggtitle("A3. Amount of chemicals in cells")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 0)+
  ylim(-7, 0)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))



#  ------------- Model comparision: Figure S3B  -------------
Val.results.pest_single <- read.csv("Val.results.pest.forplot.csv")
Val.results.pest_mixture <- read.csv("Val.results.pest.forplot_mixture.csv")
Val.results.hus <- read.csv("Val.results.hus.forplot.csv")
Val.results.sch <- read.csv("Val.results.sch.forplot.csv")
Val.results.tan <- read.csv("Val.results.tan.forplot.csv")
Val.results.pfas <- read.csv("Val.results.pfas.forplot.csv")
Val.results.unilever <- read.csv("Val.results.unilever.forplot.csv")
Val.results.Blanchette <- read.csv("Val.results.Blanchette.forplot.csv")

diff.data_single <- aggregate(.~Chemical, data = Val.results.pest_single[,c(1:11)], mean)
Val.results.pest_single <- diff.data_single[, c(1:6)]
Val.results.pest_single[c(21:40),] <- diff.data_single[, c(1,7:11)]

diff.data_mixture <- aggregate(.~Chemical, data = Val.results.pest_mixture[,c(1:11)], mean)
Val.results.pest_mixture <- diff.data_mixture[, c(1:6)]
Val.results.pest_mixture[c(21:40),] <- diff.data_mixture[, c(1,7:11)]
Val.results.pest <- rbind(Val.results.pest_single, Val.results.pest_mixture)
colnames(Val.results.pest)[2] <- "Obs"
Val.results.pest_CAS <- read.csv("Val.results.pest.forplot.csv")$CAS[1:20]
Val.results.pest$ref <- rep(c("Valdiviezo et al.(10 uM - Single)", "Valdiviezo et al.(1 uM - Single)",
                              "Valdiviezo et al.(10 uM- Mixture)", "Valdiviezo et al.(1 uM- Mixture)"), each = 20)
Val.results.pest$IOC <- "N"
Val.results.pest$IOC[which(Val.results.pest$Chemical == "2,4-Dinitrophenol")] <- "A"
Val.results.pest$CAS <- rep(Val.results.pest_CAS, 4)


diff.data_pfas <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,2)], mean)
diff.data_pfas[, c(3:6)] <- Val.results.pfas[c(1:14), c(3:6)]
diff.data_pfas$Obs.10uM <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,7)], mean)$Obs.10uM
diff.data_pfas[, c(8:11)] <- Val.results.pfas[c(1:14), c(8:11)]
Val.results.pfas <- diff.data_pfas[, c(1:6)]
Val.results.pfas[c(15:28),] <- diff.data_pfas[, c(1,7:11)]
colnames(Val.results.pfas)[2] <- "Obs"
Val.results.pfas$ref <- rep(c("PFAS (100 uM)", "PFAS (10 uM)"), each = 14)
Val.results.pfas$IOC <- "A"
Val.results.pfas$CAS <- read.csv("Val.results.pfas.forplot.csv")$CAS[1:28]


Val.results.tan$ref <- "Tanneberger et al.(Fish RTgill)"
Val.results.sch$ref <- rep(c("Schug et al.(Fish RTgut - High dose)", "Schug et al.(Fish RTgut - Low dose)"), each = 9)

comb.data <- rbind(Val.results.pest, Val.results.hus, Val.results.tan, Val.results.sch, Val.results.pfas, Val.results.Blanchette)
comb.data$MassBalanceIssue <- 0
comb.data <- rbind(comb.data, Val.results.unilever)

comb.data$'Cell type' <- c(rep("Cell-free", 80), rep("Human", 16), rep("Fish", 45),rep("Cell-free", 28), rep("Cell-free", 90), rep("Cell-free", 10))
comb.data$'Data source' <- c(rep("Our data", 80), rep("Literature", 16), rep("Literature", 45),rep("Our data", 28), rep("Literature", 90), rep("Our data", 10))
comb.data$conc <- c(rep(10,20),rep(1,20),rep(10,20),rep(1,20), 
                    c(3680,863,89.4,316,411,263,1100,215,2440,439,2.12,49.2,17.7,316,184,242,
                      0.835586139,0.049152126,14037.07518,98787.4465,39.34624697,0.075103267,4.1281912,712.8514056,2.883625129,310.5310531,93.3125972,23.00025556,4749.731472,9.110787172,1913.591359,10.12658228,14359.72461,33.3562586, 1509.328358,47.83950617,2.399314482,219.1887676,78.2369146,182.9931973,17656.06596,195.4161641,1000,
                      6.090133983,165.3439153,11.34429949,5.200208008,0.50942435,5.042864347,0.449640288,1.057529611,0.03869969,1827.040195,5291.005291,680.6579694,104.0041602,61.13092206,75.6429652	,26.97841727,21.15059222,7.73993808), rep(100,14),rep(10,14),rep(10,10),rep(10,90))
che.prop.list <- read.csv("che.prop.list.csv")
plot.data <- data.frame(Chemical = rep(comb.data$Chemical, 4),
                        CAS = rep(comb.data$CAS, 4),
                        Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 269),
                        Obs.val = rep(comb.data$Obs, 4),
                        pre.val = c(comb.data$X.free_Fischer,
                                    comb.data$X.free_VIVD,
                                    comb.data$X.free_IVMBM,
                                    comb.data$X.free_VCBA),
                        ref = rep(comb.data$ref, 4),
                        `Cell type` = rep(comb.data$`Cell type`, 4),
                        `Data source` = rep(comb.data$`Data source`, 4),
                        MassBalanceIssue = rep(comb.data$MassBalanceIssue, 4),
                        IOC = rep(comb.data$IOC, 4),
                        conc = rep(comb.data$conc, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}

#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Cell.type <- factor(plot.data$Cell.type, levels = c("Human", "Fish", "Cell-free"))
plot.data$Data.source <- factor(plot.data$Data.source, levels = c("Our data", "Literature"))
plot.data$MassBalanceIssue <- factor(plot.data$MassBalanceIssue)
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))

plot.data <- plot.data[-which(plot.data$IOC == "N"),]
plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3B1 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(color = MassBalanceIssue,shape = Data.source), alpha=0.8) +
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+8/5, y = 1-8/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 1-8/5, y = -7+8/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 1-8/5, y = -7+8/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 1-8/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental free percentage)") +
  ylab("log10(Predicted free percentage)")+
  ggtitle("B1. Fraction of free chemical concentration in media")+# (all datapoints)
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 1)+
  ylim(-7, 1)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("black", "black"))+ #"#FB684B", "#058ED9"
  scale_shape_manual(values = c(21, 16))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))

#prepare data (SOT poster data)
Val.results <- read.csv("Val.results.forplot_medconc.csv")

plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 27),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$IOC == "N"),]
plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)"),]


corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))


MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3B2 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species),color = "black", alpha=0.8) +#, shape=Cell.type
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -5+5.5/5, y = 0.5-5.5/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0.5-5.5/5, y = -5+5.5/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0.5-5.5/5, y = -5, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in media)") +
  ylab("log10(Predicted chemical amount in media)")+
  ggtitle("B2. Amount of chemicals in media")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-5, 0.5)+
  ylim(-5, 0.5)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))





Val.results <- read.csv("Val.results.forplot_cellconc.csv")
plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        CAS = rep(Val.results$CAS, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 43),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        amount = rep(Val.results$X, 4),
                        conc = rep(Val.results$Conc, 4),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4),
                        IOC = rep(Val.results$IOC, 4))
plot.data$H37 <- NA
plot.data$MW <- NA
plot.data$MP <- NA
plot.data$log.KOW <- NA
plot.data$log.KAW <- NA
plot.data$Solubility <- NA
plot.data$pKa <- NA
for(i in 1:length(plot.data$Chemical)){
  
  plot.data$H37[i] <- che.prop.list$H37[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MW[i] <- che.prop.list$MW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$MP[i] <- che.prop.list$MP[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KOW[i] <- che.prop.list$log.KOW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$log.KAW[i] <- che.prop.list$log.KAW[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$Solubility[i] <- che.prop.list$Solubility[which(che.prop.list$CAS == plot.data$CAS[i])]
  plot.data$pKa[i] <- che.prop.list$pKa[which(che.prop.list$CAS == plot.data$CAS[i])]
  
}
#plot.data <- plot.data[-which(plot.data$pKa < 1 | plot.data$pKa > 14),]
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
#plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
plot.data$IOC.lab <- "Neutral"
plot.data$IOC.lab[which(plot.data$IOC != "N")] <- "Non-Neutral"
plot.data$IOC.lab <- factor(plot.data$IOC.lab, levels = c("Non-Neutral", "Neutral"))

plot.data <- plot.data[-which(plot.data$IOC == "N"),]
plot.data <- plot.data[-which(plot.data$Model == "Zaldivar-Comenges et al. (2017)"),]

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

num <- ddply(plot.data, .(Model), summarize, num.point = length(diff), num.che = length(unique(CAS)))

Fig_S3B3 <-
  ggplot(plot.data, aes(x = logit.Obs.val, y = logit.pre.val))+
  geom_point(aes(shape = Species), color = "black", alpha=0.8) + #aes(color = Species, shape=Cell.type)
  scale_shape_manual(values = c(15, 22))+
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/4, label = paste0("rho ==", round(spearman, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = corr,
    mapping = aes(x = -7+6/5, y = 0-6/8, label = paste0("r ==", round(pearson, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  
  geom_text(
    data    = ME,
    mapping = aes(x = 0-6/5, y = -7+6/8, label = paste0("ME ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = MAD,
    mapping = aes(x = 0-6/5, y = -7+6/4, label = paste0("MAE ==", round(mean, 2))), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  geom_text(
    data    = num,
    mapping = aes(x = 0-6/5, y = -7, label = paste0("n ==", num.point,"/", num.che)), parse = TRUE,
    color = "black", size = 3.5
  )+
  
  xlab("log10(Experimental chemical amount in cells)") +
  ylab("log10(Predicted chemical amount in cells)")+
  ggtitle("B3. Amount of chemicals in cells")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(~Model)+
  xlim(-7, 0)+
  ylim(-7, 0)+
  #scale_x_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(lim = c(9E-7, 1),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))



jpeg(width = 17, height = 8.6, "Figure.S3.jpeg", units = "in", res=300)
Fig_S3A1 + Fig_S3B1 + Fig_S3A2 + Fig_S3B2 + Fig_S3A3 + Fig_S3B3 + plot_layout(height = c(1/3,1/3,1/3), widths = c(4/7, 3/7))
dev.off()

pdf(width = 17, height = 8.6, "Figure.S3.pdf")
Fig_S3A1 + Fig_S3B1 + Fig_S3A2 + Fig_S3B2 + Fig_S3A3 + Fig_S3B3 + plot_layout(height = c(1/3,1/3,1/3), widths = c(4/7, 3/7))
dev.off()
#  ------------- Sensitivity analysis: Data for Figure 4; Figure S6-9 -------------
library(ggplot2)
library(scales)
library(patchwork)
library(car)
library(plyr)
library(compositions)
library(tidyr)


Fischer.results <- read.csv("Fischer_results.csv")
Fischer.results$Model <- "Fischer"
VIVD.results <- read.csv("VIVD_results.csv")
VIVD.results$Model <- "VIVD"
VCBA_results <- read.csv("vcba_results.csv")
VCBA_results$Model <- "VCBA"
IVMBM_results <- read.csv("IVMBM_results.csv")
IVMBM_results$Model <- "IVMBM"

All_results <- rbind(Fischer.results, VIVD.results, VCBA_results, IVMBM_results)
All_results <- All_results[,-1]

All_results$F_Q_cells <- NA
All_results$F_Q_cells[which(All_results$Well == 96)] <- All_results$Q_cells[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_cells[which(All_results$Well == 384)] <- All_results$Q_cells[which(All_results$Well == 384)]/4e-5

All_results$F_Q_media <- NA
All_results$F_Q_media[which(All_results$Well == 96)] <- All_results$Q_media[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_media[which(All_results$Well == 384)] <- All_results$Q_media[which(All_results$Well == 384)]/4e-5

All_results$F_Q_media_free <- NA
All_results$F_Q_media_free[which(All_results$Well == 96)] <- All_results$Q_media_free[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_media_free[which(All_results$Well == 384)] <- All_results$Q_media_free[which(All_results$Well == 384)]/4e-5

All_results$F_C_media_free <- NA
All_results$F_C_media_free[which(All_results$Well == 96)] <- (All_results$Q_media_free[which(All_results$Well == 96)]/150)/0.000001
All_results$F_C_media_free[which(All_results$Well == 384)] <- (All_results$Q_media_free[which(All_results$Well == 384)]/40)/0.000001

All_results$FBS[which(All_results$FBS == 2)] <- "2%FBS"
All_results$FBS[which(All_results$FBS == 20)] <- "20%FBS"
All_results$Well[which(All_results$Well == 96)] <- "96well"
All_results$Well[which(All_results$Well == 384)] <- "384well"
All_results$Group <- paste(All_results$Well, All_results$FBS, sep = "_")

All_results$Model <- factor(All_results$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"))
All_results$Match.lab <- paste(All_results$Chemical, All_results$Group, sep=".")
default_data <- All_results[which(All_results$Cell.type == "Default" & All_results$Model == "IVMBM"),]
default_data$Match.lab <- paste(default_data$Chemical, default_data$Group, sep=".")
All_results$default_F_Q_cells <- NA
All_results$default_F_Q_media_free <- NA
All_results$default_F_C_media_free <- NA
for (i in 1:24128) {
  
  All_results$default_F_Q_cells[i] <- default_data$F_Q_cells[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_Q_media_free[i] <- default_data$F_Q_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_C_media_free[i] <- default_data$F_C_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
}

All_results$Fold_F_Q_cells <- log10(All_results$F_Q_cells/All_results$default_F_Q_cells)
All_results$Fold_F_Q_media_free <- log10(All_results$F_Q_media_free/All_results$default_F_Q_media_free)
All_results$Fold_F_C_media_free <- log10(All_results$F_C_media_free/All_results$default_F_C_media_free)


#Eta2 cal (ALL Chemical using 3 model)
library(lsr)
All_results_3model <- All_results[-c(which(All_results$Model == "VCBA")),]
model.F_C_media_free <- aov(F_C_media_free ~ Chemical + Cell.type + FBS + Well + Model, data = All_results_3model)
Eta2.model.F_C_media_free <- data.frame(etaSquared(model.F_C_media_free))
Eta2.model.F_C_media_free["Other",] <- NA
Eta2.model.F_C_media_free["Other",1] <- 1-sum(Eta2.model.F_C_media_free$eta.sq[1:5])

model.F_Q_cells <- aov(F_Q_cells ~ Chemical + Cell.type + FBS + Well + Model, data = All_results_3model)
Eta2.model.F_Q_cells <- data.frame(etaSquared(model.F_Q_cells))
Eta2.model.F_Q_cells["Other",] <- NA
Eta2.model.F_Q_cells["Other",1] <- 1-sum(Eta2.model.F_Q_cells$eta.sq[1:5])

Eta2.data <- data.frame(Variable = c(rownames(Eta2.model.F_Q_cells), rownames(Eta2.model.F_C_media_free)),
                        Value = c(Eta2.model.F_Q_cells$eta.sq, Eta2.model.F_C_media_free$eta.sq)*100,
                        Endpoint = rep(c("Chemical fraction in cell", "Free chemical fraction in media"), each = 6))

Eta2.data$Variable <- factor(Eta2.data$Variable, levels = c("Chemical","Model","FBS","Cell.type","Well","Other"))
Eta2.data_3model <- Eta2.data

#Eta2 cal (Neutral chemicals using 4 model)
library(lsr)
All_results_4model <- All_results[c(which(All_results$IOC == "N")),]
model.F_C_media_free <- aov(F_C_media_free ~ Chemical + Cell.type + FBS + Well + Model, data = All_results_4model)
Eta2.model.F_C_media_free <- data.frame(etaSquared(model.F_C_media_free))
Eta2.model.F_C_media_free["Other",] <- NA
Eta2.model.F_C_media_free["Other",1] <- 1-sum(Eta2.model.F_C_media_free$eta.sq[1:5])

model.F_Q_cells <- aov(F_Q_cells ~ Chemical + Cell.type + FBS + Well + Model, data = All_results_4model)
Eta2.model.F_Q_cells <- data.frame(etaSquared(model.F_Q_cells))
Eta2.model.F_Q_cells["Other",] <- NA
Eta2.model.F_Q_cells["Other",1] <- 1-sum(Eta2.model.F_Q_cells$eta.sq[1:5])

Eta2.data <- data.frame(Variable = c(rownames(Eta2.model.F_Q_cells), rownames(Eta2.model.F_C_media_free)),
                        Value = c(Eta2.model.F_Q_cells$eta.sq, Eta2.model.F_C_media_free$eta.sq)*100,
                        Endpoint = rep(c("Chemical fraction in cell", "Free chemical fraction in media"), each = 6))

Eta2.data$Variable <- factor(Eta2.data$Variable, levels = c("Chemical","Model","FBS","Cell.type","Well","Other"))
Eta2.data_4model <- Eta2.data

Eta2.data_allresults <- rbind(Eta2.data_3model, Eta2.data_4model)
Eta2.data_allresults$che.dataset <- rep(c("All chemicals", "Neutral chemicals"), each=12)
write.csv(Eta2.data_allresults, "Figure 4_sensitivity pie data.csv")



All_results$Model <- factor(All_results$Model,levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                            labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges et al. (2017)"))
All_results <- All_results[-which(All_results$Model == "Zaldivar-Comenges et al. (2017)" &
                                    All_results$IOC != "N"),]
jpeg(width = 14, height = 16, "Figure S6.jpeg", units = "in", res=300)
pdf(width = 14, height = 16, "Figure S6.pdf")
ggplot(All_results) + 
  geom_boxplot(aes(x=F_C_media_free, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  ylab("Chemical")+
  xlab("Free chemical concentration in media / Nominal concentration") +#Free chemical amount in media / Total chemical amount
  facet_grid(~Model)+
  scale_x_log10(lim = c(1E-6, 1.1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        #legend.position = c(0.02, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=12))
dev.off()

jpeg(width = 14, height = 16, "Figure S7.jpeg", units = "in", res=300)
pdf(width = 14, height = 16, "Figure S7.pdf")
ggplot(All_results) + 
  geom_boxplot(aes(x=F_Q_cells, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  ylab("Chemical")+
  xlab("Chemical amount in cell / Total chemical amount") +
  facet_grid(~Model)+
  scale_x_log10(lim = c(1E-9, 1.1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        #legend.position = c(0.02, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=12))
dev.off()

jpeg(width = 14, height = 16, "Figure S8.jpeg", units = "in", res=300)
pdf(width = 14, height = 16, "Figure S8.pdf")
ggplot(All_results) + 
  geom_boxplot(aes(x=Fold_F_C_media_free, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  geom_vline(xintercept = 0, linetype = "dashed", color="darkgray", size = 0.6)+
  ylab("Chemical")+
  xlab("Log10(Fraction of Free chemical concentration in media / Fraction of Free chemical concentration in media (Default cell by IVMBM))") +
  facet_grid(~Model)+
  xlim(-7, 4)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        #legend.position = c(0.02, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=12))
dev.off()

jpeg(width = 14, height = 16, "Figure S9.jpeg", units = "in", res=300)
pdf(width = 14, height = 16, "Figure S9.pdf")
ggplot(All_results) + 
  geom_boxplot(aes(x=Fold_F_Q_cells, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  geom_vline(xintercept = 0, linetype = "dashed", color="darkgray", size = 0.6)+
  ylab("Chemical")+
  xlab("Log10(Fraction of chemical in cell / Fraction of chemical in Default cell (IVMBM))") +
  facet_grid(~Model)+
  xlim(-7, 4)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        #legend.position = c(0.02, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=12))
dev.off()




# ---------------- Correction factor: Figure S11 ------------------
library(httk)
che.prop.list$fub_plasma <- NA
for (i in 1:116) {
  che.prop.list$fub_plasma[i] <- try(parameterize_steadystate(chem.cas = che.prop.list$CAS[i])$Funbound.plasma)
}
che.prop.list$fub_plasma <- as.numeric(che.prop.list$fub_plasma)



All_results <- rbind(Fischer.results, VIVD.results, VCBA_results, IVMBM_results)
All_results <- All_results[,-1]

All_results$F_Q_cells <- NA
All_results$F_Q_cells[which(All_results$Well == 96)] <- All_results$Q_cells[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_cells[which(All_results$Well == 384)] <- All_results$Q_cells[which(All_results$Well == 384)]/4e-5

All_results$F_Q_media <- NA
All_results$F_Q_media[which(All_results$Well == 96)] <- All_results$Q_media[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_media[which(All_results$Well == 384)] <- All_results$Q_media[which(All_results$Well == 384)]/4e-5

All_results$F_Q_media_free <- NA
All_results$F_Q_media_free[which(All_results$Well == 96)] <- All_results$Q_media_free[which(All_results$Well == 96)]/1.5e-4
All_results$F_Q_media_free[which(All_results$Well == 384)] <- All_results$Q_media_free[which(All_results$Well == 384)]/4e-5

All_results$F_C_media_free <- NA
All_results$F_C_media_free[which(All_results$Well == 96)] <- (All_results$Q_media_free[which(All_results$Well == 96)]/150)/0.000001
All_results$F_C_media_free[which(All_results$Well == 384)] <- (All_results$Q_media_free[which(All_results$Well == 384)]/40)/0.000001

All_results$FBS[which(All_results$FBS == 2)] <- "2%FBS"
All_results$FBS[which(All_results$FBS == 20)] <- "20%FBS"
All_results$Well[which(All_results$Well == 96)] <- "96well"
All_results$Well[which(All_results$Well == 384)] <- "384well"
All_results$Group <- paste(All_results$Well, All_results$FBS, sep = "_")

All_results$Model <- factor(All_results$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"))
All_results$Match.lab <- paste(All_results$Chemical, All_results$Group, sep=".")
default_data <- All_results[which(All_results$Cell.type == "Default" & All_results$Model == "IVMBM"),]
default_data$Match.lab <- paste(default_data$Chemical, default_data$Group, sep=".")
All_results$default_F_Q_cells <- NA
All_results$default_F_Q_media_free <- NA
All_results$default_F_C_media_free <- NA
for (i in 1:24128) {
  
  All_results$default_F_Q_cells[i] <- default_data$F_Q_cells[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_Q_media_free[i] <- default_data$F_Q_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_C_media_free[i] <- default_data$F_C_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
}

All_results$Fold_F_Q_cells <- log10(All_results$F_Q_cells/All_results$default_F_Q_cells)
All_results$Fold_F_Q_media_free <- log10(All_results$F_Q_media_free/All_results$default_F_Q_media_free)
All_results$Fold_F_C_media_free <- log10(All_results$F_C_media_free/All_results$default_F_C_media_free)

All_results$fub_plasma <- rep(che.prop.list$fub_plasma, 13*4*4)

Correction_fub_data <- All_results[-which(is.na(All_results$fub_plasma)),]
Correction_fub_data <- Correction_fub_data[which(Correction_fub_data$Cell.type == "Default" | 
                                                   Correction_fub_data$Cell.type == "MCF7" |
                                                   Correction_fub_data$Cell.type == "HCT116" |
                                                   Correction_fub_data$Cell.type == "HEK293T" |
                                                   Correction_fub_data$Cell.type == "HepG2" |
                                                   Correction_fub_data$Cell.type == "Me_180" |
                                                   Correction_fub_data$Cell.type == "SH-SY5Y"  |
                                                   Correction_fub_data$Cell.type == "HEK293H"),]
Correction_fub_data <- Correction_fub_data[-which(Correction_fub_data$Model == "VCBA" & 
                                                    Correction_fub_data$IOC != "N"),]
Correction_fub_data$F_correction <- Correction_fub_data$F_free_media / Correction_fub_data$fub_plasma
Correction_fub_data$Model <- factor(Correction_fub_data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                                    labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges \n et al. (2017)"))
Correction_fub_data$Chemical.dataset <- "Sensitivity analysis (n=68)"

F_correction_15che <-read.csv("F_correction_15che.csv")
F_correction_15che$IOC <- "N"
F_correction_15che <- pivot_longer(data = F_correction_15che,
                                   cols = c(2:5),
                                   names_to = "Model",
                                   values_to = "F_correction")
F_correction_15che$Model <- factor(F_correction_15che$Model, levels = c("Fischer", "IVBMB", "VIVD", "VCBA"),
                                   labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldivar-Comenges \n et al. (2017)"))
F_correction_15che$Chemical.dataset <- "Model application (n=15)"

F_correction_hisplot <- rbind(Correction_fub_data[, c("F_correction", "Model", "Chemical.dataset", "IOC")],
                              F_correction_15che[, c("F_correction", "Model", "Chemical.dataset", "IOC")])
library(scales)
ggplot(F_correction_hisplot, aes(F_correction, fill = IOC)) +
  geom_histogram(color="black",binwidth = log10(3.162278), position = "stack")+ #binwidth = 500
  scale_fill_manual(values=c("#F04248", "#2E72B2", "#6CB654"))+
  xlab("Ratio of fub_media to fub_plasma") +
  ggh4x::facet_grid2(Model~Chemical.dataset, scales = "free_y", independent = "y")+
  scale_x_log10(lim = c(1E-1, 1e+4),
                breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_rect(size = 1.1),
        #legend.position = c(0.02, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", face = "bold", size=10))
ggsave("Figure S11.jpeg", width = 8.5, height = 7)
ggsave("Figure S11.pdf", width = 8.5, height = 7)


