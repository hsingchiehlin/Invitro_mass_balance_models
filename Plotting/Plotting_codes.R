library(ggplot2)
library(scales)
library(patchwork)
library(car)
library(plyr)
library(compositions)
library(ggh4x)
#library(ggpubr)

# ---------- Figure S1: MAE&ME/ Figure S3 -------------------
#prepare data (Arimitage + Alan data)
Val.results.pest_single <- read.csv("Val.results.pest.forplot.csv")
Val.results.pest_mixture <- read.csv("Val.results.pest.forplot_mixture.csv")
Val.results.hus <- read.csv("Val.results.hus.forplot.csv")
Val.results.sch <- read.csv("Val.results.sch.forplot.csv")
Val.results.tan <- read.csv("Val.results.tan.forplot.csv")
Val.results.pfas <- read.csv("Val.results.pfas.forplot.csv")
Val.results.unilever <- read.csv("Val.results.unilever.forplot.csv")

diff.data_single <- aggregate(.~Chemical, data = Val.results.pest_single, mean)
Val.results.pest_single <- diff.data_single[, c(1:6)]
Val.results.pest_single[c(21:40),] <- diff.data_single[, c(1,7:11)]

diff.data_mixture <- aggregate(.~Chemical, data = Val.results.pest_mixture, mean)
Val.results.pest_mixture <- diff.data_mixture[, c(1:6)]
Val.results.pest_mixture[c(21:40),] <- diff.data_mixture[, c(1,7:11)]

diff.data_pfas <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,2)], mean)
diff.data_pfas[, c(3:6)] <- Val.results.pfas[c(1:14), c(3:6)]
diff.data_pfas$Obs.10uM <- aggregate(.~Chemical, data = Val.results.pfas[,c(1,7)], mean)$Obs.10uM
diff.data_pfas[, c(8:11)] <- Val.results.pfas[c(1:14), c(8:11)]
Val.results.pfas <- diff.data_pfas[, c(1:6)]
Val.results.pfas[c(15:28),] <- diff.data_pfas[, c(1,7:11)]
colnames(Val.results.pfas)[2] <- "Obs"
Val.results.pfas$ref <- rep(c("PFAS (100 uM)", "PFAS (10 uM)"), each = 14)

Val.results.pest <- rbind(Val.results.pest_single, Val.results.pest_mixture)
colnames(Val.results.pest)[2] <- "Obs"
Val.results.pest$ref <- rep(c("Valdiviezo et al.(10 uM - Single)", "Valdiviezo et al.(1 uM - Single)",
                              "Valdiviezo et al.(10 uM- Mixture)", "Valdiviezo et al.(1 uM- Mixture)"), each = 20)

Val.results.tan$ref <- "Tanneberger et al.(Fish RTgill)"
Val.results.sch$ref <- rep(c("Schug et al.(Fish RTgut - High dose)", "Schug et al.(Fish RTgut - Low dose)"), each = 9)

comb.data <- rbind(Val.results.pest, Val.results.hus, Val.results.tan, Val.results.sch, Val.results.pfas)
comb.data$MassBalanceIssue <- 0
comb.data <- rbind(comb.data, Val.results.unilever)

comb.data$'Cell type' <- c(rep("Cell-free", 80), rep("Human", 16), rep("Fish", 45),rep("Cell-free", 28), rep("Cell-free", 90))
comb.data$'Data source' <- c(rep("Our data", 80), rep("Literature", 16), rep("Literature", 45),rep("Our data", 28), rep("Literature", 90))


plot.data <- data.frame(Chemical = rep(comb.data$Chemical, 4),
                        Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 259),
                        Obs.val = rep(comb.data$Obs, 4),
                        pre.val = c(comb.data$X.free_Fischer,
                                    comb.data$X.free_VIVD,
                                    comb.data$X.free_IVMBM,
                                    comb.data$X.free_VCBA),
                        ref = rep(comb.data$ref, 4),
                        `Cell type` = rep(comb.data$`Cell type`, 4),
                        `Data source` = rep(comb.data$`Data source`, 4),
                        MassBalanceIssue = rep(comb.data$MassBalanceIssue, 4))


plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Cell.type <- factor(plot.data$Cell.type, levels = c("Human", "Fish", "Cell-free"))
plot.data$Data.source <- factor(plot.data$Data.source, levels = c("Our data", "Literature"))
plot.data$MassBalanceIssue <- factor(plot.data$MassBalanceIssue)

plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))

Fig_S1A <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(position="identity", alpha=0.5, color = "#0084EB", fill = "#0084EB") + 
  geom_text(
    data    = MAD,
    mapping = aes(x = -4, y = 50, label = paste("MAE = ", round(mean, 2), sep = "")),
    color = "black"
  )+
  geom_text(
    data    = ME,
    mapping = aes(x = -4, y = 40, label = paste("ME = ", round(mean, 2), sep = "")),
    color = "black"
  )+
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



Fig_S3A <-
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


Fig_S3B <-
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

jpeg(width = 10, height = 5.8, "Figure S3.jpeg", units = "in", res=300)
pdf(width = 10, height = 5.8, "Figure S3.pdf")
Fig_S3A + Fig_S3B + plot_layout(height = c(1/2,1/2))
dev.off()

#prepare data (SOT poster data)
Val.results <- read.csv("Val.results.forplot_medconc.csv")

plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 27),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4))

plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zald?var Comenges et al. (2017)"))
corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))


MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))


Fig_S1B <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(position="identity", alpha=0.5, color = "#0084EB", fill = "#0084EB") + 
  geom_text(
    data    = MAD,
    mapping = aes(x = -4, y = 5, label = paste("MAE = ", round(mean, 2), sep = "")),
    color = "black"
  )+
  geom_text(
    data    = ME,
    mapping = aes(x = -4, y = 4.2, label = paste("ME = ", round(mean, 2), sep = "")),
    color = "black"
  )+
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


Val.results <- read.csv("Val.results.forplot_cellconc.csv")
plot.data <- data.frame(Chemical = rep(Val.results$Chemical, 4),
                        Model = rep(c("Fischer","IVMBM","VCBA","VIVD"), each = 43),
                        Obs.val = rep(Val.results$Obs, 4),
                        pre.val = c(Val.results$Fischer,
                                    Val.results$IVMBM,
                                    Val.results$VCBA,
                                    Val.results$VIVD),
                        Cell.type = rep(Val.results$Cell, 4),
                        Species = rep(Val.results$Species, 4))
plot.data$logit.Obs.val <- log10(plot.data$Obs.val)
plot.data$logit.pre.val <- log10(plot.data$pre.val)
plot.data$diff <- plot.data$logit.pre.val-plot.data$logit.Obs.val
plot.data$abs.diff <- abs(plot.data$diff)
plot.data <- plot.data[-c(which(is.na(plot.data$diff))),]
plot.data$Model <- factor(plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                          labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zald?var Comenges et al. (2017)"))

corr<- ddply(plot.data, .(Model), summarize, pearson = cor(x=logit.Obs.val, y=logit.pre.val),
             spearman = cor(x=Obs.val, y=pre.val, method = "spearman"))

MAD <- ddply(plot.data, .(Model), summarize, mean = mean(abs.diff))

ME <- ddply(plot.data, .(Model), summarize, mean = mean(diff))


Fig_S1C <-
  ggplot(aes(diff), data = plot.data) + 
  geom_histogram(position="identity", alpha=0.5, color = "#0084EB", fill = "#0084EB") + 
  geom_text(
    data    = MAD,
    mapping = aes(x = -4, y = 8, label = paste("MAE = ", round(mean, 2), sep = "")),
    color = "black"
  )+
  geom_text(
    data    = ME,
    mapping = aes(x = -4, y = 6.8, label = paste("ME = ", round(mean, 2), sep = "")),
    color = "black"
  )+
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


pdf(width = 10, height = 8, "Figure S1.pdf")
Fig_S1A + Fig_S1B + Fig_S1C + plot_layout(height = c(1/3,1/3,1/3))
dev.off()

jpeg(width = 10, height = 8, "Figure.S1.jpeg", units = "in", res=300)
Fig_S1A + Fig_S1B + Fig_S1C + plot_layout(height = c(1/3,1/3,1/3))
dev.off()


# ------------ Figure 2 & Figure S2: Chemical prop. ------------
che_prop <- read.csv("che_prop.csv")
plot1 <- ggplot(che_prop, aes(x = MW, y = log_KOW, color = Dataset)) + 
  geom_point(aes(color = Dataset), size = 3) + 
  #geom_point(shape = 1, color = "black", size = 3) + 
  scale_color_manual(values = c("#3B7BE3", "#EE5044"))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1),
        legend.position = c(0.8, 0.2),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

dens1 <- ggplot(che_prop, aes(x = MW, fill = Dataset)) + 
  geom_density(alpha = 0.4,size = 0.4) + 
  scale_fill_manual(values = c("#3B7BE3", "#EE5044"))+
  theme_void() + 
  theme(legend.position = "none")

dens2 <- ggplot(che_prop, aes(x = log_KOW, fill = Dataset)) + 
  geom_density(alpha = 0.4,size = 0.4) + 
  scale_fill_manual(values = c("#3B7BE3", "#EE5044"))+
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()

dens1 + plot_spacer() + plot1 + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
ggsave("Figure 2.jpeg", width = 7.79, height = 6.92)
ggsave("Figure 2.pdf", width = 7.79, height = 6.92)


che_prop$log_Solubility <- log10(che_prop$Solubility)
library(GGally)
ggpairs(che_prop,                 # Data frame
        columns = c(4:7,10),
        aes(color = Dataset,  # Color by group (cat. variable)
            alpha = 0.5),
        upper = list(continuous = "blank"))+     # Transparency
  scale_colour_manual(values = c("#3B7BE3", "#EE5044"))+
  scale_fill_manual(values = c("#3B7BE3", "#EE5044"))
ggsave("Figure S2.jpeg", width = 11.1, height = 9.39)
ggsave("Figure S2.pdf", width = 11.1, height = 9.39)

# ---------- Figure S4-S5:error trend -------------
Error_med_prec.results <- read.csv("Error_med_prec.csv")

Error_med_prec.results$Diff_Fischer <- log10(Error_med_prec.results$Pred_perc_media_Fischer) - log10(Error_med_prec.results$Obs_perc_media)
Error_med_prec.results$Diff_VIVD <- log10(Error_med_prec.results$Pred_perc_media_VIVD) - log10(Error_med_prec.results$Obs_perc_media)
Error_med_prec.results$Diff_IVMBM <- log10(Error_med_prec.results$Pred_perc_media_IVMBM) - log10(Error_med_prec.results$Obs_perc_media)
Error_med_prec.results$Diff_VCBA <- log10(Error_med_prec.results$Pred_perc_media_VCBA) - log10(Error_med_prec.results$Obs_perc_media)

Error_med_prec.results$ABS_Diff_Fischer <- abs(Error_med_prec.results$Diff_Fischer)
Error_med_prec.results$ABS_Diff_VIVD <- abs(Error_med_prec.results$Diff_VIVD)
Error_med_prec.results$ABS_Diff_IVMBM <- abs(Error_med_prec.results$Diff_IVMBM)
Error_med_prec.results$ABS_Diff_VCBA <- abs(Error_med_prec.results$Diff_VCBA)


Corr.plot.data <- data.frame(Chemical = rep(Error_med_prec.results$Chemical, 4),
                             Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 286),
                             Diff = stack(Error_med_prec.results[,19:22])$values,
                             Type = rep(Error_med_prec.results$Type, 4))

Corr.plot.data <- data.frame(Chemical = rep(Corr.plot.data$Chemical, 7),
                             Model = rep(Corr.plot.data$Model, 7),
                             Diff = rep(Corr.plot.data$Diff, 7),
                             Type = rep(Corr.plot.data$Type, 7),
                             Che.prop.value = c(rep(Error_med_prec.results$log10_H37, 4),
                                                rep(Error_med_prec.results$MP, 4),
                                                rep(Error_med_prec.results$logkow, 4),
                                                rep(Error_med_prec.results$logkaw, 4),
                                                rep(Error_med_prec.results$log10_Solubility, 4),
                                                rep(Error_med_prec.results$MW, 4),
                                                rep(log10(Error_med_prec.results$Conc), 4)),
                             Che.prop.name = rep(c("log10_H37", "MP", "logkow", "logkaw", "log10_Solubility", "MW", "log10_Conc"), each = 1144))

Corr.plot.data <- Corr.plot.data[-which(Corr.plot.data$Chemical == "cyclosporine A"),]
Corr.plot.data$Model <- factor(Corr.plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                               labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))

jpeg(width = 15, height = 9.1, "Figure S4.jpeg", units = "in", res=300)
pdf(width = 15, height = 9.1, "Figure S4.pdf")
ggplot(Corr.plot.data)+
  geom_smooth(aes(x = Che.prop.value, y = Diff), method = "lm", colour = "#7BD389", fill = "#7BD38920") + 
  geom_point(aes(x = Che.prop.value, y = Diff, shape = Type), color = "#004777") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3)+
  xlab("Chemical property value") +
  ylab("Log10(Predicted % of chemical in media / Measured % of chemical in media")+
  facet_grid(Model~Che.prop.name, scales = "free")+
  scale_shape_manual(values = c(16, 21))+
  ylim(-6, 3)+
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
write.csv(corr.results.diff_media, "corr.results.diff_media_1123.csv")



Error_cell_prec.results <- read.csv("Error_cell_prec.csv")

Corr.plot.data <- data.frame(Chemical = rep(Error_cell_prec.results$Chemical, 4),
                             Model = rep(c("Fischer", "VIVD", "IVMBM", "VCBA"), each = 43),
                             Diff = stack(Error_cell_prec.results[,7:10])$values)

Corr.plot.data <- data.frame(Chemical = rep(Corr.plot.data$Chemical, 7),
                             Model = rep(Corr.plot.data$Model, 7),
                             Diff = rep(Corr.plot.data$Diff, 7),
                             Che.prop.value = c(rep(Error_cell_prec.results$log10_H37, 4),
                                                rep(Error_cell_prec.results$MP, 4),
                                                rep(Error_cell_prec.results$logkow, 4),
                                                rep(Error_cell_prec.results$logkaw, 4),
                                                rep(Error_cell_prec.results$log10_Solubility, 4),
                                                rep(Error_cell_prec.results$MW, 4),
                                                rep(log10(Error_cell_prec.results$Conc), 4)),
                             Che.prop.name = rep(c("log10_H37", "MP", "logkow", "logkaw", "log10_Solubility", "MW", "log10_Conc"), each = 172))

Corr.plot.data <- Corr.plot.data[-which(Corr.plot.data$Chemical == "cyclosporine A"),]
Corr.plot.data$Model <- factor(Corr.plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                               labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zald?var Comenges et al. (2017)"))

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

corr.results.diff_cell <- ddply(Corr.plot.data,.(Che.prop.name, Model),lm_eqn)
colnames(corr.results.diff_cell)[3:11] <- c("Intercept", "Intercept.se", "Intercept.pval", 
                                            "slope", "slope.se", "slope.pval", 
                                            "r.squared", "adj.r.squared", "r.squared.pval")
write.csv(corr.results.diff_cell, "corr.results.diff_cell.csv")




# ---------- Figure S6-10: Sensitivity analysis -------------
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
All_results$Conc <- "1 uM"

All_results$F_Q_cells <- NA
All_results$F_Q_cells[which(All_results$Well == 96 & All_results$Conc == "1 uM")] <- All_results$Q_cells[which(All_results$Well == 96 & All_results$Conc == "1 uM")]/1.5e-4
All_results$F_Q_cells[which(All_results$Well == 384 & All_results$Conc == "1 uM")] <- All_results$Q_cells[which(All_results$Well == 384 & All_results$Conc == "1 uM")]/4e-5

All_results$F_Q_media <- NA
All_results$F_Q_media[which(All_results$Well == 96 & All_results$Conc == "1 uM")] <- All_results$Q_media[which(All_results$Well == 96 & All_results$Conc == "1 uM")]/1.5e-4
All_results$F_Q_media[which(All_results$Well == 384 & All_results$Conc == "1 uM")] <- All_results$Q_media[which(All_results$Well == 384 & All_results$Conc == "1 uM")]/4e-5

All_results$F_Q_media_free <- NA
All_results$F_Q_media_free[which(All_results$Well == 96 & All_results$Conc == "1 uM")] <- All_results$Q_media_free[which(All_results$Well == 96 & All_results$Conc == "1 uM")]/1.5e-4
All_results$F_Q_media_free[which(All_results$Well == 384 & All_results$Conc == "1 uM")] <- All_results$Q_media_free[which(All_results$Well == 384 & All_results$Conc == "1 uM")]/4e-5

All_results$F_C_media_free <- NA
All_results$F_C_media_free[which(All_results$Well == 96 & All_results$Conc == "1 uM")] <- (All_results$Q_media_free[which(All_results$Well == 96 & All_results$Conc == "1 uM")]/150)/0.000001
All_results$F_C_media_free[which(All_results$Well == 384 & All_results$Conc == "1 uM")] <- (All_results$Q_media_free[which(All_results$Well == 384 & All_results$Conc == "1 uM")]/40)/0.000001

All_results$FBS[which(All_results$FBS == 2)] <- "2%FBS"
All_results$FBS[which(All_results$FBS == 20)] <- "20%FBS"
All_results$Well[which(All_results$Well == 96)] <- "96well"
All_results$Well[which(All_results$Well == 384)] <- "384well"
All_results$Group <- paste(All_results$Well, All_results$FBS, sep = "_")

All_results$Match.lab <- paste(All_results$Chemical, All_results$Group, All_results$Conc,  sep=".")
default_data <- All_results[which(All_results$Cell.type == "Default" & All_results$Model == "IVMBM"),]
default_data$Match.lab <- paste(default_data$Chemical, default_data$Group, default_data$Conc, sep=".")
All_results$default_F_Q_cells <- NA
All_results$default_F_Q_media_free <- NA
All_results$default_F_C_media_free <- NA
for (i in 1:22256) {
  
  All_results$default_F_Q_cells[i] <- default_data$F_Q_cells[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_Q_media_free[i] <- default_data$F_Q_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
  All_results$default_F_C_media_free[i] <- default_data$F_C_media_free[which(default_data$Match.lab == All_results$Match.lab[i])]
}

All_results$Fold_F_Q_cells <- log10(All_results$F_Q_cells/All_results$default_F_Q_cells)
All_results$Fold_F_Q_media_free <- log10(All_results$F_Q_media_free/All_results$default_F_Q_media_free)
All_results$Fold_F_C_media_free <- log10(All_results$F_C_media_free/All_results$default_F_C_media_free)


All_results$Model <- factor(All_results$Model,levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                            labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
jpeg(width = 14, height = 16, "Figure S6.jpeg", units = "in", res=300)
pdf(width = 14, height = 16, "Figure S6.pdf")
ggplot(All_results[which(All_results$Conc == "1 uM"),]) + 
  geom_boxplot(aes(x=F_C_media_free, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  ylab("Chemical")+
  xlab("Free chemical concentration in media / Nominal concentration") +#Free chemical amount in media / Total chemical amount
  facet_grid(~Model)+
  scale_x_log10(lim = c(1E-7, 1.1),
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
ggplot(All_results[which(All_results$Conc == "1 uM"),]) + 
  geom_boxplot(aes(x=F_Q_cells, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  ylab("Chemical")+
  xlab("Chemical amount in cell / Total chemical amount") +
  facet_grid(~Model)+
  scale_x_log10(lim = c(1E-7, 1.1),
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
ggplot(All_results[which(All_results$Conc == "1 uM"),]) + 
  geom_boxplot(aes(x=Fold_F_C_media_free, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  geom_vline(xintercept = 0, linetype = "dashed", color="darkgray", size = 0.6)+
  ylab("Chemical")+
  xlab("Log10(Fraction of Free chemical concentration in media / Fraction of Free chemical concentration in media (Default cell by IVMBM))") +
  #xlab("Log10(Percentage of Free chemical in media / Percentage of Free chemical in media (Default cell by IVMBM))") +
  facet_grid(~Model)+
  xlim(-6, 4)+
  #scale_x_log10(lim = c(1E-5, 1e+5),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
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
ggplot(All_results[which(All_results$Conc == "1 uM"),]) + 
  geom_boxplot(aes(x=Fold_F_Q_cells, y=Chemical, color=Group, fill=Group), position=position_dodge(width=.7))+
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_fill_manual(values=c("#ABCBED90", "#276EB990", "#FDAFBC90", "#FB234790"))+
  geom_vline(xintercept = 0, linetype = "dashed", color="darkgray", size = 0.6)+
  ylab("Chemical")+
  xlab("Log10(Fraction of chemical in cell / Fraction of chemical in Default cell (IVMBM))") +
  facet_grid(~Model)+
  xlim(-6, 4)+
  #scale_x_log10(lim = c(1E-5, 1e+5),
  #breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #labels = trans_format("log10", math_format(10^.x))) +
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




Scatter.plot.data <- data.frame(Chemical = rep(All_results$Chemical, 2),
                                CAS = rep(All_results$CAS, 2),
                                Cell.type = rep(All_results$Cell.type, 2),
                                FBS = rep(All_results$FBS, 2),
                                Well = rep(All_results$Well, 2),
                                Model = rep(All_results$Model, 2),
                                Group = rep(All_results$Group, 2),
                                Conc = rep(All_results$Conc, 2),
                                Def.value = c(All_results$default_F_C_media_free, All_results$default_F_Q_cells),
                                Oth.value = c(All_results$F_C_media_free, All_results$F_Q_cells),
                                Label = rep(c("Fraction of free conc. in media","Chemical fraction in cell"), each=22256))
Scatter.plot.data$Model <- factor(Scatter.plot.data$Model, levels = c("Fischer", "IVMBM", "VIVD", "VCBA"),
                                  labels = c("Fischer et al. (2017)", "Armitage et al. (2021)", "Fisher et al. (2019)", "Zaldívar Comenges et al. (2017)"))
Scatter.plot.data$Label <- factor(Scatter.plot.data$Label, levels = c("Fraction of free conc. in media","Chemical fraction in cell"))

Fig_S10 <- 
  ggplot(Scatter.plot.data, aes(x = Def.value, y = Oth.value))+
  geom_point(aes(color = Group)) +
  xlab("Percentage of chemical in default system (using IVMBM)") +
  ylab("Percentage of chemical in other systems")+
  #ggtitle("A")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  facet_grid(Label~Model, switch = "y")+
  #facet_nested(~Label+Model, nest_line = element_line(linetype = 1), ncol=4)+
  scale_x_log10(lim = c(1E-7, 1.1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1E-7, 1.1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#ABCBED", "#276EB9", "#FDAFBC", "#FB2347"))+
  scale_shape_manual(values=c(16, 21))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=10, hjust=0.5, vjust=0.5), 
        panel.border = element_rect(size = 1.1),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size=10),
        #strip.switch.pad.grid = unit(1, "cm"),
        legend.position = c(0.12, 0.09),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        #axis.title.y = element_text(vjust = -15,color = "black", face = "bold"),
        axis.title = element_text(color = "black"))
jpeg(width = 10, height = 5, "Figure S10.jpeg", units = "in", res=300)
pdf(width = 10, height = 5, "Figure S1.pdf")
Fig_S10
dev.off()


# Eta2 calculation
library(lsr)
model.F_C_media_free <- aov(F_C_media_free ~ Chemical + Cell.type + FBS + Well + Model, data = All_results[which(All_results$Conc == "1 uM"),])
Eta2.model.F_C_media_free <- data.frame(etaSquared(model.F_C_media_free, type = 1))

model.F_Q_cells <- aov(F_Q_cells ~ Chemical + Cell.type + FBS + Well + Model, data = All_results[which(All_results$Conc == "1 uM"),])
Eta2.model.F_Q_cells <- data.frame(etaSquared(model.F_Q_cells, type = 1))

Eta2.data <- data.frame(Variable = c(rownames(Eta2.model.F_Q_cells), "Residuals", rownames(Eta2.model.F_C_media_free), "Residuals"),
                        Value = c(Eta2.model.F_Q_cells$eta.sq, 1-sum(Eta2.model.F_Q_cells$eta.sq),
                                  Eta2.model.F_C_media_free$eta.sq, 1-sum(Eta2.model.F_C_media_free$eta.sq))*100,
                        Endpoint = rep(c("Chemical fraction in cell", "Fraction of free conc. in media"), each = 6))
Eta2.data$Variable <- rep(c("Chemical (n=107)", "Cell type (n=13)","FBS (2% vs. 20%)", "Labware (96-well vs. 384-well)","Model (n=4)", "Residuals"), 2)
Eta2.data$Variable <- factor(Eta2.data$Variable, levels = c("Chemical (n=107)","Model (n=4)","FBS (2% vs. 20%)",
                                                            "Cell type (n=13)","Labware (96-well vs. 384-well)", "Residuals"))
Eta2.data$Endpoint <- factor(Eta2.data$Endpoint, levels = c("Fraction of free conc. in media","Chemical fraction in cell"))



