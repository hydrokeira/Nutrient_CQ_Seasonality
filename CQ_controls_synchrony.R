require(dplyr)
require(tidyr)
require(ggplot2)

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

CQ<-read.csv("CQ_raw_slopes.csv")
colnames(CQ)[3]<-"chemical"

monthly_CQ_wide<-CQ %>%
  dplyr::select(Stream_Name, chemical, slope) %>%
  pivot_wider(names_from = chemical, values_from = slope)

CQ_pairs <- CQ %>%
  inner_join(CQ,
             by = "Stream_Name",
             suffix = c("_1", "_2")) %>%
  filter(chemical_1 < chemical_2) %>%
  transmute(
    Stream_Name,
    solute_pair = paste(chemical_1, chemical_2, sep = "-"),
    Solute1 = chemical_1,
    Solute2 = chemical_2,
    slope1  = slope_1,
    slope2  = slope_2
  )

table(CQ_pairs$solute_pair)

results_all<-read.csv("Synchrony_Shift_AllSolutePairs.csv")

CQ_pairs<-left_join(CQ_pairs, results_all[2:13])

impacted_class<-read.csv("ImpactedClass_Synchrony.csv")

p1<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward shift", "forward shift"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(slope1, slope2))+
  geom_hline(yintercept = 0, col="grey")+
  geom_vline(xintercept = 0, col="grey")+
  geom_abline(slope = 1, col="#0a4e9b", lwd=1)+
  geom_abline(slope = -1, col="#0096ce", lwd=1)+
  geom_polygon(data = tri_top, aes(x, y), fill = "grey", alpha = 0.3) +
  geom_polygon(data = tri_bottom, aes(x, y), fill = "grey", alpha = 0.3) +
  geom_polygon(data = tri_right, aes(x, y), fill = "grey", alpha = 0.3) +
  geom_polygon(data = tri_left, aes(x, y), fill = "grey", alpha = 0.3) +
  geom_point(aes(col=impacted_class, shape=impacted_class), size=3, alpha=0.8)+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.spacing = unit(0.8, "lines"),
        strip.background = element_blank(),text = element_text(size = 20), legend.position = "bottom")+
  lims(x=c(-.7, .7), y=c(-.7, .7))+
  scale_color_manual(values = c("less-impacted"="#8DAA91", 
                                "agricultural"="#E6B422", "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, "urban"=17, "agricultural & urban"=15))+
  facet_grid(rows = vars(solute_pair), cols = vars(shift_direction), scales="free")+
  labs(x="CQ Slope Solute 1", y="CQ Slope Solute 2", col="", shape="")

p1

pdf("CQ_Shift_LandUse.pdf", width = 8, height = 6.5)

p1

dev.off()

p2<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                               "agricultural & urban"))) %>%
  ggplot(aes(impacted_class, residual, col=impacted_class, shape = impacted_class))+
  geom_hline(yintercept = 0, col="#0a4e9b", lwd=1)+
  geom_violin(fill=NA)+geom_jitter()+facet_wrap(~solute_pair, nrow=4)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "null", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("less-impacted"="#8DAA91", "agricultural"="#E6B422", 
                                "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, 
                                "urban"=17, "agricultural & urban"=15))+
  labs(x="", y="Distance to 1:1 Line (S2-S1/sqrt(2))")+ylim(-1,1)

p2

p3<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward shift", "forward shift"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(impacted_class, residual, col=impacted_class, shape = impacted_class))+
  geom_hline(yintercept = 0, col="#0a4e9b", lwd=1)+
  geom_violin(fill=NA)+geom_jitter()+facet_wrap(~shift_direction, nrow=4)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "null", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("less-impacted"="#8DAA91", "agricultural"="#E6B422", 
                                "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, 
                                "urban"=17, "agricultural & urban"=15))+
  labs(x="", y="Distance to 1:1 Line (S2-S1/sqrt(2))")

p3  

p4<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward shift", "forward shift"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(solute_pair, residual))+
  geom_hline(yintercept = 0)+
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.15, ymax = Inf,
           fill = "grey", alpha = 0.5) +
  # Bottom shaded region (y <= -0.15)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = -0.15,
           fill = "grey", alpha = 0.5) +
  geom_violin(fill=NA, scale = "width")+geom_jitter(aes(col=impacted_class, shape = impacted_class))+
  facet_wrap(~shift_direction, nrow=4)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "bottom", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("less-impacted"="#8DAA91", "agricultural"="#E6B422", 
                                "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, 
                                "urban"=17, "agricultural & urban"=15))+
  labs(x="", y="Distance to 1:1 Line (S2-S1/sqrt(2))", col="", shape="")

p4

p45<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward shift", "forward shift"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(impacted_class, residual))+
  geom_hline(yintercept = 0)+
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.15, ymax = Inf,
           fill = "grey", alpha = 0.5) +
  # Bottom shaded region (y <= -0.15)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = -0.15,
           fill = "grey", alpha = 0.5) +
  geom_violin(fill=NA, scale = "width")+geom_jitter(aes(col=solute_pair, shape = solute_pair))+
  facet_wrap(~shift_direction, nrow=4)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "bottom", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  # scale_color_manual(values = c("less-impacted"="#8DAA91", "agricultural"="#E6B422", 
  #                               "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  # scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, 
  #                               "urban"=17, "agricultural & urban"=15))+
  labs(x="", y="Distance to 1:1 Line (S2-S1/sqrt(2))", col="", shape="")

p45

p5<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward shift", "forward shift"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(impacted_class, residual))+
  geom_hline(yintercept = 0)+
  geom_violin(fill=NA, scale = "width")+geom_jitter()+
  facet_grid(cols = vars(shift_direction), rows = vars(solute_pair))+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "bottom", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="", y="Distance to 1:1 Line (S2-S1/sqrt(2))")

p5

p6<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  filter(!is.na(shift)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(shift = factor(shift, levels=c("non-shape similar", "0", "1", "2", "3", "4", "5", "6"))) %>%
  #filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    shift=="non-shape similar"~"non-shape similar",
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward lag", "forward lag", "non-shape similar"))) %>%
  ggplot(aes(shift, residual))+
  geom_hline(yintercept = 0, col="#0a4e9b", lwd=1)+
  geom_violin(fill=NA, scale = "width", aes(group = shift))+geom_jitter(aes(col=shift_direction))+
  facet_wrap(~solute_pair, nrow=1)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="", y="Residual off 1:1 Line", col="")+
  scale_color_manual(values = c("in-phase"="#0a4e9b","out-of-phase"="#0096ce","forward lag" = "#E7A1B0",
                               "backward lag" = "#F6C28B", "non-shape similar"="grey"))+
  ylim(-1,1)

p6

p7<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward lag", "forward lag"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(slope1, slope2, col=shift_direction))+
  geom_hline(yintercept = 0, col="grey")+
  geom_vline(xintercept = 0, col="grey")+
  geom_abline(slope = 1, col="#0a4e9b", lwd=1)+
  geom_abline(slope = -1, col="#0096ce", lwd=1)+
  geom_point(alpha=0.5)+
  facet_wrap(~solute_pair)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "right", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="CQ Slope Solute 1", y="CQ Slope Solute 2", col="")+
  scale_color_manual(values = c("in-phase"="#0a4e9b","out-of-phase"="#0096ce","forward lag" = "#E7A1B0",
                                "backward lag" = "#F6C28B"))

p7


p8<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  filter(!is.na(shift)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  #filter(sync_char=="synchronous") %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(shift = factor(shift, levels=c("non-shape similar", "0", "1", "2", "3", "4", "5", "6"))) %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward lag", "forward lag"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(slope1, slope2, col=shift))+
  geom_hline(yintercept = 0, col="grey")+
  geom_vline(xintercept = 0, col="grey")+
  geom_abline(slope = 1, col="#0a4e9b", lwd=1, alpha=0.5)+
  geom_abline(slope = -1, col="#0096ce", lwd=1, alpha=0.5)+
  geom_point()+
  facet_wrap(~solute_pair)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="CQ Slope Solute 1", y="CQ Slope Solute 2", col="Lag\n(Months)")+
  ylim(-1.5,1.5)+xlim(-1.5, 1.5)+scale_color_manual(values=c("0"="#0a4e9b", "4"="#FF8E81", "1"="#FFC8B4", "3"="#F0A597", "5"="#EB7773",
                                                              "2"="#FCCDC3", "6"="#0096ce", "non-shape similar"="grey"))

p8

pdf("DualCQ_Residual.pdf", width = 12, height = 8)

ggarrange(p8, p6, align = "v", nrow = 2)

dev.off()

p9<-CQ_pairs %>%
  left_join(impacted_class[c(2,5)]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase",
                                                            "backward lag", "forward lag"))) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("less-impacted", "urban", "agricultural",
                                                        "agricultural & urban"))) %>%
  ggplot(aes(slope1, slope2, col=impacted_class, shape = impacted_class))+
  geom_hline(yintercept = 0, col="grey")+
  geom_vline(xintercept = 0, col="grey")+
  geom_abline(slope = 1, col="#0a4e9b", lwd=1)+
  geom_abline(slope = -1, col="#0096ce", lwd=1)+
  geom_point(size=2, alpha=0.5)+
  facet_wrap(~solute_pair)+
  theme_classic()+theme(text = element_text(size = 20), legend.position = "bottom", 
                        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="CQ Slope Solute 1", y="CQ Slope Solute 2", col="", shape="")+
  scale_color_manual(values = c("less-impacted"="#8DAA91", "agricultural"="#E6B422", 
                                "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, 
                                "urban"=17, "agricultural & urban"=15))

p9

pdf("CQ_SolutePair.pdf", width = 11, height = 15)
ggarrange(p8, p7, p9, nrow = 3)
dev.off()

pdf("Residuals_Solutes.pdf", width = 4, height = 9.6)

p2

dev.off()

pdf("Residuals_SolutePairs_Gradient.pdf", width = 8, height = 12)

ggarrange(p6, p2, align = "h")

dev.off()

pdf("Residuals_Shifts.pdf", width = 4, height = 12)

p3

dev.off()

pdf("Residuals_Shifts_2x4.pdf", width = 8, height = 12)

ggarrange(p4, p3, ncol = 2, align = "h")

dev.off()

pdf("Residuals_Shifts_Solute_colLULC.pdf", width = 4, height = 12)

p4

dev.off()

pdf("Residuals_Shifts_LULC_colSolute.pdf", width = 4, height = 12)

p45

dev.off()

pdf("Residuals_Shifts_LULC_Solute_Grid.pdf", width = 12, height = 8)

p5

dev.off()

#### old code not using####

p1<-monthly_CQ_sync %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  filter(Solute1=="DSi"&Solute2=="N") %>%
  ggplot(aes(DSi, N))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=shift_direction, shape=shift_direction), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  scale_color_manual(values = c("in-phase"="black","out-of-phase"="black","forward shift"="forestgreen",
                                "backward shift"="darkorange"))+
  scale_shape_manual(values = c("in-phase"=19, "out-of-phase"=1,"forward shift"=19,
                                "backward shift"=19))+
  theme(text = element_text(size = 20), legend.position = "null")+
  labs(x="DSi CQ", y="N CQ")+lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))

p1

p2<-monthly_CQ_sync %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  filter(Solute1=="DSi"&Solute2=="P") %>%
  ggplot(aes(DSi, P))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=shift_direction, shape=shift_direction), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  scale_color_manual(values = c("in-phase"="black","out-of-phase"="black","forward shift"="forestgreen",
                                "backward shift"="darkorange"))+
  scale_shape_manual(values = c("in-phase"=19, "out-of-phase"=1,"forward shift"=19,
                                "backward shift"=19))+
  theme(text = element_text(size = 20), legend.position = "null")+
  labs(x="DSi CQ", y="P CQ")+lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))

p2

p3<-monthly_CQ_sync %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  filter(Solute1=="N"&Solute2=="P") %>%
  ggplot(aes(N, P))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=shift_direction, shape=shift_direction), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  scale_color_manual(values = c("in-phase"="black","out-of-phase"="black","forward shift"="forestgreen",
                                "backward shift"="darkorange"))+
  scale_shape_manual(values = c("in-phase"=19, "out-of-phase"=1,"forward shift"=19,
                                "backward shift"=19))+
  theme(text = element_text(size = 20))+
  labs(x="N CQ", y="P CQ", col="Shift Direction", shape="Shift Direction")+
  lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))

p3

pdf("CQ_ShiftDirection_Synchrony.pdf", width = 15, height = 4.5)

ggarrange(p1, p2, p3, nrow = 1, widths = c(0.5, 0.5, 0.7))

dev.off()

land_use<-read.csv("Impacted_Class.csv")

k1<-monthly_CQ_wide %>%
  left_join(land_use[,c(2,5)]) %>%
  ggplot(aes(DSi, N))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=impacted_class), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  labs(x="DSi CQ", y="N CQ")+lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))+
  scale_color_manual(values = c("non-impacted"="#8DAA91", "ag"="#E6B422", "urban"="#4682B4", "ag & urban"="#A67B5B"))

k2<-monthly_CQ_wide %>%
  left_join(land_use[,c(2,5)]) %>%
  ggplot(aes(DSi, P))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=impacted_class), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null")+
  labs(x="DSi CQ", y="P CQ")+lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))+
  scale_color_manual(values = c("non-impacted"="#8DAA91", "ag"="#E6B422", "urban"="#4682B4", "ag & urban"="#A67B5B"))


k3<-monthly_CQ_wide %>%
  left_join(land_use[,c(2,5)]) %>%
  mutate(impacted_class=factor(impacted_class, levels=c("non-impacted", "urban", "ag", "ag & urban"))) %>%
  ggplot(aes(N, P))+geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+geom_point(aes(col=impacted_class), 
                                        size=2, alpha=0.7)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="N CQ", y="P CQ", col="")+lims(x=c(-1.5, 1.5), y=c(-1.5, 1.5))+
  scale_color_manual(values = c("non-impacted"="#8DAA91", "ag"="#E6B422", "urban"="#4682B4", "ag & urban"="#A67B5B"))

pdf("CQ_ShiftDirection_LandUse_Synchrony.pdf", width = 15, height = 8)

ggarrange(p1, p2, p3, k1, k2, k3, nrow = 2, ncol=3, widths = c(0.5, 0.5, 0.7))

dev.off()


CQ_land<-monthly_CQ_wide %>%
  left_join(land_use[,c(2:5)])

