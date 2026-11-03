require(dplyr)
require(tidyr)
require(ggplot2)

setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

results_all<-read.csv("Synchrony_Shift_AllSolutePairs.csv")

chem<-read.csv("WRTDS_Outputs_Clean_02022026.csv")

land_use<-read.csv("ImpactedClass_Synchrony.csv")

CQ<-read.csv("CQ_raw_slopes.csv")
colnames(CQ)[3]<-"chemical"

cvc_cvq_pairs <- CQ %>%
  inner_join(CQ,
             by = "Stream_Name",
             suffix = c("_1", "_2")) %>%
  filter(chemical_1 < chemical_2) %>%
  transmute(
    Stream_Name,
    solute_pair = paste(chemical_1, chemical_2, sep = "-"),
    Solute1 = chemical_1,
    Solute2 = chemical_2,
    cvccvq1  = cv_cq_1,
    cvccvq2  = cv_cq_2
  )

monthly_CQ <- chem %>%
  dplyr::group_by(Stream_Name, chemical) %>%
  dplyr::summarise(
    slope = coef(lm(log10(FNConc_mgL) ~ log10(Discharge_cms)))[2],  # slope = b in log(C) ~ log(Q)
    intercept = coef(lm(log10(FNConc_mgL) ~ log10(Discharge_cms)))[1],  # optional
    n = n(),  # number of points used
    r2 = summary(lm(log10(FNConc_mgL) ~ log10(Discharge_cms)))$r.squared,  # optional fit metric
    .groups = "drop"
  )

monthly_CQ_clean <- monthly_CQ %>%
  dplyr::group_by(Stream_Name, chemical) %>%
  dplyr::summarise(slope = mean(slope, na.rm = TRUE), .groups = "drop")

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

pdf("Residual_CVC_CVQ.pdf", width = 10, height = 5)

p1<-cvc_cvq_pairs %>%
  left_join(results_all) %>%
  left_join(CQ_pairs) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  #filter(sync_char=="synchronous") %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  ggplot(aes(cvccvq1, cvccvq2, col=abs(residual), shape = sync_char))+geom_point()+facet_wrap(~solute_pair)+xlim(0,2.5)+ylim(0,2.5)+
  geom_hline(yintercept = 0.5, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "null", 
        legend.text = element_text(angle = 45, hjust = 1))+
  scale_color_binned(breaks = c(0.1, 0.25, 0.5), palette = c("#0a4e9b", "lightskyblue2", "palegreen3", "darkgreen"))+
  #scale_color_gradientn(colours = c("grey80", "skyblue", "deepskyblue", "dodgerblue4"),
  #                      breaks=c(0,0.1,0.25,0.5))+
  labs(x="CV(C)/CV(Q) Solute 1", y="CV(C)/CV(Q) Solute 2", col="Abs(Residual off 1:1 Line)")+
  scale_shape_manual(values = c("synchronous"=16, "asynchronous"=5))
p1

cvc_cvq_pairs %>%
  left_join(results_all) %>%
  left_join(CQ_pairs) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  mutate(residual=(slope2-slope1)/sqrt(2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward lag",
    corr.best_shift > 6~"backward lag"
  )) %>%
  ggplot(aes(cvccvq1, cvccvq2, col=shift_direction))+geom_point()+facet_wrap(~solute_pair)+xlim(0,2.5)+ylim(0,2.5)+
  geom_hline(yintercept = 0.5, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("in-phase" = "black",
                                    "out-of-phase" = "white",     # white fill
                                    "forward lag" = "#E7A1B0",
                                    "backward lag" = "#F6C28B"))+
  #scale_color_binned(breaks = c(0.1, 0.25, 0.5), palette = c("grey85", "skyblue", "dodgerblue2", "darkblue"))+
  #scale_color_gradientn(colours = c("grey80", "skyblue", "deepskyblue", "dodgerblue4"),
  #                      breaks=c(0,0.1,0.25,0.5))+
  labs(x="CV(C)/CV(Q) Solute 1", y="CV(C)/CV(Q) Solute 2", col="Abs(Residual off 1:1 Line)")

dev.off()
                 
p2<-cvc_cvq_pairs %>%
  left_join(results_all) %>%
  left_join(CQ_pairs) %>%
  mutate(sync_char = ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  #filter(sync_char == "synchronous") %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(
    shift_direction = case_when(
      corr.best_shift == 0 ~ "in-phase",
      corr.best_shift == 6 ~ "out-of-phase",
      corr.best_shift < 6 ~ "forward lag",
      corr.best_shift > 6 ~ "backward lag"
    ),
    residual_bin = cut(abs(residual),
                       breaks = c(0, 0.1, 0.25, 0.5, Inf),
                       labels = c("0-0.1", "0.1-0.25", "0.25-0.5", "> 0.5")),
    cvc_cvq_bin=case_when(
      cvccvq1 < 0.5 & cvccvq2 < 0.5 ~ "static-static",
      cvccvq1 < 0.5 & cvccvq2 > 0.5 ~ "static-dynamic",
      cvccvq1 > 0.5 & cvccvq2 < 0.5 ~ "dynamic-static",
      cvccvq1 > 0.5 & cvccvq2 > 0.5 ~ "dynamic-dynamic"
      ),
    cvc_cvq_bin=factor(cvc_cvq_bin, levels = c("dynamic-dynamic","dynamic-static",
                                               "static-dynamic", "static-static"))) %>%
  group_by(solute_pair, residual_bin, cvc_cvq_bin) %>%
  tally() %>%
  ggplot(aes(y=cvc_cvq_bin, x=n, fill=residual_bin))+geom_bar(stat = "identity")+
  facet_wrap(~solute_pair)+
  scale_fill_manual(values = c("#0a4e9b", "lightskyblue2", "palegreen3", "darkgreen"))+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(y="", x="Number of Sites", fill="Abs(Residual) off 1:1 Line")

p2       

p3<-CQ %>%
  left_join(results_all[2:13]) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(shift = factor(shift, levels=c("non-shape similar", "0", "1", "2", "3", "4", "5", "6"))) %>%
  #filter(sync_char=="synchronous") %>%
  ggplot(aes(cv_cq, slope, col=shift))+geom_point()+facet_wrap(~chemical)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Shift")+
  scale_color_manual(values=c("0"="#0a4e9b", "4"="#FF8E81", "1"="#FFC8B4", "3"="#F0A597", "5"="#EB7773",
                              "2"="#FCCDC3", "6"="#0096ce", "non-shape similar"="grey"))

p3

p4<-CQ %>%
  left_join(results_all[2:13]) %>%
  left_join(CQ_pairs) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  #filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(residual_bin = cut(abs(residual),
                            breaks = c(0, 0.1, 0.25, 0.5, Inf),
                            labels = c("0-0.1", "0.1-0.25", "0.25-0.5", "> 0.5"))) %>%
  ggplot(aes(cv_cq, slope, col=residual_bin, shape=sync_class))+geom_point(size=2.5, alpha=0.7)+
  facet_wrap(~interaction(solute_pair, chemical), nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Residual")+
  #scale_color_manual(values = c("grey85", "skyblue", "dodgerblue2", "darkblue"))
  scale_color_manual(values = c("#0a4e9b", "lightskyblue2", "palegreen3", "darkgreen"))+
  scale_shape_manual(values = c("synchronous"=16, "asynchronous"=5))
p4

p5<-CQ %>%
  left_join(results_all[2:13]) %>%
  left_join(CQ_pairs) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  left_join(land_use[,c(2,5)]) %>%
  #filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(residual_bin = cut(abs(residual),
                            breaks = c(0, 0.1, 0.25, 0.5, Inf),
                            labels = c("0-0.1", "0.1-0.25", "0.25-0.5", "> 0.5"))) %>%
  ggplot(aes(cv_cq, slope, col=impacted_class, shape=impacted_class))+geom_point(size=2.5, alpha=0.7)+
  facet_wrap(~interaction(solute_pair, chemical), nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Residual")+
  scale_color_manual(values = c("less-impacted"="#8DAA91", 
                                "agricultural"="#E6B422", "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, "urban"=17, "agricultural & urban"=15))
p5

# p4<-monthly_CQ %>%
#   left_join(results_all) %>%
#   left_join(monthly_CQ_clean) %>%
#   left_join(CQ_pairs) %>%
#   mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
#          solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
#   filter(sync_char=="synchronous") %>%
#   filter(str_detect(solute_pair, fixed(chemical))) %>%
#   mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
#   mutate(residual_bin = cut(abs(residual),
#                             breaks = c(0, 0.1, 0.25, 0.5, Inf),
#                             labels = c("0-0.1", "0.1-0.25", "0.25-0.5", "> 0.5"))) %>%
#   ggplot(aes(cvc_cvq, slope, col=shift))+geom_point(size=2, alpha=0.7)+facet_wrap(~interaction(solute_pair, chemical), nrow=3)+xlim(0,2.5)+
#   geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
#   theme_classic()+
#   theme(text = element_text(size = 20), legend.position = "bottom")+
#   labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Shift")+
#   scale_color_viridis_c(option = "G")
# p4

p5<-CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2,5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(residual_bin = cut(abs(residual),
                            breaks = c(0, 0.1, 0.25, 0.5, Inf),
                            labels = c("0-0.1", "0.1-0.25", "0.25-0.5", "> 0.5"))) %>%
  ggplot(aes(cvc_cvq, slope, col=impacted_class, shape=impacted_class))+geom_point(size=2.5, alpha=0.7)+
  facet_wrap(~interaction(solute_pair, chemical), nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="", shape="")+
  scale_color_manual(values = c("less-impacted"="#8DAA91", 
                              "agricultural"="#E6B422", "urban"="#756BB1", "agricultural & urban"="#A67B5B"))+
  scale_shape_manual(values = c("less-impacted"=18, "agricultural"=16, "urban"=17, "agricultural & urban"=15))
p5

cq_resid_min<-read.csv("Minimum_Resid_BioControl.csv")
colnames(cq_resid_min)[3]<-"chemical"

p5<-CQ %>%
  left_join(cq_resid_min[,c(2:7)]) %>%
  mutate(
    res_bin = case_when(
      mean_res > -0.5              ~ "> -0.5",
      mean_res <= -0.5 & mean_res > -0.75 ~ "-0.5 to -0.75",
      mean_res <= -0.75 & mean_res > -1   ~ "-0.75 to -1",
      mean_res <= -1              ~ "< -1"
    )
  ) %>%
  mutate(
    res_bin = factor(
      res_bin,
      levels = c("> -0.5", "-0.5 to -0.75", "-0.75 to -1", "< -1")
    )
  ) %>%
  ggplot(aes(cv_cq, slope, col=res_bin, size=r2))+geom_point(alpha=0.7)+
  facet_wrap(~chemical, nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Residual of\nMinimum\nConcentration")+
  scale_color_manual(
    values = c(
      "> -0.5" = "grey",
      "-0.5 to -0.75" = "lightgreen",
      "-0.75 to -1" = "forestgreen",
      "< -1" = "darkgreen"
    )
  )

p5

p6<-CQ %>%
  left_join(land_use[,c(2:5)]) %>%
  ggplot(aes(cv_cq, slope, col=impacted_class, size=r2))+geom_point(alpha=0.7)+
  facet_wrap(~chemical, nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Impacted\nClass")+
  scale_color_manual(values = c("less-impacted"="#8DAA91", 
                                "agricultural"="#E6B422", "urban"="#756BB1", "agricultural & urban"="#A67B5B"))

p6

ggarrange(p6, p5)

p7<-CQ %>%
  left_join(cq_resid_min[,c(2:7)]) %>%
  ggplot(aes(cv_cq, slope, col=as.factor(month)))+geom_point(size=2.5, alpha=0.7)+
  facet_wrap(~chemical, nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Month of\nMinimum")

p7

ggarrange(p5, p7, align = "h")

cq_resid_max<-read.csv("Max_Resid_AnthoControl.csv")
colnames(cq_resid_max)[3]<-"chemical"

p6<-CQ %>%
  left_join(cq_resid_max[,c(2:7)]) %>%
  mutate(
    res_bin = case_when(
      mean_res < 0.5              ~ "< 0.5",
      mean_res >= 0.5 & mean_res < 0.75 ~ "0.5 to 0.75",
      mean_res >= 0.75 & mean_res < 1   ~ "0.75 to 1",
      mean_res >= 1              ~ "> 1"
    )
  ) %>%
  mutate(
    res_bin = factor(
      res_bin,
      levels = c("< 0.5", "0.5 to 0.75", "0.75 to 1", "> 1")
    )
  ) %>%
  ggplot(aes(cv_cq, slope, col=res_bin, size=r2))+geom_point(alpha=0.7)+
  facet_wrap(~chemical, nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "right")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Residual of\nMaximum\nConcentration")+
  scale_color_manual(
    values = c(
      "< 0.5" = "grey",
      "0.5 to 0.75" = "plum",
      "0.75 to 1" = "mediumorchid",
      "> 1" = "darkmagenta"
    )
  )

p6

pdf("CQ_CVCQ_Residual.pdf", width = 13, height = 6)
ggarrange(p6, p5, widths = c(0.56, 0.5))

dev.off()

p8<-CQ %>%
  left_join(cq_resid_max[,c(2:7)]) %>%
  ggplot(aes(cv_cq, slope, col=as.factor(month)))+geom_point(size=2.5, alpha=0.7)+
  facet_wrap(~chemical, nrow=3)+xlim(0,2.5)+
  geom_hline(yintercept = 0, lty="dashed")+geom_vline(xintercept = 0.5, lty="dashed")+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  labs(x="CV(C)/CV(Q)", y="CQ Slope", col="Month of\nMaximum")

p8

p5


cq_resid_min[,c(2:7)] %>%
  left_join(CQ[2:9]) %>%
  left_join(results_all) %>%
  left_join(CQ_pairs) %>%
  mutate(
    res_bin = case_when(
      mean_res > -0.5              ~ "> -0.5",
      mean_res <= -0.5 & mean_res > -0.75 ~ "-0.5 to -0.75",
      mean_res <= -0.75 & mean_res > -1   ~ "-0.75 to -1",
      mean_res <= -1              ~ "< -1"
    )
  ) %>%
  mutate(
    res_bin = factor(
      res_bin,
      levels = c("> -0.5", "-0.5 to -0.75", "-0.75 to -1", "< -1")
    )
  ) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(shift = factor(shift, levels=c("non-shape similar", "0", "1", "2", "3", "4", "5", "6"))) %>%
  #filter(chemical %in% c("N", "DSi")) %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  ggplot(aes(shift, residual, col=res_bin, size=r2))+geom_boxplot(aes(group = shift), outliers = F)+
  geom_jitter()+
  facet_wrap(~interaction(solute_pair, chemical), nrow=3)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  scale_color_manual(
    values = c(
      "> -0.5" = "grey",
      "-0.5 to -0.75" = "lightgreen",
      "-0.75 to -1" = "forestgreen",
      "< -1" = "darkgreen"
    )
  )

pdf("Bio_drawdown_lag_prop.pdf", width = 8, height = 12)

cq_resid_min[,c(2:7)] %>%
  left_join(CQ[2:9]) %>%
  left_join(results_all) %>%
  left_join(CQ_pairs) %>%
  mutate(
    res_bin = case_when(
      mean_res > -0.5              ~ "> -0.5",
      mean_res <= -0.5 & mean_res > -0.75 ~ "-0.5 to -0.75",
      mean_res <= -0.75 & mean_res > -1   ~ "-0.75 to -1",
      mean_res <= -1              ~ "< -1"
    )
  ) %>%
  mutate(
    res_bin = factor(
      res_bin,
      levels = c("> -0.5", "-0.5 to -0.75", "-0.75 to -1", "< -1")
    )
  ) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  mutate(shift = case_when(
    sync_char=="asynchronous"~"non-shape similar",
    .default = as.character(shift)
  )) %>%
  mutate(shift = factor(shift, levels=c("non-shape similar", "0", "1", "2", "3", "4", "5", "6"))) %>%
  #filter(chemical %in% c("N", "DSi")) %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  ggplot(aes(shift, fill=res_bin))+
  geom_bar(position = "fill")+
  facet_wrap(~interaction(solute_pair, chemical), nrow=3)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  scale_fill_manual(
    values = c(
      "> -0.5" = "grey",
      "-0.5 to -0.75" = "lightgreen",
      "-0.75 to -1" = "forestgreen",
      "< -1" = "darkgreen"
    )
  )+scale_y_continuous(trans = "reverse")

dev.off()

pdf("Bio_drawdown_month_prop.pdf", width = 8, height = 12)

cq_resid_min[,c(2:7)] %>%
  mutate(
    res_bin = case_when(
      mean_res > -0.5              ~ "> -0.5",
      mean_res <= -0.5 & mean_res > -0.75 ~ "-0.5 to -0.75",
      mean_res <= -0.75 & mean_res > -1   ~ "-0.75 to -1",
      mean_res <= -1              ~ "< -1"
    )
  ) %>%
  mutate(
    res_bin = factor(
      res_bin,
      levels = c("> -0.5", "-0.5 to -0.75", "-0.75 to -1", "< -1")
    )
  ) %>%
  ggplot(aes(month, fill=res_bin))+
  geom_bar(position = "fill")+
  facet_wrap(~chemical, nrow=3)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+
  scale_fill_manual(
    values = c(
      "> -0.5" = "grey",
      "-0.5 to -0.75" = "lightgreen",
      "-0.75 to -1" = "forestgreen",
      "< -1" = "darkgreen"
    )
  )+scale_y_continuous(trans = "reverse")

dev.off()



pdf("CVC_CQ_Residual_LULC.pdf", width = 18, height = 12)

ggarrange(p4, p5, nrow = 1)

dev.off()

pdf("Residual_CVC_CVQ_Raw.pdf", width = 13, height = 8)

ggarrange(p1, p2, nrow = 2, heights = c(0.55, 0.5), align = "v")

dev.off()

pdf("Residual_CVC_CVQ_V3.pdf", width = 10, height = 10)

ggarrange(p1, p3, nrow = 2, heights = c(0.5, 0.5))

dev.off()

monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2,5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2))
  

p1<-monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2,5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(cq_cvc_cvq_bin=case_when(
    #slope > -0.1 & slope < 0.1 ~ "chemostatic",
    slope < 0 & cvc_cvq < 0.5 ~ "static-diluting",
    slope < 0 & cvc_cvq > 0.5 ~ "dynamic-diluting",
    slope > 0 & cvc_cvq < 0.5 ~ "static-mobilizing",
    slope > 0 & cvc_cvq > 0.5 ~ "dynamic-mobilizing"
  )) %>%
  ggplot(aes(y=interaction(solute_pair, chemical), x=residual, 
             group = interaction(solute_pair, chemical, cq_cvc_cvq_bin), 
             fill=cq_cvc_cvq_bin))+
  geom_vline(xintercept=0)+geom_density_ridges(alpha=0.4, scale = 1)+theme_classic()+
  theme(legend.position = "null")+theme(text = element_text(size = 20))+
  scale_y_discrete(labels=c("DSi (DSi-N)","DSi (DSi-P)","N (DSi-N)","N (N-P)",
                            "P (DSi-P)","P (N-P)"))+
  labs(x="Residual from 1:1")

k1<-ggarrange(p1, p2, nrow = 1)

p2<-monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2,5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(cq_cvc_cvq_bin=case_when(
    #slope > -0.1 & slope < 0.1 ~ "chemostatic",
    slope < 0 & cvc_cvq < 0.5 ~ "static-diluting",
    slope < 0 & cvc_cvq > 0.5 ~ "dynamic-diluting",
    slope > 0 & cvc_cvq < 0.5 ~ "static-mobilizing",
    slope > 0 & cvc_cvq > 0.5 ~ "dynamic-mobilizing"
  )) %>%
  ggplot(aes(y=interaction(solute_pair, chemical), x=shift, 
             group = interaction(solute_pair, chemical, cq_cvc_cvq_bin), 
             fill=cq_cvc_cvq_bin))+
  geom_density_ridges(alpha=0.4, scale = 1)+theme_classic()+
  theme(text = element_text(size = 20))+
  scale_y_discrete(labels=c("DSi (DSi-N)","DSi (DSi-P)","N (DSi-N)","N (N-P)",
                            "P (DSi-P)","P (N-P)"))+
  labs(x="Lag")

p2

k2<-ggarrange(p1, p2)
k2

p3<-monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2:5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(cq_cvc_cvq_bin=case_when(
    #slope > -0.1 & slope < 0.1 ~ "chemostatic",
    slope < 0 & cvc_cvq < 0.5 ~ "static-diluting",
    slope < 0 & cvc_cvq > 0.5 ~ "dynamic-diluting",
    slope > 0 & cvc_cvq < 0.5 ~ "static-mobilizing",
    slope > 0 & cvc_cvq > 0.5 ~ "dynamic-mobilizing"
  )) %>%
  ggplot(aes(y=interaction(solute_pair, chemical), x=Cropland, 
             group = interaction(solute_pair, chemical, cq_cvc_cvq_bin), 
             fill=cq_cvc_cvq_bin))+
  geom_density_ridges(alpha=0.4, scale = 1)+theme_classic()+
  theme(legend.position = "null")+theme(text = element_text(size = 20))+
  scale_y_discrete(labels=c("DSi (DSi-N)","DSi (DSi-P)","N (DSi-N)","N (N-P)",
                            "P (DSi-P)","P (N-P)"))

k3<-ggarrange(p1, p2)

p4<-monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2:5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(cq_cvc_cvq_bin=case_when(
    slope < 0 & cvc_cvq < 0.5 ~ "static-diluting",
    slope < 0 & cvc_cvq > 0.5 ~ "dynamic-diluting",
    slope > 0 & cvc_cvq < 0.5 ~ "static-mobilizing",
    slope > 0 & cvc_cvq > 0.5 ~ "dynamic-mobilizing"
  )) %>%
  ggplot(aes(y=interaction(solute_pair, chemical), x=log(Impervious+0.001), 
             group = interaction(solute_pair, chemical, cq_cvc_cvq_bin), 
             fill=cq_cvc_cvq_bin))+
  geom_density_ridges(alpha=0.4, scale = 1)+theme_classic()+
  labs(x="Log(Impervious)")+theme(text = element_text(size = 20))+
  scale_y_discrete(labels=c("DSi (DSi-N)","DSi (DSi-P)","N (DSi-N)","N (N-P)",
                            "P (DSi-P)","P (N-P)"))
p4
pdf("RidgePlots.pdf", width = 15, height = 14)

ggarrange(p1, p2, p3, p4, widths = c(0.5, 0.65))

dev.off()

pdf("RidgePlots_test.pdf", width = 15, height = 7)
k1
k2
k3
dev.off()

ridge_plot_df<-monthly_CQ %>%
  left_join(results_all) %>%
  left_join(monthly_CQ_clean) %>%
  left_join(CQ_pairs) %>%
  left_join(land_use[,c(2:5)]) %>%
  #filter(impacted_class!="less-impacted") %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous"),
         solute_pair = paste(Solute1, Solute2, sep = "-")) %>%
  filter(sync_char=="synchronous") %>%
  filter(str_detect(solute_pair, fixed(chemical))) %>%
  mutate(residual = (slope2 - slope1) / sqrt(2)) %>%
  mutate(cq_cvc_cvq_bin=case_when(
    slope < 0 & cvc_cvq < 0.5 ~ "static-diluting",
    slope < 0 & cvc_cvq > 0.5 ~ "dynamic-diluting",
    slope > 0 & cvc_cvq < 0.5 ~ "static-mobilizing",
    slope > 0 & cvc_cvq > 0.5 ~ "dynamic-mobilizing"
  ))

write.csv(ridge_plot_df, "Data_for_RidgePlot.csv")
