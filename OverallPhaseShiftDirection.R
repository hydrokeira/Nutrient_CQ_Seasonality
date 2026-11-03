setwd("/Users/keirajohnson/Library/CloudStorage/Box-Box/Keira_Johnson/SiSyn/NutrientSynchrony")

results_all<-read.csv("Synchrony_Shift_AllSolutePairs.csv")

plot_b<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "low shape similarity", "high shape similarity")) %>%
  mutate(sync_char=factor(sync_char, levels=c("low shape similarity", "high shape similarity"))) %>%
  ggplot(aes(y=1))+geom_bar(aes(fill=sync_char), position = "fill")+
  theme_classic()+labs(x="", y="", fill="")+
  scale_fill_manual(values = c("grey", "black"))+
  theme(text = element_text(size = 20), legend.position = "top", axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  facet_wrap(~solute_pair, scales="free_y", nrow=3)+
  scale_x_continuous(labels = scales::percent)

plot_b

results_sum<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "low shape similarity", "high shape similarity")) %>%
  group_by(solute_pair, sync_class) %>%
  tally() %>%
  pivot_wider(names_from = sync_class, values_from = n)

tab <- results_sum %>%
  column_to_rownames("solute_pair") %>%
  as.matrix()

chisq.test(tab)

prop.table(tab, margin = 1)

chisq.test(tab)$stdres

success <- tab[, "synchronous"]
total <- rowSums(tab)

pairwise.prop.test(success, total, p.adjust.method = "BH")

results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  group_by(solute_pair, sync_char) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(solute_pair) %>%
  mutate(prop = n / sum(n))

plot_c<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  ggplot(aes(shift))+geom_bar(fill="black")+
  theme_classic()+labs(x="Shift (Months)", y="Number of Sites", fill="")+
  #ggtitle("Distribution of Phase Offset for Synchronous Sites")+
  scale_x_continuous(labels = seq(0,6,1), breaks = seq(0,6,1))+
  theme(text = element_text(size = 20))+
  facet_wrap(~solute_pair)

plot_c

prop_shift_summary<-results_all %>%
  mutate(solute_pair=paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char=ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char=="synchronous") %>%
  group_by(shift, solute_pair) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(solute_pair) %>%
  mutate(prop = n / sum(n))

results_all<-results_all %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  ))

results_all_summary<-results_all %>%
  mutate(solute_pair = paste0(Solute1, "-", Solute2)) %>%
  mutate(shift_direction=case_when(
    corr.best_shift==0~"in-phase",
    corr.best_shift==6~"out-of-phase",
    corr.best_shift < 6~"forward shift",
    corr.best_shift > 6~"backward shift",
  )) %>%
  group_by(solute_pair, shift_direction) %>%
  tally() %>%
  pivot_wider(names_from = shift_direction, values_from = n) %>%
  select(solute_pair, `backward shift`, `forward shift`)

tab <- results_all_summary %>%
  column_to_rownames("solute_pair") %>%
  as.matrix()

chisq.test(tab)

prop.table(tab, margin = 1)

chisq.test(tab)$stdres

success <- tab[, "forward shift"]
total <- rowSums(tab)

pairwise.prop.test(success, total, p.adjust.method = "BH")

p1<-results_all %>%
  mutate(solute_pair = paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char = ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char == "synchronous") %>%
  mutate(shift_direction = case_when(
    corr.best_shift == 0 ~ "in-phase",
    corr.best_shift == 6 ~ "out-of-phase",
    corr.best_shift < 6  ~ "forward shift",
    corr.best_shift > 6  ~ "backward shift"
  )) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("in-phase", "out-of-phase", "backward shift", "forward shift"))) %>%
  group_by(solute_pair, shift_direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(solute_pair) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "", y = prop, fill = shift_direction)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # outlines every slice
  coord_polar(theta = "y") +
  facet_wrap(~solute_pair) +
  scale_fill_manual(values = c(
    "in-phase" = "black",
    "out-of-phase" = "white",     # white fill
    "forward shift" = "#E7A1B0",
    "backward shift" = "#F6C28B"
  )) +
  theme_void() +
  theme(text = element_text(size = 20), legend.position = "bottom")+labs(fill="")

p1
p2<-results_all %>%
  mutate(solute_pair = paste0(Solute1, "-", Solute2)) %>%
  mutate(sync_char = ifelse(corr.synchrony < 0.5, "asynchronous", "synchronous")) %>%
  filter(sync_char == "synchronous") %>%
  mutate(shift_direction = case_when(
    corr.best_shift == 0 ~ "in-phase",
    corr.best_shift == 6 ~ "out-of-phase",
    corr.best_shift < 6  ~ "forward shift",
    corr.best_shift > 6  ~ "backward shift"
  )) %>%
  filter(shift_direction %in% c("forward shift", "backward shift")) %>%
  mutate(shift_direction=factor(shift_direction, levels = c("backward shift", "forward shift"))) %>%
  ggplot(aes(x=shift, fill=shift_direction))+geom_bar(position = "fill")+
  facet_wrap(~solute_pair)+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "bottom")+labs(fill="")+
  scale_fill_manual(values = c(
    "in-phase" = "black",
    "out-of-phase" = "white",     # white fill
    "forward shift" = "#E7A1B0",
    "backward shift" = "#F6C28B"
  ))+
  labs(x="Shift (Months)", y="Proportion")+
  scale_x_continuous(
    breaks = 0:6, labels = 0:6,
    limits = c(0, 6)
  )

p2

pdf("Prop_Synchrony_Lag.pdf", width = 10, height = 8)

ggarrange(plot_c, p2, nrow = 2, align = "v", heights = c(0.4, 0.45))

dev.off()

pdf("Overall_Prop_Synchrony.pdf", width = 4, height = 8)

plot_b

dev.off()
