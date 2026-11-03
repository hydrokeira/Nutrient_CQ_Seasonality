library(dplyr)
library(ggplot2)

set.seed(123)

# ==========================================================
# COMMON Si:N TIME SERIES
# ==========================================================

dat <- expand.grid(
  Year = 2000:2024,
  Month = 1:12
) %>%
  arrange(Year, Month) %>%
  mutate(
    
    season = sin(2 * pi * (Month - 3) / 12),
    
    summer_effect = case_when(
      Month %in% 6:8 ~ 1,
      Month %in% c(5, 9) ~ 0.5,
      TRUE ~ 0.1
    ),
    
    # SAME Si:N in all three scenarios
    # Seasonal pattern + summer-specific temporal change
    log_Si_N =
      2 -
      0.06 * (Year - 2000) * summer_effect +
      0.25 * season +
      rnorm(n(), 0, 0.05),
    
    Si_N = exp(log_Si_N)
  )

l1<-ggplot(dat, aes(Month, Si_N))+geom_line(aes(col=Year, group = Year))+
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Si:N Change",
    x = "Month",
    y = "Si:N"
  )

x3<-ggarrange(l1, l1, l1, l1, nrow = 4)

dat_Q_controlled <- dat %>%
  mutate(
    Q = exp(
      7 -
        2 * log_Si_N +
        rnorm(n(), 0, 0.15)
    )
  )


dat_Q_controlled <- dat %>%
  mutate(
    Q = exp(
      7 -
        2 * log_Si_N +
        rnorm(n(), 0, 0.15)
    )
  )

p1<-ggplot(dat_Q_controlled, aes(Month, Q))+geom_line(aes(col=Year, group = Year))+
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Non-linear discharge change",
    x = "Month",
    y = "Discharge"
  )

dat_Q_controlled_shifted <- dat %>%
  mutate(
    Q = exp(
      7 -
        2 * log_Si_N +
        0.05 * (Year - 2000) + 
        rnorm(n(), 0, 0.15)
    )
  )

p2<-ggplot(dat_Q_controlled_shifted, aes(Month, Q))+geom_line(aes(col=Year, group = Year))+
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Increasingly non-linear discharge change",
    x = "Month",
    y = "Discharge"
  )

dat_constant_Q <- dat %>%
  mutate(
    Q = exp(
      5 +
        0.8 * season +
        rnorm(n(), 0, 0.1)
    )
  )

p3<-ggplot(dat_constant_Q, aes(Month, Q))+geom_line(aes(col=Year, group = Year))+
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "No discharge change",
    x = "Month",
    y = "Discharge"
  )

dat_Q_decoupled <- dat %>%
  mutate(
    Q = exp(
      5 +
        0.6 * season +
        0.03 * (Year - 2000) +
        rnorm(n(), 0, 0.1)
    )
  )

p4<-ggplot(dat_Q_decoupled, aes(Month, Q))+geom_line(aes(col=Year, group = Year))+
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Constant discharge change",
    x = "Month",
    y = "Discharge"
  )

x1<-ggarrange(p3, p4, p1, p2, nrow = 4)

x1

k1<-ggplot(
  dat_Q_controlled_shifted,
  aes(log(Q), log(Si_N), colour = Year)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group=Year, col=Year),
              method = "lm",
              se = FALSE) +
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Discharge change coupled with Si:N",
    x = "log(Q)",
    y = "log(Si:N)"
  )

k2<-ggplot(
  dat_Q_controlled,
  aes(log(Q), log(Si_N), colour = Year)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group=Year, col=Year),
              method = "lm",
              se = FALSE) +
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Discharge change coupled with Si:N",
    x = "log(Q)",
    y = "log(Si:N)"
  )

k3<-ggplot(
  dat_constant_Q,
  aes(log(Q), log(Si_N), colour = Year)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group=Year, col=Year),
              method = "lm",
              se = FALSE) +
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "No discharge change",
    x = "log(Q)",
    y = "log(Si:N)"
  )

k4<-ggplot(
  dat_Q_decoupled,
  aes(log(Q), log(Si_N), colour = Year)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group=Year, col=Year),
              method = "lm",
              se = FALSE) +
  scale_color_viridis_c(option = "G") +
  theme_classic() +
  theme(legend.position = "null", text=element_text(size=15))+
  labs(
    title = "Discharge change decoupled from Si:N",
    x = "log(Q)",
    y = "log(Si:N)"
  )


dat_Q_decoupled %>%
  group_by(Year) %>%
  summarise(
    n = sum(complete.cases(Si_N, Q)),
    SQ_slope = if(n >= 3) {
      coef(lm(log(Si_N) ~ log(Q)))[2]
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  summarise(sd(SQ_slope))

x2<-ggarrange(k3, k4, k2,k1, nrow = 4)

pdf("Example_Change_CQ_Slope.pdf", width = 15, height = 10)
ggarrange(x3, x1, x2, nrow = 1)
dev.off()
