# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

n <- read.csv("C:/Users/diamo/Desktop/loire_n.csv")
n$tn = n$TKN+n$N.NO3
plot(n$Annee, n$tn)
cor.test(x = n$Annee, y = n$TKN, method = "kendall")
lm(TKN~Annee, data = n)


q <- readRDS("C:/Users/diamo/Desktop/dampierre_discharge_daily") %>%
  filter(date >= ymd("1993-01-01")) %>%
  mutate(Q = ifelse(discharge.daily < 0, NA, discharge.daily),
         Date = as.POSIXct(date)) %>%
  select(Q, Date)

q_s <- q %>%
  mutate(year = year(Date)) %>%
  group_by(year) %>%
  summarize(l = sum(Q<200))
library(hydrostats)
# Hydrostatistics
df_q_stat <- q %>%
  group_by(year(Date)) %>%
  nest() %>%
  mutate(ls = map(data, low.spells, threshold = 200)) %>%
  unnest(ls) %>%
  rename(year = `year(Date)`)

mean(df_q_stat$avg.low.spell.duration)
df_lowflows <- df_q_stat %>%
  select(year = `year(Date)`, lf_dur = max.low.duration)
saveRDS(df_lowflows, "Data/Discharge/low_flow_duration")
df_lowflows <- readRDS("Data/Discharge/low_flow_duration")
df_q_stat2 <- df_q %>%
  mutate(threshold = ifelse(Q < 200, 1, 0)) %>%
  group_by(year(Date)) %>%
  summarize(run = rle(df_q_stat$threshold))

qa <- dat %>%
  select(date, GPP, ER) %>%
  mutate(date = mdy(date),
         year = year(date),
         month = month(date)) %>%
  filter(between(month, 4, 10)) %>%
  group_by(year) %>%
  summarize(GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  left_join(df_q_stat %>%
              select(-data)) %>%
  left_join(q_s)

# normalize gpp
qad <- qa %>%
  transmute(gppn = GPP / (avg.low.spell.duration),
            gppn2 = GPP / (max.low.duration),
            gpp3 = GPP /l,
            period = if_else(year < 2014, 1, 2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize_all(.funs = median, na.rm = T)
plot(qad$gppn2)

  
plot(qad$gppn)
df_met_l %>%
  mutate(year = year(date)) %>%
  group_by(year, key) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE)) %>%
  left_join(df_lowflows) %>%
  # left_join(df_q_stat %>%
  #             select(year = "year(Date)", ls = avg.low.spell.duration)) %>%
  mutate(rel = mean / lf_dur,
         relmed = median /lf_dur,
         pd = if_else(year<2013, 0, 1)) %>%
  group_by(pd, key) %>% 
  summarize(mean = mean(rel, na.rm = TRUE), 
            median = median(relmed, na.rm = TRUE))


library(trend)
mk <- dat %>%
  select(date, GPP, ER) %>%
  mutate(date = mdy(date)) %>%
  left_join(q2)
mk <- mk[-(9486:9487),]
mk <- na_kalman(mk)

mk2 <- mk %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(year, month) %>%
  summarize_all(.funs = median, na.rm = T)


partial.mk.test(mk$GPP, log(mk$q))


gpp_ts <- ts(mk$GPP, frequency = 365)
stl_gpp = stl(gpp_ts, "periodic")
seasonal_stl_gpp   <- stl_gpp$time.series[,1]
trend_stl_gpp     <- stl_gpp$time.series[,2]
random_stl_gpp  <- stl_gpp$time.series[,3]
plot(trend_stl_gpp)
sens.slope(trend_stl_gpp)
sens.slope(gpp_ts)
q_ts <- ts(log(mk$q), frequency = 365)
stl_q <- stl(q_ts, "periodic")
seasonal_stl_q   <- stl_q$time.series[,1]
trend_stl_q     <- stl_q$time.series[,2]
random_stl_q  <- stl_q$time.series[,3]
plot(trend_stl_q)
sens.slope(q_ts)
partial.mk.test(trend_stl_gpp, trend_stl_q)
plot(gpp_ts)




gpp_ts2 <- ts(mk2$GPP, frequency = 12)
stl_gpp2 = stl(gpp_ts, "periodic")
seasonal_stl_gpp2   <- stl_gpp2$time.series[,1]
trend_stl_gpp2    <- stl_gpp2$time.series[,2]
random_stl_gpp2  <- stl_gpp2$time.series[,3]
plot(trend_stl_gpp)
sens.slope(trend_stl_gpp2)
sens.slope(gpp_ts2)
mk.test(gpp_ts2)

q_ts2 <- ts(log(mk2$q), frequency = 12)
stl_q2 <- stl(q_ts2, "periodic")
seasonal_stl_q2   <- stl_q2$time.series[,1]
trend_stl_q2     <- stl_q2$time.series[,2]
random_stl_q2  <- stl_q2$time.series[,3]
plot(trend_stl_q2)
sens.slope(q_ts2)
mk.test(q_ts2)
partial.mk.test(gpp_ts2, q_ts2)
partial.cor.trend.test(gpp_ts2, q_ts2)

partial.mk.test(trend_stl_gpp2, trend_stl_q2)
plot(gpp_ts)

