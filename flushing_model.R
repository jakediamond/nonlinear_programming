library(deSolve)
library(zoo)
library(tidyverse)
set.seed(42)

# Set up times for modeling
simulation_time = 30*365 #days
dt = 1
times <- seq(1, simulation_time, by = dt)

# Light function
light <- data.frame(times = times,
                    import =  0.5*(540+440*sin(2*pi*times/365-1.4)))
inputPAR <- approxfun(light, rule = 2)

# Flushing function
flush_fun <- data.frame(times = times,
                        import =  c(floor(runif(365*3, 1, 9))/10,
                                    floor(runif(365*3, 1, 9))/10,
                                    floor(runif(365*3, 1, 9))/10,
                                    floor(runif(365*3, 1, 9))/10,
                                    floor(runif(365*3, 1, 8))/10,
                                    floor(runif(365*3, 1, 6))/10,
                                    round(runif(365*3, 1, 3))/10,
                                    round(runif(365*3, 1, 2))/10,
                                    runif(365*3, 1, 2)/10,
                                    runif(365*3, 1, 2)/10))
inputFlush <- approxfun(flush_fun, rule = 2)
#-------------------------------------------------
# Model
loire_toy_model <- function(time,states,parms){
  # set.seed(42)
  #Calculate the environment
  PAR = inputPAR(time)
  #unpack states
  BIOMASS = states[1]
  #unpack parms
  maxGPP = parms[1]
  kspar = parms[2]
  autoResp = parms[3]
  heteroResp = parms[4]
  allocInput = parms[5]
  # Internal calculation of new flushing rate based on time
  newFlush = inputFlush(time)
  #calculate fluxes
  gpp = maxGPP * (PAR/(PAR + kspar))
  flushing = newFlush * BIOMASS
  ar = autoResp * gpp
  hr = heteroResp * BIOMASS
  alloc = allocInput
  #calculate derivatives
  dBIOMASS = gpp + alloc - flushing - hr - ar
  # Calculate some other fluxes
  er = -hr - ar
  nep = gpp + er
  #return the list of derivatives plus any other variables that you are interested in
  return(list(c(dBIOMASS),                          # the rate of change
              c(GPP = gpp, ER = er, NEP = nep, flush = newFlush, PAR = PAR)))
}
# Parameters
parms <- c(
  maxGPP = 10,
  kspar = 120,
  autoResp = 0.5,
  heteroResp = 0.2,
  allocInput = 2)
# Initial conditions
yini <- c(
  BIOMASS = 1)

# Model running and analysis ----------------------------------------------
#uses the R deSolve function 
out = ode(y = yini, times = times, func = loire_toy_model,
          parms = parms, method = "lsoda")
# Get the data in dataframe
df <- as.tibble(out) %>%
  mutate_all(as.numeric) %>%
  rename(GPP = GPP.maxGPP, ER = ER.heteroResp, NEP = NEP.maxGPP,
         flushrate = flush, PAR = PAR)
# Quick look at biomass
ggplot(data = df,
       aes(x = time, y = NEP)) +
  geom_point()

# Calculate annual ar1 for NEP
year_ar1 <- df %>%
  mutate(year = ntile(time, 30)) %>%
  group_by(year) %>%
  mutate(day = ntile(time, 365)) %>%
  filter(between(day, 120, 260)) %>%
  # summarize_all(~sd(.))
  summarize_all(~acf(., lag.max= 3, plot=FALSE)$acf[2])

# Plot
p_nep_ar1 <- ggplot(data = left_join(year_ar1,
                                     mutate(flush_fun, year = ntile(times, 30)) %>%
                                       group_by(year) %>%
                                       summarize(flush = mean(import))),
                    aes(x = year, y = NEP,
                        color = flush,
                        group = 1)) +
  stat_smooth(se = FALSE) +
  geom_point() +
  geom_line() +
  theme_bw(base_size = 6) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_color_viridis_c(name = "flushing rate") +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2, 0.75),
        legend.direction = "horizontal") +
  guides(colourbar = guide_legend(title.position = "top")) +
  xlab("year") +
  ylab("ar(1) of summer NEP")
p_nep_ar1

ggsave(plot = p_nep_ar1,
        filename = "//LY-LHQ-SRV/jake.diamond/Loire_DO/Figures/Middle_Loire/Figure4_final.svg",
        dpi = 300,
        width = 9,
        height = 4,
        units = "cm")
library(earlywarnings)
x = as.ts(df$NEP)
generic_ews(x, winsize = 3.33)

library(lubridate)
library(imputeTS)
df_wq <- readRDS("//LY-LHQ-SRV/jake.diamond/Loire_DO/Data/Loire_DO/middle_loire_wq")
p <- df_wq %>%
  filter(between(year,1979,2006),
         solute == "CHLA") %>%
  arrange(date) %>%
  group_by(year, month, solute) %>%
  distinct(date, .keep_all = TRUE) %>%
  summarize(mean = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(solute) %>%
  transmute(date = ymd(paste(year, month, "01", sep = "-")),
            value = mean) %>%
  na.trim() %>%
  na_kalman()
generic_ews(p$value, winsize = 25, logtransform = TRUE, detrending = "gaussian",
            bandwidth = 2)
sensitivity_ews(p$value, logtransform = TRUE, detrending = "gaussian")
