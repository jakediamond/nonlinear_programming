library(streamMetabolizer)
library(tidyverse)

bayes_name <- mm_name(type='bayes', pool_K600 = "none", 
                      err_proc_iid = FALSE, err_obs_iid = TRUE)
bayes_name <- mm_name(type="mle", GPP_fun='satlight', ER_fun='constant')
bayes_specs <- specs(bayes_name)
bayes_specs


df_sm$solar.time <- lubridate::force_tz(df_sm$datetime, 'Etc/GMT+0')

lubridate::tz(df_sm$solar.time)
df_sm$solar.time <- streamMetabolizer::calc_solar_time(df_sm$solar.time, longitude=0)
df_sm$datetime <- NULL
df_sm <- df_sm[-c(1:4, 174:193), ]
mm <- metab(bayes_specs, data=df_sm)

predict_metab(mm)
plot_metab_preds(mm)
get_params(mm)
plot_DO_preds(mm)



c<- get_fit(mm)$inst
traceplot(get_mcmc(mm))
get_fit(mm)$warnings
get_fit(mm)$errors
select(predict_metab(mm), warnings, errors)
mcmc <- get_mcmc(mm)
rstan::traceplot(mcmc, pars='K600_daily', nrow=3)
get_fit(mm)$overall %>%
  select(ends_with('Rhat'))
get_fit(mm)$inst$err_proc_iid_mean


dat <- data_metab(num_days='3', res='15')
head(data_metab)
