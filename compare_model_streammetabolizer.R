library(streamMetabolizer)


bayes_name <- mm_name(type='bayes', pool_K600 = "none", err_obs_iid = FALSE)
bayes_name <- mm_name(type="mle")
bayes_specs <- specs(bayes_name)
bayes_specs


df_sm$solar.time <- lubridate::force_tz(df_sm$datetime, 'Etc/GMT+0')

lubridate::tz(df_sm$solar.time)
df_sm$solar.time <- streamMetabolizer::calc_solar_time(df_sm$solar.time, longitude=0)
df_sm$datetime <- NULL
# df_use <- df_sm %>% dplyr::filter(row_number() %% 23 != 1)
mm <- metab(bayes_specs, data=df_sm)

predict_metab(mm)
plot_metab_preds(mm)
get_params(mm)
plot_DO_preds(mm)
get_fit(mm)
traceplot(get_mcmc(mm))
get_fit(mm)$warnings
get_fit(mm)$errors
select(predict_metab(mm), warnings, errors)
mcmc <- get_mcmc(mm)
rstan::traceplot(mcmc, pars='K600_daily', nrow=3)
get_fit(mm)$overall %>%
  select(ends_with('Rhat'))
