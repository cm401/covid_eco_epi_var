### Modelling the interaction of SARS-CoV-2 transmission intensity, excess deaths and the economy ###

library(here)

if(FALSE)
{
  # initialize packages, data & source functions from files
  source(paste0(here(),"/R/setup.R"))
  
  # specify and run model estimation      
  model_1 <- brm_VAR(covid_data_weekly %>% filter(VALIDATION_TYPE=="TRAINING") %>% 
                       group_by(CountryName) %>% 
                       #mutate(fd_NPI=(StringencyIndex-dplyr::lag(StringencyIndex,1))/100,
                       mutate(fd_NPI=(ContainmentHealthIndex-dplyr::lag(ContainmentHealthIndex,1))/100,
                              fd_NPI_l1 = dplyr::lag(fd_NPI,1),
                              fd_NPI_l2 = dplyr::lag(fd_NPI,2),
                              fd_NPI_l3 = dplyr::lag(fd_NPI,3),
                              fd_NPI_l4 = dplyr::lag(fd_NPI,4),
                              fd_NPI_l5 = dplyr::lag(fd_NPI,5)) %>% ungroup(), 
                     response_variables  = c( "log_R","eco_ed","fd_NPI","fd_GDP","fd_transit"),
                     exogenous_variables = c( ),
                     p                   = 1,
                     priors              = c(prior(normal(0,1),class=b,resp=fdNPI),
                                             prior(normal(0,1),class=b,resp=fdGDP),
                                             prior(normal(0,1),class=b,resp=ecoed),
                                             prior(normal(0,1),class=b,resp=fdtransit),
                                             prior(normal(0,1),class=b,resp=logR),
                                             prior(cauchy(0,2),class=sd,resp=fdNPI),
                                             prior(cauchy(0,2),class=sd,resp=fdGDP),
                                             prior(cauchy(0,2),class=sd,resp=ecoed),
                                             prior(cauchy(0,2),class=sd,resp=fdtransit),
                                             prior(cauchy(0,2),class=sd,resp=logR),
                                             prior(lkj(2), class = rescor)),
                     panel_id            = "CountryName",
                     group_by_variant    = FALSE,
                     iterations          = 4000,
                     no_population_intercept = TRUE,
                     control             = list(adapt_delta = 0.9, step_size = 0.01, max_treedepth = 10),
                     interaction_term    = "vaccination_doses_per_person + WT_variant + Alpha_variant + Delta_variant + Omicron_variant",
  )
  
  model_2 <- brm_VAR(covid_data_weekly %>% filter(VALIDATION_TYPE=="TRAINING"), 
                   response_variables  = c( "log_R","eco_ed","fd_GDP","fd_transit"),
                   exogenous_variables = c( "fd_C1", "lvl_C1L_lag", "fd_C2", "lvl_C2L_lag", "fd_C4", "lvl_C4L_lag", "fd_C5", "lvl_C5L_lag", "fd_C8", "lvl_C8L_lag", 
                                            "fd_H2", "lvl_H2L_lag", "fd_H6", "lvl_H6L_lag", "fd_H8", "lvl_H8L_lag", "fd_E1", "lvl_E1L_lag" ),
                   p                   = 1,
                   priors              = c(prior(normal(0,1),class=b,resp=fdGDP),
                                           prior(normal(0,1),class=b,resp=ecoed),
                                           prior(normal(0,1),class=b,resp=fdtransit),
                                           prior(normal(0,1),class=b,resp=logR),
                                           prior(cauchy(0,2),class=sd,resp=fdGDP),
                                           prior(cauchy(0,2),class=sd,resp=ecoed),
                                           prior(cauchy(0,2),class=sd,resp=fdtransit),
                                           prior(cauchy(0,2),class=sd,resp=logR),
                                           prior(lkj(2), class = rescor)),
                   panel_id            = "CountryName",
                   group_by_variant    = FALSE,
                   iterations          = 4000,
                   no_population_intercept = TRUE,
                   control             = list(adapt_delta = 0.9, step_size = 0.01, max_treedepth = 10),
                   interaction_term    = "vaccination_doses_per_person + WT_variant + Alpha_variant + Delta_variant + Omicron_variant",
  )
  
  model_3 <- brm_VAR(covid_data_weekly %>% filter(VALIDATION_TYPE=="TRAINING"), 
                   response_variables  = c( "log_R","eco_ed","fd_GDP","fd_transit"),
                   exogenous_variables = c( "fd_C1", "lvl_C1L_lag", "fd_C2", "lvl_C2L_lag", "fd_C4", "lvl_C4L_lag", "fd_C5", "lvl_C5L_lag", "fd_C8", "lvl_C8L_lag", 
                                            "fd_H2", "lvl_H2L_lag", "fd_H6", "lvl_H6L_lag", "fd_H8", "lvl_H8L_lag", "fd_E1", "lvl_E1L_lag" ),
                   p                   = 1,
                   priors              = c(prior(normal(0,1),class=b,resp=fdGDP),
                                           prior(normal(0,1),class=b,resp=ecoed),
                                           prior(normal(0,1),class=b,resp=fdtransit),
                                           prior(normal(0,1),class=b,resp=logR),
                                           prior(cauchy(0,2),class=sd,resp=fdGDP),
                                           prior(cauchy(0,2),class=sd,resp=ecoed),
                                           prior(cauchy(0,2),class=sd,resp=fdtransit),
                                           prior(cauchy(0,2),class=sd,resp=logR)),
                   panel_id            = "CountryName",
                   group_by_variant    = FALSE,
                   iterations          = 4000,
                   no_population_intercept = TRUE,
                   control             = list(adapt_delta = 0.9, step_size = 0.01, max_treedepth = 10),
                   interaction_term    = "vaccination_doses_per_person + WT_variant + Alpha_variant + Delta_variant + Omicron_variant",
  )
  
  which_model_to_use <- "model_2"
  
  if(which_model_to_use=="model_1") model <- model_1
  if(which_model_to_use=="model_2") model <- model_2
  if(which_model_to_use=="model_3") model <- model_3
  
  # plot NPI coefficients
  plot_npi_coefficients(model, resp_names = c("log~ED", "Delta~GDP","Delta~Transit", "log~R"))
  
  # plot VAR coefficients
  plot_var_coefficients(model, panel_label = "A.")
  
  # plot generalized impulse response function
  if(which_model_to_use=="model_2") 
  {
    irf_calced <- compute_irf_mcmc(model$models,1,lag=12,facet_scales = "free", plot_labels = c("log~R","log~ED","Delta~GDP","Delta~Transit"),level_order = c("logR","ecoed","fdGDP","fdtransit"))
    irf_calced$oirf_plot + scale_x_continuous(breaks = seq(4, 12, by = 4)) + 
      theme(strip.background=element_rect(fill="#CCFFFF"), strip.text=element_text(color="black")) + 
      ggtitle("B.")
  }
  
  if(which_model_to_use=="model_1") 
  {
    irf_calced <- compute_irf_mcmc(model$models,1,lag=12,facet_scales = "free", plot_labels = c("log~R","log~ED","Delta~NPI","Delta~GDP","Delta~Transit"),level_order = c("logR","ecoed","fdNPI","fdGDP","fdtransit"))
    irf_calced$oirf_plot + scale_x_continuous(breaks = seq(4, 12, by = 4)) + 
      theme(strip.background=element_rect(fill="#CCFFFF"), strip.text=element_text(color="black")) + 
      ggtitle("B.")
  }
  
  # mcmc diagnostic plots
  mcmc_rhat_hist(rhat(model$models)) + mcmc_neff_hist(neff_ratio(model$models), size = 3)
}

