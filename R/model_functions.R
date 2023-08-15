### Model functions ###


#' Title  Create & Run Bayesian multi-level VAR model
#'
#' @param data                      Covid input data tibble with all data items for model / estimation
#' @param start_date      
#' @param end_date 
#' @param response_variables        Response Variables for Vector Autoregressive structure
#' @param exogenous_variables       NPI etc
#' @param p                         order of VAR 
#' @param priors                    Piors for model
#' @param panel_id                  Idenitifier in dataset for grouping
#' @param interaction_term          Any term which may have an interaction (e.g. vaccination & dominant variant)
#' @param group_by_variant          either run model as full time period or separately for each variant
#' @param variants                  Variant to run for if running for each variant separately
#' @param no_population_intercept   exclude population intercept 
#' @param save_to_file              target location of where to save output to
#' @param iterations                how many mcmc iterations to run for
#' @param control                   stan parameters
#'
#' @return
#' @export
#'
#' @examples
brm_VAR <- function(
    data,      
    start_date              = as.Date("2020-01-01"),
    end_date                = as.Date("2022-01-01"),
    response_variables      = c( "log_R","fd_GDP","fd_transit"),
    exogenous_variables     = c( "fd_C1", "fd_C2", "fd_C4", "fd_C5", "fd_C8", "fd_H2", "fd_H6", "fd_H8", "fd_E1" ),
    p                       = 2,
    priors                  = c(prior(normal(0,1),class=b,resp=fdGDP),
                                prior(normal(0,1),class=b,resp=fdtransit),
                                prior(normal(0,1),class=b,resp=logR),
                                prior(lkj(2), class = cor)), 
    panel_id                = "CountryName",
    interaction_term        = "",
    group_by_variant        = FALSE,
    variants                = c("1_wildtype", "2_alpha", "3_delta"),
    no_population_intercept = FALSE,
    save_to_file            = NA,      #provide file name here
    iterations              = 2000,
    control                 = NULL     #list(adapt_delta = 0.9, stepsize = 0.01, max_treedepth = 15)
)
{
  data <- data %>% filter(Date>=start_date & Date<end_date) #enforce start & end date
  
  VAR_components <- c()
  
  for(i in 1:p)
  {
    append <- unlist(lapply(response_variables,paste0, paste0("_l",i)))
    VAR_components[i] <- paste(append,collapse="+")
  }
  
  if(str_length(interaction_term))
    VAR_components[i+1] <- interaction_term
  
  VAR_components <- paste(VAR_components,collapse="+")
  
  if(no_population_intercept)
    VAR_components <- paste0( "0 + ", VAR_components)
  
  if(str_length(panel_id))
  {
    bf_string <- paste0("mvbind(",paste(response_variables,collapse=","),") ~ ", 
                        VAR_components,
                        if(length(exogenous_variables)) "+", 
                        if(length(exogenous_variables)) paste(exogenous_variables,collapse="+"),
                        " + (1|",panel_id,")")
  } else {
    bf_string <- paste0("mvbind(",paste(response_variables,collapse=","),") ~ ", 
                        VAR_components,
                        if(length(exogenous_variables)) "+", 
                        if(length(exogenous_variables)) paste(exogenous_variables,collapse="+"))
  }
  
  brm_mv <- bf(formula(bf_string)) + set_rescor(TRUE)  
  
  if(group_by_variant)
  {
    models           <- list()
    stan_code        <- list()
    model_draws      <- as_tibble()
    model_hypothesis <- as_tibble()
    
    for(variant_in in variants)
    {
      models[[variant_in]]    <- brm(brm_mv, data= data %>% filter(variant==variant_in), prior=priors, iter = iterations, control = control )
      stan_code[[variant_in]] <- stancode(models[[variant_in]])
      
      models[[variant_in]]    <- add_criterion(models[[variant_in]],"waic")
      models[[variant_in]]    <- add_criterion(models[[variant_in]],"loo")
      
      tmp         <- bind_rows(lapply(as_draws(models[[variant_in]]),as_tibble))
      tmp$variant <- variant_in
      model_draws <- bind_rows(model_draws,tmp)
      
      hh_tmp_greater0  <- as_tibble(hypothesis(models[[variant_in]], paste(rownames(fixef(models[[variant_in]]))[-1], "> 0"),class='b')$hypothesis)
      hh_tmp           <- as_tibble(hypothesis(models[[variant_in]], paste(rownames(fixef(models[[variant_in]]))[-1], "< 0"),class='b')$hypothesis) %>%
        rename(`Star<0`=Star)
      hh_tmp$'Star>0'  <- hh_tmp_greater0$Star
      hh_tmp$Variant   <- variant_in
      model_hypothesis <- bind_rows(model_hypothesis,hh_tmp)
    }
    
    model_draws <- model_draws %>% pivot_longer(-variant, names_to = "variable", values_to = "coef")
    model_draws$resp <- sapply(str_split(model_draws$variable,"_"),function(x)if( x[2] %in% str_replace(response_variables,"_","") )x[2] else NA)
    model_draws$NPI  <- sapply(str_split(model_draws$variable,"_"),function(x)if( length(x)>3 & is_character(x[4]) & (str_length(x[4])==2|str_length(x[4])==3) & x[4]!="GDP") x[4] else NA )
    
  } else {
    models    <- brm(brm_mv, data= data, prior=priors, iter = iterations, control = control ) 
    stan_code <- stancode(models)
    
    models    <- add_criterion(models,"waic")
    models    <- add_criterion(models,"loo")
    
    model_draws <- bind_rows(lapply(as_draws(models),as_tibble)) %>% mutate(period = "full") %>%
      pivot_longer(-period, names_to = "variable", values_to = "coef")
    model_draws$resp <- sapply(str_split(model_draws$variable,"_"),function(x)if( x[2] %in% str_replace(response_variables,"_","") )x[2] else NA)
    model_draws$NPI  <- sapply(str_split(model_draws$variable,"_"),function(x)if( length(x)>3 & is_character(x[4]) & (str_length(x[4])==2|str_length(x[4])==3) & x[4]!="GDP") x[4] else NA )
    
    hh_tmp_greater0  <- as_tibble(hypothesis(models, paste(rownames(fixef(models))[-1], "> 0"),class='b')$hypothesis)
    hh_tmp           <- as_tibble(hypothesis(models, paste(rownames(fixef(models))[-1], "< 0"),class='b')$hypothesis) %>%
      rename(`Star<0`=Star)
    hh_tmp$'Star>0'  <- hh_tmp_greater0$Star
    hh_tmp$Variant   <- "full"
    model_hypothesis <- hh_tmp
  }
  
  output <- list(models           = models,
                 stan_code        = stan_code,
                 model_hypothesis = model_hypothesis,
                 model_draws      = model_draws,
                 data_used        = data)
  
  if(!is.na(save_to_file)) saveRDS(output,paste0(here(),"/data/model_outputs/",save_to_file,".rds"))
  
  return(output)
}

# predict y(t+1) given estimated model and new data at y(t)
brm_predict <- function(brm_obj, new_data_in)
{
  ndi           <- new_data_in
  colnames(ndi) <- str_replace_all( colnames(new_data_in), "_", "" )
  resp          <- brm_obj$ranef$resp
  
  pred <- predict(brm_obj,newdata = new_data_in )
  
  forecast_error_change <- list()
  prediction_plot       <- list()
  rmse_calc             <- list()
  
  for(resp_in in resp)
  {
    tt   <- as_tibble(list(pred    = pred[,,resp_in][,1], 
                           resp    = ndi %>% pull(!!sym(resp_in)),
                           resp_l1 = ndi %>% pull(!!sym(paste0(resp_in,"l1")))))
    
    prediction_plot[[resp_in]]       <- tt %>% filter(!is.na(pred)&!is.na(resp)) %>% ggplot(aes(x=pred,y=resp)) + geom_point() + geom_point(aes(x=resp_l1,col="red")) + geom_smooth()
    rmse_calc[[resp_in]]             <- tt %>% mutate(resp_abs_err = abs(resp - pred), l1_abs_error = abs(resp - resp_l1)) %>%
      filter(!is.na(resp_abs_err)&!is.na(l1_abs_error) ) %>% 
      summarise(RMSE_pred = sqrt(mean(resp_abs_err^2)),
                RMSE_l1 = sqrt(mean(l1_abs_error^2)),
                MSE_pred = mean(resp_abs_err^2))
    forecast_error_change[[resp_in]] <- as.double((rmse_calc[[resp_in]]$RMSE_pred / rmse_calc[[resp_in]]$RMSE_l1 - 1)*100)
  }
  
  return(list(prediction_plot       = prediction_plot,
              rmse_calc             = rmse_calc,
              forecast_error_change = forecast_error_change,
              predictions           = pred ))
}

brm_NPI_sensitivity <- function(
    data,      
    start_date              = as.Date("2020-01-01"),
    end_date                = as.Date("2022-01-01"),
    response_variables      = c( "log_R","eco_ed","fd_GDP","fd_transit"),
    exogenous_variables     = c( "fd_C1", "lvl_C1L_lag", "fd_C2", "lvl_C2L_lag", "fd_C4", "lvl_C4L_lag", "fd_C5", "lvl_C5L_lag", "fd_C8", "lvl_C8L_lag", 
                                 "fd_H2", "lvl_H2L_lag", "fd_H6", "lvl_H6L_lag", "fd_H8", "lvl_H8L_lag", "fd_E1", "lvl_E1L_lag" ),
    p                       = 1,
    priors                  = c(prior(normal(0,1),class=b,resp=fdGDP),
                                prior(normal(0,1),class=b,resp=ecoed),
                                prior(normal(0,1),class=b,resp=fdtransit),
                                prior(normal(0,1),class=b,resp=logR),
                                prior(cauchy(0,2),class=sd,resp=fdGDP),
                                prior(cauchy(0,2),class=sd,resp=ecoed),
                                prior(cauchy(0,2),class=sd,resp=fdtransit),
                                prior(cauchy(0,2),class=sd,resp=logR)),
    panel_id                = "CountryName",
    iterations              = 4000,
    no_population_intercept = TRUE,
    control                 = list(adapt_delta = 0.9, step_size = 0.01, max_treedepth = 10),
    interaction_term        = "vaccination_doses_per_person + WT_variant + Alpha_variant + Delta_variant + Omicron_variant",
    save_to_file            = NA      #provide file name here
)
{
  data <- data %>% filter(Date>=start_date & Date<end_date) #enforce start & end date
  
  VAR_components <- c()
  VAR_eqns       <- c()
  
  for(i in 1:p)
    VAR_components <- c(VAR_components, unlist(lapply(response_variables,paste0, paste0("_l",i))))
  
  for(i in 1:length(response_variables))
  {
    VAR_eqns[i] <- paste(response_variables[i], "~" )  
    VAR_comp_in <- VAR_components[i]
    
    if(no_population_intercept)
    {
      VAR_comp_in <- paste0( "0 + ", paste(VAR_comp_in,collapse="+"))
    } else {
      VAR_comp_in <- paste(VAR_comp_in,collapse="+")
    }
    
    VAR_eqns[i] <- paste(VAR_eqns[i],VAR_comp_in, "+", paste(exogenous_variables,collapse="+"), "+" ,interaction_term,
                         " + (1|",panel_id,")")
  }
  
  brm_mv <- mvbf(VAR_eqns[1],VAR_eqns[2],VAR_eqns[3],VAR_eqns[4],rescor = FALSE)  
  
  models    <- brm(brm_mv, data= data, prior=priors, iter = iterations, control = control ) 
  stan_code <- stancode(models)
  
  models    <- add_criterion(models,"waic")
  models    <- add_criterion(models,"loo")
  
  model_draws <- bind_rows(lapply(as_draws(models),as_tibble)) %>% mutate(period = "full") %>%
    pivot_longer(-period, names_to = "variable", values_to = "coef")
  model_draws$resp <- sapply(str_split(model_draws$variable,"_"),function(x)if( x[2] %in% str_replace(response_variables,"_","") )x[2] else NA)
  model_draws$NPI  <- sapply(str_split(model_draws$variable,"_"),function(x)if( length(x)>3 & is_character(x[4]) & (str_length(x[4])==2|str_length(x[4])==3) & x[4]!="GDP") x[4] else NA )
  
  hypothesis_model <- as_tibble(hypothesis(models, paste(rownames(fixef(models)), "= 0"),class='b',alpha=0.05)$hypothesis)
  
  output <- list(models           = models,
                 stan_code        = stan_code,
                 model_hypothesis = hypothesis_model,
                 model_draws      = model_draws,
                 data_used        = data)
  
  if(!is.na(save_to_file)) saveRDS(output,paste0(here(),"/data/model_outputs/",save_to_file,".rds"))
  
  return(output)
}

brm_extract_Amat <- function(brm_obj,n_resp,p,get="Estimate")   # "Estimate", "l-95% CI" or "u-95% CI"
{
  coefs         <- summary(brm_obj)$fixed           # Get coefficients to construct A_i matrices
  coef_tbl      <- as_tibble(coefs)
  coef_tbl$term <- rownames(coefs)
  
  A_mat         <- array(NA,dim=c(n_resp,n_resp,p))
  
  for(i in 1:p)
  {
    tmp <- coef_tbl %>% filter(str_detect(term, paste0("_l",i)))
    tmp$response <- sapply(str_split(tmp$term,"_"),function(x)x[1] )
    tmp$variable <- sapply(str_split(tmp$term,"_"),function(x)x[3] )
    A_mat[,,i]   <- tmp %>% dplyr::select(!!sym(get), response, variable) %>% 
      pivot_wider(names_from = variable,values_from = !!sym(get)) %>%
      dplyr::select(-response) %>% as.matrix() %>% t()
  }
  
  return(A_mat)
}

# helper function to compute PCA loadings for NPIs
compute_NPI_PCA <- function(covid_data,NPI_variables,NPI_schema)
{
  pca.data           <- covid_data %>% dplyr::select(NPI_variables)
  colnames(pca.data) <- NPI_schema %>% filter(data_label %in% NPI_variables) %>% pull(Codebook)
  pca.res            <- prcomp(na.omit(pca.data),scale=TRUE)
  NPI_PCA1           <- as.matrix(pca.data) %*% as.vector(pca.res$rotation[,1])
  
  return(list(PCA1=NPI_PCA1,pca.res=pca.res))
}

compute_irf_mcmc <- function(brm_obj,p,lags=6,plot_labels=NA,facet_scales = "fixed",level_order = c("logR","ecoed","fdGDP","fdtransit"))
{
  rescov     <- VarCorr(brm_obj)                    # Extract Variance and Correlation Components 
  cov        <- as_tibble(rescov$residual__$cov)    # Get residual covariance terms
  cov_matrix <- cov %>% dplyr::select(colnames(cov)[str_starts(colnames(cov),"Estimate")]) %>% as.matrix()
  
  resp_names <- str_replace( colnames(cov_matrix), "Estimate.", "" )
  n_resp     <- length(resp_names)                  # number of response variables
  
  P_matrix   <- t(chol(cov_matrix))                 # Cholesky decomposition (lower triangular) of covariance matrix
  
  A <- brm_extract_Amat(brm_obj,n_resp,p,get="Estimate")
  
  A_mcmc <- array(NA,dim=c(n_resp,n_resp,p,1000))
  for(i in 1:n_resp)
  {
    for(j in 1:p)
      A_mcmc[i,,j,] <- t(mvrnorm(n=1000,mu=as.vector(A[i,,j]),Sigma=as.matrix(cov_matrix)))
  }
  
  oirf_tbl   <- as_tibble(x=NULL)
  girf_tbl   <- as_tibble(x=NULL)
  
  for(k in 1:1000)
  {
    A_mat    <- array(A_mcmc[,,,k],dim=c(n_resp,n_resp,p))
    
    phi      <- array(NA,dim=c(n_resp,n_resp,lags+1)) # Construction of the phi matrices up to lags 
    phi[,,1] <- diag(n_resp)  # this is phi 0
    phi[,,2] <- phi[,,1] %*% A_mat[,,1]
    
    for(i in 3:(lags+1))
    {
      phi[,,i] <- phi[,,i-1] %*% A_mat[,,1]
      if(p>1)
        for(j in 2:p)
          if( i > j )
            phi[,,i] <- phi[,,i] + phi[,,i-j] %*% A_mat[,,j]  
    }
    
    sigma_jj     <- diag(cov_matrix) 
    
    oirf_calced   <- array(NA,dim=c(n_resp,n_resp,lags+1))  # Construct orthogonal impulse response
    girf_calced   <- array(NA,dim=c(n_resp,n_resp,lags+1))
    
    for(i in 1:(lags+1))
    {
      oirf_calced[,,i]   <- phi[,,i] %*% P_matrix
      girf_calced[,,i]   <- phi[,,i] %*% diag(sigma_jj^(-1/2)) %*% cov_matrix
    }
    
    colnames(oirf_calced) <- colnames(girf_calced) <- resp_names
    rownames(oirf_calced) <- rownames(girf_calced) <- resp_names
    
    for(i in 1:(lags+1))                                         # put oirf into a tibble so that we can plot easily
    {
      tmp          <- as_tibble(oirf_calced[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      oirf_tbl     <- bind_rows(oirf_tbl,tmp) 
      
      tmp          <- as_tibble(girf_calced[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      girf_tbl     <- bind_rows(girf_tbl,tmp) 
    }
  }
  
  if(length(resp_names)==3)
  {
    tt <- girf_tbl %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),pctle="pct_0.5")
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),pctle="pct_0.025"))
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==4)
  {
    tt <- girf_tbl %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),fdtransit=median(!!sym(resp_names[4])),pctle="pct_0.5")
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),fdtransit=quantile(!!sym(resp_names[4]),0.025),pctle="pct_0.025"))
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),fdtransit=quantile(!!sym(resp_names[4]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==5)
  {
    tt <- girf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=median(!!sym(resp_names[3])),logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[4])),fdtransit=median(!!sym(resp_names[5])),pctle="pct_0.5")
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=quantile(!!sym(resp_names[3]),0.025),logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[4]),0.025),fdtransit=quantile(!!sym(resp_names[5]),0.025),pctle="pct_0.025"))
    tt <- bind_rows(tt, girf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=quantile(!!sym(resp_names[3]),0.975),logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[4]),0.975),fdtransit=quantile(!!sym(resp_names[5]),0.975),pctle="pct_0.975"))
  }
  
  
  girf      <- tt %>% pivot_longer(-c(lag,variable,pctle),names_to = "response", values_to = "gir") 
  
  girf$variable <- factor(as.character(girf$variable),levels = level_order)
  girf$response <- factor(as.character(girf$response),levels = level_order)
  
  if(length(plot_labels)==1)
  {
    levels(girf$variable) <- c("log~R","log~ED","Delta~GDP","Delta~Transit")
    levels(girf$response) <- c("log~R","log~ED","Delta~GDP","Delta~Transit")
  } else {
    levels(girf$variable) <- plot_labels
    levels(girf$response) <- plot_labels
  }
  
  girf_plot <- girf %>% group_by(response,variable,pctle) %>% ggplot(aes(x=lag,y=gir,col=pctle,linetype=pctle)) + geom_line() + 
    facet_grid(response~variable, labeller = label_parsed, scales = facet_scales ) + 
    theme(legend.position = "none" ) +
    scale_color_manual(values=c(pct_0.025="red",pct_0.5="black",pct_0.975="red")) +
    scale_linetype_manual(values = c(pct_0.025="dashed",pct_0.5="solid",pct_0.975="dashed")) +
    ylab("Generalised Impulse Response") + xlab("Lags")
  
  if(length(resp_names)==3)
  {
    tto <- oirf_tbl %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),pctle="pct_0.5")
    tto <- bind_rows(tto, oirf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),pctle="pct_0.025"))
    tto <- bind_rows(tto, oirf_tbl%>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==4)
  {
    tto <- oirf_tbl %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),fdtransit=median(!!sym(resp_names[4])),pctle="pct_0.5")
    tto <- bind_rows(tto, oirf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),fdtransit=quantile(!!sym(resp_names[4]),0.025),pctle="pct_0.025"))
    tto <- bind_rows(tto, oirf_tbl %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),fdtransit=quantile(!!sym(resp_names[4]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==5)
  {
    tto <- oirf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=median(!!sym(resp_names[3])),logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[4])),fdtransit=median(!!sym(resp_names[5])),pctle="pct_0.5")
    tto <- bind_rows(tto, oirf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=quantile(!!sym(resp_names[3]),0.025),logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[4]),0.025),fdtransit=quantile(!!sym(resp_names[5]),0.025),pctle="pct_0.025"))
    tto <- bind_rows(tto, oirf_tbl %>% group_by(variable,lag) %>% summarise(fdNPI=quantile(!!sym(resp_names[3]),0.975),logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[4]),0.975),fdtransit=quantile(!!sym(resp_names[5]),0.975),pctle="pct_0.975"))
  }
  
  oirf      <- tto %>% pivot_longer(-c(lag,variable,pctle),names_to = "response", values_to = "oir") 
  
  oirf$variable <- factor(as.character(oirf$variable),levels = level_order)
  oirf$response <- factor(as.character(oirf$response),levels = level_order)
  
  if(length(plot_labels)==1)
  {
    levels(oirf$variable) <- c("log~R","log~ED","Delta~GDP","Delta~Transit")
    levels(oirf$response) <- c("log~R","log~ED","Delta~GDP","Delta~Transit")
  } else {
    levels(oirf$variable) <- plot_labels
    levels(oirf$response) <- plot_labels
  }
  
  oirf_plot <- oirf %>% group_by(response,variable,pctle) %>% ggplot(aes(x=lag,y=oir,col=pctle,linetype=pctle)) + geom_line() + 
    facet_grid(response~variable, labeller = label_parsed, scales = facet_scales ) + 
    scale_color_manual(values=c(pct_0.025="red",pct_0.5="black",pct_0.975="red")) +
    scale_linetype_manual(values = c(pct_0.025="dashed",pct_0.5="solid",pct_0.975="dashed")) + 
    ylab("Orthogonolised Impulse Response") + xlab("Lags")
  
  return(list(oirf_plot=oirf_plot,oirf_tbl=oirf,
              girf_plot=girf_plot,girf_tbl=girf))
}