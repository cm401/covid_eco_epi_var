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
                        VAR_components, "+", 
                        paste(exogenous_variables,collapse="+"),
                        " + (1|",panel_id,")")
  } else {
    bf_string <- paste0("mvbind(",paste(response_variables,collapse=","),") ~ ", 
                        VAR_components, "+", 
                        paste(exogenous_variables,collapse="+"))
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

# helper function to get VAR coefficients
brm_extract_Amat_from_draws <- function(draw,n_resp,p,resp_names)  
{
  df       <- draw
  df$head  <- 0
  df       <- df %>% dplyr::select(head,colnames(draw))
  coef_tbl <- df %>%
    gather(key = var_name, value = value, 2:ncol(df)) %>% 
    spread(key = names(df)[1],value = 'value')
  
  colnames(coef_tbl) <-c("term","value") 
  
  A_mat         <- array(NA,dim=c(n_resp,n_resp,p))
  
  for(i in 1:p)
  {
    tmp <- coef_tbl %>% filter(str_detect(term, paste0("_l",i)))
    tmp$response <- sapply(str_split(tmp$term,"_"),function(x)x[2] )
    tmp$variable <- sapply(str_split(tmp$term,"_"),function(x)paste0(x[3],x[4]) )
    A_mat[,,i]   <- tmp %>% dplyr::select(value, response, variable) %>% 
      pivot_wider(names_from = variable,values_from = value) %>% 
      dplyr::select(c("response",resp_names)) %>%
      arrange(match(response, resp_names)) %>%
      dplyr::select(-response) %>% as.matrix() %>% t()
  }
  
  return(A_mat)
}

# Compute Impulse Response Functions
compute_irf_from_draws <- function(brm_obj,p,lags=6,plot_labels=NA)
{
  rescov     <- VarCorr(brm_obj)                    # Extract Variance and Correlation Components 
  cov        <- as_tibble(rescov$residual__$cov)    # Get residual covariance terms
  cov_matrix <- cov %>% dplyr::select(colnames(cov)[str_starts(colnames(cov),"Estimate")]) %>% as.matrix()
  cov_mat_u  <- cov %>% dplyr::select(colnames(cov)[str_starts(colnames(cov),"Q97.5")]) %>% as.matrix()
  cov_mat_l  <- cov %>% dplyr::select(colnames(cov)[str_starts(colnames(cov),"Q2.5")]) %>% as.matrix()
  
  resp_names <- str_replace( colnames(cov_matrix), "Estimate.", "" )
  n_resp     <- length(resp_names)                  # number of response variables
  
  P_matrix   <- t(chol(cov_matrix))                 # Cholesky decomposition (lower triangular) of covariance matrix
  P_mat_u    <- t(chol(cov_mat_u))                 
  P_mat_l    <- t(chol(cov_mat_l)) 
  
  draws      <- bind_rows(lapply(as_draws(brm_obj),as_tibble))   
  
  oirf_tbl   <- as_tibble()
  girf_tbl   <- as_tibble()
  
  samples    <- round(runif(1000,1,8000)) 
  
  for(k in samples)#1:dim(draws)[1])
  {
    draw <- draws[k,]
    A_mat <- brm_extract_Amat_from_draws(draw,n_resp,p,resp_names)
    
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
    sigma_jj_u   <- diag(cov_mat_u) 
    sigma_jj_l   <- diag(cov_mat_l) 
    
    oirf_calced   <- array(NA,dim=c(n_resp,n_resp,lags+1))  # Construct orthogonal impulse response
    girf_calced   <- array(NA,dim=c(n_resp,n_resp,lags+1))
    oirf_calced_u <- array(NA,dim=c(n_resp,n_resp,lags+1))  # Construct orthogonal impulse response
    girf_calced_u <- array(NA,dim=c(n_resp,n_resp,lags+1))
    oirf_calced_l <- array(NA,dim=c(n_resp,n_resp,lags+1))  # Construct orthogonal impulse response
    girf_calced_l <- array(NA,dim=c(n_resp,n_resp,lags+1))
    
    for(i in 1:(lags+1))
    {
      oirf_calced[,,i]   <- phi[,,i] %*% P_matrix
      girf_calced[,,i]   <- phi[,,i] %*% diag(sigma_jj^(-1/2)) %*% cov_matrix
      oirf_calced_u[,,i] <- phi[,,i] %*% P_mat_u 
      girf_calced_u[,,i] <- phi[,,i] %*% diag(sigma_jj_u^(-1/2)) %*% cov_mat_u 
      oirf_calced_l[,,i] <- phi[,,i] %*% P_mat_l 
      girf_calced_l[,,i] <- phi[,,i] %*% diag(sigma_jj_l^(-1/2)) %*% cov_mat_l
    }
    
    colnames(oirf_calced) <- colnames(girf_calced) <- colnames(oirf_calced_u) <- colnames(girf_calced_u) <- colnames(oirf_calced_l) <- colnames(girf_calced_l) <- resp_names
    rownames(oirf_calced) <- rownames(girf_calced) <- rownames(oirf_calced_u) <- rownames(girf_calced_u) <- rownames(oirf_calced_l) <- rownames(girf_calced_l) <- resp_names
    
    for(i in 1:(lags+1))                                         # put oirf into a tibble so that we can plot easily
    {
      tmp          <- as_tibble(oirf_calced[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.5"
      oirf_tbl     <- bind_rows(oirf_tbl,tmp) 
      
      tmp          <- as_tibble(oirf_calced_u[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.975"
      oirf_tbl     <- bind_rows(oirf_tbl,tmp) 
      
      tmp          <- as_tibble(oirf_calced_l[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.025"
      oirf_tbl     <- bind_rows(oirf_tbl,tmp) 
      
      tmp          <- as_tibble(girf_calced[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.5"
      girf_tbl     <- bind_rows(girf_tbl,tmp) 
      
      tmp          <- as_tibble(girf_calced_u[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.975"
      girf_tbl     <- bind_rows(girf_tbl,tmp) 
      
      tmp          <- as_tibble(girf_calced_l[,,i])
      tmp$variable <- resp_names
      tmp$lag      <- i
      tmp$draw     <- k
      tmp$pctle    <- "pct_0.025"
      girf_tbl     <- bind_rows(girf_tbl,tmp)
    }
  }
  
  if(length(resp_names)==3)
  {
    tt <- girf_tbl %>% filter(pctle=="pct_0.5") %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),pctle="pct_0.5")
    tt <- bind_rows(tt, girf_tbl %>% filter(pctle=="pct_0.025") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),pctle="pct_0.025"))
    tt <- bind_rows(tt, girf_tbl %>% filter(pctle=="pct_0.975") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==4)
  {
    tt <- girf_tbl %>% filter(pctle=="pct_0.5") %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),fdtransit=median(!!sym(resp_names[4])),pctle="pct_0.5")
    tt <- bind_rows(tt, girf_tbl %>% filter(pctle=="pct_0.025") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),fdtransit=quantile(!!sym(resp_names[4]),0.025),pctle="pct_0.025"))
    tt <- bind_rows(tt, girf_tbl %>% filter(pctle=="pct_0.975") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),fdtransit=quantile(!!sym(resp_names[4]),0.975),pctle="pct_0.975"))
  }
  
  girf      <- tt %>% pivot_longer(-c(lag,variable,pctle),names_to = "response", values_to = "gir") 
  
  girf$variable <- as.factor(as.character(girf$variable))
  girf$response <- as.factor(as.character(girf$response))
  
  if(length(plot_labels)==1)
  {
    levels(girf$variable) <- c("log~ED", "Delta~GDP","Delta~Transit", "log~R")
    levels(girf$response) <- c("log~ED", "Delta~GDP","Delta~Transit", "log~R")
  } else {
    levels(girf$variable) <- plot_labels
    levels(girf$response) <- plot_labels
  }
  
  girf_plot <- girf %>% group_by(response,variable,pctle) %>% ggplot(aes(x=lag,y=gir,col=pctle,linetype=pctle)) + geom_line() + 
    facet_grid(response~variable, labeller = label_parsed) + 
    theme(legend.position = "none" ) +
    scale_color_manual(values=c(pct_0.025="red",pct_0.5="black",pct_0.975="red")) +
    scale_linetype_manual(values = c(pct_0.025="dashed",pct_0.5="solid",pct_0.975="dashed")) +
    ylab("Generalised Impulse Response") + xlab("Lags")
  
  if(length(resp_names)==3)
  {
    tto <- oirf_tbl %>% filter(pctle=="pct_0.5") %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),pctle="pct_0.5")
    tto <- bind_rows(tto, oirf_tbl %>% filter(pctle=="pct_0.025") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),pctle="pct_0.025"))
    tto <- bind_rows(tto, oirf_tbl %>% filter(pctle=="pct_0.975") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),pctle="pct_0.975"))
  } else if(length(resp_names)==4)
  {
    tto <- oirf_tbl %>% filter(pctle=="pct_0.5") %>% group_by(variable,lag) %>% summarise(logR=median(!!sym(resp_names[1])),ecoed=median(!!sym(resp_names[2])),fdGDP=median(!!sym(resp_names[3])),fdtransit=median(!!sym(resp_names[4])),pctle="pct_0.5")
    tto <- bind_rows(tto, oirf_tbl %>% filter(pctle=="pct_0.025") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.025),ecoed=quantile(!!sym(resp_names[2]),0.025),fdGDP=quantile(!!sym(resp_names[3]),0.025),fdtransit=quantile(!!sym(resp_names[4]),0.025),pctle="pct_0.025"))
    tto <- bind_rows(tto, oirf_tbl %>% filter(pctle=="pct_0.975") %>% group_by(variable,lag) %>% summarise(logR=quantile(!!sym(resp_names[1]),0.975),ecoed=quantile(!!sym(resp_names[2]),0.975),fdGDP=quantile(!!sym(resp_names[3]),0.975),fdtransit=quantile(!!sym(resp_names[4]),0.975),pctle="pct_0.975"))
  }
  
  oirf      <- tto %>% pivot_longer(-c(lag,variable,pctle),names_to = "response", values_to = "oir") 
  
  oirf$variable <- as.factor(as.character(oirf$variable))
  oirf$response <- as.factor(as.character(oirf$response))
  
  if(length(plot_labels)==1)
  {
    levels(oirf$variable) <- c("Excess~Deaths", "Delta~GDP","Delta~Transit", "log~R")
    levels(oirf$response) <- c("Excess~Deaths", "Delta~GDP","Delta~Transit", "log~R")
  } else {
    levels(oirf$variable) <- plot_labels
    levels(oirf$response) <- plot_labels
  }
  
  oirf_plot <- oirf %>% group_by(response,variable,pctle) %>% ggplot(aes(x=lag,y=oir,col=pctle,linetype=pctle)) + geom_line() + 
    facet_grid(response~variable, labeller = label_parsed) + 
    scale_color_manual(values=c(pct_0.025="red",pct_0.5="black",pct_0.975="red")) +
    scale_linetype_manual(values = c(pct_0.025="dashed",pct_0.5="solid",pct_0.975="dashed")) + 
    ylab("Orthogonolised Impulse Response") + xlab("Lags")
  
  return(list(oirf_plot=oirf_plot,oirf_tbl=oirf,
              girf_plot=girf_plot,girf_tbl=girf))
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