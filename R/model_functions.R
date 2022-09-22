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

# helper function to compute PCA loadings for NPIs
compute_NPI_PCA <- function(covid_data,NPI_variables,NPI_schema)
{
  pca.data           <- covid_data %>% dplyr::select(NPI_variables)
  colnames(pca.data) <- NPI_schema %>% filter(data_label %in% NPI_variables) %>% pull(Codebook)
  pca.res            <- prcomp(na.omit(pca.data),scale=TRUE)
  NPI_PCA1           <- as.matrix(pca.data) %*% as.vector(pca.res$rotation[,1])
  
  return(list(PCA1=NPI_PCA1,pca.res=pca.res))
}