### Plot functions ###

plot_npi_coefficients <- function(model, resp_names = c("log~ED", "Delta~GDP","Delta~Transit", "log~R") )
{
  hypothesis_model <- as_tibble(hypothesis(model$models, paste(rownames(fixef(model$models)), "= 0"),class='b',alpha=0.05)$hypothesis)
  
  npi_results_plt  <- hypothesis_model %>% # need to create Latex table from this
    mutate( Hypothesis = str_remove(Hypothesis,"[(]"),
            Hypothesis = str_remove(Hypothesis,"[)] = 0")) %>%                    # use subtables
    dplyr::select(-c(Evid.Ratio, Post.Prob,Star)) %>%
    dplyr::rename(variable=Hypothesis) %>%
    filter(!str_ends(variable,"Intercept")&!str_ends(variable,"l1")&!str_detect(variable,"SARS_C")&!str_detect(variable,"variant")&!str_detect(variable,"doses")) %>%
    mutate(Level = if_else(str_detect(variable,"L_")==TRUE,"NPI~Level","NPI~Change"),
           NPI   = stringr::str_split(variable, "_") %>% map_chr(., 3),
           resp  = stringr::str_split(variable, "_") %>% map_chr(., 1)) 
  
  npi_results_plt$resp         <- as.factor(npi_results_plt$resp)
  levels(npi_results_plt$resp) <- resp_names
  
  plot <- npi_results_plt %>%
    mutate(NPI=replace(NPI,NPI=="C1","Schools Closing"),
           NPI=replace(NPI,NPI=="C1L","Schools Closing"),
           NPI=replace(NPI,NPI=="C2","Workplace closure"),
           NPI=replace(NPI,NPI=="C2L","Workplace closure"),
           NPI=replace(NPI,NPI=="C4","Restrictions on gatherings"),
           NPI=replace(NPI,NPI=="C4L","Restrictions on gatherings"),
           NPI=replace(NPI,NPI=="C5","Close public transport"),
           NPI=replace(NPI,NPI=="C5L","Close public transport"),
           NPI=replace(NPI,NPI=="C8","International travel controls"),
           NPI=replace(NPI,NPI=="C8L","International travel controls"),
           NPI=replace(NPI,NPI=="E1","Income support"),
           NPI=replace(NPI,NPI=="E1L","Income support"),
           NPI=replace(NPI,NPI=="H2","Testing policy"),
           NPI=replace(NPI,NPI=="H2L","Testing policy"),
           NPI=replace(NPI,NPI=="H6","Facial Coverings"),
           NPI=replace(NPI,NPI=="H6L","Facial Coverings"),
           NPI=replace(NPI,NPI=="H8","Protection of elderly people"),
           NPI=replace(NPI,NPI=="H8L","Protection of elderly people")) %>%
    mutate(col=case_when(CI.Upper < 0 ~ "red",
                         CI.Lower > 0 ~ "blue",
                         CI.Upper > 0 & CI.Lower <0 ~ "black")) %>%
    ggplot(aes(y=NPI,x=Estimate,xmin = CI.Lower, xmax = CI.Upper,col=col)) + geom_pointinterval() + 
    geom_vline(aes(xintercept=0),linetype="dashed") +
    theme(strip.background=element_rect(fill="#CCFFFF"),
          strip.text=element_text(color="black"),
          legend.position = "none") +
    scale_color_manual(values = c("red" = "red",
                                  "blue"="blue",
                                  "black"="black")) +
    facet_grid(resp~Level,scales = "free", labeller = label_parsed) + 
    ylab("Non Pharmaceutical Intervention") + xlab("")  # 1000 x 800
  
  return( plot )
}

plot_var_coefficients <- function(model, panel_label = "A.")
{
  hypothesis_model <- as_tibble(hypothesis(model$models, paste(rownames(fixef(model$models)), "= 0"),class='b',alpha=0.05)$hypothesis)
  
  var_results_plt  <- hypothesis_model %>% # need to create Latex table from this
    mutate( Hypothesis = str_remove(Hypothesis,"[(]"),
            Hypothesis = str_remove(Hypothesis,"[)] = 0")) %>%                    # use subtables
    dplyr::select(-c(Evid.Ratio, Post.Prob,Star)) %>%
    dplyr::rename(variable=Hypothesis) %>%
    filter(!str_ends(variable,"Intercept")& str_ends(variable,"l1")) %>%
    mutate(Level  = if_else(str_detect(variable,"L_")==TRUE,"NPI~Level","NPI~Change"),
           resp_y = paste0(stringr::str_split(variable, "_") %>% map_chr(., 2),stringr::str_split(variable, "_") %>% map_chr(., 3)),
           resp   = stringr::str_split(variable, "_") %>% map_chr(., 1)) 
  
  var_results_plt$resp <- as.factor(var_results_plt$resp)
  levels(var_results_plt$resp) <- c("log~ED", "Delta~GDP","Delta~Transit", "log~R")
  var_results_plt$resp_y <- as.factor(var_results_plt$resp_y)
  levels(var_results_plt$resp_y) <- c("log~ED", "Delta~GDP","Delta~Transit", "log~R")
  
  var_results_plt <- var_results_plt %>% mutate(col=case_when(CI.Upper < 0 ~ "red",
                                                              CI.Lower > 0 ~ "blue",
                                                              CI.Upper > 0 & CI.Lower <0 ~ "black"),
                                                label  = paste( "$ ", str_replace(resp,"~"," "), ", ", str_replace(resp_y,"~"," "), " \\_lag 1 $")) 
  
  plot <- var_results_plt %>% 
    ggplot(aes(y=reorder(label,-desc(Estimate)),x=Estimate,xmin = CI.Lower, xmax = CI.Upper,col=col)) + geom_pointinterval() + 
    geom_vline(aes(xintercept=0),linetype="dashed") +
    theme(legend.position = "none" ) +
    scale_color_manual(values = c("red" = "red",
                                  "blue"="blue",
                                  "black"="black")) +
    scale_y_discrete(labels = TeX(var_results_plt$label)) + 
    ylab("") + xlab("") + ggtitle(panel_label)
  
  return(plot)
}