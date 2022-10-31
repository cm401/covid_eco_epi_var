### Data functions ### 

# Create weekly dataset from daily data
create_weekly_dataset <- function(file_name='/data/covid_eu_data.rds')
{
  covid_data <- readRDS(paste0(here(),file_name))
  NPI_schema <- read_csv(paste0(here(),"/data/oxford_BSG.csv"))
  
  # create face 2 face sector group 
  face2face_group <- covid_data %>% dplyr::select(CountryCode,f2f_yes) %>% 
    filter(!is.na(f2f_yes)) %>% unique() %>% group_by(CountryCode) %>% 
    summarise(face2face_mean = mean(f2f_yes)) %>% 
    mutate(face2face_group = cat_fn_face_face(face2face_mean))
  
  # join face to face sector group data
  covid_data <- covid_data %>% left_join(face2face_group,by=c("CountryCode")) %>% 
    group_by(CountryName) %>% arrange(Date) %>%
    fill(starts_with("SARS_CoV_2_")|ends_with("_variant")) %>%
    fill(population,total_vaccinations) %>%
    ungroup()
  
  NPI_lag    <- 5 # 5 day lag for NPI level data
  
  covid_data <- covid_data %>% group_by(CountryName) %>% fill(f2f_yes) %>% 
    fill(estimated_daily_excess_deaths_per_100k) %>%
    mutate(
      lvl_C1L_lag = dplyr::lag(C1_combined_numeric,NPI_lag),
      lvl_C2L_lag = dplyr::lag(C2_combined_numeric,NPI_lag),
      lvl_C3L_lag = dplyr::lag(C3_combined_numeric,NPI_lag),
      lvl_C4L_lag = dplyr::lag(C4_combined_numeric,NPI_lag),
      lvl_C5L_lag = dplyr::lag(C5_combined_numeric,NPI_lag),
      lvl_C6L_lag = dplyr::lag(C6_combined_numeric,NPI_lag),
      lvl_C7L_lag = dplyr::lag(C7_combined_numeric,NPI_lag),
      lvl_C8L_lag = dplyr::lag(C8_combined_numeric,NPI_lag),
      lvl_H1L_lag = dplyr::lag(H1_combined_numeric,NPI_lag),
      lvl_H2L_lag = dplyr::lag(H2_combined_numeric,NPI_lag),
      lvl_H3L_lag = dplyr::lag(H3_combined_numeric,NPI_lag),
      lvl_H6L_lag = dplyr::lag(H6_combined_numeric,NPI_lag),
      lvl_H8L_lag = dplyr::lag(H8_combined_numeric,NPI_lag),
      lvl_E1L_lag = dplyr::lag(E1_combined_numeric,NPI_lag),
      lvl_E2L_lag = dplyr::lag(E2_combined_numeric,NPI_lag)
    ) %>% ungroup()
  
  # add PCA of NPIs
  NPI_variables   <- c("C1_combined_numeric", "C2_combined_numeric", "C3_combined_numeric", "C4_combined_numeric",     "C5_combined_numeric", "C6_combined_numeric", "C7_combined_numeric", "C8_combined_numeric",
                       "E1_combined_numeric", "E2_combined_numeric", "H2_combined_numeric", "H6_combined_numeric" )
  NPI_C_variables <- c("C1_combined_numeric", "C2_combined_numeric", "C3_combined_numeric", "C4_combined_numeric",     "C5_combined_numeric", "C6_combined_numeric", "C7_combined_numeric", "C8_combined_numeric" )
  NPI_E_variables <- c("E1_combined_numeric", "E2_combined_numeric" )
  NPI_H_variables <- c("H1_combined_numeric", "H2_combined_numeric", "H3_combined_numeric", "H6_combined_numeric", "H7_combined_numeric","H8_combined_numeric" )
  
  covid_data$PCA1_C   <- compute_NPI_PCA(covid_data,NPI_variables = NPI_C_variables,NPI_schema)$PCA1
  covid_data$PCA1_E   <- compute_NPI_PCA(covid_data,NPI_variables = NPI_E_variables,NPI_schema)$PCA1
  covid_data$PCA1_H   <- compute_NPI_PCA(covid_data,NPI_variables = NPI_H_variables,NPI_schema)$PCA1
  covid_data$PCA1_NPI <- compute_NPI_PCA(covid_data,NPI_variables = NPI_variables,NPI_schema)$PCA1
  
  #filter down to weekly data (using OECD GDP data to determine dates)
  covid_data_weekly <- covid_data %>% group_by(CountryName) %>% 
    fill(IC_median_ensamble_si_2_rt,f2f_yes,Region_v3) %>% filter(!is.na(`Tracker (counterfactual)`)) %>%
    mutate( 
      centre_R      = R - 1, 
      log_R         = log(lshtm_rt_median),
      fd_log_R      = log_R - dplyr::lag(log_R,1),
      fd_GDP        = (`Tracker (level)` - dplyr::lag(`Tracker (level)`,1) - prepandemic_trend_growth)/10,
      fd_GDP_l1     = dplyr::lag(fd_GDP,1),
      fd_GDP_l2     = dplyr::lag(fd_GDP,2),
      fd_GDP_l3     = dplyr::lag(fd_GDP,3),
      fd_GDP_l4     = dplyr::lag(fd_GDP,4),
      fd_GDP_l5     = dplyr::lag(fd_GDP,5),
      eco_ed_raw    = log(estimated_daily_excess_deaths_per_100k+1),
      eco_ed        = dplyr::lead(eco_ed_raw,3), # economist excess deaths
      eco_ed_l1     = dplyr::lead(eco_ed_raw,2), # economist excess deaths
      eco_ed_l2     = dplyr::lead(eco_ed_raw,1), # economist excess deaths
      eco_ed_l3     = eco_ed_raw, # economist excess deaths
      eco_ed_l4     = dplyr::lag(eco_ed_raw,1), # economist excess deaths
      eco_ed_l5     = dplyr::lag(eco_ed_raw,2), # economist excess deaths
      log_R_l1      = dplyr::lag(log_R,1),
      log_R_l2      = dplyr::lag(log_R,2),
      log_R_l3      = dplyr::lag(log_R,3),
      log_R_l4      = dplyr::lag(log_R,4),
      log_R_l5      = dplyr::lag(log_R,5),
      pct_eq_px     = PX_LAST/dplyr::lag(PX_LAST,1)-1,
      fd_cds_px     = CDS_PX_LAST - dplyr::lag(CDS_PX_LAST,1),
      fd_ensoe      = day_mean_quantity - dplyr::lag(day_mean_quantity,1),
      fd_NPI        = dplyr::lag(StringencyIndex - dplyr::lag(StringencyIndex,1),1),
      fd_GRI        = dplyr::lag(GovernmentResponseIndex - dplyr::lag(GovernmentResponseIndex,1),1),
      fd_CHI        = dplyr::lag(ContainmentHealthIndex - dplyr::lag(ContainmentHealthIndex,1),1),
      fd_ESI        = dplyr::lag(EconomicSupportIndex - dplyr::lag(EconomicSupportIndex,1),1),
      fd_PCN        = dplyr::lag(PCA1_NPI - dplyr::lag(PCA1_NPI,1),1),
      fd_PCC        = dplyr::lag(PCA1_C - dplyr::lag(PCA1_C,1),1),
      fd_PCE        = dplyr::lag(PCA1_E - dplyr::lag(PCA1_E,1),1),
      fd_PCH        = dplyr::lag(PCA1_H - dplyr::lag(PCA1_H,1),1),
      fd_C1         = dplyr::lag(C1_combined_numeric - dplyr::lag(C1_combined_numeric,1),1),
      fd_C2         = dplyr::lag(C2_combined_numeric - dplyr::lag(C2_combined_numeric,1),1),
      fd_C3         = dplyr::lag(C3_combined_numeric - dplyr::lag(C3_combined_numeric,1),1),
      fd_C4         = dplyr::lag(C4_combined_numeric - dplyr::lag(C4_combined_numeric,1),1),
      fd_C5         = dplyr::lag(C5_combined_numeric - dplyr::lag(C5_combined_numeric,1),1),
      fd_C6         = dplyr::lag(C6_combined_numeric - dplyr::lag(C6_combined_numeric,1),1),
      fd_C7         = dplyr::lag(C7_combined_numeric - dplyr::lag(C7_combined_numeric,1),1),
      fd_C8         = dplyr::lag(C8_combined_numeric - dplyr::lag(C8_combined_numeric,1),1),
      fd_E1         = dplyr::lag(E1_combined_numeric - dplyr::lag(E1_combined_numeric,1),1),
      fd_E2         = dplyr::lag(E2_combined_numeric - dplyr::lag(E2_combined_numeric,1),1),
      fd_H1         = dplyr::lag(H1_combined_numeric - dplyr::lag(H1_combined_numeric,1),1),
      fd_H2         = dplyr::lag(H2_combined_numeric - dplyr::lag(H2_combined_numeric,1),1),
      fd_H3         = dplyr::lag(H3_combined_numeric - dplyr::lag(H3_combined_numeric,1),1),
      fd_H6         = dplyr::lag(H6_combined_numeric - dplyr::lag(H6_combined_numeric,1),1),
      fd_H8         = dplyr::lag(H8_combined_numeric - dplyr::lag(H8_combined_numeric,1),1),
      fd_workplaces = dplyr::lag(workplaces - dplyr::lag(workplaces,1),1),
      fd_transit    = dplyr::lag((`transit stations` - dplyr::lag(`transit stations`,1))/100,1),
      fd_transit_l1 = dplyr::lag(fd_transit,1),
      fd_transit_l2 = dplyr::lag(fd_transit,2),
      fd_transit_l3 = dplyr::lag(fd_transit,3),
      fd_transit_l4 = dplyr::lag(fd_transit,4),
      fd_transit_l5 = dplyr::lag(fd_transit,5),
      vaccination_doses_per_person = replace_na(total_vaccinations,0) / population
    ) %>%
    ungroup() %>% rename(GDP_weekly_counterfactual = `Tracker (counterfactual)`) 
  
  # assign which part of the data is training vs testing data
  covid_data_weekly$VALIDATION_TYPE <- "TRAINING"
  for( comp in unique(covid_data_weekly$variant) )
    if(!is.na(comp))
    {
      n <- length(covid_data_weekly$variant==comp)
      covid_data_weekly$VALIDATION_TYPE[ runif(n) < 0.1 ] <- "TESTING"
    }
  
  return(covid_data_weekly)
}