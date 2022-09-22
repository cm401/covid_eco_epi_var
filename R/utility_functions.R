### Utility functions ###

# Text formatting function
format_fn <- function(x) {
  formatC(x, digits = 3, format = "f", drop0trailing = FALSE)
}

# Convert week of year to date
year_week_to_date <- function(year_week)
{
  if(is.character(year_week)&length(year_week)==1)
  {
    split <- unlist(lapply(str_split(year_week,"-")[[1]],as.double))
    return(lubridate::ymd( paste0(split[1],"-01-01") ) + lubridate::weeks( split[2] - 1 ))  
  } else {
    split <- lapply(str_split(year_week,"-"),as.double)
    
    out_date <- c()
    
    for(comp in split)
      out_date <- c( out_date, as.character(lubridate::ymd( paste0(comp[1],"-01-01") ) + lubridate::weeks( comp[2] - 1 )))
    
    return(out_date)
  }
}

# create groups for face 2 face sectors
cat_fn_face_face <- function(x)
{
  if(length(x)>1)
  {
    res = c()
    
    for(z in x)
    {
      if(z>0.45) res = c(res,"f2f_high")
      else if(z>0.40) res = c(res,"f2f_medium_high")
      else if(z>0.35) res = c(res,"f2f_medium_low")
      else res = c(res,"f2f_low")
    }
  } else {
    res = c()
    if(x<0.35) res = c(res,"f2f_low")
    if(x<0.40) res = c(res,"f2f_medium_low")
    if(x<0.45) res = c(res,"f2f_medium_high")
    res = c(res,"f2f_high")
  }
  
  return(res)
}

# create SARS-CoV-2 time periods
get_covid_wave_dates <- function()
{
  return(list(
    "1_first_wave"   = seq(as.Date("2020-01-06"),as.Date("2020-07-01"),by="days"),
    "2_summer_2020"  = seq(as.Date("2020-07-02"),as.Date("2020-09-08"),by="days"),
    "3_alpha_wave"   = seq(as.Date("2020-09-09"),as.Date("2021-03-03"),by="days"),
    "4_delta_wave"   = seq(as.Date("2021-03-04"),as.Date("2021-07-06"),by="days"),
    "5_H2_2021"      = seq(as.Date("2021-07-07"),as.Date("2021-11-23"),by="days"),
    "6_omicron_wave" = seq(as.Date("2021-11-24"),as.Date("2022-02-01"),by="days")
  ))
}

# mechanically assign SARS-CoV-2 wave by date
get_covid_variant_dates <- function()
{
  return(list(
    "1_wildtype" = seq(as.Date("2020-01-06"),as.Date("2020-12-15"),by="days"),
    "2_alpha"    = seq(as.Date("2020-12-16"),as.Date("2021-05-14"),by="days"),
    "3_delta"    = seq(as.Date("2021-05-15"),as.Date("2021-12-12"),by="days"),
    "4_omicron"  = seq(as.Date("2021-12-13"),Sys.Date(),by="days")
  ))
}

# create SARS-CoV-2 time periods by quarters / year
get_covid_alternat_dates <- function()
{
  return(list(
    "1_phase"  = seq(as.Date("2020-03-01"),as.Date("2020-08-31"),by="days"),
    "2_phase"  = seq(as.Date("2020-09-01"),as.Date("2021-02-28"),by="days"),
    "3_phase"  = seq(as.Date("2021-03-01"),as.Date("2021-08-31"),by="days"),
    "4_phase"  = seq(as.Date("2021-09-01"),Sys.Date(),by="days")
  ))
}

# add dates to the covid data set
add_dates_to_covid_data <- function(data_in)
{
  wave_dates     <- get_covid_wave_dates()
  variant_dates  <- get_covid_variant_dates()
  alternat_dates <- get_covid_alternat_dates()
  
  data_in$wave <- NA
  for(wave in names(wave_dates))
    data_in[data_in$Date %in% wave_dates[[wave]], ]$wave <- wave 
  
  data_in$variant <- NA
  for(variant in names(variant_dates))
    data_in[data_in$Date %in% variant_dates[[variant]], ]$variant <- variant 
  
  seasons <- list(spring=c(3,4,5), summer=c(6,7,8),autumn=c(9,10,11),winter=c(12,1,2))
  
  data_in$season <- NA
  for(season in names(seasons))
    data_in[month(data_in$Date) %in% seasons[[season]], ]$season <- season 
  
  data_in$alternate_date <- NA
  for(ad in names(alternat_dates))
    data_in[data_in$Date %in% alternat_dates[[ad]], ]$alternate_date <- ad
  
  return(data_in)
}

# create pdata structure for panel models
create_pdata          <- function(x)
{ 
  return(pdata.frame(x,index=c('CountryName','Date'), drop.index = TRUE, row.names = TRUE)) 
}
