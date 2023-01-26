library(pacman)

# load packages
pacman::p_load(tidyverse,
               tidytext,
               kableExtra,
               ggdist,
               ggrepel,
               viridis,
               lubridate,
               countrycode,
               patchwork,
               latex2exp,
               bayesplot,
               brms,
               here,
               zoo,
               MASS)

# initialize here functionality 
here()

# load functions
source( paste0(here(),"/R/utility_functions.R" ) )
source( paste0(here(),"/R/data_functions.R" ) )
source( paste0(here(),"/R/model_functions.R" ) )
source( paste0(here(),"/R/plot_functions.R" ) )

# load data
covid_data_weekly <- create_weekly_dataset()

# set up some initial settings
set.seed(1234)
options(mc.cores = 4, brms.backend = "cmdstanr")
theme_set(theme_light())
