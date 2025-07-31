library(tidyverse)

# follows the model and policy experiments of Hsieh and Moretti (2019, AEJ: Macro).
# see tables 4 and 5 of the paper

params <- 
  c(0.65, # alpha
    0.25, # eta 
    0.4, # beta
    0.05, # rate R
    1-0.65-0.25, # land exponent
    0.52 # tax parameter
  )

brain_hubs <- 
  c(
    5600,
    7400,
    7360
  )

rust_belts <- 
  c(
    2640, 
    7610, 
    5280, 
    4800, 
    9320, 
    1320, 
    4040, 
    6960, 
    6840, 
    3680, 
    3200, 
    1280, 
    4320, 
    8680, 
    2360, 
    960, 
    2000, 
    8400, 
    870, 
    80, 
    7040, 
    8160, 
    3720, 
    8320, 
    3000, 
    9140, 
    3920, 
    1680, 
    2040, 
    2160, 
    6280, 
    280, 
    3740, 
    7800, 
    1280, 
    3620, 
    6880, 
    2760, 
    8160
  )


south_msas <- 
  c(
    40,
    120,
    220,
    320,
    450,
    50,
    520,
    600,
    640,
    720,
    760,
    840,
    920,
    1000,
    1240,
    1440,
    1520,
    1540,
    1560,
    1760,
    1800,
    1880,
    1920,
    1950, 
    2020,
    2320,
    2560,
    2580,
    2680,
    2700,
    2750,
    2880,
    2900,
    2920,
    3120,
    3160, 
    3180, 
    3290, 
    3360, 
    3440, 
    3560, 
    3600, 
    3660, 
    3810, 
    3840, 
    3880, 
    3960, 
    3980, 
    4280, 
    4400,
    4420,
    4600,
    4640,
    4680,
    4880,
    4900,
    4920,
    5160,
    5200,
    5240,
    5360, 
    5560, 
    5720, 
    5790, 
    5800, 
    5880, 
    5960, 
    6080, 
    6640, 
    6760, 
    6800, 
    7240, 
    7510, 
    7520, 
    7680,
    8240,
    8280,
    8560,
    8600,
    8640,
    8800,
    8840,
    8960,
    9080,
    9160,
    9200
  )

oth_large_msas <- 
  c(
    1123, 
    1600, 
    1640, 
    1840, 
    2080, 
    3480, 
    3760, 
    4120, 
    4480, 
    5080, 
    5120, 
    6160, 
    6200, 
    6440, 
    6780, 
    6920, 
    7160, 
    7320, 
    7600
  )

westcoast_hubs <- 
  c(1600,
    1920,
    3360,
    4480,
    5600,
    7360,
    7400
  )

msa_list <- 
  read_csv("data/raw/hm_replication_files/hsieh_moretti_msa_labels.csv") %>%
  select( msa ) %>% 
  unlist() %>%
  as.numeric()

base_data <- 
  haven::read_dta("data/raw/hm_replication_files/data3.dta")

hm_elasticity_experiment <- 
  function(
    theta, dataf, parameters, treatmsas, experiment = "brains"
  ){
    
    if( !(experiment %in% c("brains", "others")) ){
      stop(" please select 'brains' or 'others' for 'experiment.' ")
    }
    
    alpha_param <- parameters[1]
    eta_param <- parameters[2]
    beta_param <- parameters[3]
    rate <- parameters[4]
    land_param <- parameters[5]
    tax_param <- parameters[6]
    theta_param <- theta
    
    treated_msa_names <- 
      read_csv("data/raw/hm_replication_files/hsieh_moretti_msa_labels.csv") %>%
      filter( msa %in% treatmsas ) %>%
      select( msa_name ) %>% 
      paste( collapse = ", ")
    
    hsieh_moretti_data <- 
      dataf %>%
      select( msa, 
              state, 
              division, 
              WRLURI, 
              unaval, 
              inverse, 
              populat_saiz, 
              emp1964, 
              emp2009, 
              wage1964, 
              wage2009, 
              logwage64, 
              logwage09, 
              logwage_condit64, 
              logwage_condit09, 
              HP64, 
              HP09 ) |> 
      group_by(state) |> 
      mutate( 
        st_avg_wr = mean(WRLURI, na.rm = T),
        st_avg_unaval = mean(unaval, na.rm = T),
        st_avg_pop = mean(populat_saiz, na.rm = T)
      ) |> 
      ungroup() |> 
      mutate( 
        treat = ifelse( msa %in% treatmsas, 1, 0),
        wage_condit_1964 = exp(logwage_condit64),
        wage_condit_2009  = (exp(logwage_condit09)/weighted.mean(exp(logwage_condit09), w = emp2009))*weighted.mean(wage_condit_1964, w = emp1964),
        logwage64 = log(wage_condit_1964),
        logwage09 = log(wage_condit_2009),
        WRLURI = ifelse( is.na(WRLURI) == T, st_avg_wr,WRLURI ),
        unaval = ifelse( is.na(unaval) == T, st_avg_unaval,unaval),
        populat_saiz = ifelse( is.na(populat_saiz) == T, st_avg_pop, populat_saiz ),
        adj_tfp_64 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp1964*wage_condit_1964^( (1-eta_param)/land_param),
        adj_tfp_09 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp2009*wage_condit_2009^( (1-eta_param)/land_param),
        tfp_term_64 =  ( (emp1964^land_param)*(wage_condit_1964^(1-eta_param)) ),
        tfp_term_09 =  ( (emp2009^land_param)*(wage_condit_2009^(1-eta_param)) ),
        amenity_pm_64 = mean(wage_condit_2009)*(1 + (.33 * ( (HP64 - mean(HP64))/mean(HP64) )) - (.51 * ( (wage_condit_1964 - mean(wage_condit_1964))/mean(wage_condit_1964) )) ), 
        amenity_pm_09 = mean(wage_condit_2009)*(1 + (.33 * ( (HP09 - mean(HP09))/mean(HP09) )) - (.51 * ( (wage_condit_2009 - mean(wage_condit_2009))/mean(wage_condit_2009) )) ), 
        log_amenity_pm_09 = log(amenity_pm_09),
        amenity_im_64 = mean(wage_condit_2009)*(1 + (.33 * ( (HP64 - mean(HP64))/mean(HP64) )) - (.51 * ( (wage_condit_1964 - mean(wage_condit_1964))/mean(wage_condit_1964) )) + ( theta_param*((emp1964-mean(emp1964))/mean(emp2009)) ) ),
        amenity_im_09 = mean(wage_condit_2009)*(1 + (.33 * ( (HP09 - mean(HP09))/mean(HP09) )) - (.51 * ( (wage_condit_2009 - mean(wage_condit_2009))/mean(wage_condit_2009) )) + ( theta_param*((emp2009-mean(emp2009))/mean(emp1964)) ) ),
        emp_share_64 = (emp1964 / sum(emp1964)),
        emp_share_09 = (emp2009 / sum(emp2009)),
        Pbar64 = sum(emp_share_64*wage_condit_1964),
        Pbar09 = sum(emp_share_09*wage_condit_2009),
        tfp_exponent = ( beta_param*inverse + theta_param)     / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        amenity_exponent = land_param / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        invelasticity_intercept = mean(inverse) - mean(-5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*WRLURI),
        mean_WRI = median(WRLURI),
        brain_WRI = mean(WRLURI[msa %in% c( 5600, 7400, 7360 )]),
        new_inver = 
          case_when(
            treat == 1 & experiment == "brains" ~ invelasticity_intercept -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
            treat == 1 & experiment == "others" ~ invelasticity_intercept -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((brain_WRI)),
            treat == 0 ~ inverse
          ),
        tfp_exponent_policy = ( beta_param*new_inver + theta_param) / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        amenity_exponent_policy = ( 1-alpha_param -eta_param)  / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        expon_diff_tfp = (tfp_exponent_policy - tfp_exponent),
        expon_diff_amenity = (amenity_exponent_policy - amenity_exponent), 
        log_wage_policy = logwage09 + expon_diff_tfp*log(tfp_term_09) - expon_diff_amenity*log_amenity_pm_09,
        wage_policy = exp(log_wage_policy),
        outout_arg_09_actual = adj_tfp_09*( (Pbar09/wage_condit_2009)^((1-eta_param)/land_param)  ),
        agg_output_arg_09_actual = sum(outout_arg_09_actual),
        output_09_actual = (eta_param/rate)^(eta_param/(1-eta_param))*(agg_output_arg_09_actual^(land_param/(1-eta_param))  )/1000000,
        emp2009_policy =  ( ( ( alpha_param^(1-eta_param)*(eta_param^eta_param)/(rate^eta_param) )/(wage_policy^(1-eta_param)) )^(1/land_param) )*adj_tfp_09,
        new_emp2009_policy = (1 - sum(emp2009_policy - emp2009)/sum(emp2009_policy) )*emp2009_policy  ,
        emp_share_09_policy = (new_emp2009_policy / sum(new_emp2009_policy) ),
        diff_empc = emp_share_09_policy - (emp1964/sum(emp1964)),
        Pbar64c = sum(emp_share_09_policy*wage_policy),
        output_arg_09_policy = adj_tfp_09*( (Pbar64c/wage_policy)^((1-eta_param)/land_param)  ),
        agg_output_arg_09_policy = sum(output_arg_09_policy),
        output_09_policy = (eta_param/rate)^(eta_param/(1-eta_param))*(  agg_output_arg_09_policy^(land_param/(1-eta_param))  )/1000000,
        output_arg_64_actual = adj_tfp_64 * ( (Pbar64/wage_condit_1964)^((1-eta_param)/land_param)  ),
        agg_output_arg_64_actual = sum(output_arg_64_actual),
        output_64_actual = (eta_param/rate)^(eta_param/(1-eta_param))*(agg_output_arg_64_actual^(land_param/(1-eta_param))  )/1000000,
        diff09 = 
          ( (output_09_policy - output_64_actual) / 
              (output_09_actual - output_64_actual) ) - 1,
        diff09_welf_adj=
          (( (output_09_policy/(Pbar64c/1e6) ) - (output_64_actual/(Pbar64/1e6) ) ) /  # hsieh renormalizes pbar by 1e6
             ( (output_09_actual/(Pbar09/1e6) ) - (output_64_actual/(Pbar64/1e6) ) )) -1
      )
    
    out_sum_tbl <- 
      hsieh_moretti_data %>%
      summarise(
        policy_effect = mean(diff09),
        welfare_effect  = mean(diff09_welf_adj),
        theta = theta_param
      )
    
    return(out_sum_tbl)
    
  }

theta_plot =
  map(
    seq(0,3,0.05),
    hm_elasticity_experiment,
    dataf = base_data,
    parameters = params,
    treatmsas = brain_hubs
  ) |> 
  bind_rows()

theta_plot |> 
  ggplot( aes( x = theta) ) +
  geom_line( aes( y = policy_effect, color = 'Output'), linewidth = 1 ) +
  geom_line( aes( y = welfare_effect, color = 'Welfare'), linewidth = 1) +
  scale_y_continuous( labels = scales::percent_format()) +
  scale_color_manual( values = c('Output' = 'black', 'Welfare' = 'blue') ) +
  labs( x = "Theta parameter", y = "Percent change", color = "" ) +
  theme_minimal()
