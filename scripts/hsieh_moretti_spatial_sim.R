library(tidyverse)

# follows the model and policy experiments of Hsieh and Moretti (2019, AEJ: Macro).
# see tables 4 and 5 of the paper

params <- 
  c(0.65, # alpha
    0.25, # eta 
    0.4, # beta
    0.05 # rate R
  )

brain_hubs <- 
  c( 5600, 7400, 7360 )

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
  c( 1123, 1600, 1640, 1840, 2080, 3480, 3760, 4120, 4480, 5080, 5120, 6160, 6200, 6440, 6780, 6920, 7160, 7320, 7600 )

westcoast_hubs <- 
  c( 1600, 1920, 3360, 4480, 5600, 7360, 7400 )

msa_list <- 
  read_csv("data/raw/hm_replication_files/hsieh_moretti_msa_labels.csv") %>%
  select( msa ) %>% 
  unlist() %>%
  as.numeric()

base_data <- 
  haven::read_dta("data/raw/hm_replication_files/data3.dta")

hm_elasticity_experiment <- 
  function(
    theta, dataf, parameters, treatmsas, output, experiment = "brains", label = T
  ){
    
    if( !(experiment %in% c("brains", "others")) ){
      stop(" please select 'brains' or 'others' for 'experiment.' ")
    }
    if( !(output %in% c("own", "nat")) ){
      stop(" please select 'own' or 'nat' for 'output.' ")
    }
    
    alpha_param <- parameters[1]
    eta_param <- parameters[2]
    beta_param <- parameters[3]
    rate <- parameters[4]
    land_param = 1 - alpha_param - eta_param
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
        price_avg_ratio_64 = sum(emp_share_64*wage_condit_1964),
        price_avg_ratio_09 = sum(emp_share_09*wage_condit_2009),
        tfp_exponent = ( beta_param*inverse + theta_param) / ( land_param + (beta_param*inverse + theta_param)*(1 - eta_param)),
        amenity_exponent = land_param / ( land_param + (beta_param*inverse + theta_param)*(1 - eta_param)),
        invelasticity_intercept = mean(inverse) - mean(-5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*WRLURI),
        mean_WRI = median(WRLURI),
        brain_WRI = mean(WRLURI[msa %in% c( 5600, 7400, 7360 )]),
        new_inver = 
          case_when(
          treat == 1 & experiment == "brains" ~ invelasticity_intercept -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
          treat == 1 & experiment == "others" ~ invelasticity_intercept -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((brain_WRI)),
          treat == 0 ~ inverse
        ),
        tfp_exponent_policy = ( beta_param*new_inver + theta_param) / ( land_param + (beta_param*new_inver + theta_param)*(1 - eta_param) ),
        amenity_exponent_policy = land_param / ( land_param + (beta_param*new_inver + theta_param)*(1 - eta_param) ),
        expon_diff_tfp = (tfp_exponent_policy - tfp_exponent),
        expon_diff_amenity = (amenity_exponent_policy - amenity_exponent), 
        log_wage_policy = log(wage_condit_2009) + expon_diff_tfp*log(tfp_term_09) - expon_diff_amenity*log_amenity_pm_09,
        wage_policy = exp(log_wage_policy),
        outout_arg_09_actual = adj_tfp_09*( (price_avg_ratio_09/wage_condit_2009)^((1-eta_param)/land_param)  ),
        agg_output_arg_09_actual = sum(outout_arg_09_actual),
        output_09_actual = (eta_param/rate)^(eta_param/(1-eta_param))*(agg_output_arg_09_actual^(land_param/(1-eta_param))  ),
        emp2009_policy =  ( ( ( alpha_param^(1-eta_param)*(eta_param^eta_param)/(rate^eta_param) )/(wage_policy^(1-eta_param)) )^(1/land_param) )*adj_tfp_09,
        new_emp2009_policy = (1 - sum(emp2009_policy - emp2009)/sum(emp2009_policy) )*emp2009_policy  ,
        emp_share_09_policy = (new_emp2009_policy / sum(new_emp2009_policy) ),
        diff_empc = emp_share_09_policy - (emp1964/sum(emp1964)),
        price_avg_ratio_64_policy = sum(emp_share_09_policy*wage_policy),
        output_arg_09_policy = adj_tfp_09*( (price_avg_ratio_64_policy/wage_policy)^((1-eta_param)/land_param)  ),
        agg_output_arg_09_policy = sum(output_arg_09_policy),
        output_09_policy = (eta_param/rate)^(eta_param/(1-eta_param))*(  agg_output_arg_09_policy^(land_param/(1-eta_param))  ),
        output_arg_64_actual = adj_tfp_64 * ( (price_avg_ratio_64/wage_condit_1964)^((1-eta_param)/land_param)  ),
        agg_output_arg_64_actual = sum(output_arg_64_actual),
        output_64_actual = (eta_param/rate)^(eta_param/(1-eta_param))*(agg_output_arg_64_actual^(land_param/(1-eta_param))  ),
        diff09 = 
          ( (output_09_policy - output_64_actual) / 
              (output_09_actual - output_64_actual) ) - 1,
        diff09_welf_adj=
          (( (output_09_policy/price_avg_ratio_64_policy ) - (output_64_actual/price_avg_ratio_64 ) ) / 
          ( (output_09_actual/price_avg_ratio_09 ) - (output_64_actual/price_avg_ratio_64 ) )) -1
      )
    
    out_sum_tbl <- 
      hsieh_moretti_data %>%
      summarise(
        policy_effect = 100*mean(diff09),
        welfare_effect  = 100*mean(diff09_welf_adj)
      )
    
    out_data <- 
      hsieh_moretti_data |> 
      filter( msa %in% treatmsas ) |> 
      select(
        -c(
          price_avg_ratio_64, 
          price_avg_ratio_64_policy, 
          price_avg_ratio_09, 
          mean_WRI, 
          agg_output_arg_64_actual, 
          agg_output_arg_09_policy, 
          agg_output_arg_09_actual, 
          output_64_actual, 
          output_09_actual, 
          diff09_welf_adj 
        )
      )
    
    if( output == "nat"){
      if(label == T){
        out_sum_tbl <- 
          out_sum_tbl |> 
          mutate( 
            model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" )
          )
      }
      
      return(out_sum_tbl)
    }else{
      if(label == T){
        out_data <- 
          out_data |> 
          mutate( 
            model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" )
          )
      }
      
      return(out_data)
    }
    
  }

hsieh_tbl_4 <- 
  bind_rows(
    map(
      0,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs,
      experiment = "brains",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "New York, San Fransisco, San Jose"),
    map(
      0,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = rust_belts,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "Rust Belt"),
    map(
      0,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = south_msas,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "South"),
    map(
      0,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = oth_large_msas,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "Other large cities")
  )

hsieh_tbl_5 <- 
  bind_rows(
    map(
      0.3,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs,
      experiment = "brains",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "New York, San Fransisco, San Jose"),
    map(
      0.3,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = rust_belts,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "Rust Belt"),
    map(
      0.3,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = south_msas,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |> 
      mutate( change_in = "South"),
    map(
      0.3,
      hm_elasticity_experiment,
      dataf = base_data, 
      parameters = params , 
      treatmsas = oth_large_msas,
      experiment = "others",
      output = 'nat',
      label = F
    )[[1]] |>
      mutate( change_in = "Other large cities")
  )

hsieh_tbl_4 |> View()
hsieh_tbl_5 |> View()

westcoast_experiment_nat <- 
  map(
    c(0,0.3),
    hm_elasticity_experiment,
    dataf = base_data, 
    parameters = params , 
    treatmsas = westcoast_hubs,
    output = 'nat'
  ) |> 
  bind_rows()

westcoast_experiment_own <- 
  map(
    c(0,0.3),
    hm_elasticity_experiment,
    dataf = base_data, 
    parameters = params , 
    treatmsas = westcoast_hubs,
    output = 'own'
  ) |> 
  bind_rows()
