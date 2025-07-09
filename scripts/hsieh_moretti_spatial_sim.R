library(tidyverse)

# follows the model and policy experiments of Hsieh and Moretti (2019, AEJ: Macro).
# see tables 4 and 5 of the paper

params =
  c(0.65, # alpha
    0.25, # eta 
    0.4, # beta
    0.05, # rate R
    1-0.65-0.25, # land exponent
    0.52 # tax parameter
  )

brain_hubs = 
  c(
    5600,
    7400,
    7360
  )

rust_belts = 
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


south_msas = 
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

oth_large_msas = 
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

base_data = 
  haven::read_dta("data/raw/hm_replication_files/data3.dta")

hm_elasticity_experiment_own <- 
  function(
    theta, dataf, parameters, treatmsas, experiment = "brains", label = T
  ){
    
    if( !(experiment %in% c("brains", "others")) ){
      stop("please select 'brains' or 'others'.")
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
    
    hsieh_moretti_states <- 
      dataf %>%
      select( msa, 
              state, 
              division, 
              WRLURI, 
              unaval, 
              elasticity, 
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
              HP09 ) %>%
      group_by( state ) %>%
      transmute( 
        msa = msa, 
        st_avg_wr = mean(WRLURI, na.rm = T),
        st_avg_unaval = mean(unaval, na.rm = T),
        st_avg_elasticity = mean(elasticity, na.rm = T),
        st_avg_pop = mean(populat_saiz, na.rm = T)
      )
    
    hsieh_moretti_data <- 
      dataf %>%
      select( msa, 
              state, 
              division, 
              WRLURI, 
              unaval, 
              elasticity, 
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
              HP09 ) %>%
      merge( hsieh_moretti_states ) %>%
      mutate( 
        treat = ifelse( msa %in% treatmsas, 1, 0),
        wage1964 = exp(logwage_condit64),
        wage2009 = exp(logwage_condit09),
        tmp64 = weighted.mean(wage1964, w = emp1964),
        tmp09 = weighted.mean(wage2009, w = emp2009),
        wage2009  = (wage2009/tmp09)*tmp64,
        logwage64 = log(wage1964),
        logwage09 = log(wage2009),
        WRLURI = ifelse( is.na(WRLURI) == T, st_avg_wr,WRLURI ),
        unaval = ifelse( is.na(unaval) == T, st_avg_unaval,unaval),
        elasticity = ifelse( is.na(elasticity) == T, st_avg_elasticity, elasticity ),
        populat_saiz = ifelse( is.na(populat_saiz) == T, st_avg_pop, populat_saiz ),
        AT64 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp1964*wage1964^( (1-eta_param)/land_param),
        AT09 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp2009*wage2009^( (1-eta_param)/land_param),
        Atfp64 =  ( (emp1964^land_param)*(wage1964^(1-eta_param)) ) ,
        Atfp09 =  ( (emp2009^land_param)*(wage2009^(1-eta_param)) ) ,
        tfp64 =  log(Atfp64) ,
        tfp09 =  log(Atfp09) ,
        ww64 = mean(wage1964),
        ww09 = mean(wage2009),
        pp64 = mean(HP64),
        pp09 = mean(HP09),
        Q64 = (.33 * ( (HP64 - pp64)/pp64 )) - (.51 * ( (wage1964 - ww64)/ww64 )),
        Q09 = (.33 * ( (HP09 - pp09)/pp09 )) - (.51 * ( (wage2009 - ww09)/ww09 )),
        QQ64 = ww09*(1+Q64), 
        QQ09 = ww09*(1+Q09), 
        logQQ09 = log(QQ09),
        B64 = (.33 * ( (HP64 - pp64)/pp64 )) - (.51 * ( (wage1964 - ww64)/ww64 )) + ( theta_param*((emp1964-mean(emp1964))/mean(emp2009)) ),
        B09 = (.33 * ( (HP09 - pp09)/pp09 )) - (.51 * ( (wage2009 - ww09)/ww09 )) + ( theta_param*((emp2009-mean(emp2009))/mean(emp1964)) ),
        BB64 = ww09*(1+B64),
        BB09 = ww09*(1+B09),
        P64 = wage1964,
        P09 = wage2009,
        e1964 = (emp1964 / sum(emp1964)),
        e2009 = (emp2009 / sum(emp2009)),
        Pbar64 = sum(e1964*wage1964),
        Pbar09 = sum(e2009*wage2009),
        P_Pbar64 = P64/Pbar64,
        P_Pbar09 = P09/Pbar09,
        c1 = ( beta_param*inverse + theta_param)     / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        c2 = land_param / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        d3 = mean(inverse) - mean(-5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*WRLURI),
        mean_WRI = median(WRLURI),
        brain_WRI = mean(WRLURI[msa %in% c( 5600, 7400, 7360 )]),
        new_inver = 
          case_when(
          treat == 1 & experiment == "brains" ~ d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
          treat == 1 & experiment == "others" ~ d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((brain_WRI)),
          treat == 0 ~ inverse
        ),
        h1 = ( beta_param*new_inver + theta_param) / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        h2 = ( 1-alpha_param -eta_param)  / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        diff1 = (h1 - c1),
        diff2 = (h2 - c2), 
        logx = logwage09 + diff1*log(Atfp09) - diff2*logQQ09,
        x = exp(logx),
        h09a = AT09*( (Pbar09/P09)^((1-eta_param)/land_param)  ),
        hh09a = sum(h09a),
        y09 = (eta_param/rate)^(eta_param/(1-eta_param))*(hh09a^(land_param/(1-eta_param))  )/1000000,
        emp09c =  ( ( ( alpha_param^(1-eta_param)*(eta_param^eta_param)/(rate^eta_param) )/(x^(1-eta_param)) )^(1/land_param) )*AT09,
        DDx = emp09c - emp2009,
        h = sum(DDx)/sum(emp09c),
        new_emp09c = (1-h)*emp09c  ,
        e2009c = (new_emp09c / sum(new_emp09c) ),
        diff_empc = e2009c - (emp1964/sum(emp1964)),
        Pbar64c = sum(e2009c*x),
        h09 = AT09*( (Pbar64c/x)^((1-eta_param)/land_param)  ),
        hh09 = sum(h09),
        yy09 = (eta_param/rate)^(eta_param/(1-eta_param))*(  hh09^(land_param/(1-eta_param))  )/1000000,
        h64a = AT64 * ( (Pbar64/P64)^((1-eta_param)/land_param)  ),
        hh64a = sum(h64a),
        y64 = (eta_param/rate)^(eta_param/(1-eta_param))*(hh64a^(land_param/(1-eta_param))  )/1000000,
        diff09 = ( (yy09 - y64) / (y09 - y64) ) - 1,
        Pbar09  = Pbar09/1000000,
        Pbar64c = Pbar64c/1000000,
        Pbar64  = Pbar64/1000000,
        diff09_welf_adj=( (yy09/(Pbar64c) ) - (y64/(Pbar64) ) ) /  ( (y09/(Pbar09) ) - (y64/(Pbar64) ) ) -1
      )
    
    out.list = 
      hsieh_moretti_data |> 
      filter( msa %in% treatmsas ) |> 
      select(
        -c(
          tmp64, 
          tmp09, 
          ww64, 
          ww09, 
          pp64, 
          pp09, 
          Pbar64, 
          Pbar64c, 
          Pbar09, 
          mean_WRI, 
          hh64a, 
          hh09, 
          hh09a, 
          y64, 
          y09, 
          diff09_welf_adj 
        )
      )
    
    if(label == T){
       out.list = 
         out.list |> 
         mutate( 
           model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" )
         )
    }
    
    return(out.list)
  }

hm_elasticity_experiment_nat <- 
  function(
    theta, dataf, parameters, treatmsas, experiment = "brains", label = T
  ){
    
    if( !(experiment %in% c("brains", "others")) ){
      stop("please select 'brains' or 'others'.")
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
    
    hsieh_moretti_states <- 
      dataf %>%
      select( msa, 
              state, 
              division, 
              WRLURI, 
              unaval, 
              elasticity, 
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
              HP09 ) %>%
      group_by( state ) %>%
      transmute( 
        msa = msa, 
        st_avg_wr = mean(WRLURI, na.rm = T),
        st_avg_unaval = mean(unaval, na.rm = T),
        st_avg_elasticity = mean(elasticity, na.rm = T),
        st_avg_pop = mean(populat_saiz, na.rm = T)
      )
    
    out_data_hsieh_moretti <- 
      dataf %>%
      select( msa, 
              state, 
              division, 
              WRLURI, 
              unaval, 
              elasticity, 
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
              HP09 ) %>%
      merge( hsieh_moretti_states ) %>%
      mutate( 
        treat = ifelse( msa %in% treatmsas, 1, 0),
        wage1964 = exp(logwage_condit64),
        wage2009 = exp(logwage_condit09),
        tmp64 = weighted.mean(wage1964, w = emp1964),
        tmp09 = weighted.mean(wage2009, w = emp2009),
        wage2009  = (wage2009/tmp09)*tmp64,
        logwage64 = log(wage1964),
        logwage09 = log(wage2009),
        WRLURI = ifelse( is.na(WRLURI) == T, st_avg_wr,WRLURI ),
        unaval = ifelse( is.na(unaval) == T, st_avg_unaval,unaval),
        elasticity = ifelse( is.na(elasticity) == T, st_avg_elasticity, elasticity ),
        populat_saiz = ifelse( is.na(populat_saiz) == T, st_avg_pop, populat_saiz ),
        AT64 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp1964*wage1964^( (1-eta_param)/land_param),
        AT09 = ( (alpha_param^(1-eta_param))*(eta_param^eta_param) / (rate^eta_param) )*emp2009*wage2009^( (1-eta_param)/land_param),
        Atfp64 =  ( (emp1964^land_param)*(wage1964^(1-eta_param)) ) ,
        Atfp09 =  ( (emp2009^land_param)*(wage2009^(1-eta_param)) ) ,
        tfp64 =  log(Atfp64) ,
        tfp09 =  log(Atfp09) ,
        ww64 = mean(wage1964),
        ww09 = mean(wage2009),
        pp64 = mean(HP64),
        pp09 = mean(HP09),
        Q64 = (.33 * ( (HP64 - pp64)/pp64 )) - (.51 * ( (wage1964 - ww64)/ww64 )),
        Q09 = (.33 * ( (HP09 - pp09)/pp09 )) - (.51 * ( (wage2009 - ww09)/ww09 )),
        QQ64 = ww09*(1+Q64), 
        QQ09 = ww09*(1+Q09), 
        logQQ09 = log(QQ09),
        B64 = (.33 * ( (HP64 - pp64)/pp64 )) - (.51 * ( (wage1964 - ww64)/ww64 )) + ( theta_param*((emp1964-mean(emp1964))/mean(emp2009)) ),
        B09 = (.33 * ( (HP09 - pp09)/pp09 )) - (.51 * ( (wage2009 - ww09)/ww09 )) + ( theta_param*((emp2009-mean(emp2009))/mean(emp1964)) ),
        BB64 = ww09*(1+B64),
        BB09 = ww09*(1+B09),
        P64 = wage1964,
        P09 = wage2009,
        e1964 = (emp1964 / sum(emp1964)),
        e2009 = (emp2009 / sum(emp2009)),
        Pbar64 = sum(e1964*wage1964),
        Pbar09 = sum(e2009*wage2009),
        P_Pbar64 = P64/Pbar64,
        P_Pbar09 = P09/Pbar09,
        c1 = ( beta_param*inverse + theta_param)     / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        c2 = land_param / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
        d3 = mean(inverse) - mean(-5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*WRLURI),
        mean_WRI = median(WRLURI),
        brain_WRI = mean(WRLURI[msa %in% c( 5600, 7400, 7360 )]),
        new_inver = 
          case_when(
            treat == 1 & experiment == "brains" ~ d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
            treat == 1 & experiment == "others" ~ d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((brain_WRI)),
            treat == 0 ~ inverse
          ),
        h1 = ( beta_param*new_inver + theta_param) / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        h2 = ( 1-alpha_param -eta_param)  / ( beta_param*new_inver + 1-alpha_param -eta_param - beta_param*new_inver*eta_param +theta_param - eta_param*theta_param),
        diff1 = (h1 - c1),
        diff2 = (h2 - c2), 
        logx = logwage09 + diff1*log(Atfp09) - diff2*logQQ09,
        x = exp(logx),
        h09a = AT09*( (Pbar09/P09)^((1-eta_param)/land_param)  ),
        hh09a = sum(h09a),
        y09 = (eta_param/rate)^(eta_param/(1-eta_param))*(hh09a^(land_param/(1-eta_param))  )/1000000,
        emp09c =  ( ( ( alpha_param^(1-eta_param)*(eta_param^eta_param)/(rate^eta_param) )/(x^(1-eta_param)) )^(1/land_param) )*AT09,
        DDx = emp09c - emp2009,
        h = sum(DDx)/sum(emp09c),
        new_emp09c = (1-h)*emp09c  ,
        e2009c = (new_emp09c / sum(new_emp09c) ),
        Pbar64c = sum(e2009c*x),
        h09 = AT09*( (Pbar64c/x)^((1-eta_param)/land_param)  ),
        hh09 = sum(h09),
        yy09 = (eta_param/rate)^(eta_param/(1-eta_param))*(  hh09^(land_param/(1-eta_param))  )/1000000,
        h64a = AT64 * ( (Pbar64/P64)^((1-eta_param)/land_param)  ),
        hh64a = sum(h64a),
        y64 = (eta_param/rate)^(eta_param/(1-eta_param))*(hh64a^(land_param/(1-eta_param))  )/1000000,
        diff09 = ( (yy09 - y64) / (y09 - y64) ) - 1,
        Pbar09  = Pbar09/1000000,
        Pbar64c = Pbar64c/1000000,
        Pbar64  = Pbar64/1000000,
        diff09_welf_adj=( (yy09/(Pbar64c) ) - (y64/(Pbar64) ) ) /  ( (y09/(Pbar09) ) - (y64/(Pbar64) ) ) -1
      ) %>%
      summarise(
                policy_effect = 100*mean(diff09),
                welfare_effect  = 100*mean(diff09_welf_adj)
                )
    
    if(label == T){
      out.list = 
        out_data_hsieh_moretti |> 
        mutate( 
          model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" )
        )
    }else{
      
      out.list = 
        out_data_hsieh_moretti
    }
    
    
    return(out.list)
  }

brainhub_experiment_nat <- 
  bind_rows(
    map(
      c(0,0.3),
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs
    )
  )

brainhub_experiment_own <- 
  bind_rows(
    map(
      c(0,0.3),
      hm_elasticity_experiment_own,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs
    )
  )

westcoast_experiment_own <- 
  map(
    c(0,0.3),
    hm_elasticity_experiment_own,
    dataf = base_data, 
    parameters = params , 
    treatmsas = westcoast_hubs
  ) |> 
  bind_rows()

westcoast_experiment_nat <- 
  map(
    c(0,0.3),
    hm_elasticity_experiment_nat,
    dataf = base_data, 
    parameters = params , 
    treatmsas = westcoast_hubs
  ) |> 
  bind_rows()

pdx_experiment_own = 
  map(
    c(0,0.3),
    hm_elasticity_experiment_own,
    dataf = base_data, 
    parameters = params , 
    treatmsas = 6440
  ) |> 
  bind_rows()

pdx_experiment_nat = 
  map(
    c(0,0.3),
    hm_elasticity_experiment_nat,
    dataf = base_data, 
    parameters = params , 
    treatmsas = 6440
  ) |> 
  bind_rows()

pnw_experiment_own = 
  map(
    c(0,0.3),
    hm_elasticity_experiment_own,
    dataf = base_data, 
    parameters = params , 
    treatmsas = c(6440, 7600)
  ) |> 
  bind_rows()

pnw_experiment_nat = 
  map(
    c(0,0.3),
    hm_elasticity_experiment_nat,
    dataf = base_data, 
    parameters = params , 
    treatmsas = c(6440, 7600)
  ) |> 
  bind_rows()

hsieh_tbl_4 = 
  bind_rows(
    map(
      0,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs,
      experiment = "brains",
      label = F
    )[[1]] |> 
      mutate( change_in = "New York, San Fransisco, San Jose"),
    map(
      0,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = rust_belts,
      experiment = "others",
      label = F
    )[[1]] |> 
      mutate( change_in = "Rust Belt"),
    map(
      0,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = south_msas,
      experiment = "others",
      label = F
    )[[1]] |> 
      mutate( change_in = "South"),
    map(
      0,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = oth_large_msas,
      experiment = "others",
      label = F
    )[[1]] |> 
      mutate( change_in = "Other large cities")
  )

hsieh_tbl_5 = 
  bind_rows(
    map(
      0.3,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = brain_hubs,
      experiment = "brains",
      label = F
    )[[1]] |> 
      mutate( change_in = "New York, San Fransisco, San Jose"),
    map(
      0.3,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = rust_belts,
      experiment = "others",
      label = F
    )[[1]] |> 
      mutate( change_in = "Rust Belt"),
    map(
      0.3,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = south_msas,
      experiment = "others",
      label = F
    )[[1]] |> 
      mutate( change_in = "South"),
    map(
      0.3,
      hm_elasticity_experiment_nat,
      dataf = base_data, 
      parameters = params , 
      treatmsas = oth_large_msas,
      experiment = "others",
      label = F
    )[[1]] |>
      mutate( change_in = "Other large cities")
  )

hsieh_tbl_4 |> View()
hsieh_tbl_5 |> View()
