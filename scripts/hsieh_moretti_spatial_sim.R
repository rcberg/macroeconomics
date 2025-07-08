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

westcoast_hubs <- 
  c(1600,
    1920,
    3360,
    4480,
    5600,
    7360,
    7400
  )

brain_hubs = 
  c(
    5600,
    7400,
    7360
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
    theta, dataf, parameters, treatmsas
  ){
    
    
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
        south = ifelse( division %in% c(5, 6, 7), 1, 0),
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
        south_WR = weighted.mean(WRLURI, w = south),
        new_inver =  ifelse( treat == 1,
                             d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
                             inverse),
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
      mutate( 
        model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" )
      ) |> 
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
          south_WR, 
          hh64a, 
          hh09, 
          hh09a, 
          y64, 
          y09, 
          diff09_welf_adj 
        )
      )
    
    return(out.list)
  }

hm_elasticity_experiment_nat <- 
  function(
    theta, dataf, parameters, treatmsas
  ){
    
    
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
        south = ifelse( division %in% c(5, 6, 7), 1, 0),
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
        south_WR = weighted.mean(WRLURI, w = south),
        new_inver =  ifelse( treat == 1,
                             d3 -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
                             inverse),
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
      summarise(output_64 = mean(y64),
                output_act_09 = mean(y09),
                output_cf_09 = mean(yy09),
                change_in = treated_msa_names,
                policy_effect = 100*mean(diff09),
                welfare_effect  = 100*mean(diff09_welf_adj),
                mean_wage = weighted.mean(wage2009, w = emp2009),
                mean_cf_wage = weighted.mean(x, w = new_emp09c)) %>%
      mutate( model = ifelse( theta_param> 0, "Imperfect Mobility","Perfect Mobility" ))
    
    out.list = 
      out_data_hsieh_moretti
    
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

bind_rows(
  westcoast_experiment[[1]][[2]],
  westcoast_experiment[[2]][[2]]
)

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
