library(tidyverse)

# follows the model and policy experiments of Hsieh and Moretti (2019, AEJ: Macro).
# see tables 4 and 5 of the paper

params =
  c(0.65, # alpha
    0.25, # eta 
    0.4, # beta
    0.05, # rate R
    0.3   # theta (inverse) parameter
  )

alpha_param <- params[1]
eta_param <- params[2]
beta_param <- params[3]
rate <- params[4]
theta_param <- params[5]
land_param = 1 - alpha_param - eta_param

treated_msas = 
  c(
    5600,
    7400,
    7360
  )

hsieh_moretti_data <- 
  haven::read_dta("data/raw/hm_replication_files/data3.dta") %>%
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
    wage_condit_1964 = exp(logwage_condit64),
    wage_condit_2009  = (exp(logwage_condit09)/weighted.mean(exp(logwage_condit09), w = emp2009))*weighted.mean(wage_condit_1964, w = emp1964),
    WRLURI = ifelse( is.na(WRLURI) == T, st_avg_wr,WRLURI ),
    unaval = ifelse( is.na(unaval) == T, st_avg_unaval,unaval),
    populat_saiz = ifelse( is.na(populat_saiz) == T, st_avg_pop, populat_saiz ),
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
    tfp_exponent = ( beta_param*inverse + theta_param) / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param),
    amenity_exponent = land_param / ( beta_param*inverse + 1-alpha_param -eta_param - beta_param*inverse*eta_param +theta_param - eta_param*theta_param)
  )

policy_df = 
  hsieh_moretti_data |> 
  mutate(
    invelasticity_intercept = mean(inverse) - mean(-5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*WRLURI),
    mean_WRI = median(WRLURI),
    new_inver = 
      case_when(
        msa %in% treated_msas ~ invelasticity_intercept -5.260*unaval + .475*(log(populat_saiz)*unaval) + .280*((mean_WRI)),
        .default =  inverse
      ),
    tfp_exponent_policy = ( beta_param*new_inver + theta_param) / ( land_param + (beta_param*new_inver + theta_param)*(1 - eta_param) ),
    amenity_exponent_policy = land_param / ( land_param + (beta_param*new_inver + theta_param)*(1 - eta_param) ),
    q_term_wage = 
      wage_condit_2009/(tfp_term_09^tfp_exponent),
    q_term_labor = 
      ( emp2009/(tfp_term_09^(amenity_exponent/land_param)) )^(-1),
    q_raw_wage = 
      (wage_condit_2009/(tfp_term_09^tfp_exponent))^(1/amenity_exponent),
    q_raw_labor = 
      ( tfp_term_09/(emp2009^(1 - alpha_param - eta_param + (beta_param*inverse + theta_param)*(1-eta_param))) )^(1/(1-eta_param)) ,
    foc_labor = 
      (tfp_term_09/(wage_condit_2009^(1-eta_param)) )^(1/land_param),
    wage_policy = 
      (tfp_term_09^tfp_exponent_policy)*(q_raw_wage^amenity_exponent_policy),
    emp_policy = 
      (tfp_term_09/(wage_policy^(1-eta_param)) )^(1/land_param)
  )

# verification
#
#policy_df |> 
#  ggplot( aes( x = q_raw_labor, y = q_raw_wage ) ) +
#  geom_point() + 
#  theme_minimal()
#
#policy_df |> 
#  lm( formula = q_raw_wage ~ q_raw_labor ) |> summary()
