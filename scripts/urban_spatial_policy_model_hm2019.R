library(tidyverse)

# follows the model and policy experiments of Hsieh and Moretti (2019, AEJ: Macro).
# see tables 4 and 5 of the paper

alpha_param <- 0.65
eta_param <- 0.25
beta_param <- 0.4
rate <- 0.05
theta_param <- 0.3
land_param = 1 - alpha_param - eta_param

tfp_constant <- 
  (alpha_param^(1-eta_param))*((eta_param/rate)^eta_param)

treated_msas <- 
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
    local_tfp_64 = emp1964*(wage_condit_1964^((1-eta_param)/land_param))/tfp_constant,
    local_tfp_09 = emp2009*(wage_condit_2009^((1-eta_param)/land_param))/tfp_constant,
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

policy_df <- 
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

# policy experiment output table
policy_df |> 
  mutate(
    emp_diff = emp_policy - emp2009,
    agg_emp_diff = sum(emp_diff),
    rescale_factor = agg_emp_diff/sum(emp_policy), # hm do this rescaling move. is it valid?
    emp_policy_prime = (1-rescale_factor)*emp_policy,
    emp_share_prime = emp_policy_prime/sum(emp_policy_prime),
    price_ratio_policy = sum(emp_share_prime*wage_policy),
    agg_output_term_policy = sum( ((tfp_term_09/tfp_constant)*((price_ratio_policy/wage_policy)^(1-eta_param)))^(1/((1-eta_param)*(1+theta_param) - alpha_param)) ),
    agg_output_policy = ((eta_param/rate)^(eta_param/(1-eta_param)))*(agg_output_term_policy^(((1-eta_param)*(1+theta_param) - alpha_param)/(1-eta_param))),
    output_term_64 = sum( ((tfp_term_64/tfp_constant)*((price_avg_ratio_64/wage_condit_1964)^(1-eta_param)))^(1/((1-eta_param)*(1+theta_param) - alpha_param)) ),
    output_term_09 = sum( ((tfp_term_09/tfp_constant)*((price_avg_ratio_09/wage_condit_2009)^(1-eta_param)))^(1/((1-eta_param)*(1+theta_param) - alpha_param)) ),
    output_64 = ((eta_param/rate)^(eta_param/(1-eta_param)))*(output_term_64^(((1-eta_param)*(1+theta_param) - alpha_param)/(1-eta_param))),
    output_09 = ((eta_param/rate)^(eta_param/(1-eta_param)))*(output_term_09^(((1-eta_param)*(1+theta_param) - alpha_param)/(1-eta_param))),
    output_diff = (agg_output_policy - output_64)/(output_09 - output_64) - 1,
    welfare_diff = ((agg_output_policy/price_ratio_policy) - (output_64/price_avg_ratio_64))/((output_09/price_avg_ratio_09) - (output_64/price_avg_ratio_64)) - 1
  ) |>
  summarize(
    output_pct_diff = mean(output_diff)*100,
    welfare_pct_diff = mean(welfare_diff)*100
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
