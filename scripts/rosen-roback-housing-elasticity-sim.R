pacman::p_load(
  tidyverse,
  tidycensus
)

msa_elasticity_df = 
  read_csv("data/raw/rosen_roback_sim_2022data.csv") |> 
  select( msafips, saiz_elasticity, saiz_inv_elasticity, unaval, log_pop_1990, log_wrl ) |> 
  filter( !is.na(saiz_inv_elasticity)) |> 
  rename( 'cbsa_fips' = 'msafips')

usa_msa_df = 
  get_acs(
    geography = 'cbsa',
    variables = 
      c(
        'B01001_001',
        'B25002_002',
        'B25003_003',
        'B25088_001',
        'B20017_001',
        'B19013_001',
        'B25064_001',
        'B25077_001'
      ),
    output = 'wide',
    year = 2022
  ) |> 
  filter(grepl("Metro Area",NAME),
         !grepl(" PR Metro Area", NAME)) |> 
  rename(
    'cbsa_fips' = 'GEOID',
    'cbsa_name' = 'NAME',
    'population' = 'B01001_001E',
    'occ_h' = 'B25002_002E',
    'renter_occ_h' = 'B25003_003E',
    'earnings' = 'B20017_001E',
    'hh_income' = 'B19013_001E',
    'monthly_rent' = 'B25064_001E',
    'monthly_cost' = 'B25088_001E',
    'house_value' = 'B25077_001E'
  ) |> 
  mutate( yearly_rent = 12*monthly_rent,
          yearly_cost = 12*monthly_cost
  ) |> 
  select(
    cbsa_fips,
    cbsa_name,
    population,
    occ_h,
    renter_occ_h,
    earnings,
    hh_income,
    monthly_rent,
    monthly_cost,
    yearly_rent,
    yearly_cost,
    house_value
  )

city_msa_df = 
  merge( usa_msa_df, msa_elasticity_df ) |> 
  mutate( city = row_number())

rent_income_regression = 
  lm(
    data = city_msa_df,
    formula = log(hh_income) ~ log(yearly_cost)
  )

rent_income_invert_regression = 
  lm(
    data = city_msa_df, 
    formula = log(yearly_cost)~ log(hh_income)
  )


###############
# PARAMETERS
###############
# log_rent = r_int + a*log_pop
# log_hh_inc = k + g*log_rent - s

g = rent_income_regression$coefficients["log(yearly_cost)"] # housing utility parameter? not sure why this works so well
a =  city_msa_df$saiz_inv_elasticity
r_int = log(city_msa_df$yearly_cost) - a*log(city_msa_df$population) # rent intercept
n = nrow(city_msa_df)  # number of cities
L_max = round(sum(log(city_msa_df$population)))
w_baseline = log(city_msa_df$hh_income) # wage distribution
s = -rent_income_regression$residuals # amenity distribution

A = matrix( data = 0, nrow = n, ncol = n)

A[,1] = rep(1,n)

for(i in 1:n){
  if(i == n){
    A[i,] = c(1,rep(-g*a[i],(n-1)))
  }else{
    A[i,(i+1)] = g*a[i]
  }
}

y_baseline = w_baseline+s-g*r_int
y_baseline[n] = w_baseline[n] + s[n] - g*r_int[n] - g*a[n]*L_max

x_baseline = solve(A,y_baseline)
L_baseline = c(x_baseline[2:n], L_max-sum(x_baseline[2:n]) )

output_baseline = 
  data.frame(
    city = 1:n,
    wage = w_baseline,
    labor = L_baseline,
    rent = r_int+a*L_baseline,
    amenity = s,
    utility = w_baseline - g*a*L_baseline + s - g*r_int
  )

output_baseline$log_output = output_baseline$wage+output_baseline$labor
print(output_baseline)

ggplot() + 
  geom_point(data = output_baseline, 
             aes(
               x = rent, 
               y = wage),
             color = '#F8766D') +
  geom_point( data = city_msa_df,
              aes(
                x = log(yearly_cost), 
                y = log(hh_income)), 
              color = '#00BFC4') + 
  theme_minimal()

treated_cbsas = c(41940,41860,37100,41740,31080,42660,38900)
median_wrluri = median(city_msa_df$log_wrl)
treated_city = city_msa_df$city[ city_msa_df$cbsa_fips %in% treated_cbsas ]
a_policy = city_msa_df$saiz_inv_elasticity
a_policy[treated_city] = 
  -5.26*city_msa_df$unaval[treated_city] + 
  0.475*city_msa_df$unaval[treated_city]*city_msa_df$log_pop_1990[treated_city] +
  0.28*median_wrluri

A_policy = matrix( data = 0, nrow = n, ncol = n)

A_policy[,1] = rep(1,n)

for(i in 1:n){
  if(i == n){
    A_policy[i,] = c(1,rep(-g*a_policy[i],(n-1)))
  }else{
    A_policy[i,(i+1)] = g*a_policy[i]
  }
}

y_policy = w_baseline+s-g*r_int
y_policy[n] = w_baseline[n] + s[n] - g*r_int[n] - g*a_policy[n]*L_max

x_policy = solve(A_policy,y_policy)
L_policy = c(x_policy[2:n], L_max-sum(x_policy[2:n]) )

output_policy = 
  data.frame(
    city = 1:n,
    wage_policy = w_baseline,
    labor_policy = L_policy,
    rent_policy = r_int+a_policy*L_policy,
    amenity_policy = s,
    utility_policy = w_baseline - g*a_policy*L_policy + s - g*r_int
  )

output_policy$log_output_policy = output_policy$wage_policy+output_policy$labor_policy
print(output_policy)

ggplot() + 
  geom_point(data = output_policy, 
             aes(
               x = rent_policy, 
               y = wage_policy),
             color = '#F8766D') +
  geom_point( data = city_msa_df,
              aes(
                x = log(yearly_cost), 
                y = log(hh_income)), 
              color = '#00BFC4') + 
  theme_minimal()

policy_comparison = 
  merge(
    output_baseline,
    output_policy
  ) |> 
  mutate( policy_diff_pop = labor_policy - labor,
          policy_diff_utility = utility_policy - utility,
          policy_diff_output = log_output_policy - log_output,
          policy_diff_rent = rent_policy - rent,
          policy_city = if_else( city %in% treated_city, 1, 0) ,
          policy_diff_output_total = sum(policy_diff_output)
  )
