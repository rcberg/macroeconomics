pacman::p_load(
  tidyverse
)

state_pce_df = 
  read_csv("data/raw/state_pce/SAPCE1__ALL_AREAS_1997_2023.csv") |> 
  mutate(
    across(`1997`:`2023`, .fns = ~.x*1000000),
    fips = as.numeric(GeoFIPS)/1000
  ) |> 
  filter( !is.na(fips), fips > 0, fips < 57) |> 
  pivot_longer(
    cols = `1997`:`2023`,
    names_to = 'year',
    values_to = 'value'
  ) |> 
  select( fips, GeoName, Description, year, value ) |> 
  pivot_wider(
    names_from = Description,
    values_from = value
  ) |> 
  janitor::clean_names() |> 
  rename('state' = 'geo_name')

state_total_income_df = 
  read_csv("data/raw/state_inc/SAINC1__ALL_AREAS_1929_2024.csv") |> 
  mutate(
    fips = as.numeric(GeoFIPS)/1000
  ) |> 
  filter( !is.na(fips), fips > 0, fips < 57, Description == "Personal income (millions of dollars)") |> 
  pivot_longer(
    cols = `1929`:`2024`,
    names_to = 'year',
    values_to = 'total_income'
  ) |> 
  select( fips, GeoName, year, total_income ) |> 
  rename( 'state' = 'GeoName')

state_disposable_income_df = 
  read_csv("data/raw/state_inc/SAINC51__ALL_AREAS_1948_2024.csv") |> 
  mutate(
    fips = as.numeric(GeoFIPS)/1000
  ) |> 
  filter( !is.na(fips), fips > 0, fips < 57, Description %in% c("Disposable personal income (millions of dollars)", "Population (persons) 1/") ) |> 
  pivot_longer(
    cols = `1948`:`2024`,
    names_to = 'year',
    values_to = 'value'
  ) |> 
  select( fips, GeoName, Description, year, value ) |> 
  pivot_wider(
    id_cols = c("fips","GeoName","year"),
    names_from = Description,
    values_from = value
  ) |> 
  janitor::clean_names() |> 
  rename( 
    'state' = 'geo_name',
    'disposable_income' = 'disposable_personal_income_millions_of_dollars',
    'population' = 'population_persons_1'
  )

state_income_df = 
  merge(state_total_income_df,
  state_disposable_income_df) |> 
  mutate( 
    across(total_income:disposable_income, .fns = ~.x*1000000),
    taxes = total_income - disposable_income,
    tax_share = taxes/total_income,
    total_inc_capita = total_income/population,
    disposable_inc_capita = disposable_income/population
 )

state_inc_consume_df = 
  select( 
    state_pce_df,
    fips,state,year,personal_consumption_expenditures,goods,services,housing_and_utilities
  ) |> 
  rename(
    'consumption_total' = 'personal_consumption_expenditures',
    'consumption_goods' = 'goods',
    'consumption_services' = 'services',
    'housing' = 'housing_and_utilities'
  ) |> 
  merge(
    state_income_df
  ) |> 
  mutate(
    housing_capita = housing/population,
    consumption_capita = consumption_total/population,
    housing_grossinc_share = housing/total_income,
    housing_dispoinc_share = housing/disposable_income,
    tax_capita = taxes/population
  )

#out_df = 
#  state_income_df |> 
#  merge( state_pce_df )
#
#write_csv(out_df, "data/export/state_consumption_income_1997_2023.csv")
