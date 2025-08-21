pacman::p_load(
  tidyverse,
  tidycensus,
  ipumsr
)

ddi = read_ipums_ddi("data/raw/usa_00088.xml")
census_2000_df = 
  read_ipums_micro(ddi) |> 
  janitor::clean_names() |> 
  mutate(
    immigrant = ifelse(bpl > 120, 1, 0)
  )

wager_df = 
  census_2000_df |> 
  filter( incwage > 0, age > 24, age < 56, uhrswork > 39 ) |> 
  mutate(
    hourly_wage = incwage/(uhrswork*wkswork1),
    log_wage = log(hourly_wage),
    wage_diff = log_wage - log(mean(hourly_wage)),
    pexp = age - educ - 6,
    educ_yrs = educ + 6,
    ind_var_1 = ifelse(ind1950/100 < 2, 1, 0),
    ind_var_2 = ifelse(ind1950 >= 2 & ind1950 < 3, 1, 0),
    ind_var_3 = ifelse(ind1950 >= 3 & ind1950 < 4, 1, 0),
    ind_var_4 = ifelse(ind1950 >= 4 & ind1950 < 5, 1, 0),
    ind_var_5 = ifelse(ind1950 >= 5 & ind1950 < 6, 1, 0),
    ind_var_6 = ifelse(ind1950 >= 6 & ind1950 < 7, 1, 0),
    ind_var_7 = ifelse(ind1950 >= 7 & ind1950 < 8, 1, 0),
    ind_var_8 = ifelse(ind1950 >= 8 & ind1950 < 9, 1, 0),
    ind_var_9 = ifelse(ind1950 >= 9, 1, 0),
    occ_var_1 = ifelse(occ1950/100 < 2, 1, 0),
    occ_var_2 = ifelse(occ1950 >= 2 & occ1950 < 3, 1, 0),
    occ_var_3 = ifelse(occ1950 >= 3 & occ1950 < 4, 1, 0),
    occ_var_4 = ifelse(occ1950 >= 4 & occ1950 < 5, 1, 0),
    occ_var_5 = ifelse(occ1950 >= 5 & occ1950 < 6, 1, 0),
    occ_var_6 = ifelse(occ1950 >= 6 & occ1950 < 7, 1, 0),
    occ_var_7 = ifelse(occ1950 >= 7 & occ1950 < 8, 1, 0),
    occ_var_8 = ifelse(occ1950 >= 8 & occ1950 < 9, 1, 0),
    occ_var_9 = ifelse(occ1950 >= 9, 1, 0),
    married = ifelse(marst < 3, 1, 0),
    separated = ifelse( marst == 3, 1, 0),
    divorced = ifelse( marst == 4, 1, 0),
    widowed = ifelse( marst == 3, 1, 0),
    vet_var = ifelse( vetstat == 2, 1, 0),
    black = ifelse( race == 2, 1, 0),
    indigenous = ifelse( race == 3, 1, 0),
    asian = ifelse( race %in% 4:6, 1, 0),
    hispanic = ifelse(hispan < 9 & hispan > 0, 1, 0),
    other_race = ifelse(race > 6, 1, 0),
    no_english = ifelse(speakeng == 1, 1, 0),
    low_english = ifelse( speakeng  == 6, 1, 0)
  )

houser_df = 
  census_2000_df |> 
  filter( pernum == 1, ownershp > 0 ) |> 
  mutate(
    renter = ifelse( ownershp == 2, 1, 0 ),
    imputed_rent = ifelse( ownershp == 1, 0.0785*valueh + costelec + costgas + costwatr + costfuel, rentgrs*12 ),
    log_rent = log(imputed_rent),
    fams_per_room = famsize/rooms,
    lotsize_sm = ifelse( acreprop == 2 | acreprop == 7, 1, 0 ),
    lotsize_lg = ifelse( acreprop == 8, 1, 0 ),
    old_build_10 = ifelse( builtyr < 4, 1, 0 ),
    old_build_20 = ifelse( builtyr == 4, 1, 0 ),
    old_build_30 = ifelse( builtyr == 5, 1, 0 ),
    old_build_40 = ifelse( builtyr == 6, 1, 0 ),
    old_build_50 = ifelse( builtyr == 7, 1, 0 ),
    old_build_60 = ifelse( builtyr == 8, 1, 0 ),
    old_build_61 = ifelse( builtyr == 9, 1, 0 ),
    has_plumbing = ifelse( plumbing > 19, 1, 0 ),
    has_kitchen = ifelse( kitchen > 2, 1, 0 ),
    commercial = ifelse( commuse == 2, 1, 0 ),
    condo_var = ifelse( condo == 2 & ownershp == 1, 1, 0 )
  ) |> 
  filter( imputed_rent > 0 )

wage_model_weights =
  wager_df |> 
  lm( 
    formula = 
      log_wage ~ 
        as.factor(educ) +
        pexp + I(pexp^2) + I(pexp^3) + I(pexp^4) + I(educ_yrs*pexp) +
        ind_var_1 + ind_var_2 + ind_var_3 + ind_var_4 + ind_var_5 + ind_var_6 + ind_var_7 + ind_var_8 + ind_var_9 + 
        occ_var_1 + occ_var_2 + occ_var_3 + occ_var_4 + occ_var_5 + occ_var_6 + occ_var_7 + occ_var_8 + occ_var_9 +
        married + separated + widowed + divorced +
        vet_var + I(vet_var*age) +
        black + hispanic + asian + indigenous + other_race + 
        I(immigrant*black) + I(immigrant*hispanic) + I(immigrant*asian) + I(immigrant*other_race) +
        no_english +
        low_english +
        #as.factor(met2013) +
        as.factor(sex):as.factor(educ) +
        as.factor(sex):pexp + as.factor(sex):I(pexp^2) + as.factor(sex):I(pexp^3) + as.factor(sex):I(pexp^4) + as.factor(sex):I(educ_yrs*pexp) +
        as.factor(sex):ind_var_1 + as.factor(sex):ind_var_2 + as.factor(sex):ind_var_3 + as.factor(sex):ind_var_4 + as.factor(sex):ind_var_5 + as.factor(sex):ind_var_6 + as.factor(sex):ind_var_7 + as.factor(sex):ind_var_8 + as.factor(sex):ind_var_9 + 
        as.factor(sex):occ_var_1 + as.factor(sex):occ_var_2 + as.factor(sex):occ_var_3 + as.factor(sex):occ_var_4 + as.factor(sex):occ_var_5 + as.factor(sex):occ_var_6 + as.factor(sex):occ_var_7 + as.factor(sex):occ_var_8 + as.factor(sex):occ_var_9 +
        as.factor(sex):married + as.factor(sex):separated + as.factor(sex):widowed + as.factor(sex):divorced +
        as.factor(sex):vet_var + as.factor(sex):I(vet_var*age) +
        as.factor(sex):black + as.factor(sex):hispanic + as.factor(sex):asian + as.factor(sex):indigenous + as.factor(sex):other_race + 
        as.factor(sex):I(immigrant*black) + as.factor(sex):I(immigrant*hispanic) + as.factor(sex):I(immigrant*asian) + as.factor(sex):I(immigrant*other_race) +
        as.factor(sex):no_english +
        as.factor(sex):low_english,
      weights = perwt
    )

wage_weights_df = 
  wage_model_weights |> 
  broom::augment( newdata = wager_df ) |> 
  mutate(  regression_weights = .fitted*perwt )

wage_model = 
  wage_weights_df |> 
  lm( 
    formula = 
      log_wage ~ 
        as.factor(educ) +
        pexp + I(pexp^2) + I(pexp^3) + I(pexp^4) + I(educ_yrs*pexp) +
        ind_var_1 + ind_var_2 + ind_var_3 + ind_var_4 + ind_var_5 + ind_var_6 + ind_var_7 + ind_var_8 + ind_var_9 + 
        occ_var_1 + occ_var_2 + occ_var_3 + occ_var_4 + occ_var_5 + occ_var_6 + occ_var_7 + occ_var_8 + occ_var_9 +
        married + separated + widowed + divorced +
        vet_var + I(vet_var*age) +
        black + hispanic + asian + indigenous + other_race + 
        I(immigrant*black) + I(immigrant*hispanic) + I(immigrant*asian) + I(immigrant*other_race) +
        yrsusa1 +
        no_english +
        low_english +
        as.factor(sex):as.factor(educ) +
        as.factor(sex):pexp + as.factor(sex):I(pexp^2) + as.factor(sex):I(pexp^3) + as.factor(sex):I(pexp^4) + as.factor(sex):I(educ_yrs*pexp) +
        as.factor(sex):ind_var_1 + as.factor(sex):ind_var_2 + as.factor(sex):ind_var_3 + as.factor(sex):ind_var_4 + as.factor(sex):ind_var_5 + as.factor(sex):ind_var_6 + as.factor(sex):ind_var_7 + as.factor(sex):ind_var_8 + as.factor(sex):ind_var_9 + 
        as.factor(sex):occ_var_1 + as.factor(sex):occ_var_2 + as.factor(sex):occ_var_3 + as.factor(sex):occ_var_4 + as.factor(sex):occ_var_5 + as.factor(sex):occ_var_6 + as.factor(sex):occ_var_7 + as.factor(sex):occ_var_8 + as.factor(sex):occ_var_9 +
        as.factor(sex):married + as.factor(sex):separated + as.factor(sex):widowed + as.factor(sex):divorced +
        as.factor(sex):vet_var + as.factor(sex):I(vet_var*age) +
        as.factor(sex):black + as.factor(sex):hispanic + as.factor(sex):asian + as.factor(sex):indigenous + as.factor(sex):other_race + 
        as.factor(sex):I(immigrant*black) + as.factor(sex):I(immigrant*hispanic) + as.factor(sex):I(immigrant*asian) + as.factor(sex):I(immigrant*other_race) +
        as.factor(sex):yrsusa1 +
        as.factor(sex):no_english +
        as.factor(sex):low_english +
        as.factor(metarea),
      weights = regression_weights
    )

housing_model_weights = 
  houser_df |> 
  filter( ownershp == 1 ) |> 
  lm(
    formula = 
    log_rent ~
      as.factor(unitsstr) +
      as.factor(rooms)*as.factor(bedrooms) +
      fams_per_room +
      lotsize_lg +
      lotsize_sm +
      old_build_10 +
      old_build_20 +
      old_build_30 +
      old_build_40 +
      old_build_50 +
      old_build_60 +
      old_build_61 +
      has_kitchen +
      has_plumbing +
      commercial +
      condo_var +
      as.factor(metarea),
    weights = hhwt
  )

house_weights_df = 
  housing_model_weights |> 
  broom::augment( newdata = houser_df ) |> 
  mutate(  regression_weights = .fitted*hhwt )

housing_model = 
  house_weights_df |> 
  lm(
    formula = 
    log_rent ~
      as.factor(unitsstr) +
      as.factor(rooms)*as.factor(bedrooms) +
      fams_per_room +
      lotsize_lg +
      lotsize_sm +
      old_build_10 +
      old_build_20 +
      old_build_30 +
      old_build_40 +
      old_build_50 +
      old_build_60 +
      old_build_61 +
      has_kitchen +
      has_plumbing +
      commercial +
      condo_var +
      renter:as.factor(unitsstr) +
      renter:as.factor(rooms)*as.factor(bedrooms) +
      I(renter*fams_per_room) +
      I(renter*lotsize_lg) +
      I(renter*lotsize_sm) +
      I(renter*old_build_10) +
      I(renter*old_build_20) +
      I(renter*old_build_30) +
      I(renter*old_build_40) +
      I(renter*old_build_50) +
      I(renter*old_build_60) +
      I(renter*old_build_61) +
      I(renter*has_kitchen) +
      I(renter*has_plumbing) +
      I(renter*commercial) +
      I(renter*condo_var) +
      as.factor(metarea),
    weights = regression_weights
  )

wage_effects = 
  broom::tidy(wage_model) |> 
  filter( grepl("metarea", term) ) |> 
  transmute(
    msa = str_sub( term, start = 19) |> as.numeric(),
    wage_differential = estimate
  )

housecost_effects = 
  broom::tidy(housing_model) |> 
  filter( grepl("metarea", term) ) |> 
  transmute(
    msa = str_sub( term, start = 19) |> as.numeric(),
    housecost_differential = estimate
  )

msa_differentials_df = 
  tibble( msa = 0,
  wage_differential = 0,
  housecost_differential = 0) |> 
  bind_rows(
    merge(wage_effects,
      housecost_effects)
    )

write_csv(msa_differentials_df, "data/export/albouy_stuart_2019_differentials.csv")
