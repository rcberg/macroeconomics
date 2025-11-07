library(tidyverse)
# CAUTIONARY NOTE: this 'sector-level' I-O is amazing to practice with, but for many
#                  reasons is BAD to actually use! see reference texts.

# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output

ordinary_commodities = 
  c( "11","21","22","23","31G","42","44RT","48TW","51","FIRE","PROF","6","7","81","G" )

scale_fct = 1e6 # all units in table are millions of USD

use_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/use_before_redef_prod_price.xlsx",
    skip = 8
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Commodity Output`, ~replace_na(.,0))) 

use_matrix_total = 
  use_df |> 
  filter( Code %in% c(ordinary_commodities,"Other","Used") ) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
  as.matrix()

use_matrix = 
  use_df |> 
  filter( Code %in% ordinary_commodities ) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
  as.matrix()

make_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/make.xlsx",
    skip = 9
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Industry Output`, ~replace_na(.,0))) 

make_matrix_total = 
  make_df |> 
  filter( !is.na(Code) ) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Noncomparable imports and rest-of-the-world adjustment /2/`) |> 
  as.matrix()

make_matrix = 
  make_matrix_total[,!grepl("/",colnames(make_matrix_total)) ]

commodity_output_col = 
  make_df |> 
  filter( Name == "Total Commodity Output") |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
  unlist()

q_hat = 
  diag(commodity_output_col)
q_hat.inv =
  solve(q_hat)

industry_total_output_col = 
  make_df$`Total Industry Output`[!is.na(make_df$Code)]

scrap_col = 
  make_df$`Scrap, used and secondhand goods /1/`[!is.na(make_df$Code)]

industry_output_col = 
  industry_total_output_col - scrap_col

nonscrap_ratio =
  industry_output_col/industry_total_output_col

nonscrap_ratio_hat = diag(nonscrap_ratio)
nonscrap_ratio_hat.inv = solve(nonscrap_ratio_hat)

g_hat = diag(industry_total_output_col)
g_hat.inv = solve(g_hat)

direct_input_coef_matrix = 
  use_matrix%*%g_hat.inv

# industry technology assumption

mkt_share_matrix = 
  make_matrix%*%q_hat.inv

transformation_matrix = 
  nonscrap_ratio_hat.inv%*%mkt_share_matrix

direct_coef_matrix = 
  direct_input_coef_matrix%*%transformation_matrix

com_total_req_matrix_uninv = 
  diag(nrow(direct_coef_matrix)) - direct_coef_matrix
com_total_req_matrix = 
  solve(com_total_req_matrix_uninv) # this seems to work!

ind_com_total_req_matrix = # industry x commodity
  transformation_matrix%*%com_total_req_matrix

wb_matrix = 
  transformation_matrix%*%direct_input_coef_matrix

reverse_direct_coefs.uninv = # industry x industry
 diag(nrow(wb_matrix)) - wb_matrix

ind_total_req_matrix = # industry x industry
  solve(reverse_direct_coefs.uninv)

# validating the total requirements tables
# see last paragraph on p.'12-14' of reference text 1

total_final_uses = 
  use_df$`Total Final Uses (GDP)`[use_df$Code %in% ordinary_commodities]

tot_req_validate_df = 
  tibble(
    check_column = 
      c(
      as.numeric(com_total_req_matrix%*%total_final_uses), 
      as.numeric(ind_com_total_req_matrix%*%total_final_uses)
    ),
    actual_column = c(commodity_output_col, industry_output_col),
    error_rate_vs_actual = (check_column/actual_column) - 1,
    output_type = 
      c(
      rep("commodity", length(commodity_output_col)), 
      rep("industry", length(industry_output_col))
     )
  )

ggplot( tot_req_validate_df, aes( x = check_column, y = actual_column) ) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs( x = "Calculated from matrix ($millions)", y = "Actual data ($millions)") +
  scale_x_continuous(labels = scales::dollar_format()) +
  scale_y_continuous(labels = scales::dollar_format()) +
  theme_minimal()

ggplot( tot_req_validate_df, aes( x = error_rate_vs_actual, fill = output_type) ) + 
  geom_histogram( bins = 20, color = NA, position = 'dodge') + 
  labs( x = "Error rate (% from actual)", y = "Number of output categories", fill = "Output type") +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal()

lm( actual_column ~ check_column, tot_req_validate_df ) |> summary()

# plots of errors compared to official BEA total requirements matrix
# feel free to comment these out if not needed

#com_comp_total_req_matrix = 
#  readxl::read_excel(
#    "data/raw/input-output-tables/sector-level/com-com_total_requirement.xlsx",
#    skip = 6
#  ) |> 
#  filter( Code %in% ordinary_commodities ) |> 
#  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
#  as.matrix()
#
#com_error_matrix =  # commodity x commodity
#  com_comp_total_req_matrix - com_total_req_matrix
#
#ind_com_comp_total_req_matrix = 
#  readxl::read_excel(
#    "data/raw/input-output-tables/sector-level/ind-com_total_requirement.xlsx",
#    skip = 6
#  ) |> 
#  filter( !is.na(Code) ) |> 
#  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
#  as.matrix()
#
#ind_com_error_matrix = 
#  ind_com_comp_total_req_matrix - ind_com_total_req_matrix
#
#ind_comp_total_req_matrix = 
#    readxl::read_excel(
#    "data/raw/input-output-tables/sector-level/ind-ind_total_requirement.xlsx",
#    skip = 6
#  ) |> 
#  filter( !is.na(Code) ) |> 
#  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
#  as.matrix()
#
#ind_error_matrix = 
#  ind_comp_total_req_matrix - ind_total_req_matrix
#
#hist(
#  com_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Commodity-by-commodity requirement table error"
#)
#text(0.02, 80, bquote(sigma == .(round(sd(as.vector(com_error_matrix)), 5 ))))
#
#hist(
#  ind_com_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Industry-by-commodity requirement table error"
#)
#text(0.12, 70, bquote(sigma == .(round(sd(as.vector(ind_com_error_matrix)), 5 ))))
#
#hist(
#  ind_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Industry-by-industry requirement table error"
#)
#text(0.0175, 42, bquote(sigma == .(round(sd(as.vector(ind_error_matrix)), 5 ))))
#
#tibble(
#  computed = as.vector(ind_com_total_req_matrix),
#  bea = as.vector(ind_com_comp_total_req_matrix)
#) |> 
#  ggplot( aes( x = computed, y = bea)) + 
#  geom_point() + 
#  geom_abline(slope = 1, intercept = , linetype = 2) + 
#  theme_minimal()
#