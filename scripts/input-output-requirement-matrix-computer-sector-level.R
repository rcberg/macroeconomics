library(tidyverse)
# CAUTIONARY NOTE: this 'sector-level' I-O is amazing to practice with, but for many
#                  reasons is BAD to actually use! see reference texts.

# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output
use_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/use_before_redef_prod_price.xlsx",
    skip = 8
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Commodity Output`, ~replace_na(.,0))) 

use_matrix = 
  use_df |> 
  filter( !is.na(Code), !startsWith(Code,"V")) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
  as.matrix()

make_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/make.xlsx",
    skip = 9
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Industry Output`, ~replace_na(.,0))) 

make_matrix = 
  make_df |> 
  filter( !is.na(Code) ) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Noncomparable imports and rest-of-the-world adjustment /2/`) |> 
  as.matrix()

commodity_output_col = 
  make_df |> 
  filter( Name == "Total Commodity Output") |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Noncomparable imports and rest-of-the-world adjustment /2/`) |> 
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

com_comp_total_req_matrix = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/com-com_total_requirement.xlsx",
    skip = 6
  ) |> 
  filter( !is.na(Code) ) |> 
  select(!c(Code,`Commodity Description`)) |> 
  as.matrix()

com_error_matrix =  # commodity x commodity
  com_comp_total_req_matrix - com_total_req_matrix

ind_com_total_req_matrix = # industry x commodity
  transformation_matrix%*%com_total_req_matrix

wb_matrix = 
  transformation_matrix%*%direct_input_coef_matrix

reverse_direct_coefs.uninv = # industry x industry
 diag(nrow(wb_matrix)) - wb_matrix

ind_total_req_matrix = # industry x industry
  solve(reverse_direct_coefs.uninv)

ind_com_comp_total_req_matrix = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/ind-com_total_requirement.xlsx",
    skip = 6
  ) |> 
  filter( !is.na(Code) ) |> 
  select(!c(Code,`Industry Description`)) |> 
  as.matrix()

ind_com_error_matrix = 
  ind_com_comp_total_req_matrix - ind_com_total_req_matrix

ind_comp_total_req_matrix = 
    readxl::read_excel(
    "data/raw/input-output-tables/sector-level/ind-ind_total_requirement.xlsx",
    skip = 6
  ) |> 
  filter( !is.na(Code) ) |> 
  select(!c(Code,`Industry Description`)) |> 
  as.matrix()

ind_error_matrix = 
  ind_comp_total_req_matrix - ind_total_req_matrix

# plots of errors compared to official BEA total requirements matrix
hist(
  com_error_matrix,
  breaks = 30,
  xlab = "Deviation from BEA",
  ylab = "Number of matrix entries",
  main = "Commodity-by-commodity requirement table error"
)
text(0.02, 80, bquote(sigma == .(round(sd(as.vector(com_error_matrix)), 5 ))))

hist(
  ind_com_error_matrix,
  breaks = 30,
  xlab = "Deviation from BEA",
  ylab = "Number of matrix entries",
  main = "Industry-by-commodity requirement table error"
)
text(0.12, 70, bquote(sigma == .(round(sd(as.vector(ind_com_error_matrix)), 5 ))))


hist(
  ind_error_matrix,
  breaks = 30,
  xlab = "Deviation from BEA",
  ylab = "Number of matrix entries",
  main = "Industry-by-industry requirement table error"
)
text(0.0175, 42, bquote(sigma == .(round(sd(as.vector(ind_error_matrix)), 5 ))))
