using LinearAlgebra
using XLSX, CSV, DataFrames
using Plots

ic_total_req = Matrix{Float64}( CSV.read("data/export/input-output-tables/ind_com_total_requirement_computed.csv", DataFrame) )
cc_total_req = Matrix{Float64}( CSV.read("data/export/input-output-tables/com_com_total_requirement_computed.csv", DataFrame) )
#bea_ic_total_req = Matrix{Float64}( XLSX.readdata("data/raw/input-output-tables/detail-level/ind-com_total_requirement.xlsx", "2017", "C6:ON407") )
#bea_cc_total_req = Matrix{Float64}( XLSX.readdata("data/raw/input-output-tables/detail-level/com-com_total_requirement.xlsx", "2017", "C6:ON407") )

final_demand = Matrix{Float64}( CSV.read("data/export/input-output-tables/final_uses_2017_adjusted.csv", DataFrame) )

# include 2 below to compare computed output with data
#ind_output_tot = XLSX.readdata("data/raw/input-output-tables/detail-level/use.xlsx", "2017","C414:ON414")
#ind_output = ind_output_tot[:,Not(279)] # "Customs Duties;" not a producing industry
com_labels = XLSX.readdata("data/raw/input-output-tables/detail-level/ind-com_total_requirement.xlsx","2017","C4:ON5")
ind_labels = XLSX.readdata("data/raw/input-output-tables/detail-level/ind-com_total_requirement.xlsx","2017","A6:B407")

total_req_eigs = eigvals(Matrix{Float64}(cc_total_req))
# BEA supplied total requirements matrix has insane eigenvalues, include these when looking at those
Plots.scatter(total_req_eigs, label = "Eigenvlaues")
#Plots.plot!(exp.(im*(0:0.001π:2π)), label = "Complex unit circle", linecolor = "red", linewidth = 2)

unit_demand_del = ones(size(ic_total_req,2))
unit_demand_req = ic_total_req*unit_demand_del
Plots.scatter(unit_demand_req, label = "Industry requrement per unit")
Plots.plot!(x -> 1.0, linecolor = "red", linewidth = 2, label = "=1")

output_computed = 
    ic_total_req*final_demand
#Plots.scatter(log.(output_computed[:]))
Plots.scatter(output_computed)
#comparison = Matrix{Float64}( hcat(transpose(ind_output), output_computed) )
#Plots.scatter( comparison[:,1], comparison[:,2])

final_uses_columns = 
    replace!(XLSX.readdata( "data/raw/input-output-tables/detail-level/use.xlsx", 
        "2017", "OP7:PH408" ), 
        missing => 0)[setdiff(1:end, (279, 399, 400, 401, 402)),:]
domestic_final_uses = 
    final_uses_columns[:,setdiff(1:end,7)]
final_uses_vec = 
    sum( domestic_final_uses, dims = 2 )

use_export_col = 
    final_uses_columns[:,7]

use_va_row = 
    replace!(XLSX.readdata( "data/raw/input-output-tables/detail-level/use.xlsx", 
        "2017", "C413:ON413" ), 
        missing => 0)[setdiff(1:end, 278)]

intermed_use_matrix = 
    replace!(XLSX.readdata( "data/raw/input-output-tables/detail-level/use.xlsx", 
        "2017", "C7:ON408" ), 
        missing => 0)[setdiff(1:end, (279, 399, 400, 401, 402)), setdiff(1:end, 278)]

use_matrix_shares = 
    copy(intermed_use_matrix)
use_va_shares = 
    use_va_row/sum(use_va_row)

import_row = 
    transpose(
        sum(replace!(XLSX.readdata( "data/raw/input-output-tables/detail-level/supply.xlsx", 
        "2017", "OP7:OQ408" ), 
        missing => 0)[setdiff(1:end, (279, 399, 400, 401, 402)),:], dims = 2)
    )

for i in axes(use_matrix_shares, 1)
    if sum(use_matrix_shares[i,:]) == 0 
        use_matrix_shares[i,:] = zeros( size(use_matrix_shares, 2) )
    else 
        use_matrix_shares[i,:] = use_matrix_shares[i,:]/sum(use_matrix_shares[i,:])
    end
end

intermed_make_matrix = 
    transpose(
        replace!(XLSX.readdata( "data/raw/input-output-tables/detail-level/supply.xlsx", 
        "2017", "C7:ON408" ), 
        missing => 0)[setdiff(1:end, (279, 399, 400, 401, 402)), setdiff(1:end, 278)]
    )

