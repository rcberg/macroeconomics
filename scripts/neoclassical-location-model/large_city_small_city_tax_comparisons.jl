using JuMP
using Ipopt
using Plots, StatsPlots

# Albouy and Stuart's parameterization
#
# 1. Model constants/calibration targets. (These I don't want changed)
const σ_D = 0.667
const σ_X = 0.667
const σ_Y = 0.667
const α = (σ_D - 1) / σ_D
const β = (σ_X - 1) / σ_X
const χ = (σ_Y - 1) / σ_Y
const i_bar = 1.0
const L_large = 1000.0
const L_small = L_large*1e-6
const N_large = 1000.0

# Targets from Albouy and Stuart's Table 1.
const target_sy = 0.36
const target_θL = 0.025 
const target_θN = 0.825
const target_ϕL = 0.233
const target_ϕN = 0.617

## These can change!
# AS (2019) replication baseline
λ_L = 0.17
λ_N = 0.7

τ_target = 0.0
share_match = "no"

include("large_city_small_city_simulator_taxsim.jl")
large_city_baseline = vcat([all_variables(large_city) value.(all_variables(large_city))][setdiff(1:size(all_variables(large_city),1), 16:18),:], [τ 0.0])
small_city_baseline = copy(meta_results)

λ_N_large = value(large_city_baseline[9,2])/(value(large_city_baseline[9,2]) + value(large_city_baseline[10,2]))
λ_L_large = value(large_city_baseline[11,2])/(value(large_city_baseline[11,2]) + value(large_city_baseline[12,2]))
# Trying out τ > 0
τ_target = 0.361
share_match = "τ"

include("large_city_small_city_simulator_taxsim.jl")
large_city_tax = [all_variables(large_city) value.(all_variables(large_city))][setdiff(1:size(all_variables(large_city),1), 16:18),:]
small_city_tax = copy(meta_results)
small_city_tax_plot = deepcopy(small_plot)

λ_N_small = value(large_city_tax[9,2])/(value(large_city_tax[9,2]) + value(large_city_tax[10,2]))
λ_L_small = value(large_city_tax[11,2])/(value(large_city_tax[11,2]) + value(large_city_tax[12,2]))

tax_diffs = (value.(large_city_tax[:,2]) .- value.(large_city_baseline[:,2]) ) ./ value.(large_city_tax[:,2])
#tax_diff_matrix = hcat(large_city_tax, tax_diffs)

comparison_matrix = hcat(large_city_tax,large_city_baseline[:,2])
comparison_matrix = hcat(comparison_matrix, tax_diffs)
bar(value.(comparison_matrix[:,4]),
    xticks = (1:length(string.(comparison_matrix[:,1])),string.(comparison_matrix[:,1])),
    labels=:none
)

bar(value.(comparison_matrix[1:20,4]),
    xticks = (1:length(string.(comparison_matrix[1:20,1])),string.(comparison_matrix[1:20,1])),
    labels=:none
)

q_compare = [collect(-0.2:0.01:0.2) small_city_baseline[1][:,2] small_city_tax[1][:,2]]
ax_compare = [collect(-0.2:0.01:0.2) (small_city_baseline[2][:,2]./(1 - target_sy)) small_city_tax[2][:,2]./(1 - target_sy)]
ay_compare = [collect(-0.2:0.01:0.2) small_city_baseline[3][:,2]./target_sy small_city_tax[3][:,2]./target_sy]

plot(q_compare[:,1], q_compare[:,2], title = "Quality of Life", label = "Baseline")
comp_plot_q = plot!(q_compare[:,1],q_compare[:,3], label = "Tax")

plot(ax_compare[:,1], ax_compare[:,2], title = "Traded productivity", label = "Baseline")
comp_plot_ax = plot!(ax_compare[:,1],ax_compare[:,3], label = "Tax")

plot(ay_compare[:,1], ay_compare[:,2], title = "Home productivity", label = "Baseline")
comp_plot_ay = plot!(ay_compare[:,1],ay_compare[:,3], label = "Tax")

dw, dh = default(:size).*1.5
using Plots.Measures
small_city_comp_plot = plot(comp_plot_q, comp_plot_ax, comp_plot_ay, layout = (3,1), size = (dw, dh * 3), leftmargin = 10mm)
small_city_comp_plot

# Un-comment to save plot.
#pth_out = pwd() * "\\reports\\plots\\as_small_city_tax_comp_plot.png"
#savefig(small_city_comp_plot, pth_out)