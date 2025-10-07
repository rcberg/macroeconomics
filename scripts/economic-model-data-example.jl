using NonlinearSolve
using CSV
using DataFrames
import Plots

include("ge-model-function.jl")

msa_csv = CSV.read("data/raw/cbsa-population-density-2023.csv", DataFrame)
msa_csv.density = msa_csv.population ./ (msa_csv.land_sq_meters/2.59e6) # sq meters -> sq miles
lnth = length(msa_csv.met2013)

# initialize
dens = copy(msa_csv.density[1])
    params = 
        (s_y  = 0.36, 
        θ_l = 0.025, 
        θ_n = 0.825, 
        ϕ_l = 0.233, 
        ϕ_n = 0.617, 
        L_tot= 1.0, 
        N_tot = dens,
        ρ = 0.05, 
        σ = 0.667, 
        τ = 0.361, 
        q = 1, 
        a_x = 1, 
        a_y = 1)

prob = NonlinearProblem(f!, u0, params)
sol = solve(prob, maxiters = 10000)
sol_transformed = sol.u .^2
sol_mat = copy(sol_transformed)

u0 = ones(15)
for i=2:lnth
    dens = copy(msa_csv.density[i])
    params = 
        (s_y  = 0.36, 
        θ_l = 0.025, 
        θ_n = 0.825, 
        ϕ_l = 0.233, 
        ϕ_n = 0.617, 
        L_tot= 1.0, 
        N_tot = dens,
        ρ = 0.05, 
        σ = 0.667, 
        τ = 0.361, 
        q = 1, 
        a_x = 1, 
        a_y = 1)
        
        prob = NonlinearProblem(f!, u0, params)
        sol = solve(prob, maxiters = 10000)
        sol_transformed = sol.u .^2
        sol_mat = hcat(sol_mat, sol_transformed)
end

Plots.histogram(sol_mat[3,:], bins = 30) # histogram of model-solved land prices