using NonlinearSolve
using CSV
using DataFrames
import Plots

include("ge-symbolic-model-function.jl")

msa_csv = CSV.read("data/raw/cbsa-population-density-2023.csv", DataFrame)
msa_csv.land_sq_mi = (msa_csv.land_sq_meters/2.59e6)
msa_csv.density = msa_csv.population ./ msa_csv.land_sq_mi # sq meters -> sq miles
lnth = length(msa_csv.met2013)
u0 = ones(15)

# initialize
pop = float(copy(msa_csv.population[1]))
land = float(copy(msa_csv.land_sq_mi[1]))
dens = pop/land
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
sol = solve(prob, maxiters = 1000)
sol_transformed = sol.u .^2
sol_mat = copy(sol_transformed)
sol_resid = copy(sol.resid)
sol_worked = if Symbol(sol.retcode) == :Success 1 else 0 end

for i=2:lnth
    pop = copy(msa_csv.population[i])
    land = copy(msa_csv.land_sq_mi[i])
    dens = pop/land
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
        q = 1.0, 
        a_x = 1.0, 
        a_y = 1.0)
        
        prob = NonlinearProblem(f!, u0, params)
        sol = solve(prob, maxiters = 1000)
        sol_transformed = sol.u .^2
        sol_mat = hcat(sol_mat, sol_transformed)
        sol_resid = hcat(sol_resid, sol.resid)
        sol_worked = vcat(sol_worked, if Symbol(sol.retcode) == :Success 1 else 0 end)
end

Plots.histogram(sol_mat[3,:], bins = 30) # histogram of model-solved land prices
sum(abs.(sol_resid))
sum(sol_worked)