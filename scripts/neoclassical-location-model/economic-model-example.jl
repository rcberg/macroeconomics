using NonlinearSolve
import Plots

#u0 = [1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 6.0, 4.0, 6.0, 4.0, 1.0, 9.0, 13.0, 100.0, 120.0, 1.0]
u0 = ones(15)
params = 
    (s_y = 0.36, 
    θ_l = 0.025, 
    θ_n = 0.825, 
    ϕ_l = 0.233, 
    ϕ_n = 0.617, 
    L_tot= 1.0, 
    N_tot = 10.0, 
    ρ = 0.1, 
    σ = 0.667, 
    τ = 0.361, 
    q = 1.0, 
    a_x = 1.0, 
    a_y = 1.0)

if params.σ == 1.0
    include("ge-symbolic-cobbd-model-function.jl")
else
    include("ge-symbolic-model-function.jl")
end

prob = NonlinearProblem(f!, u0, params)
sol = solve(prob, LevenbergMarquardt(), maxiters = 10000)
sol_transformed = sol.u .^2
Plots.plot(sol.resid)
#println(sol_transformed)
Plots.plot(sol_transformed)

# checks:
cmc_check = sol_transformed[13] + sol_transformed[14] ≈ sol_transformed[15] # capital market clearing
nmc_check = sol_transformed[9] + sol_transformed[10] ≈ params[7] # labor market clearing
lmc_check = sol_transformed[11] + sol_transformed[12] ≈ params[6] # land market clearing
hmc_check = sol_transformed[8] ≈ sol_transformed[6]*params[7] # home market clearing
expenditure_function_optim = sol_transformed[1]*( (1-params[1])^params[9] + (params[1]^params[9])*sol_transformed[4]^(1-params[9]) )^(1/(1-params[9]))
income = (sol_transformed[2] + sol_transformed[3]*params[6]/params[7] + params[8]*sol_transformed[15]/params[7])
tax = params[10]*(income)
cons_mc_check = sol_transformed[5] + sol_transformed[4]*sol_transformed[6] ≈ sol_transformed[2] + sol_transformed[3]*params[6]/params[7] + params[8]*sol_transformed[15]/params[7]
share_x = sol_transformed[5]/(sol_transformed[5] + sol_transformed[4]*sol_transformed[6]) # traded good expenditure share
share_y = sol_transformed[4]*sol_transformed[6]/(sol_transformed[5] + sol_transformed[4]*sol_transformed[6]) # home good expenditure share
ces_share_y = (params[1]^params[9])*(sol_transformed[4]^(1-params[9]))/( (1-params[1])^params[9] + (params[1]^params[9])*(sol_transformed[4]^(1-params[9])) )
ces_share_x = ((1-params[1])^params[9])/( (1-params[1])^params[9] + (params[1]^params[9])*(sol_transformed[4]^(1-params[9])) )
share_check_x = share_x ≈ ces_share_x # do expenditure shares match CES?
share_check_y = share_y ≈ ces_share_y


println("Capital market clears? ", cmc_check)
println("Labor market clears? ", nmc_check)
println("Land market clears? ", lmc_check)
println("Home market clears? ", hmc_check)
println("Income = expenditure? ", cons_mc_check)
println("Spending share on x correct? ", share_check_x)
println("Spending share on y correct? ", share_check_y)