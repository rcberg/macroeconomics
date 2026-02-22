using Symbolics
using NonlinearSolve, LinearAlgebra
using Plots

include("general_equilibrium_function_1hh2good3fct.jl")

# the parameters and guesses need to be supplied by name in a tuple, like so
# these are just sample values for L = N_T = 1.0 
#parameter_vec = (
#    γ_L = 0.0134,
#    γ_N = 0.9284,
#    ρ_L = 0.3668,
#    ρ_N = 0.5773,
#    η_X = 0.6949,
#    η_Y = 0.1,
#    q = 1.0,
#    a_x = 1.0,
#    a_y = 1.0,
#    L = 1.0,
#    σ_D = 1.667,
#    σ_X = 1.667,
#    σ_Y = 1.667,
#    N_T = 1.0, 
#    ι = 1.0) 
#
#guess_vector = (
#    x = 1.0,
#    y = 1.0,
#    p = 1.0,
#    w = 1.0,
#    r = 1.0,
#    λ = 1.0,
#    L_x = 0.5*parameter_vec.L,
#    L_y = 0.5*parameter_vec.L,
#    N = (1 - parameter_vec.η_X - parameter_vec.η_Y)*parameter_vec.N_T,
#    N_x = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
#    N_y = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
#    K = 1.0,
#    K_x = 0.5,
#    K_y = 0.5 )
#
var_labs =  ["x" "y" "p" "w" "r" "λ" "L_x" "L_y" "N" "N_x" "N_y" "K" "K_x" "K_y"]

solution_container = []
convergence_container = []
did_solve = []
final_guess = []

# this exercise demonstrates symmetric scaling-up of the factors.
scale_max = 100.0
for i in 1.0:1.0:scale_max
    scale = i
    
    parameter_vec = (
        γ_L = 0.0134,
        γ_N = 0.9284,
        ρ_L = 0.3668,
        ρ_N = 0.5773,
        η_X = 0.6949,
        η_Y = 0.1,
        q = 1.0,
        a_x = 1.0,
        a_y = 1.0,
        L = scale,
        σ_D = 1.667,
        σ_X = 1.667,
        σ_Y = 1.667,
        N_T = scale, 
        ι = 1.0 )
    
    guess_vector = (
        x = scale,
        y = scale,
        p = 1.0,
        w = 1.0,
        r = 1.0,
        λ = 1.0,
        L_x = 0.5*parameter_vec.L,
        L_y = 0.5*parameter_vec.L,
        N = (1 - parameter_vec.η_X - parameter_vec.η_Y)*parameter_vec.N_T,
        N_x = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
        N_y = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
        K = scale,
        K_x = 0.5*scale,
        K_y = 0.5*scale )
    
    convergence_output, solution_output, out_success = economy_model( parameter_vec, guess_vector )

    push!(solution_container, solution_output)
    push!(convergence_container, convergence_output)
    push!(did_solve, out_success)
    if scale == scale_max
        push!(final_guess, guess_vector)
    end
end

solution_matrix = reduce(hcat, solution_container)
convergence_matrix = reduce(hcat, convergence_container)
Plots.plot(solution_matrix', xlabel = "Scale value", ylabel = "Solution value", label = var_labs)
# can check that prices obey homogeneity of degree one
Plots.plot(solution_matrix[3:5,:]', xlabel = "Scale value", ylabel = "Solution value", label = [var_labs[3] var_labs[4] var_labs[5]])

# look at the FOC residuals. kind of interesting patterns.
Plots.plot(convergence_matrix', xlabel = "Scale value", ylabel="FOC error", legend=:none)
Plots.plot(convergence_matrix[1:13,:]', xlabel = "Scale value", ylabel="FOC error", legend=:none)