using Symbolics
using NonlinearSolve, LinearAlgebra
using Plots

include("general_equilibrium_function_1hh2good3fct.jl")

var_labs =  ["x" "y" "p" "w" "r" "λ" "L_x" "L_y" "N" "N_x" "N_y" "K" "K_x" "K_y"]
foc_labs = ["walras_law" "l_mkt_clearing" "n_mkt_clearing" "k_mkt_clearing" "π_x" "π_y" "π_total" "λ" "∂x/∂N_x - w" "p*∂y/∂N_y - w" "∂x/∂L_x - r" "p*∂y/∂L_y - r" "∂x/∂K_x - ι" "p*∂y/∂K_y - ι"]
scale_max = 100.0

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

# this function demonstrates symmetric scaling-up of the factors.
function symmetric_scale( scalemax )
    solution_container = []
    convergence_container = []
    did_solve = []

    for i in 1.0:scalemax
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
    end
    out_solution = reduce(hcat, solution_container)
    out_convergence = reduce(hcat, convergence_container)

    return out_solution, out_convergence
end

symm_solution_matrix, symm_convergence_matrix = symmetric_scale(scale_max)

Plots.plot(symm_solution_matrix', xlabel = "Scale value", ylabel = "Solution value", label = var_labs)
# look at the FOC residuals. kind of interesting patterns.
Plots.plot(symm_convergence_matrix', xlabel = "Scale value", ylabel="FOC error", legend=:none)

for j in 1:(size(symm_solution_matrix)[1])
    p = Plots.plot(symm_solution_matrix[j,:], 
            xlabel = "Scale value", 
            ylabel = "Solution value", 
            xticks = (1:12:scale_max,1:12:scale_max), 
            label = foc_labs[j])
    display(p)
    sleep(0.5)
end

for j in 1:(size(symm_solution_matrix)[1])
    p = Plots.plot(symm_convergence_matrix[j,:], 
            xlabel = "Scale value", 
            ylabel = "FOC error", 
            xticks = (1:12:scale_max,1:12:scale_max), 
            label = foc_labs[j])
    display(p)
    sleep(0.5)
end

# the following function executes an asymmetric scaling-down of L (or N_T) from L = N_T = 100.0.
# as the scale approaches ~1/10th of scale_max, solutions start becoming numerically unstable.
# this instability is worse for N_T than L, which i believe is due to the (N_T - N)^α utility term.

scaledown_choice = "L"

function asymmetric_scale( choice, scalemax )
    solution_container = []
    convergence_container = []
    did_solve = []
    last_good_guess = []

    if choice == "N_T" J = 94.0 else J = scalemax end # numerical instability overcomes N_T at >94

    for i in 1.0:J
        scale = scalemax + 1 - i
        
        # if you want to mess with any of the other "fixed" parameters, can add an "elseif..." below
        if scaledown_choice == "L"
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
                N_T = scalemax, 
                ι = 1.0 )
        else 
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
                L = scalemax,
                σ_D = 1.667,
                σ_X = 1.667,
                σ_Y = 1.667,
                N_T = scale, 
                ι = 1.0 )
        end
            if i == 1.0
                guess_vector = (
                    x = (parameter_vec.γ_N*parameter_vec.N_T + parameter_vec.γ_L*parameter_vec.L)/(parameter_vec.γ_L + parameter_vec.γ_N),
                    y = (parameter_vec.ρ_N*parameter_vec.N_T + parameter_vec.ρ_L*parameter_vec.L)/(parameter_vec.ρ_L + parameter_vec.ρ_N),
                    p = 1.0,
                    w = 1.0,
                    r = 1.0,
                    λ = 1.0,
                    L_x = 0.5*parameter_vec.L,
                    L_y = 0.5*parameter_vec.L,
                    N = (1 - parameter_vec.η_X - parameter_vec.η_Y)*parameter_vec.N_T,
                    N_x = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
                    N_y = (1 - parameter_vec.η_X - parameter_vec.η_Y)*0.5*parameter_vec.N_T,
                    K = (scale + scalemax) / 2,
                    K_x = 0.5*(scale + scalemax) / 2,
                    K_y = 0.5*(scale + scalemax) / 2 )
            else
                guess_vector = (
                    x = copy(last_good_guess[1]),
                    y = copy(last_good_guess[2]),
                    p = copy(last_good_guess[3]),
                    w = copy(last_good_guess[4]),
                    r = copy(last_good_guess[5]),
                    λ = copy(last_good_guess[6]),
                    L_x = copy(last_good_guess[7]),
                    L_y = copy(last_good_guess[8]),
                    N = copy(last_good_guess[9]),
                    N_x = copy(last_good_guess[10]),
                    N_y = copy(last_good_guess[11]),
                    K = copy(last_good_guess[12]),
                    K_x = copy(last_good_guess[13]),
                    K_y = copy(last_good_guess[14]) )
            end
        
        convergence_output, solution_output, out_success = economy_model( parameter_vec, guess_vector )
    
        push!(solution_container, solution_output)
        push!(convergence_container, convergence_output)
        push!(did_solve, out_success)
        if out_success == 1 || i == 1.0
            last_good_guess = copy(solution_output)
        end
    end
    
    which_succeeded = findall(did_solve .== 1.0)
    out_solution = reduce(hcat, solution_container)
    out_convergence = reduce(hcat, convergence_container)

    return which_succeeded, out_solution, out_convergence
end

# as previously mentioned, numerical stability breaks down eventually around scale ~= 1/10 * scale_max
# use this "well_behaved_flag" to restrict domain to the "well-behaved" region
idx_success, asym_solution_matrix, asym_convergence_matrix = asymmetric_scale(scaledown_choice, scale_max)

scale_success = (100 + 1) .- idx_success
well_behaved_flag = 1

if well_behaved_flag == 1
    if scaledown_choice == "L"
        domain = 1:91
    else
        domain = 1:87
    end
else
    domain = copy(idx_success)
end

Plots.plot(asym_solution_matrix[:,domain]', 
    xticks = (1:12:length(domain),scale_success[1:12:length(domain)]), 
    xlabel = "Scale "*scaledown_choice*" value", 
    ylabel = "Solution value", 
    label = var_labs)
# look at the FOC residuals. kind of interesting patterns.
Plots.plot(asym_convergence_matrix[:,domain]', 
     xticks = (1:12:length(domain),scale_success[1:12:length(domain)]), 
     xlabel = "Scale "*scaledown_choice*" value", 
     ylabel="FOC error", 
     label = var_labs)

# for inspecting individual solutions/FOC residuals
for j in 1:(size(asym_solution_matrix)[1])
    p = Plots.plot(asym_solution_matrix[j,domain], 
            xlabel = "Scale "*scaledown_choice*" value", 
            ylabel = "Solution value", 
            xticks = (1:12:length(domain),scale_success[1:12:length(domain)]), 
            label = var_labs[j],
            dpi = 300)
    display(p) # required to show plot inside a loop
    sleep(0.5)
end

for j in 1:(size(asym_solution_matrix)[1])
    p = Plots.plot(asym_convergence_matrix[j,domain], 
            xlabel = "Scale "*scaledown_choice*" value", 
            ylabel = "FOC error", 
            xticks = (1:12:length(domain),scale_success[1:12:length(domain)]), 
            label = foc_labs[j])
    display(p)
    sleep(0.5)
end