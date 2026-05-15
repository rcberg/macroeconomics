using JuMP, Ipopt
using Symbolics
using Plots

# load up our functions
# note: they're separated so you can tweak the "_func_compiler.jl" with new functions, if desired
include("symbolic_function_city_simulator_func_compiler.jl")
include("symbolic_function_city_simulator_builder.jl")
#include("stuart_initial_val_fun.jl")
# loads the compiled cost/expenditure functions; outputting `fns`
fns = compile_step()

# takes `fns` from `compile_step()` and builds the JuMP model
model_lg = build_model(fns;
    Q = 1.0, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
    τ = 0.0, ι_bar = 1.0, 
    N_fixed = 1000.0
)

# solve the model output and acquire the solutions
optimize!(model_lg)
solution_lg = NamedTuple{Tuple(Symbol.(name.(all_variables(model_lg))))}(value.(all_variables(model_lg)))

support = 0.81:0.01:1.2

# using "warm start"
### Q ###
results_q = []

model_sm = build_model(fns;
        Q = 0.8, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_q, (Q = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support

    global model_sm, prev_model
    
    model_sm = build_model(fns;
        Q = i, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)))
    optimize!(model_sm)
    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_q, (Q = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
end

plot_x = [row.Q for row in results_q] .- 1.0
plot_y = [row.N for row in results_q] ./ [row.N for row in results_q][plot_x .== 0.0] .- 1.0
Plots.plot(plot_x, plot_y, label = "Q",
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-1.5, 2), yticks = -1.5:0.5:2.0)

### AX ###
results_ax = []

model_sm = build_model(fns;
        Q = 1.0, A_X = 0.8, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_ax, (AX = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support
    global model_sm, prev_model

    model_sm = build_model(fns;
        Q = 1.0, A_X = i, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
    set_optimizer_attribute(model_sm, "max_iter", 5000)
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)))
    optimize!(model_sm)
    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_ax, (AX = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
end

plot_x = [row.AX for row in results_ax] .- 1.0
plot_y = ([row.N for row in results_ax] ./ [row.N for row in results_ax][plot_x .== 0.0] .- 1.0) ./ (1-0.36)
Plots.plot(plot_x, plot_y, label = "AX",
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-1.5, 1.5), yticks = -1.5:0.5:1.5
            )

### AY ###
results_ay = []

model_sm = build_model(fns;
        Q = 1.0, A_X = 1.0, A_Y = 0.8, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
#set_optimizer_attribute(model_sm, "max_iter", 5000)
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_ay, (AY = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support
    global model_sm, prev_model

    model_sm = build_model(fns;
        Q = 1.0, A_X = 1.0, A_Y = i, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
        u_bar = solution_lg.u_bar_v
        )
    set_optimizer_attribute(model_sm, "max_iter", 5000)
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)))
    optimize!(model_sm)
    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_ay, (AY = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
end

plot_x = [row.AY for row in results_ay] .- 1.0
plot_y = ( [row.N for row in results_ay] ./ [row.N for row in results_ay][plot_x .== 0.0] .- 1.0 ) ./0.36
Plots.plot(plot_x, plot_y, label = "AY", 
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-2.0, 2.0), yticks = -2.0:0.5:2.0 )

            # using stuart supplied initial values
#syms = [:N_total, :w, :p, :r, :x, :y, :X, :Y, :NX, :NY, :LX, :LY, :KX, :KY, :K]
#    
#### Q ###
#results_q = []
#
#model_sm = build_model(fns;
#        Q = 0.8, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
#        u_bar = solution_lg.u_bar_v
#        )
#set_optimizer_attribute(model_sm, "print_level", 0)
#optimize!(model_sm)
#push!(results_q, (Q = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#for i in support
#
#    global model_sm
#
#    ivals = stuart_init_val_fn(syms; variable = "Q", val = i)
#    
#    model_sm = build_model(fns;
#        Q = i, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
#        u_bar = solution_lg.u_bar_v
#        )
#    set_optimizer_attribute(model_sm, "print_level", 0)
#    set_start_value.([model_sm[i] for i in syms], ivals)
#    optimize!(model_sm)
#    
#    push!(results_q, (Q = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#end
#