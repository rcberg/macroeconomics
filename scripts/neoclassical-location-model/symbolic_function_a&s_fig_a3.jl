using JuMP, Ipopt
using Symbolics
using Plots

# load up our functions
# note: they're separated so you can tweak the "_func_compiler.jl" with new functions, if desired
include("symbolic_function_city_simulator_func_compiler.jl")
include("symbolic_function_city_simulator_builder.jl")

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

# Q has a couple versions as i tried matching it to a&s better.
# supplied values work well for Q

### Q ###

# using default initial values
#support = 0.81:0.01:1.2
#results_q = []
#
#model_sm = build_model(fns;
#        Q = 0.8, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#set_optimizer_attribute(model_sm, "max_iter", 500)
#set_optimizer_attribute(model_sm, "print_level", 0)
#optimize!(model_sm)
#push!(results_q, (Q = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#for i in support
#    global model_sm
#    model_sm = build_model(fns;
#        Q = i, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#    set_optimizer_attribute(model_sm, "max_iter", 500)
#    set_optimizer_attribute(model_sm, "print_level", 0)
#    optimize!(model_sm)
#    push!(results_q, (Q = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#end

# using supplied initial values from stuart
init_syms = [:N_total, :w, :p, :r, :x, :y, :X, :Y, :NX, :NY, :LX, :LY, :KX, :KY, :K]
v1 = [0.0001, 0.65, 0.65, 0.005, 0.55,  0.55, 0.00001,  0.00001, 0.00001,  0.00001, 0.0004,  0.0004, 0.000001,  0.000001, 0.000001]
v2 = [0.0001, 0.7, 0.7, 0.05, 0.4,  0.4, 0.0003,  0.0003, 0.0003,  0.0003, 0.0004,  0.0004, 0.00001,  0.00001, 0.00001  ]
v3 = [0.001, 0.8, 0.8, 0.07, 0.4,  0.4, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001  ]

support = 0.81:0.01:1.2
results_q = []

model_sm = build_model(fns;
        Q = 0.8, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        u_bar = solution_lg.u_bar_v
        )
set_optimizer_attribute(model_sm, "max_iter", 5000)
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_q, (Q = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
for i in support

    global model_sm
    
    if 0.8 <= i <= 0.84
        iv = v1
    elseif 0.841 <= i <= 1.0
        iv = v2
    else
        iv = v3
    end
    
    model_sm = build_model(fns;
        Q = i, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        u_bar = solution_lg.u_bar_v
        )
    set_optimizer_attribute(model_sm, "max_iter", 5000)
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.([model_sm[j] for j in init_syms], iv)
    optimize!(model_sm)
    push!(results_q, (Q = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
end

### using "warm start"
#support = 0.81:0.01:1.2
#results_q = []
#
#model_sm = build_model(fns;
#        Q = 0.8, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#set_optimizer_attribute(model_sm, "max_iter", 5000)
#set_optimizer_attribute(model_sm, "print_level", 0)
#optimize!(model_sm)
#push!(results_q, (Q = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#prev_model = model_sm
#for i in support
#    
#    model_sm = build_model(fns;
#        Q = i, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#    set_optimizer_attribute(model_sm, "max_iter", 5000)
#    set_optimizer_attribute(model_sm, "print_level", 0)
#    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)))
#    optimize!(model_sm)
#    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
#        prev_model = model_sm
#    end
#    push!(results_q, (Q = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#end

plot_x = [row.Q for row in results_q] .- 1.0
plot_y = [row.N for row in results_q] ./ [row.N for row in results_q][plot_x .== 0.0] .- 1.0
Plots.plot(plot_x, plot_y, label = "Q")

### AX ###

# using supplied initial values from stuart
#init_syms = [:N_total, :w, :p, :r, :x, :y, :X, :Y, :NX, :NY, :LX, :LY, :KX, :KY, :K]
#v1 = [0.0001, 0.5, 0.5, 0.01, 0.5,  0.5, 0.0001,  0.0001, 0.0001,  0.0001, 0.0004,  0.0004, 0.00001,  0.00001, 0.00001  ]
#v2 = [0.0001, 0.7, 0.7, 0.05, 0.5,  0.5, 0.0001,  0.0001, 0.0001,  0.0001, 0.0004,  0.0004, 0.0001,  0.0001, 0.0001]
#v3 = [0.001, 0.8, 1.2, 0.1, 0.45,  0.45, 0.0008,  0.0008, 0.0008,  0.0008, 0.0005,  0.0005, 0.0001,  0.0001, 0.0001 ] 
#v4 = [0.001, 0.8, 1.3, 0.2, 0.6,  0.3, 0.001,  0.001, 0.0008,  0.0008, 0.0005,  0.0005, 0.0001,  0.0001, 0.0001  ]
#
#support = 0.81:0.01:1.2
#results_ax = []
#
#model_sm = build_model(fns;
#        Q = 1.0, A_X = 0.8, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#set_optimizer_attribute(model_sm, "max_iter", 5000)
#set_optimizer_attribute(model_sm, "print_level", 0)
#optimize!(model_sm)
#push!(results_ax, (AX = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#for i in support
#    
#    if 0.8 <= i <= 0.95
#        iv = v1
#    elseif 0.951 <= i <= 1.15
#        iv = v2
#    elseif 1.151 <= i <= 1.19
#        iv = v3
#    else
#        iv = v4
#    end
#    
#    model_sm = build_model(fns;
#        Q = 1.0, A_X = i, A_Y = 1.0, L_0 = 1000.0, 
#        τ = 0.0, ι_bar = 1.0, 
#        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
#        u_bar = solution_lg.u_bar_v
#        )
#    set_optimizer_attribute(model_sm, "max_iter", 5000)
#    set_optimizer_attribute(model_sm, "print_level", 0)
#    set_start_value.([model_sm[j] for j in init_syms], iv)
#    optimize!(model_sm)
#    push!(results_ax, (AX = i, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
#end

# using "warm start"
support = 0.81:0.01:1.2
results_ax = []

model_sm = build_model(fns;
        Q = 1.0, A_X = 0.8, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        u_bar = solution_lg.u_bar_v
        )
set_optimizer_attribute(model_sm, "max_iter", 5000)
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_ax, (AX = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support
    global model_sm

    model_sm = build_model(fns;
        Q = 1.0, A_X = i, A_Y = 1.0, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
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
Plots.plot(plot_x, plot_y, label = "AX")

### AY ###

# using "warm start"
support = 0.81:0.01:1.2
results_ay = []

model_sm = build_model(fns;
        Q = 1.0, A_X = 1.0, A_Y = 0.8, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        u_bar = solution_lg.u_bar_v
        )
set_optimizer_attribute(model_sm, "max_iter", 5000)
set_optimizer_attribute(model_sm, "print_level", 0)
optimize!(model_sm)
push!(results_ay, (AY = 0.8, N = value(model_sm[:N_total]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support
    global model_sm

    model_sm = build_model(fns;
        Q = 1.0, A_X = 1.0, A_Y = i, L_0 = 1000.0, 
        τ = 0.0, ι_bar = 1.0, 
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
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
plot_y = ( [row.N for row in results_ay] ./ [row.N for row in results_ay][plot_x .== 0.0] .- 1.0 )./0.36
Plots.plot(plot_x, plot_y, label = "AY")

###stuart initial value bank
### #Q: 0.8 <= val <= 0.84
### #[0.0001, 0.65, 0.65, 0.005, 0.55,  0.55, 0.00001,  0.00001, 0.00001,  0.00001, 0.0004,  0.0004, 0.000001,  0.000001, 0.000001]
### #Q: 0.841 <= val <= 1.0
### #[0.0001, 0.7, 0.7, 0.05, 0.4,  0.4, 0.0003,  0.0003, 0.0003,  0.0003, 0.0004,  0.0004, 0.00001,  0.00001, 0.00001  ]
### #Q: 1.001 <= val <= 1.2
### [0.001, 0.8, 0.8, 0.07, 0.4,  0.4, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001,  0.0001, 0.0001  ]
### #AX: 0.8 <= val <= 0.95
### [0.0001, 0.5, 0.5, 0.01, 0.5,  0.5, 0.0001,  0.0001, 0.0001,  0.0001, 0.0004,  0.0004, 0.00001,  0.00001, 0.00001  ]
### #AX: 0.951 <= val <= 1.15
### [0.0001, 0.7, 0.7, 0.05, 0.5,  0.5, 0.0001,  0.0001, 0.0001,  0.0001, 0.0004,  0.0004, 0.0001,  0.0001, 0.0001]
### #AX: 1.151 <= val <= 1.19
### [0.001, 0.8, 1.2, 0.1, 0.45,  0.45, 0.0008,  0.0008, 0.0008,  0.0008, 0.0005,  0.0005, 0.0001,  0.0001, 0.0001 ] 
### #AX: val == 1.2
### [0.001, 0.8, 1.3, 0.2, 0.6,  0.3, 0.001,  0.001, 0.0008,  0.0008, 0.0005,  0.0005, 0.0001,  0.0001, 0.0001  ]
### #AY: 0.8 <= val <= 0.92
### [0.0001, 0.7, 0.7, 0.005, 0.55,  0.55, 0.0001,  0.0001, 0.0001,  0.0001, 0.0004,  0.0004, 0.00001,  0.00001, 0.00001]
### #AY: 0.921 <= val <= 1.03
### [0.0006, 0.6, 0.8, 0.05, 0.45,  0.45, 0.0003,  0.0003, 0.0003,  0.0003, 0.0004,  0.0004, 0.0001,  0.00001, 0.0001]
### #AY: 1.04 <= val
### [0.0001, 0.6, 0.8, 0.05, 0.45,  0.45, 0.0003,  0.0003, 0.0003,  0.0003, 0.0004,  0.0004, 0.0001,  0.0001, 0.0001  ]
