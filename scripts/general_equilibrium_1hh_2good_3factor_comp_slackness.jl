using Symbolics
using NonlinearSolve, LinearAlgebra
using Plots

# these are going to do heavy lifting making this problem converge
# we need to enforce constraints via transformation of variables
#
# sigmoid will constrain the solution to be nonnegative
# logit will transform the solution back to the correct form
function logit(val, limit)
    if val isa Number
        ϵ = 1e-12
        x_safe = clamp(val, ϵ, limit - ϵ)
        return log(x_safe / (limit - x_safe) )
    else
        return log(val / (limit - val))
    end
end

sigmoid(val, limit) = limit / (1.0 + exp(-val))

@variables v[1:15]
@variables q a_x a_y L σ_D σ_X σ_Y N_T
@variables γ_L γ_N ρ_L ρ_N η_X η_Y η_N

idx_u = (
    x = 1,
    y = 2,
    p = 3,
    w = 4,
    r = 5,
    ι = 6,
    λ = 7,
    L_x = 8,
    L_y = 9,
    N = 10,
    N_x = 11,
    N_y = 12,
    K = 13,
    K_x = 14,
    K_y = 15
)

idx_p = (
    γ_L = 1, 
    γ_N = 2, 
    ρ_L = 3, 
    ρ_N = 4, 
    η_X = 5, 
    η_Y = 6, 
    q = 7, 
    a_x = 8, 
    a_y = 9, 
    L = 10, 
    σ_D = 11, 
    σ_X = 12, 
    σ_Y = 13, 
    N_T = 14
)

d = (
    x = exp(v[idx_u.x]),
    y = exp(v[idx_u.y]),
    p = exp(v[idx_u.p]),
    w  = exp(v[idx_u.w]),
    r  = exp(v[idx_u.r]),
    ι  = exp(v[6]),
    λ = exp(v[idx_u.λ]),
    L_x = sigmoid(v[idx_u.L_x], L),
    L_y = sigmoid(v[idx_u.L_y], L),
    N = exp(v[idx_u.N]),
    N_x = exp(v[idx_u.N_x]),
    N_y = exp(v[idx_u.N_y]),
    K = exp(v[idx_u.K]),
    K_x = exp(v[idx_u.K_x]),
    K_y = exp(v[idx_u.K_y])
    )

α = (σ_D - 1)/σ_D
β = (σ_X - 1)/σ_X
χ = (σ_Y - 1)/σ_Y

u_xy = q*( η_X*(d.x)^α + η_Y*(d.y)^α + (1 - η_X - η_Y)*(N_T - d.N) )^(1/α)
pf_x = a_x*( γ_L*(d.L_x)^β + γ_N*(d.N_x)^β + (1 - γ_L - γ_N)*(d.K_x)^β )^(1/β)
pf_y = a_y*( ρ_L*(d.L_y)^χ + ρ_N*(d.N_y)^χ + (1 - ρ_L - ρ_N)*(d.K_y)^χ )^(1/χ)

Λ_u = u_xy + d.λ * ( d.w*d.N + d.r*L + d.ι*d.K - d.x - d.p*d.y )
Π_x = pf_x - ( d.w*d.N_x + d.r*d.L_x + d.ι*d.K_x )
Π_y = (d.p)*pf_y - ( d.w*d.N_y + d.r*d.L_y + d.ι*d.K_y )

#h_vars = [ d.x, d.y, d.N, d.λ ]
#x_vars = [ d.L_x, d.N_x, d.K_x ]
#y_vars = [ d.L_y, d.N_y, d.K_y ]

h_focs = Symbolics.gradient(Λ_u, [v[idx_u.x], v[idx_u.y], v[idx_u.N], v[idx_u.λ]] )
firm_x_focs = Symbolics.gradient(Π_x, [v[idx_u.L_x], v[idx_u.N_x], v[idx_u.K_x]] )
firm_y_focs = Symbolics.gradient(Π_y, [v[idx_u.L_y], v[idx_u.N_y], v[idx_u.K_y]] )

k_clearing = sqrt( d.ι^2.0 + (d.K - d.K_x - d.K_y)^2.0 ) - (d.ι + d.K - d.K_x - d.K_y)
l_clearing = sqrt( d.r^2.0 + (L - d.L_x - d.L_y)^2.0 ) - (d.r + L - d.L_x - d.L_y)
n_clearing = sqrt( d.w^2.0 + (d.N - d.N_x - d.N_y)^2.0 ) - (d.w + d.N - d.N_x - d.N_y)

x_clearing = d.x - pf_x
y_clearing = sqrt(d.p^2.0 + (d.y - pf_y)^2.0) - (d.p + d.y - pf_y)

system_eqs = [firm_x_focs; firm_y_focs; h_focs; x_clearing; y_clearing; k_clearing; l_clearing; n_clearing]

params = [γ_L, γ_N, ρ_L, ρ_N, η_X, η_Y, q, a_x, a_y, L, σ_D, σ_X, σ_Y, N_T]
full_sys = Symbolics.expand(system_eqs)

j_matrix = Symbolics.jacobian(full_sys, v)
#det_j = Symbolics.det(j_matrix) #apparently this takes >1.3 trillion terms to compute and sum. oops.7
# reminder
# parameters: [γ_L, γ_N, ρ_L, ρ_N, η_X, η_Y, q, a_x, a_y, L, σ_D, σ_X, σ_Y, N_T]
# variables:  [x,   y,   p,   w,   r,   ι,   λ, L_x, L_y, N, N_x, N_y, K,   K_x, K_y]
p_vals = [ 0.0134, 0.9284, 0.3668, 0.5773, 0.6949, 0.1, 1.0, 1.0, 1.0, 1.0, 1.667, 1.667, 1.667, 1.0]
guessvec = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5*p_vals[idx_p.L], 0.5*p_vals[idx_p.L], p_vals[idx_p.N_T], 0.5*p_vals[idx_p.N_T], 0.5*p_vals[idx_p.N_T], 1.0, 0.5, 0.5 ]

u0 = [ log.(guessvec[1:7]); logit.(guessvec[8:9], p_vals[idx_p.L]); log.(guessvec[10:15]) ]
f_ge! = build_function(full_sys, v, params, expression=false)[2]
j_ip! = build_function(j_matrix, v, params, expression=false)[2]

#LBs = ones(length(vars)).*1e-6
#UBs = ones(length(vars)).*Inf

ge_system = NonlinearFunction(f_ge!; jac = j_ip!)
problem = NonlinearProblem( ge_system, u0, p_vals )

sol = solve(problem, LevenbergMarquardt())

raw_solution = sol.u
if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.StalledSuccess
    println("**Solved successfully**")
    println("")
else
    println("Solver failed with code: $(sol.retcode).")
    println("")
end
solution = [ exp.(raw_solution[1:7]); sigmoid.(raw_solution[8:9], p_vals[idx_p.L]); exp.(raw_solution[10:15]) ]

alfa = (-1 + p_vals[idx_p.σ_D])/p_vals[idx_p.σ_D]
beta = (-1 + p_vals[idx_p.σ_X])/p_vals[idx_p.σ_X]
chi =  (-1 + p_vals[idx_p.σ_Y])/p_vals[idx_p.σ_Y]

profit_x = solution[idx_u.x] -             solution[idx_u.w]*solution[idx_u.N_x] - solution[idx_u.r]*solution[idx_u.L_x] -  solution[idx_u.ι]*solution[idx_u.K_x]
profit_y = solution[idx_u.p]*solution[idx_u.y] - solution[idx_u.w]*solution[idx_u.N_y] - solution[idx_u.r]*solution[idx_u.L_y] -  solution[idx_u.ι]*solution[idx_u.K_y]
profit_total = profit_x + profit_y

# gradient check
# these should all be tiny on the order of << 1e-10
println("(HH:) income - expenses = $(solution[idx_u.w]*solution[idx_u.N] + solution[idx_u.r]*p_vals[idx_p.L] + solution[idx_u.ι]*solution[idx_u.K] - solution[idx_u.x] - solution[idx_u.p]*solution[idx_u.y])")
println("(L:)  demand - supply = $(p_vals[idx_p.L] - solution[idx_u.L_x] - solution[idx_u.L_y])") # land market clearing
println("(N:)  demand - supply = $(solution[idx_u.N] - solution[idx_u.N_x] - solution[idx_u.N_y])") # labor market clearing 
println("(K:)  demand - supply = $(solution[idx_u.K] - solution[idx_u.K_x] - solution[idx_u.K_y])") # capital market clearing
println("good X profit = $profit_x") # firm X zero-profit condition
#println("good X revenue = $(solution[idx_u.x])") # these were needed when i didn't have the zpc nailed-down
println("good Y profit = $profit_y") # firm Y zero-profit condition
#println("good Y revenue = $( solution[idx_u.p]*solution[idx_u.y] )") # these were needed when i didn't have the zpc nailed-down
#println("total industry revenue = $( solution[idx_u.x] + solution[idx_u.p]*solution[idx_u.y])")
println("total industry profit = $profit_total") # total zero profit condition
println("λ - ∂U/∂x = $(    solution[idx_u.λ] -                                              p_vals[idx_p.η_X]  * (solution[idx_u.x]^(  alfa - 1)) * p_vals[idx_p.q]   * ( p_vals[idx_p.η_X]*(solution[idx_u.x])^alfa   + p_vals[idx_p.η_Y]*(solution[idx_u.y])^alfa   + (1 - p_vals[idx_p.η_X] - p_vals[idx_p.η_Y])*(p_vals[idx_p.N_T] - solution[idx_u.N]) )^((1/alfa) - 1) ) ") # consumption foc
println("w - ∂X/∂N_x = $(  solution[idx_u.w] -                                              p_vals[idx_p.γ_N]  * (solution[idx_u.N_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) )") # good x wage foc
println("w - p*∂Y/∂N_y = $(solution[idx_u.w] - solution[idx_u.p] *                          p_vals[idx_p.ρ_N]  * (solution[idx_u.N_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) )") # good y wage foc
println("r - ∂X/∂L_x = $(  solution[idx_u.r] -                                              p_vals[idx_p.γ_L]  * (solution[idx_u.L_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) )") # good x land foc
println("r - p*∂Y/∂L_y = $(solution[idx_u.r] - solution[idx_u.p] *                          p_vals[idx_p.ρ_L]  * (solution[idx_u.L_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) )") # good y land foc
println("ι - ∂X/∂K_x = $(  solution[idx_u.ι] -                     (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N]) * (solution[idx_u.K_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) )") # good x capt foc
println("ι - p*∂Y/∂K_y = $(solution[idx_u.ι] - solution[idx_u.p] * (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N]) * (solution[idx_u.K_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) )") # good y capt foc
println("")

# show the solution to the system
println(join(string.(keys(idx_u)) .* " = " .* string.(solution),"\n"))
Plots.bar([i for i in String.(keys(d))], solution, legend=:none)