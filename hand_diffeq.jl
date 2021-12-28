using DifferentialEquations, Flux, DiffEqFlux, Plots, JLD2
include("./src/5_circuits.jl")
include("./src/params.jl")
include("./src/visualization.jl")
pyplot()

u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Simulation interval and intermediary points
t_end = 5
dt = 0.05
tspan = (0.0, t_end)
tsteps = 0.0:dt:t_end

# define the external Force
F_ext = [1, 0]
u_des = 1 .+ 1 *sin.(2*tsteps)

# LV equation parameters
p = zeros(39)
prob = ODEProblem(diffeq_hand_exo!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=dt)
plot(sol)

anim = Plots.Animation()
draw_exo_traj(sol, p, u_des)
gif(anim, "./res/hand_new.gif", fps=20)


"""all losses in one func"""
function loss(p)
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    u0 += 0.1 .* randn(6)
    prob = ODEProblem(diffeq_hand_exo!, u0, tspan, p)
    sol = solve(prob, Tsit5(), p=p, saveat=tsteps)
    loss  = sum(abs2, sol[1, :])# .- u_des)
    loss += sum(abs2, sol[2, :] .- u_des)
    loss += sum(abs2, sol[3, :] .- u_des)
    loss += apply_length_safety(sol, p)

    return loss, sol
end

p_opt = copy(p)

callback = function (p_run, l, pred)
    @save "test1.jld2" p_run
    # Tell sciml_train to not halt the optimization. If return true, then
    # optimization stops.
    cb_hand(p_run, l, pred)
    return false
end


anim = Plots.Animation()
l_all = []
l_max = loss(p)[1]
i_all = []

iter = 0
max_iter = 100
result_ode = DiffEqFlux.sciml_train(loss, p, ADAM(0.005), cb=callback, maxiters=max_iter, save_best=true)
gif(anim, "./data/results/hand_training.gif", fps=20)

p_opt = result_ode.minimizer

prob = ODEProblem(diffeq_hand_exo!, u0, tspan, p_opt)
sol = solve(prob, Tsit5(), saveat=dt)
plot(sol)

anim = Plots.Animation()

@load "test1.jld2"
prob = ODEProblem(diffeq_hand_exo!, u0, tspan, p_run)
sol = solve(prob, Tsit5(), saveat=dt)
draw_exo_traj(sol, p_run, u_des)
gif(anim, "./data/results/hand_opt.gif", fps=20)
