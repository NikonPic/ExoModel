using Plots

"""draw many points as one body"""
function draw_body!(fig, body, label; lw=3)
    x = [body_i[1] for body_i in body]
    y = [body_i[2] for body_i in body]
    plot!(fig, x, y, label=label, lw=lw)
end

"""draw many points as one body"""
function draw_body_force!(fig, body, label; lw=3)
    x = [body_i[1] for body_i in body]
    y = [body_i[2] for body_i in body]
    plot!(fig, x, y, label=label, lw=lw, arrow=0.4, color="red")
end

"""plot a single line between a and v"""
function plot_line!(fig, A, B, label; lw=3)
    length = round(norm(A, B), digits=2)
    plot!(fig, [A[1], B[1]], [A[2], B[2]], label=string(label, ":", length), lw=lw)
end

"""
draw the complete anmiation
"""
function draw_exo_traj(sol, p, u_des)
    u, t = sol.u, sol.t

    tend = t[end]

    # go over all timesteps
    for (i, ui) in enumerate(u)
        alpha, beta, gamma = ui[1], ui[2], ui[3]
        draw_current_image(u, t, i, alpha, beta, gamma, tend, p, u_des)
    end
end

"""
draw the ull image for one timestep
"""
function draw_current_image(u, t, i, alpha, beta, gamma, t_end, p, u_des)
    Tₒ₁, Tₒ₂, Tₒ₃, _, _ = get_transformations(alpha, beta, gamma)
    all_points = forward_circuits(alpha, beta, gamma, p)
    body_0, body_1, body_2, body_3, body_stab = get_draw_bodies(all_points, Tₒ₁, Tₒ₂, Tₒ₃, p)

    f_attach = (all_points[5][4], all_points[5][4] + [5*sin(2*t[i]), 0])

    plotsize = (800, 960)
    fig = plot(;size = plotsize, layout=grid(3, 1, heights=[0.8, 0.1, 0.1]))

    # draw the body
    draw_body!(fig[1], body_0, "body_0")
    draw_body!(fig[1], body_1, "body_1")
    draw_body!(fig[1], body_2, "body_2")
    draw_body!(fig[1], body_3, "body_3")
    draw_body!(fig[1], body_stab, "exo")
    draw_body_force!(fig[1], f_attach, "external force")
    xlims!(-17.5, 22.5)
    ylims!(-20, 20)

    # draw angles over t
    _alpha = [_u[1] for _u in u][1:i]
    _beta = [_u[2] for _u in u][1:i]
    _gamma = [_u[3] for _u in u][1:i]
    d_alpha = [_u[4] for _u in u][1:i]
    d_beta = [_u[5] for _u in u][1:i]
    d_gamma = [_u[6] for _u in u][1:i]
    t_i = t[1:i]
    u_des_i = u_des[1:i]

    plot!(fig[2], t_i, _alpha, label="alpha")
    plot!(fig[2], t_i, _beta, label="beta")
    plot!(fig[2], t_i, _gamma, label="gamma")
    plot!(fig[2], t_i, u_des_i, label="u_des", line=:dot, color=:green)
    xlims!(fig[2], 0, t_end)
    ylims!(fig[2], -2, 2)
    ylabel!(fig[2], "rad")
    xlabel!(fig[2], "t")

    plot!(fig[3], t_i, d_alpha, label="d_alpha")
    plot!(fig[3], t_i, d_beta, label="d_beta")
    plot!(fig[3], t_i, d_gamma, label="d_gamma")
    xlims!(fig[3], 0, t_end)
    ylims!(fig[3], -2, 2)
    ylabel!(fig[3], "rad/s")
    xlabel!(fig[3], "t")

    frame(anim)
end

function get_draw_bodies(all_points, Tₒ₁, Tₒ₂, Tₒ₃, p)
    l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = safety_p(p)

    # extract points
    A_all, B_all, C_all, D_all, E_all = all_points

    # define body_0
    body_0 = (E_all[5], [E_all[5][1], 0], [E_all[1][1], 0], E_all[1], [E_all[1][1], 0], C_all[1])

    # body 1
    D5ₒ_bot =  Tₒ₁ * [a11, 0]
    D1ₒ_bot =  Tₒ₁ * [a12, 0]

    body_1 = (C_all[1], D5ₒ_bot, D_all[5], D5ₒ_bot, D1ₒ_bot, D_all[1], D1ₒ_bot, B_all[1])

    # body 2
    B2ₒ_bot = Tₒ₂ * [a2 / 2, 0] + B_all[1]
    body_2 = (B_all[1], B2ₒ_bot, B_all[2], B2ₒ_bot, A_all[1])

    # body 3
    A2ₒ_bot = Tₒ₃ * [a3 / 2, 0] + A_all[1]
    endpoint = Tₒ₃ * [a3, 0] + A_all[1]
    body_3 = (A_all[1], A2ₒ_bot, A_all[2], A2ₒ_bot, endpoint)

    # fachwerk
    body_stab = (E_all[5], E_all[4], E_all[3], E_all[2], E_all[1], E_all[2], C_all[2], E_all[2], E_all[3], D_all[2], D_all[1], D_all[2], A_all[3], A_all[4], A_all[3], A_all[2])

    return body_0, body_1, body_2, body_3, body_stab
end



function cb_hand(p, l, pred)
    alpha, beta, gamma = 0, 0, 0
    Tₒ₁, Tₒ₂, Tₒ₃, _, _ = get_transformations(alpha, beta, gamma)
    all_points = forward_circuits(alpha, beta, gamma, p)
    body_0, body_1, body_2, body_3, body_stab = get_draw_bodies(all_points, Tₒ₁, Tₒ₂, Tₒ₃, p)

    plotsize = (800, 960)
    fig = plot(;size = plotsize, layout=grid(3, 1, heights=[0.8, 0.1, 0.1]))

    # draw the body
    draw_body!(fig[1], body_0, "body_0")
    draw_body!(fig[1], body_1, "body_1")
    draw_body!(fig[1], body_2, "body_2")
    draw_body!(fig[1], body_3, "body_3")
    draw_body!(fig[1], body_stab, "exo")
    title!(fig[1], string(round(l, digits=2)))

    xlims!(-17.5, 22.5)
    ylims!(-20, 20)

    plot!(fig[2], tsteps, pred[1,:], label="alpha")
    plot!(fig[2], tsteps, pred[2,:], label="beta")
    plot!(fig[2], tsteps, pred[3,:], label="gamma")
    plot!(fig[2], tsteps, u_des, label="u_des", line=:dot, color=:blue)

    xlims!(fig[2], 0, t_end)
    ylims!(fig[2], -2, 2)
    ylabel!(fig[2], "rad")
    xlabel!(fig[2], "t")

    append!(l_all, l)
    append!(i_all, global iter += 1)

    plot!(fig[3], i_all, l_all, ylim=[50, l_max], xlim=[1, max_iter], lab="loss")
    xlabel!(fig[3], "iterations")
    frame(anim)
    gui()
end
