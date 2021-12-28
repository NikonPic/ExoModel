using LinearAlgebra
include("./mass_matrix.jl")


"""
special characters to select
ₒ₁₂₃₄

      B     D    E
   A           C
---0----0-----0-|
  PIP  DIP   MCP
"""

function get_transformations(alpha, beta, gamma)
    # manage the separate rotations
    T₁ₒ = rotx(alpha)
    T₂₁ = rotx(beta)
    T₃₂ = rotx(gamma)
    T₂ₒ = T₂₁ * T₁ₒ
    T₃ₒ = T₃₂ * T₂ₒ

    # for clean notation save transpose
    Tₒ₁, T₁₂, T₂₃, Tₒ₂, Tₒ₃ = T₁ₒ', T₂₁', T₃₂', T₂ₒ', T₃ₒ'

    return Tₒ₁, Tₒ₂, Tₒ₃, T₂₃, T₂₁
end

"""manage all points of the forward circuits"""
function forward_circuits(alpha, beta, gamma, p)
    l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = safety_p(p)

    # get the rotation matrices
    Tₒ₁, Tₒ₂, Tₒ₃, T₂₃, T₂₁ = get_transformations(alpha, beta, gamma)

    # define the vector between coo sys 2 and 0 (from 0 to 2 in 0)
    diff_20₁ = [a1, 0]
    diff_20ₒ = Tₒ₁ * diff_20₁

    # define all relevant coordinates of the links / Forward kinematics
    A1₂, A2₂, A3₂, A4₂      = get_4link_A(l_s9, l_s10, a2, a3, b2, b3, T₂₃)
    A1ₒ, A2ₒ, A3ₒ, A4ₒ      = f2t0(A1₂, Tₒ₂, diff_20ₒ), f2t0(A2₂, Tₒ₂, diff_20ₒ), f2t0(A3₂, Tₒ₂, diff_20ₒ), f2t0(A4₂, Tₒ₂, diff_20ₒ)
    B1₂, B2₂, B3₂, B4₂, B5₂ = get_5link_B(l_s7, l_s8, a1, a12, b12, T₂₁, A3₂, A4₂)
    B1ₒ, B2ₒ, B3ₒ, B4ₒ, B5ₒ = f2t0(B1₂, Tₒ₂, diff_20ₒ), f2t0(B2₂, Tₒ₂, diff_20ₒ), f2t0(B3₂, Tₒ₂, diff_20ₒ), f2t0(B4₂, Tₒ₂, diff_20ₒ), f2t0(B5₂, Tₒ₂, diff_20ₒ)
    C1ₒ, C2ₒ, C3ₒ, C4ₒ      = get_4link_C(l_s4, l_s5, a02, a11, b02, b11, Tₒ₁)
    D1ₒ, D2ₒ, D3ₒ, D4ₒ, D5ₒ = get_5link_D(l_s3, l_s6, B4ₒ, B5ₒ, C2ₒ, C3ₒ)
    E1ₒ, E2ₒ, E3ₒ, E4ₒ, E5ₒ = get_5link_E(l_s1, l_s2, a01, b01, C3ₒ, C4ₒ, D3ₒ)

    # concatenate
    A_all = (A1ₒ, A2ₒ, A3ₒ, A4ₒ)
    B_all = (B1ₒ, B2ₒ, B3ₒ, B4ₒ, B5ₒ)
    C_all = (C1ₒ, C2ₒ, C3ₒ, C4ₒ)
    D_all = (D1ₒ, D2ₒ, D3ₒ, D4ₒ, D5ₒ)
    E_all = (E1ₒ, E2ₒ, E3ₒ, E4ₒ, E5ₒ)

    all_points = (A_all, B_all, C_all, D_all, E_all)

    return all_points
end


"""calculate the stab directions"""
function stab_dircetions(all_points)
    A_all, B_all, C_all, D_all, E_all = all_points

    # get all stab directions
    S1_dir  = get_dir(E_all[5], E_all[4])
    S2_dir  = get_dir(E_all[4], E_all[3])
    S3_dir  = get_dir(E_all[3], E_all[2])
    S4_dir  = get_dir(E_all[2], C_all[4])
    S5_dir  = get_dir(C_all[3], C_all[2])
    S6_dir  = get_dir(D_all[3], D_all[2])
    S7_dir  = get_dir(D_all[2], D_all[1])
    S8_dir  = get_dir(B_all[4], B_all[3])
    S9_dir  = get_dir(B_all[3], B_all[2])
    S10_dir = get_dir(A_all[3], A_all[2])

    # accumulate to one field
    S_dirs = (S1_dir, S2_dir, S3_dir, S4_dir, S5_dir, S6_dir, S7_dir, S8_dir, S9_dir, S10_dir)

    return S_dirs
end

"""calculate all relevant positions of the hand using forward kinematics"""
function diffeq_hand_exo!(du, u, p, t)
    # get the state
    alpha, beta, gamma, d_alpha, d_beta, d_gamma = u

    # apply forward kinematics
    all_points = forward_circuits(alpha, beta, gamma, p)
    S_dirs = stab_dircetions(all_points)
    
    # external Force
    F_ext = [80 * sin(2*t), 0]

    # derive force distribution
    F_stab = extract_symbolic_stab_forces(S_dirs, F_ext)

    # S_dirs with force
    S_dirs_f = F_stab .* S_dirs

    # calculate the Mass matrix, h and Q_ext
    M_all, h_all, Q_all = get_M_h_Q_single(alpha, beta, gamma, d_alpha, d_beta, d_gamma, p, S_dirs_f)
    ddq = get_ddq(M_all, h_all, Q_all)

    # define derivatives
    du[1], du[2], du[3] = u[4], u[5], u[6]
    du[4], du[5], du[6] = ddq[1], ddq[2], ddq[3]
end


"""calculate all relevant positions of the hand using forward kinematics"""
function diffeq_hand!(du, u, p, t)
    # get the state
    alpha, beta, gamma, d_alpha, d_beta, d_gamma = u

    # no external forces
    S_dirs_f = [zeros(2,1) for _ in range(0, 100; length=10)]

    # calculate the Mass matrix, h and Q_ext
    M_all, h_all, Q_all = get_M_h_Q_single(alpha, beta, gamma, d_alpha, d_beta, d_gamma, p, S_dirs_f)
    ddq = get_ddq(M_all, h_all, Q_all)

    # define derivatives
    du[1], du[2], du[3] = u[4], u[5], u[6]
    du[4], du[5], du[6] = ddq[1], ddq[2], ddq[3]
end


"""
hardcoded matrix of all stab eqautions
"""
function get_stab_forces(S_dirs, F_ext)
    # redefine parameters
    S1_dir, S2_dir, S3_dir, S4_dir, S5_dir, S6_dir, S7_dir, S8_dir, S9_dir, S10_dir = S_dirs
    F_ext_x, F_ext_y = F_ext
    
    # calculate stab forces -> using point eqautions
    M_stab = [
        S1_dir[1] -S2_dir[1] 0 0 0 0 0 0 0 0;
        S1_dir[2] -S2_dir[2] 0 0 0 0 0 0 0 0;
        0 S2_dir[1] S3_dir[1] 0 0 -S6_dir[1] 0 0 0 0;
        0 S2_dir[2] S3_dir[2] 0 0 -S6_dir[2] 0 0 0 0;
        0 0 -S3_dir[1] -S4_dir[1] -S5_dir[1] 0 0 0 0 0;
        0 0 -S3_dir[2] -S4_dir[2] -S5_dir[2] 0 0 0 0 0;
        0 0 0 0 0 S6_dir[1] -S7_dir[1] -S8_dir[1] 0 0;
        0 0 0 0 0 S6_dir[2] -S7_dir[2] -S8_dir[2] 0 0;
        0 0 0 0 0 0 0 S8_dir[1] -S9_dir[1] -S10_dir[1];
        0 0 0 0 0 0 0 S8_dir[2] -S9_dir[2] -S10_dir[2]
    ]
    b_stab = [-F_ext_x; -F_ext_y; 0; 0; 0; 0; 0; 0; 0; 0]

    F_stab = inv(M_stab) * b_stab
    return F_stab
end


"""calculate the normed direction"""
function get_dir(A, B)
    return (B - A) / (norm(B - A))
end

"""get rot matrix"""
function rotx(angle)
    return [
        cos(angle) -sin(angle);
        sin(angle) cos(angle);
    ]
end


"""get the position of Pₒ from P₂"""
function f2t0(P₂, Tₒ₂, diff_20ₒ)
    Pₒ = Tₒ₂ * P₂ + diff_20ₒ
    return Pₒ
end


"""seafty for stab lengths"""
function apply_length_safety(sol, p)

    ah_all = 0

    for sol_i in sol.u
        a, b, c = sol_i[1:3]
        all_points = forward_circuits(a, b, c, p)
        ah_all += get_ah_loss(all_points, p)
    end 

    return ah_all / length(sol)
end


"""calculate the safety for each length option"""
function get_ah_loss(all_points, p)
    l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = safety_p(p)
    A_all, B_all, C_all, D_all, E_all = all_points

    h_a = calculate_h(A_all[2], l_s10, A_all[4], l_s9)
    h_b = calculate_h(B_all[3], l_s8, B_all[5], l_s7)
    h_c = calculate_h(C_all[2], l_s5, C_all[4], l_s4)
    h_d = calculate_h(D_all[2], l_s6, D_all[4], l_s3)
    h_e = calculate_h(E_all[3], l_s2, E_all[5], l_s1)

    return h_a + h_b + h_c + h_d + h_e

end

"""compute the security loss for a single linkage"""
function calculate_h(pₒ, rₒ, p₁, r₁; mode=-1)
    d2 = (p₁ - pₒ)' * (p₁ - pₒ)
    d = sqrt(d2)
    a = (rₒ^2 - r₁^2 + d^2) / (2 * d)
    h2 = rₒ^2 - a^2
    h = sqrt(h2)

    # define diff_center
    p₂ = pₒ + (a / d) * (p₁ - pₒ)

    x_inter = p₂[1] + mode * (h / d) * (p₁[2] - pₒ[2])
    y_inter = p₂[2] - mode * (h / d) * (p₁[1] - pₒ[1])

    p_mid = [x_inter, y_inter] 
    theta = get_angle(p_mid, pₒ, p₁)

    pₒ₁ = pₒ + 0.5* (p₁ - pₒ)

    diff_center = (pₒ₁ - p₂)' * (pₒ₁ - p₂)

    return (10 / h) + 0.2 * diff_center + (10 / theta)
end

function get_angle(p1, p2, p3)
    a = [p1[1] - p2[1], p1[2] - p2[2]]  
    b = [p1[1] - p3[1], p1[2] - p3[2]]
    
    theta = acos((a'b) / (norm(a)*norm(b)))
    return theta
end

"""
determine the midpoint given two points and lengths
https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
"""
function get_two_point_center(pₒ, rₒ, p₁, r₁, mode=-1)
    # predefine distance d, a, h
    d = sqrt((p₁ - pₒ)' * (p₁ - pₒ))
    a = (rₒ^2 - r₁^2 + d^2) / (2 * d)
    h = sqrt(rₒ^2 - a^2)

    # define p2
    p₂ = pₒ + (a / d) * (p₁ - pₒ)

    # define intersection coordinates
    x_inter = p₂[1] + mode * (h / d) * (p₁[2] - pₒ[2])
    y_inter = p₂[2] - mode * (h / d) * (p₁[1] - pₒ[1])

    # return p3
    return [x_inter, y_inter]
end


"""
Calculate the 4link: A
    A3
 l_a2=l_s10  l_a3=l_s9
A2       A4

    A1 (=PIP)    
"""
function get_4link_A(l_s9, l_s10, a₂, a₃, b₂, b₃, T₂₃)
    # calculate the edge positions
    A1₂ = [a₂, 0]
    A2₂ = [a₂, 0] + T₂₃ * [a₃ / 2, b₃]
    A4₂ = [a₂ / 2, b₂]

    # determine position of top point
    A3₂ = get_two_point_center(A4₂, l_s9, A2₂, l_s10)
    
    # return all four coordinates
    return A1₂, A2₂, A3₂, A4₂
end



"""
Calculate the 5link: B
        B4
     l_b3=l_s8  
    B3        l_b4=l_s7
               
B2                  B5

    B1 (=DIP)    
"""
function get_5link_B(l_s7, l_s8, a₁, a₁₂, b₁₂, T₂₁, A3₂, A4₂)    
    # calculate the edge positions
    B1₂ = [0, 0]
    B2₂ = A4₂
    B3₂ = A3₂
    B5₂ = T₂₁ * [a₁₂ - a₁, b₁₂]

    # determine position of top point
    B4₂ = get_two_point_center(B5₂, l_s7, B3₂, l_s8)
    
    # return all five coordinates
    return B1₂, B2₂, B3₂, B4₂, B5₂
end


"""
Calculate the 4link: C
         C3
 l_c2=l_s5   l_c3 = l_s4
C2                C4

         C1 (=PIP)    
"""
function get_4link_C(l_s4, l_s5, aₒ₂, a₁₁, bₒ₂, b₁₁, Tₒ₁)
    # calculate the edge positions
    C1ₒ = [0, 0]
    C2ₒ = Tₒ₁ * [a₁₁, b₁₁]
    C4ₒ = [-aₒ₂, bₒ₂]

    # determine position of top point
    C3ₒ = get_two_point_center(C4ₒ, l_s4, C2ₒ, l_s5)
    
    # return all four coordinates
    return C1ₒ, C2ₒ, C3ₒ, C4ₒ
end



"""
Calculate the 5link: D
            D3
   l_d2=l_s6   l_d3=l_s3      
D2                  D4
              D5
    D1    
"""
function get_5link_D(l_s3, l_s6, B4ₒ, B5ₒ, C2ₒ, C3ₒ)
    # calculate the edge positions
    D1ₒ = B5ₒ
    D2ₒ = B4ₒ
    D4ₒ = C3ₒ
    D5ₒ = C2ₒ

    # determine position of top point
    D3ₒ = get_two_point_center(D4ₒ, l_s3, D2ₒ, l_s6)
    
    # return all four coordinates
    return D1ₒ, D2ₒ, D3ₒ, D4ₒ, D5ₒ
end


"""
Calculate the 5link: E
        E4
   l_e4=l_s2    l_e5=l_s1   
E3             E5
    E2    
        E1    
"""
function get_5link_E(l_s1, l_s2, aₒ₁, bₒ₁, C3ₒ, C4ₒ, D3ₒ)   
    # calculate the edge positions
    E1ₒ = C4ₒ
    E2ₒ = C3ₒ
    E3ₒ = D3ₒ
    E5ₒ = [-aₒ₁, bₒ₁]

    # determine position of top point
    E4ₒ = get_two_point_center(E5ₒ, l_s1, E3ₒ, l_s2)
    
    # return all four coordinates
    return E1ₒ, E2ₒ, E3ₒ, E4ₒ, E5ₒ
end
