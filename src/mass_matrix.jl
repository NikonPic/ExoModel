"""
Contains all results of sympy extracted solutions
"""


"""
M h and q as tuple return
"""
function get_M_h_Q_single(alpha, beta, gamma, d_alpha, d_beta, d_gamma, p, S_dirs)
    l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = safety_p(p)
    # get all relevant forces
    S5x, S5y, S7x, S7y, S9x, S9y, S10x, S10y = S_dirs[5][1], S_dirs[5][2], S_dirs[7][1], S_dirs[7][2], S_dirs[9][1], S_dirs[9][2], S_dirs[10][1],S_dirs[10][2]

    # Mass matrix
    M1 = 0.5*a1.^2 .* m1 + 1.0*a1.^2 .* m2 + 1.0*a1.^2 .* m3 + 1.0*a1.*a2.*m2.*cos(beta) + 2.0*a1.*a2.*m3.*cos(beta) + 2.0*a1.*a3.*m3.*cos(beta + gamma) + 0.0833333333333333*a1.*m1 + 0.5*a2.^2 .*m2 + 1.0*a2.^2 .*m3 + 2.0*a2.*a3.*m3.*cos(gamma) + 0.0833333333333333*a2.*m2 + 1.25*a3.^2 .*m3 + 0.0833333333333333*a3.*m3 + 0.25*m1.*r1.^2 + 0.25*m2.*r2.^2 + 0.25*m3.*r3.^2
    M2 = 0.5*a1.*a2.*m2.*cos(beta) + 1.0*a1.*a2.*m3.*cos(beta) + 1.0*a1.*a3.*m3.*cos(beta + gamma) + 0.5*a2.^2 .*m2 + 1.0*a2.^2 .*m3 + 2.0*a2.*a3.*m3.*cos(gamma) + 0.0833333333333333*a2.*m2 + 1.25*a3.^2 .*m3 + 0.0833333333333333*a3.*m3 + 0.25*m2.*r2.^2 + 0.25*m3.*r3.^2
    M3 = m3.*(1.0*a1.*a3.*cos(beta + gamma) + 1.0*a2.*a3.*cos(gamma) + 1.25*a3.^2 + 0.083*a3 + 0.25*r3.^2)

    M4 = 0.5*a1.*a2.*m2.*cos(beta) + 1.0*a1.*a2.*m3.*cos(beta) + 1.0*a1.*a3.*m3.*cos(beta + gamma) + 0.5*a2.^2 .*m2 + 1.0*a2.^2 .*m3 + 2.0*a2.*a3.*m3.*cos(gamma) + 0.0833333333333333*a2.*m2 + 1.25*a3.^2 .*m3 + 0.0833333333333333*a3.*m3 + 0.25*m2.*r2.^2 + 0.25*m3.*r3.^2
    M5 = 0.5*a2.^2 .*m2 + 1.0*a2.^2 .*m3 + 2.0*a2.*a3.*m3.*cos(gamma) + 0.0833333333333333*a2.*m2 + 1.25*a3.^2 .*m3 + 0.0833333333333333*a3.*m3 + 0.25*m2.*r2.^2 + 0.25*m3.*r3.^2
    M6 = m3.*(1.0*a2.*a3.*cos(gamma) + 1.25*a3.^2 + 0.0833333333333333*a3 + 0.25*r3.^2)

    M7 = m3.*(1.0*a1.*a3.*cos(beta + gamma) + 1.0*a2.*a3.*cos(gamma) + 1.25*a3.^2 + 0.0833333333333333*a3 + 0.25*r3.^2)
    M8 = m3.*(1.0*a2.*a3.*cos(gamma) + 1.25*a3.^2 + 0.0833333333333333*a3 + 0.25*r3.^2)
    M9 = m3.*(1.25*a3.^2 + 0.0833333333333333*a3 + 0.25*r3.^2)

    # h vector
    h1 = -1.0*a1.*a2.*d_alpha.*d_beta.*m2.*sin(beta) - 2.0*a1.*a2.*d_alpha.*d_beta.*m3.*sin(beta) - 0.5*a1.*a2.*d_beta.^2 .*m2.*sin(beta) - 1.0*a1.*a2.*d_beta.^2 .*m3.*sin(beta) - 2.0*a1.*a3.*d_alpha.*d_beta.*m3.*sin(beta + gamma) - 2.0*a1.*a3.*d_alpha.*d_gamma.*m3.*sin(beta + gamma) - 1.0*a1.*a3.*d_beta.^2 .*m3.*sin(beta + gamma) - 2.0*a1.*a3.*d_beta.*d_gamma.*m3.*sin(beta + gamma) - 1.0*a1.*a3.*d_gamma.^2 .*m3.*sin(beta + gamma) - 0.5*a1.*g.*m1.*cos(alpha) - 1.0*a1.*g.*m2.*cos(alpha) - 1.0*a1.*g.*m3.*cos(alpha) - 2.0*a2.*a3.*d_alpha.*d_gamma.*m3.*sin(gamma) - 2.0*a2.*a3.*d_beta.*d_gamma.*m3.*sin(gamma) - 1.0*a2.*a3.*d_gamma.^2 .*m3.*sin(gamma) - 0.5*a2.*g.*m2.*cos(alpha + beta) - 1.0*a2.*g.*m3.*cos(alpha + beta) - 1.0*a3.*g.*m3.*cos(alpha + beta + gamma) + 1.0*alpha.*c1 - 1.0*alpha_0.*c1
    h2 = 0.5*a1.*a2.*d_alpha.^2 .*m2.*sin(beta) + 1.0*a1.*a2.*d_alpha.^2 .*m3.*sin(beta) + 1.0*a1.*a3.*d_alpha.^2 .*m3.*sin(beta + gamma) - 2.0*a2.*a3.*d_alpha.*d_gamma.*m3.*sin(gamma) - 2.0*a2.*a3.*d_beta.*d_gamma.*m3.*sin(gamma) - 1.0*a2.*a3.*d_gamma.^2 .*m3.*sin(gamma) - 0.5*a2.*g.*m2.*cos(alpha + beta) - 1.0*a2.*g.*m3.*cos(alpha + beta) - 1.0*a3.*g.*m3.*cos(alpha + beta + gamma) + 1.0*beta.*c2 - 1.0*beta_0.*c2
    h3 = 1.0*a1.*a3.*d_alpha.^2 .*m3.*sin(beta + gamma) + 1.0*a2.*a3.*d_alpha.^2 .*m3.*sin(gamma) + 2.0*a2.*a3.*d_alpha.*d_beta.*m3.*sin(gamma) + 1.0*a2.*a3.*d_beta.^2 .*m3.*sin(gamma) - 1.0*a3.*g.*m3.*cos(alpha + beta + gamma) + 1.0*c3.*gamma - 1.0*c3.*gamma_0

    # Q_ext vector
    Q1 = -d1.*d_alpha - sqrt(a11.^2 + b11.^2).*(S5x.*sin(alpha - atan(b11./a11)) + S5y.*cos(alpha - atan(b11./a11))) - sqrt(a12.^2 + b12.^2).*(S7x.*sin(alpha - atan(b12./a12)) + S7y.*cos(alpha - atan(b12./a12)))
    Q2 = -d2.*d_beta - sqrt(a2.^2 + 4*b2.^2).*(S9x.*sin(alpha + beta - atan(2*b2./a2)) + S9y.*cos(alpha + beta - atan(2*b2./a2)))/2
    Q3 = -d3.*d_gamma - sqrt(a3.^2 + 4*b3.^2).*(S10x.*sin(alpha + beta + gamma - atan(2*b3./a3)) + S10y.*cos(alpha + beta + gamma - atan(2*b3./a3)))/2

    # concatenate
    M_all = (M1, M2, M3, M4, M5, M6, M7, M8, M9)
    h_all = (h1, h2, h3)
    Q_all = (Q1, Q2, Q3)

    return M_all, h_all, Q_all
end


"""
substitution for the inverse mass matrix solution
"""
function get_ddq(M_all, h_all, Q_all)
    M1, M2, M3, M4, M5, M6, M7, M8, M9 = M_all
    h1, h2, h3 = h_all
    Q1, Q2, Q3 = Q_all


    ddq1 = ((Q1 - h1).*(M5.*M9 - M6.*M8) - (Q2 - h2).*(M4.*M9 - M6.*M7) + (Q3 - h3).*(M4.*M8 - M5.*M7))./(M1.*M5.*M9 - M1.*M6.*M8 - M2.*M4.*M9 + M2.*M6.*M7 + M3.*M4.*M8 - M3.*M5.*M7)
    ddq2 = (-(Q1 - h1).*(M2.*M9 - M3.*M8) + (Q2 - h2).*(M1.*M9 - M3.*M7) - (Q3 - h3).*(M1.*M8 - M2.*M7))./(M1.*M5.*M9 - M1.*M6.*M8 - M2.*M4.*M9 + M2.*M6.*M7 + M3.*M4.*M8 - M3.*M5.*M7)
    ddq3 = ((Q1 - h1).*(M2.*M6 - M3.*M5) - (Q2 - h2).*(M1.*M6 - M3.*M4) + (Q3 - h3).*(M1.*M5 - M2.*M4))./(M1.*M5.*M9 - M1.*M6.*M8 - M2.*M4.*M9 + M2.*M6.*M7 + M3.*M4.*M8 - M3.*M5.*M7)

    ddq_all = (ddq1, ddq2, ddq3)

    return ddq_all
end


"""
solution of the stab forces
"""
function extract_symbolic_stab_forces(S_dirs, F_ext)
    # extract input
    S1_dir, S2_dir, S3_dir, S4_dir, S5_dir, S6_dir, S7_dir, S8_dir, S9_dir, S10_dir = S_dirs
    F_ext_x, F_ext_y = F_ext
    S1x, S1y = S1_dir
    S2x, S2y = S2_dir
    S3x, S3y = S3_dir
    S4x, S4y = S4_dir
    S5x, S5y = S5_dir
    S6x, S6y = S6_dir
    S7x, S7y = S7_dir
    S8x, S8y = S8_dir
    S9x, S9y = S9_dir
    S10x, S10y = S10_dir

    # symbolic solution
    F_1  =  (F_ext_y*S2x - F_ext_x*S2y)/(S1x*S2y - S2x*S1y)
    F_2  =  (F_ext_y*S1x - F_ext_x*S1y)/(S1x*S2y - S2x*S1y)
    F_3  = -((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S6y - S6x*S2y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y))
    F_4  =  ((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S6y - S6x*S2y)*(S3x*S5y - S5x*S3y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S4x*S5y - S5x*S4y))
    F_5  = -((F_ext_y*S1x - F_ext_x*S1y)*(S3x*S4y - S4x*S3y)*(S2x*S6y - S6x*S2y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S4x*S5y - S5x*S4y))
    F_6  = -((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S3y - S3x*S2y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y))
    F_7  = -((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S3y - S3x*S2y)*(S6x*S8y - S8x*S6y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S7x*S8y - S8x*S7y))
    F_8  =  ((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S3y - S3x*S2y)*(S6x*S7y - S7x*S6y))/((S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S7x*S8y - S8x*S7y))
    F_9  =  ((F_ext_y*S1x - F_ext_x*S1y)*(S10y*S8x - S10x*S8y)*(S2x*S3y - S3x*S2y)*(S6x*S7y - S7x*S6y))/((S10y*S9x - S10x*S9y)*(S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S7x*S8y - S8x*S7y))
    F_10 = -((F_ext_y*S1x - F_ext_x*S1y)*(S2x*S3y - S3x*S2y)*(S6x*S7y - S7x*S6y)*(S8x*S9y - S9x*S8y))/((S10y*S9x - S10x*S9y)*(S1x*S2y - S2x*S1y)*(S3x*S6y - S6x*S3y)*(S7x*S8y - S8x*S7y))
    # gather solution
    F_stab = (F_1, F_2, F_3, F_4, F_5, F_6, F_7, F_8, F_9, F_10)

    return F_stab
end