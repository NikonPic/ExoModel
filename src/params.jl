"""bulky function which manages all relevant parameters"""
function deliver_all(p, p_stat)
    # l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10 = p
    # a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = p_stat

    m1, m2, m3 = p
    a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = p_stat


    return l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0
end

function deliver_p_stat()
    p, p_stat = init_params()
    return p_stat
end

"""init"""
function init_params()
    # lengths
    l_s1  = 8.0
    l_s2  = 8.0
    l_s3  = 6.0
    l_s4  = 4.0
    l_s5  = 4.0
    l_s6  = 6.0
    l_s7  = 6.0
    l_s8  = 6.0
    l_s9  = 4.0
    l_s10 = 4.0

    # offset lengths
    a01 = 2.
    a02 = 1.
    a1  = 5.
    a11 = 2.
    a12 = 3.
    a2  = 5.
    a3  = 5.
    b01 = 2.
    b02 = 1.
    b11 = 1.
    b12 = 1.
    b2  = 1.
    b3  = 1.

    # mass
    m1 = 0.1
    m2 = 0.05
    m3 = 0.03

    # radius finger elements
    r1 = 1.
    r2 = 1.
    r3 = 1.

    # gravity
    g = 9.81

    # stiffnes
    c1 = 50.
    c2 = 50.
    c3 = 100.

    # zero angles
    alpha_0, beta_0, gamma_0 = 0.5, 0.5, 0.8

    # damping 
    d1 = 100.
    d2 = 100.
    d3 = 100.

    p = [l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0]
    return p
end


"""
this function handles as a sefty of all parameters and ensures smooth behaviour
"""
function safety_p(p; fac=20)
    l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0 = p

    p_init = init_params()
    l_s1i, l_s2i, l_s3i, l_s4i, l_s5i, l_s6i, l_s7i, l_s8i, l_s9i, l_s10i, a01i, a02i, a1i, a11i, a12i, a2i, a3i, b01i, b02i, b11i, b12i, b2i, b3i, m1i, m2i, m3i, gi, c1i, c2i, c3i, d1i, d2i, d3i, r1i, r2i, r3i, alpha_0i, beta_0i, gamma_0i = p_init

    # lengths
    l_s1  = l_s1i  + fac * tanh(l_s1)
    l_s2  = l_s2i  + fac * tanh(l_s2)
    l_s3  = l_s3i  + fac * tanh(l_s3)
    l_s4  = l_s4i  + fac * tanh(l_s4)
    l_s5  = l_s5i  + fac * tanh(l_s5)
    l_s6  = l_s6i  + fac * tanh(l_s6)
    l_s7  = l_s7i  + fac * tanh(l_s7)
    l_s8  = l_s8i  + fac * tanh(l_s8)
    l_s9  = l_s9i  + fac * tanh(l_s9) 
    l_s10 = l_s10i + fac * tanh(l_s10)

    # offset lengths
    a01 = a01i
    a02 = a02i
    a1  = a1i
    a11 = a11i
    a12 = a12i
    a2  = a2i
    a3  = a3i
    b01 = b01i
    b02 = b02i
    b11 = b11i
    b12 = b12i
    b2  = b2i
    b3  = b3i

    # mass
    m1 = m1i
    m2 = m2i
    m3 = m3i

    # radius finger elements
    r1 = r1i
    r2 = r2i
    r3 = r3i

    # gravity
    g = gi

    # stiffnes
    c1 = c1i
    c2 = c2i
    c3 = c3i

    # zero angles
    alpha_0 = alpha_0i
    beta_0  = beta_0i
    gamma_0 = gamma_0i

    # damping 
    d1 = d1i
    d2 = d2i
    d3 = d3i

    p = [l_s1, l_s2, l_s3, l_s4, l_s5, l_s6, l_s7, l_s8, l_s9, l_s10, a01, a02, a1, a11, a12, a2, a3, b01, b02, b11, b12, b2, b3, m1, m2, m3, g, c1, c2, c3, d1, d2, d3, r1, r2, r3, alpha_0, beta_0, gamma_0]
    return p
end