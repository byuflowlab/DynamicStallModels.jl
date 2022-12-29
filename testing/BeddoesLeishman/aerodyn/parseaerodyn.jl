

function parseaerodyn(entries, items)
    ni = length(items)
    ne = Int(length(entries)/ni)
    mat = zeros(ne, ni)
    
    for i = 1:ne
        idx = (i-1)*ni+1:(i*ni)
        mat[i,:] = entries[idx]
    end
    return mat
end

function getintermediatestates(states, U, a, c, dcndalpha, A1, A2, b1, b2)
    Ka = states[4]
    Kpa = states[10]
    Kq = states[5]
    Kpq = states[11]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)/2  
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 
    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b #Checked that uses slope. 
    Tq = 3*k_q*TI/4 #Equation 1.10b
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q

    alpha = states[1]
    # q = states[3]
    qf = states[29]

    X1 = states[6]
    X2 = states[7]
    X3 = states[8] #EQ 1.16a 
    X4 = states[9] #EQ 1.16a
    alphae = (alpha - alpha0) - X1 - X2

    dcndalpha_circ = dcndalpha/beta
    Ncirc_q = dcndalpha_circ*qf/2 - X3 -X4
    Ncirc_aq = dcndalpha_circ*alphae 

    fpp = states[19]
    fterm = (1 + 2*sqrt(fpp))/3

    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) 

    Cvn = states[23]

    Cpot = dcndalpha_circ*alphae*alpha

    return Cfsn, Cvn, Nnoncirc_aq, Ncirc_q, Ncirc_aq, fpp, Nnoncirc_a, Nnoncirc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot
end

function extractintermediatestates(states, Uvec, c, airfoil; a=343.0)
    n, m = size(states)

    Cfsn = zeros(n)
    Cvn = zeros(n)
    Nnc_aq = zeros(n)
    Nc_q = zeros(n)
    Nc_aq = zeros(n)
    fpp = zeros(n)
    Nnc_a = zeros(n)
    Nnc_q = zeros(n)
    Talpha = zeros(n)
    M = zeros(n)
    k_alpha = zeros(n)
    TI = zeros(n)
    alphae = zeros(n)
    k_q = zeros(n)
    Cpot = zeros(n)

    for i = 1:n
        Cfsn[i], Cvn[i], Nnc_aq[i], Nc_q[i], Nc_aq[i], fpp[i], Nnc_a[i], Nnc_q[i], Talpha[i], M[i], k_alpha[i], TI[i], alphae[i], k_q[i], Cpot[i] = getintermediatestates(states[i,:], Uvec[i], a, c, airfoil.dcndalpha, airfoil.A[1], airfoil.A[2], airfoil.b[1], airfoil.b[2])
    end

    return Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot
end