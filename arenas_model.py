## Main script to run the Arenas epidemic model of [1,2]

import numpy as np
from helper_functions import *

# Define one markov step
def markov_step(x, M):
    '''
    Computes one step in the markovian model. We assume that x_{t+1} = M x_t, where M = M(x_t, t).

    Inputs:
    `x`: state variables (S,E,A,I,H,R,D)
    `M`: list of non-zero elements of the transition matrix M
    '''
    # state variables
    S,E,A,I,H,R,D = x
    # interaction terms
    M_SS, M_ES, M_EE, M_AE, M_AA, M_IA, M_II, M_HI, M_HH, M_RI, M_RH, M_RR, M_DH, M_DD = M

    return np.array([
        M_SS * S,
        M_ES * S + M_EE * E,
        M_AE * E + M_AA * A,
        M_IA * A + M_II * I,
        M_HI * I + M_HH * H,
        M_RI * I + M_RH * H + M_RR * R,
        M_DH * H + M_DD * D
    ])

def iterate_model(x0, T, params, ext_params):
    '''
    Solves the markovian model for `T` time steps (days) for the initial conditions `x0` and the set of parameters `params` and `ext_params`.

    Inputs:
    `x0`: list with the initial compartiment densities (S0, E0, A0, I0, H0, R0, D0)
    `params`': list of parameters in the same order than in Arenas report [2]: (β, kg, ηg, αg, ν, μg, γg, ωg, ψg, χg, n_ig, R_ij, C_gh, ξ, pg, σ, κ0, ϕ)
    `ext_params`: list of pre-computed quantities necessary for the model: (zg * kg, f(n_i_eff/s_i), n_ig_eff)

    Note: `ext_params` are not included in params for efficiency reasons.
    `ext_params` includes the normalization factor times the number of contacts `zk_g`, the effective density population vector `f(x_i)` and the effective population matrix `n_ig_eff`.
    For the bayes approach, we wouldn't want to compute these quantities for every simulation.
    '''

    ## READING ##

    # Read parameters
    β = params[0]
    kg = params[1]#[1]
    ηg = params[2]
    αg = params[3]#[1]
    ν  = params[4]
    μg = params[5]#[1]
    γg = params[6]#[1]
    ωg = params[7]
    ψg = params[8]
    χg = params[9]
    n_ig = params[10]#[0][0]
    R_ij = params[11]
    C_gh = params[12]
    ξ  = params[13]
    pg = params[14]#[1]
    σ  = params[15]
    κ0 = params[16]
    ϕ  = params[17]
    tc = params[18]
    tf = params[19] # containtment end
    # Recurrent computation
    one_minus_pg = 1 - pg

    # Read external parameters related to population (and number of contacts) They don't change at all
    zk_g, f_i, n_ig_eff = ext_params

    ## PRECOMPUTING
    # Compute probability of infection
    # >1-D
    ρ_t_eff = get_ρ_ig_eff(x0[2]+ν*x0[3] , n_ig, n_ig_eff, R_ij, C_gh, pg, one_minus_pg)
    Q_t = Q_ig( zk_g, f_i, ρ_t_eff )
    P_t = P_ig(β, Q_t)
    Π_t = Π_ig( P_t, R_ij, pg, one_minus_pg)

    # 1-D treatment
    # Π_t = Π_1D( x0[2]+ ν*x0[3], β, kg)

    # Compute interaction terms
    M_SS = 1 - Π_t
    M_ES = Π_t
    M_EE = 1 - ηg
    M_AE = ηg
    M_AA = 1 - αg
    M_IA = αg
    M_II = 1 - μg
    M_HI = μg * γg
    M_HH = ωg * (1 - ψg) + (1 - ωg)*(1 - χg)
    M_RI = μg * (1 - γg)
    M_RH = (1 - ωg) * χg
    M_RR = 1
    M_DH = ωg * ψg
    M_DD = 1

    M = [M_SS, M_ES, M_EE, M_AE, M_AA, M_IA, M_II, M_HI, M_HH, M_RI, M_RH, M_RR, M_DH, M_DD]

    ## PREALLOCATION
    flow = np.zeros( [T+1, *x0.shape] )
    flow[0,:] = x0

    x_old = x0

    ## MODEL DYNAMICS
    for t in range(T):

        # Take markov step
        x_new = markov_step(x_old, M)
        # Update flow vector
        flow[t+1,:] = x_new

        # Containtment
        if t+1 == tc:
            # mid-agers (1-D treatment)
            kg[1] = (1-κ0)*kg[1] + κ0*(σ-1)
            pg = (1-κ0)*pg

            # In the model, the mobility matrix R_ij remains constant even after containment
            ρ_t_eff = get_ρ_ig_eff( x_new[2]+ν*x_new[3] , n_ig, n_ig_eff, R_ij, C_gh, pg, one_minus_pg)
            Q_t = Q_ig( zk_g, f_i, ρ_t_eff )
            P_t = P_ig(β, Q_t)
            Π_t = Π_ig( P_t, R_ij, pg, one_minus_pg)

            ## Contained people (susceptible + recovered)
            C_tc = ( get_n_i( (x_new[0]+x_new[5]) * n_ig ) / get_n_i(n_ig) )**σ
            C_tc = C_tc.reshape(len(C_tc), 1 ) # This reshape is only for dimension coherency

            # 1-D treatment
            # Π_t = Π_1D( x_new[2]+ν*x_new[3], β, kg )

            # 1-D treatment
            # C_tc = ( x_new[0]+x_new[5] )**σ

            # update dynamic interaction terms
            M[0] = (1 - Π_t)*(1 - (1 - ϕ)*κ0*C_tc)
            M[1] = Π_t*(1 - (1 - ϕ)*κ0*C_tc)
        else:
            # compute new prob. of infection
            ρ_t_eff = get_ρ_ig_eff( x_new[2]+ν*x_new[3] , n_ig, n_ig_eff, R_ij, C_gh, pg, one_minus_pg)
            Q_t = Q_ig( zk_g, f_i, ρ_t_eff )
            P_t = P_ig(β, Q_t)
            Π_t = Π_ig( P_t, R_ij, pg, one_minus_pg)

            # 1-D treatment
            # Π_t = Π_1D( x_new[2]+ν*x_new[3], β, kg )

            # update dynamic interaction terms
            M[0] = (1 - Π_t)
            M[1] = Π_t

        # end of containtment
        if t+1 == tc+tf:
            # mid-agers (1-D treatment)
            kg = (kg - κ0*(σ-1))/(1 - κ0)
            pg = pg/(1 - κ0)

            ρ_t_eff = get_ρ_ig_eff( x_new[2]+ν*x_new[3] , n_ig, n_ig_eff, R_ij, C_gh, pg, one_minus_pg)
            Q_t = Q_ig( zk_g, f_i, ρ_t_eff )
            P_t = P_ig(β, Q_t)
            Π_t = Π_ig( P_t, R_ij, pg, one_minus_pg)

            # 1-D treatment
            # Π_t = Π_1D( x_new[2]+ν*x_new[3], β, kg )

            # update dynamic interaction terms
            M[0] = (1 - Π_t)*(1 + (1 - ϕ)*κ0*C_tc)
            M[1] = Π_t*(1 + (1 - ϕ)*κ0*C_tc)

        x_old = x_new

    return flow
