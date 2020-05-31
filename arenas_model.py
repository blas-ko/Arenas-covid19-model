## Main script to run the Arenas epidemic model of [1,2] aggregating spatial and age couplings.

import numpy as np
# from helper_functions import *

def runtest():
    '''
    Runs the test.py.
    '''
    from test import main
    main()

# Define one markov step
def markov_step(x, M):
    '''
    Computes one step in the markovian model. We assume that x_{t+1} = M x_t, where M = M(x_t, t).

    Inputs:
    `x`: state variables (S,E,A,I,H,R,D)
    `M`: list of non-zero elements of the transition matrix M
    '''
    # state variables
    S,E,A,I,H,Rᴵ,Rᴴ,D = x
    # interaction terms
    M_SS, M_ES, M_EE, M_AE, M_AA, M_IA, M_II, M_HI, M_HH, M_RᴵI, M_RᴵRᴵ, M_RᴴH, M_RᴴRᴴ, M_DH, M_DD = M

    return np.array([
        M_SS * S,                 # S
        M_ES * S + M_EE * E,      # E
        M_AE * E + M_AA * A,      # A
        M_IA * A + M_II * I,      # I
        M_HI * I + M_HH * H,      # H
        M_RᴵI * I + M_RᴵRᴵ * Rᴵ,  # Rᴵ
        M_RᴴH * H + M_RᴴRᴴ * Rᴴ,  # Rᴴ
        M_DH * H + M_DD * D       # D
    ])

def iterate_model(x0, T, params):
    '''
    Solves the markovian model for `T` time steps (days) for the initial conditions `x0` and the set of parameters `params` and `ext_params`.

    Inputs:
    `x0`: list with the initial compartiment densities (S0, E0, A0, I0, H0, R0, D0)
    `params`': list of parameters in the same order than in Arenas report [2]: (β, kg, ηg, αg, ν, μg, γg, ωg, ψg, χg, n_ig, σ, κ0, ϕ, tc, tf, κf)

    Output:
    `flow`: 7-dimensional time series. Each dimension corresponds to S(t), E(t), A(t), I(t), H(t), R(t), D(t) respectively.
    '''

    ## READING ##

    # Read parameters (1-D treatment. In the general treatment, suffix `g` indicates an NG-sized vector)
    β = params[0]
    k = params[1]
    η = params[2]
    α = params[3]
    ν = params[4]
    μ = params[5]
    γ = params[6]
    ω = params[7]
    ψ = params[8]
    χᴵ = params[9]
    χᴴ = params[10]
    n = params[11] # population
    # containtment params
    σ  = params[12]
    κ0 = params[13]
    ϕ  = params[14]
    tc = params[15]
    tf = params[16]
    κf = params[17]

    # Compute infection probability
    Π_t = Π_1D( x0[2]+ ν*x0[3], β, k)

    # Compute interaction terms
    M_SS   = 1 - Π_t
    M_ES   = Π_t
    M_EE   = 1 - η
    M_AE   = η
    M_AA   = 1 - α
    M_IA   = α
    M_II   = γ * (1 - μ) + (1 - γ) * (1 - χᴵ)
    M_HI   = γ * μ
    M_HH   = ω * (1 - ψ) + (1 - ω) * (1 - χᴴ)
    M_RᴵI  = (1 - γ) * χᴵ
    M_RᴵRᴵ = 1
    M_RᴴH  = (1 - ω) * χᴴ
    M_RᴴRᴴ = 1
    M_DH   = ω * ψ
    M_DD   = 1

    ## Non-zero interactions for transition-like matrix
    M = [M_SS, M_ES, M_EE, M_AE, M_AA, M_IA, M_II, M_HI, M_HH, M_RᴵI, M_RᴵRᴵ, M_RᴴH, M_RᴴRᴴ, M_DH, M_DD]

    ## PREALLOCATION
    flow = np.zeros( [T+1, *x0.shape] )
    flow[0,:] = x0

    # Case when containtment happens before initial conditions
    # This is only approximate, as we don't know how many people where susceptible at containtment date
    if tc < 0:
        # Lower avg. number of contacts as a function of containtment
        k = (1-κ0)*k + κ0*(σ-1)

        # Compute infection probability
        Π_t = Π_1D( x0[2]+ν*x0[3], β, k )

        ## Contained people (susceptible + recovered)
        C_tc = ( x0[0]+x0[5] )**σ

        # update dynamic interaction terms
        M[0] = (1 - Π_t)*(1 - (1 - ϕ)*κ0*C_tc)
        M[1] = Π_t*(1 - (1 - ϕ)*κ0*C_tc)

    x_old = x0

    ## MODEL DYNAMICS
    for t in range(T):

        # Take markov step
        x_new = markov_step(x_old, M)
        # Update flow vector
        flow[t+1,:] = x_new

        # Containtment
        if t+1 == tc:
            # Lower avg. number of contacts as a function of containtment
            k = (1-κ0)*k + κ0*(σ-1)

            # Compute infection probability
            Π_t = Π_1D( x_new[2]+ν*x_new[3], β, k )

            ## Contained people (susceptible + recovered)
            C_tc = ( x_new[0]+x_new[5] )**σ

            # update dynamic interaction terms
            M[0] = (1 - Π_t)*(1 - (1 - ϕ)*κ0*C_tc)
            M[1] = Π_t*(1 - (1 - ϕ)*κ0*C_tc)
        else:

            # Update probability of transmission
            Π_t = Π_1D( x_new[2]+ν*x_new[3], β, k )

            # update dynamic interaction terms
            M[0] = (1 - Π_t)
            M[1] = Π_t

        # end of containtment
        if t+1 == tc+tf:
            # mid-agers (1-D treatment)
            k = ( k - κ0*(σ-1) ) / (1 - κ0)
            k = (1-κf)*k + κf*(σ-1)  # new normality containment

            # Compute infection probability
            Π_t = Π_1D( x_new[2]+ν*x_new[3], β, k )

            # update dynamic interaction terms
            M[0] = (1 - Π_t)*(1 + (1 - ϕ)*κf*C_tc)
            M[1] = Π_t*(1 + (1 - ϕ)*κf*C_tc)

        x_old = x_new

    return flow

## Helper functions
# Probability of infection for 1D treatment of the model
def Π_1D(ρ, β, k):
    '''
    Returns the probability of infection per patch per age strata per day considering the effective mobility patterns.
    This function is designed for the aggregate model, where NP = NG = 1.
    The output is a scalar.
    '''
    return 1 - (1 - β)**(k*ρ)
