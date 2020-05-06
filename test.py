# Simple test for arenas_model.py with artificial data
# Run this script to see the model in action with artificial parameters and initial conditions.

import numpy as np
import matplotlib.pyplot as plt
from arenas_model import *

if __name__ == '__main__':

    # Set model dimensions. The model has NC compartiments, NP patches and NG age stratas
    NP = 10
    NG = 3
    NC = 7

    ##--SET ARTIFICIAL POPULATIONS, AREAS, AND MOBILITIES--##

    # Artificial population distribution
    N = 100_000
    n_ig = np.random.rand( NP, NG ) * N # number of people

    # Artificial area of regions
    s_i = np.random.rand(NP) * 1000 # km^2

    # Artificial column stochastic mobility matrix
    R_ij = np.random.rand( NP,NP )
    R_ij = R_ij / R_ij.sum(axis=0)


    ##--SET MODEL PARAMETERS--##

    # Arenas parameters
    import arenas_params as ap
    params = [
        ap.β, ap.kg, ap.η, ap.αg, ap.ν, ap.μg, ap.γg, ap.ωg, ap.ψg, ap.χg, n_ig, R_ij, ap.Cgh, ap.ξ, ap.pg, ap.σ, ap.κ0, ap.ϕ, ap.tc, ap.tf
    ]

    ## Containtment parameters
    tc = 10
    params[-2] = tc
    tf = 10
    params[-1] = tf

    ##--PRECOMPUTE ADDITIONAL PARAMETERS--##

    from ext_params import get_ext_params

    # get some necessary parameters for computation
    pg = params[14]
    ξ  = params[13]
    kg = params[1]
    one_minus_pg = 1 - pg # recurrent term

    ## Set external fixed parameters
    ext_params = get_ext_params( n_ig, s_i, R_ij, pg, one_minus_pg, ξ, kg )


    ##--SET ARTIFICIAL INITIAL CONDITIONS--##
    # I'll put a small amount of asymptomatic cases. The rest will be all susceptible.

    # fraction of asymptomatics
    A_frac = 0.01
    # initial conditions (NCxNPxNG matrix)
    x0 = np.zeros( [NC, NP, NG] )
    x0[0] = 1 - A_frac # susceptibles
    x0[2] = A_frac # asymptomatics

    ###-------------###
    ###--RUN MODEL--###
    ###-------------###

    # Set number of markovian steps
    T = 40
    # Obtain model's flow
    flow = iterate_model(x0, T, params, ext_params)

    ###--SOME PLOTS--###

    # 1.- Plot specific patch for specific age strata (still missing its respective population)
    i = 0 # patch
    g = 1 # age strata

    S_ig =  [flow[t][0][i][g] for t in range(flow.shape[0])]
    E_ig =  [flow[t][1][i][g] for t in range(flow.shape[0])]
    A_ig =  [flow[t][2][i][g] for t in range(flow.shape[0])]
    I_ig =  [flow[t][3][i][g] for t in range(flow.shape[0])]
    H_ig =  [flow[t][4][i][g] for t in range(flow.shape[0])]
    R_ig =  [flow[t][5][i][g] for t in range(flow.shape[0])]
    D_ig =  [flow[t][6][i][g] for t in range(flow.shape[0])]

    plt.plot(S_ig, label='S')
    plt.plot(E_ig, label='E')
    plt.plot(A_ig, label='A')
    plt.plot(I_ig, label='I')
    plt.plot(H_ig, label='H')
    plt.plot(R_ig, label='R')
    plt.plot(D_ig, label='D')

    # containment dates
    plt.axvline([tc], label='containtment', c='black', lw=1)
    plt.axvline([tc+tf], label='release', c='black', lw=1, ls='--')

    plt.title('Individual city and age strata')
    plt.legend()
    plt.show()

    # Plot the aggregate dynamics (still missing their respective populations; this is just a sum of densities)

    S =  [flow[t][0].sum() for t in range(flow.shape[0])]
    E =  [flow[t][1].sum() for t in range(flow.shape[0])]
    A =  [flow[t][2].sum() for t in range(flow.shape[0])]
    I =  [flow[t][3].sum() for t in range(flow.shape[0])]
    H =  [flow[t][4].sum() for t in range(flow.shape[0])]
    R =  [flow[t][5].sum() for t in range(flow.shape[0])]
    D =  [flow[t][6].sum() for t in range(flow.shape[0])]

    plt.plot(S, label='S')
    plt.plot(E, label='E')
    plt.plot(A, label='A')
    plt.plot(I, label='I')
    plt.plot(H, label='H')
    plt.plot(R, label='R')
    plt.plot(D, label='D')

    # containment dates
    plt.axvline([tc], label='containtment', c='black', lw=1)
    plt.axvline([tc+tf], label='release', c='black', lw=1, ls='--')

    plt.title('Aggregate of all cases')
    plt.legend()
    plt.show()
