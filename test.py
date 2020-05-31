# Simple test for arenas_model.py with artificial data
# Run this script to see the model in action with artificial parameters and initial conditions.

import numpy as np
import matplotlib.pyplot as plt
from arenas_model import *
from  data_handling.parameters import get_params_arenas


def main():
    # Set number of compartiments (S,E,A,I,H,Rᴵ,Rᴴ,D -> 8)
    NC = 8

    ##--SET ARTIFICIAL POPULATIONS, AREAS, AND MOBILITIES--##

    # Artificial population
    N = 100_000

    ##--SET MODEL PARAMETERS--##

    # Arenas parameters (χ is repeated for the new model with one more compartiment)
    params = get_params_arenas()

    ## Containtment parameters
    tc = 10
    params[-2] = tc
    tf = 10
    params[-1] = tf

    ##--SET ARTIFICIAL INITIAL CONDITIONS--##
    # I'll put a small amount of asymptomatic cases. The rest will be all susceptible.

    # fraction of asymptomatics
    A_frac = 0.01
    # initial conditions (NC vector)
    x0 = np.zeros( NC )
    x0[0] = 1 - A_frac # susceptibles
    x0[2] = A_frac # asymptomatics

    ###-------------###
    ###--RUN MODEL--###
    ###-------------###

    # Set number of days for simulation
    T = 40
    # Obtain model's flow
    flow = iterate_model(x0, T, params)
    # Transform densities to absolute population
    flow = flow * N

    ###--SOME PLOTS--###

    # Plot the aggregate dynamics

    S  =  [flow[t][0] for t in range(flow.shape[0])]
    E  =  [flow[t][1] for t in range(flow.shape[0])]
    A  =  [flow[t][2] for t in range(flow.shape[0])]
    I  =  [flow[t][3] for t in range(flow.shape[0])]
    H  =  [flow[t][4] for t in range(flow.shape[0])]
    Rᴵ =  [flow[t][5] for t in range(flow.shape[0])]
    Rᴴ =  [flow[t][6] for t in range(flow.shape[0])]
    D  =  [flow[t][7] for t in range(flow.shape[0])]

    plt.plot(S, label='S')
    plt.plot(E, label='E')
    plt.plot(A, label='A')
    plt.plot(I, label='I')
    plt.plot(H, label='H')
    plt.plot(Rᴵ, label='Rᴵ')
    plt.plot(Rᴴ, label='Rᴴ')
    plt.plot(D, label='D')

    # containment dates
    plt.axvline([tc], label='containtment', c='black', lw=1)
    plt.axvline([tc+tf], label='release', c='black', lw=1, ls='--')

    plt.title('Aggregate of all cases')
    plt.xlabel('Days since t0')
    plt.xlabel('Number of people')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
