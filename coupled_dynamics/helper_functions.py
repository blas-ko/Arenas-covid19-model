## Helper functions for the Arenas et al model [1,2]

import numpy as np

# NG : cardinality of the age strata
# NP : number of patches (regions)

def f(pop_dens_i, ξ):
    '''
    Returns the influence of population density, where
    `pop_dens_i`: effective population in patch i / surface in patch i in km^2
    '''
    return 2 - np.exp(-ξ*pop_dens_i)

# populations
def get_n_g(n_ig):
    '''
    Returns the aggregate population in all patches by age strata.

    Input:
    `n_ig`: NPxNG matrix with the population of patch i and age strata g
    '''
    try: # NP, NG != 1
        return n_ig.sum( axis=0 )
    except:
        return n_ig

def get_n_i(n_ig):
    '''
    Returns the population in every patch aggregating their age stratas.

    Input:
    `n_ig`: NPxNG matrix with the population of patch i and age strata g
    '''
    try: # NP, NG != 1
        return n_ig.sum( axis=1 )
    except:
        try:
            return n_ig.sum( axis=0 )
        except:
            return n_ig

def get_n_i_eff(n_ig, R_ij, pg):
    '''
    Returns the effective population matrix per patch considering the mobility patterns of the model aggregating their age strata.
    '''
    # if pg is a scalar, make it an array so the inner product dimension makes sense
    try:
        pg.shape[0]
    except:
        pg = np.array([pg])

    return np.dot( np.dot( R_ij, n_ig ), pg ) + np.dot( n_ig, 1 - pg )

def get_n_ig_eff(n_ig, R_ij, pg):
    '''
    Returns the effective population matrix per patch per age strata considering the mobility patterns of the model.
    '''
    return (1 - pg) * n_ig  +  pg *  np.dot( R_ij, n_ig )

# densities
def get_ρ_ig_eff(ρ_ig, n_ig, n_ig_eff, R_ij, C_gh, pg, one_minus_pg):
    '''
    Returns the effective compartiment population density matrix per patch per age strata considering the mobility patterns and the contact-by-age matrix of the model, where
    `ρ_ig`: population density matrix of a given compartiment.
    '''
    return np.dot( get_ρ_ig_no_coupling(ρ_ig, n_ig, R_ij, pg, one_minus_pg) / n_ig_eff , np.transpose(C_gh) )

def get_ρ_ig_no_coupling(ρ_ig, n_ig, R_ij, pg, one_minus_pg):
    '''
    Returns the effective compartiment population density matrix per patch per age strata considering the mobility patterns but not the contact-by-age matrix of the model, where
    `ρ_ig`: population density matrix of a given compartiment.
    '''
    nρ_ig = n_ig * ρ_ig
    return one_minus_pg * nρ_ig  +  pg * np.dot( R_ij, nρ_ig  )


def Q_ig(zk_g, f_i, ρ_ig_eff):
    '''
    Returns the expected number of contacts per patch per age strata per day.
    The output is an NPxNG matrix.
    '''
    return zk_g * f_i * ρ_ig_eff

def P_ig( β, Q_t ):
    '''
    Returns the probability of infection per patch per age strata per day only considering the mobility patterns in the average number of contacts.
    The output is an NPxNG matrix.
    '''
    return 1 - (1 - β)**Q_t

def Π_ig( P_t, R_ij, pg, one_minus_pg ):
    '''
    Returns the probability of infection per patch per age strata per day considering the effective mobility patterns.
    The output is an NPxNG matrix.
    '''
    return one_minus_pg * P_t  +  pg *  np.dot( R_ij, P_t )


## Use this when running the aggregate one-dimensional model
# 1D treatment
def Π_1D(ρ, β, kg):
    '''
    Returns the probability of infection per patch per age strata per day considering the effective mobility patterns.
    This function is designed for the aggregate model, where NP = NG = 1.
    The output is a scalar.
    '''
    return 1 - (1 - β)**(kg*ρ)
