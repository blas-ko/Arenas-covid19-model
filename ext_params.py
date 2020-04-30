import numpy as np
from helper_functions import *
# from arenas_params import pg, ξ, kg

def get_ext_params( n_ig, s_i, R_ij, pg, one_minus_pg, ξ, kg ):

    # number of patches
    try:
        NP = np.array(n_ig).shape[0]
    except:
        NP = 1

    # effective population given the mobility parameters
    n_i_eff = get_n_i_eff(n_ig, R_ij, pg)
    n_ig_eff = get_n_ig_eff(n_ig, R_ij, pg)
    n_g = get_n_g(n_ig)

    # patch related effective density
    f_i = f(n_i_eff/s_i, ξ)
    if NP != 1:
        f_i = f_i.reshape(NP,1)

    # age related normalization factor
    z_g = n_g / np.dot( np.transpose(f_i), n_ig_eff )
    # precoputation of age related fixed params (number of contacts could enter the bayesian formalism later...)
    zk_g = z_g * kg

    ##  EXTERNAL FIXED PARAMETERS
    return [zk_g, f_i, n_ig_eff]
