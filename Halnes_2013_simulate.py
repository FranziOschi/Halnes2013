import numpy as np
import matplotlib.pyplot as plt
from Parameters import p
from Astro import Astro_multi_compartment

def simulate(params = p):

    # simulate system
    astro = Astro_multi_compartment(params = params)

    # extract simulation results
    K = astro.K
    Na = astro.Na
    Cl = astro.Cl
    K_o = astro.K_o
    Na_o = astro.Na_o
    Cl_o = astro.Cl_o

    return K, Na, Cl, K_o, Na_o, Cl_o

# run the simlation and sve the output
K, Na, Cl, K_o, Na_o, Cl_o = simulate()

np.savez('output.npz', K = K, Na = Na, Cl = Cl,
         K_o = K_o, Na_o = Na_o, Cl_o = Cl_o)

