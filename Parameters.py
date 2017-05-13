# generate dictionary for parameters
p = {}

# simulation duration and timestep
p['time'] = 600. #seconds
p['dt'] = .1 #seconds

# initial concentrations
p['Na_0']    =  15. + 0.189 #mole/meter**3
p['Na_o_0']          = 145. - 0.378 #mole/meter**3
p['K_0']          = 100. - 0.041 #mole/meter**3
p['K_o_0']          = 3. + 0.082 #mole/meter**3
p['Cl_0']          = 5. + 0.145 #mole/meter**3
p['Cl_o_0']          = 134. - 0.29 #mole/meter**3
p['V_0'] = (-85. - 1.4) * 1e-3 #volt

# physical constants
p['F'] = 96500.0 # coulomb / mole
p['R'] =  8.315 # joule / (mole*kelvin)
p['T'] = 298. # kelvin
p['psifac'] = ((p['R']*p['T'])/p['F'])

# geometrical parameters
p['N'] = 100
p['length'] = 3 * 1e-4 #m
p['h'] = p['length']/p['N'] #length of each compartment
p['a_i'] = 0.4 #Tissue volume fraction being astrocytes
p['a_o'] = 0.2 #(Tissue volume fraction being ECS
p['O_m'] = p['a_i']/(5.00 * 1e-8) #meter, Astrocytic membrane area per tissue volume

# membrane parameters
p['gK'] = 16.96 #S/m**2
p['gNa'] = 1. #S/m**2
p['gCl'] = 0.5 #S/m**2
p['C_m'] = 1.0 *1e-2 #farad/meter**2

# diffusion coefficients
p['D_Na'] = 1.33 * 1e-9 #meter**2 * second **-1
p['D_Cl'] = 2.03 * 1e-9 #meter**2 * second **-1
p['D_K'] = 1.96 * 1e-9 #meter**2 * second **-1

# tortuosity
p['lamb_intra'] = 3.2
p['lamb_extra'] = 1.6

# P (NKA)
p['K_mN_NKA'] = 10. #mole/meter**3
p['K_mK'] = 1.5 #mole/meter**3
p['P_max'] = 1.115*1e-6 #mol/(meter**2 * second)

# valence
p['z_Cl'] = -1.
p['z_Na'] = 1.
p['z_K'] = 1.

# influx/outflux input zone
p['k_dec'] = 2.9 * 1e-8 #Output factor (m/s)
p['Imax_const'] = 5.5 * 1e-7 # Input amplitude K/Na-(exchange) (mol/(s m^2))