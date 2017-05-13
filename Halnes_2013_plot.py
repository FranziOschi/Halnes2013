import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from Parameters import *

params = {#'backend': 'Agg',
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'axes.linewidth' : 1.,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'axes.spines.right' : False,
          'axes.spines.top' : False,
          'xtick.major.size' : 0,
          'ytick.major.size' : 0,
          'font.family': 'sans-serif',
          'mathtext.rm': 'sans',
          'lines.linewidth': 2.0,
          'savefig.dpi' : 300,
          'figure.subplot.wspace' : 0.8,
          'figure.subplot.hspace' : 0.6,
          'legend.fontsize' : 14
          }
rcParams.update(params)

data = np.load('output.npz')
K = data['K']
Na = data['Na']
Cl = data['Cl']
K_o = data['K_o']
Na_o = data['Na_o']
Cl_o = data['Cl_o']

# put all variables from p into the local namespace
locals().update(p)

############################# Figure 4 ##############################################

# membrane potential
X_oZ_o = -((O_m*C_m*V_0)/(a_o)) - F*(K_o_0 + Na_o_0 - Cl_o_0)
V_400 = -((a_o)/(C_m * O_m))*((F*(K_o[int(400/dt),:] + Na_o[int(400/dt),:] - Cl_o[int(400/dt),:])) + X_oZ_o)
V_x0 = -((a_o)/(C_m * O_m))*((F*(K_o[:,0] + Na_o[:,0] - Cl_o[:,0])) + X_oZ_o)


fig = plt.figure(figsize=(8,10))
plt.subplots_adjust(wspace = 0.4, hspace = 0.4)
ax1 = plt.subplot(4,2,1)
I_max = np.zeros(int(time/dt))
I_max[int(100./dt):int(400./dt)] = Imax_const
plt.plot(np.arange(0,time,dt), (I_max)*1e6, 'k', label = r'$\mathsf{j^{in}}$')
plt.plot(np.arange(0,time,dt), -(k_dec*(K_o[:,0] - K_o_0))*1e6, 'k--', label = r'$\mathsf{j^{out}}$')
plt.xticks(np.arange(0,time+1,200.), ['0', '200', '400', '600'])
plt.ylabel(r'$\mathsf{j\/[\mu mol\/m^{-2}\/s^-1]}$')
ax1.legend(bbox_to_anchor=(1.2, 1.05))

plt.subplot(4,2,2)
I_max = np.zeros(N)
I_max[0:int(N*0.1)] = Imax_const
plt.plot(np.linspace(0,length,N)*1e3, (I_max)*1e6, 'k', label = r'$\mathsf{j^{in}}$')
plt.plot(np.linspace(0,length,N)*1e3, -(k_dec*(K_o[int(400/dt),:] - K_o_0))*1e6, 'k--', label = r'$\mathsf{j^{out}}$')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)

ax2 = plt.subplot(4,2,3)
plt.plot(np.arange(0,time,dt), K_o[:,0] - K_o_0, label = r'$\mathsf{K^+}$', color = '#1f77b4')
plt.plot(np.arange(0,time,dt), Na_o[:,0] - Na_o_0, label = r'$\mathsf{Na^+}$', color = '#ff7f0e')
plt.plot(np.arange(0,time,dt), Cl_o[:,0] - Cl_o_0, label = r'$\mathsf{Cl^-}$', color = '#2ca02c')
plt.plot(np.arange(0,time,dt), (K_o[:,0] + Na_o[:,0] - Cl_o[:,0]) + (-((O_m*C_m*V_0)/(a_o*F)) - (K_o_0 + Na_o_0 - Cl_o_0)), label = r'$\mathsf{e^+}$', color = '#d62728')
plt.ylabel(r'$\mathsf{\Delta [k]_{E}\/[mM]}$')
plt.xticks(np.arange(0,time+1,200.), ['0', '200', '400', '600'])
ax2.legend(bbox_to_anchor=(1.2, 0.6))

plt.subplot(4,2,4)
plt.plot(np.linspace(0,length,N)*1e3, K_o[int(400/dt),:] - K_o_0, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, Na_o[int(400/dt),:] - Na_o_0, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, Cl_o[int(400/dt),:] - Cl_o_0, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, (K_o[int(400/dt),:] + Na_o[int(400/dt),:] - Cl_o[int(400/dt),:]) + (-((O_m*C_m*V_0)/(a_o*F)) - (K_o_0 + Na_o_0 - Cl_o_0)), color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)

plt.subplot(4,2,5)
plt.plot(np.arange(0,time,dt), K[:,0] - K_0, color = '#1f77b4')
plt.plot(np.arange(0,time,dt), Na[:,0] - Na_0, color = '#ff7f0e')
plt.plot(np.arange(0,time,dt), Cl[:,0] - Cl_0, color = '#2ca02c')
plt.plot(np.arange(0,time,dt), (K[:,0] + Na[:,0] - Cl[:,0]) + (((O_m*C_m*V_0)/(a_i*F)) - (K_0 + Na_0 - Cl_0)), color = '#d62728')
plt.ylabel(r'$\mathsf{\Delta [k]_{I}\/[mM]}$')
plt.xticks(np.arange(0,time+1,200.), ['0', '200', '400', '600'])

plt.subplot(4,2,6)
plt.plot(np.linspace(0,length,N)*1e3, K[int(400/dt),:] - K_0, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, Na[int(400/dt),:] - Na_0, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, Cl[int(400/dt),:] - Cl_0, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, (K[int(400/dt),:] + Na[int(400/dt),:] - Cl[int(400/dt),:]) + (((O_m*C_m*V_0)/(a_i*F)) - (K_0 + Na_0 - Cl_0)), color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)

plt.subplot(4,2,7)
plt.plot(np.arange(0,time,dt), V_x0*1000., 'k')
plt.xlabel('Time [sec]')
plt.ylabel(r'$\mathsf{v_M\/[mV]}$')
plt.xticks(np.arange(0,time+1,200.), ['0', '200', '400', '600'])

plt.subplot(4,2,8)
plt.plot(np.linspace(0,length,N)*1e3, V_400*1000., 'k')
plt.xlabel('x [mm]')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)

plt.savefig('Halnes_2013_Fig4.png', dpi = 300)
plt.show()


############################# Figure 5 ##############################################

# finite difference method for calculation of pde: matrix for first derivative
I_1a = np.zeros((N,N))
for i in range(0,len(I_1a)-1):
    I_1a[i,i] = -1.
    I_1a[i,i+1] = 1.
I_1a /= (h)
I_1a[-1,-1] = -1./(-h)
I_1a[-1,-2] = 1./(-h)

# membrane potential
X_oZ_o = -((O_m*C_m*V_0)/(a_o)) - F*(K_o_0 + Na_o_0 - Cl_o_0)
V = -((a_o)/(C_m * O_m))*((F*(K_o[int(400/dt),:] + Na_o[int(400/dt),:] - Cl_o[int(400/dt),:])) + X_oZ_o)

# transmemrane currents
P = P_max * ((Na[int(400/dt),:]**1.5)/(Na[int(400/dt),:]**1.5 + K_mN_NKA**1.5)) * ((K_o[int(400/dt),:])/(K_o[int(400/dt),:] + K_mK))
E_K = psifac*np.log(K_o[int(400/dt),:]/K[int(400/dt),:])
E_Cl = psifac*np.log(Cl[int(400/dt),:]/Cl_o[int(400/dt),:])
E_K_0 = psifac*np.log(K_o_0/K_0)
E_Na = psifac*np.log(Na_o[int(400/dt),:]/Na[int(400/dt),:])
delta_V_mil = (V - E_K)*1000. #in mV
E_K_0_mil = E_K_0*1000.
V_mil = V*1000.
f_Kir = np.sqrt(K_o[int(400/dt),:]/K_o_0) * ((1. + np.exp(18.4/42.4))/(1. + np.exp((delta_V_mil + 18.5)/(42.5)))) * \
        ((1. + np.exp(-(118.6 + E_K_0_mil)/(44.1)))/(1. + np.exp(-(118.6 + V_mil)/(44.1))))
P = P_max * ((Na[int(400/dt),:]**1.5)/(Na[int(400/dt),:]**1.5 + K_mN_NKA**1.5)) * ((K_o[int(400/dt),:])/(K_o[int(400/dt),:] + K_mK))

# transmembrane flux densities
J_K_m = ((gK*f_Kir)/F) * (V - E_K) - 2.*P
J_Na_m = (gNa/F) * (V - E_Na) + 3.*P
J_Cl_m = (-gCl/F) * (V - E_Cl)

# diffusive flux
J_KiD = -(D_K/(lamb_intra**2)) * I_1a.dot(K[int(400/dt),:])
J_NaiD = -(D_Na/(lamb_intra**2)) * I_1a.dot(Na[int(400/dt),:])
J_CliD = -(D_Cl/(lamb_intra**2)) * I_1a.dot(Cl[int(400/dt),:])
J_KoD = -(D_K/(lamb_extra**2)) * I_1a.dot(K_o[int(400/dt),:])
J_NaoD = -(D_Na/(lamb_extra**2)) * I_1a.dot(Na_o[int(400/dt),:])
J_CloD = -(D_Cl/(lamb_extra**2)) * I_1a.dot(Cl_o[int(400/dt),:])

# intra- and extracellular resistivity
r_o = (psifac*(lamb_extra**2))/(F*(D_Na*Na_o[int(400/dt),:]+D_K*K_o[int(400/dt),:]+D_Cl*Cl_o[int(400/dt),:]))
r_i = (psifac*(lamb_intra**2))/(F*(D_Na*Na[int(400/dt),:]+D_K*K[int(400/dt),:]+D_Cl*Cl[int(400/dt),:]))

# current densities due to diffusion
i_odiff = F*(z_K*J_KoD + z_Na*J_NaoD + z_Cl*J_CloD)
i_idiff = F*(z_K*J_KiD + z_Na*J_NaiD + z_Cl*J_CliD)

# calculate intra- and extracellular membrane voltage
dVidx = (I_1a.dot(V) + ((r_o*a_i*i_idiff)/a_o) + (r_o*i_odiff))*((1. + ((r_o*a_i)/(r_i*a_o)))**-1)
dVodx = (-I_1a.dot(V) + ((r_i*a_o*i_odiff)/a_i) + (r_i*i_idiff))*((1. + ((r_i*a_o)/(r_o*a_i)))**-1)

# field flux
J_KiV = -((D_K*z_K)/((lamb_intra**2) * psifac)) * (K[int(400/dt),:]*dVidx)
J_NaiV = -((D_Na*z_Na)/((lamb_intra**2) * psifac)) * (Na[int(400/dt),:]*dVidx)
J_CliV = -((D_Cl*z_Cl)/((lamb_intra**2) * psifac)) * (Cl[int(400/dt),:]*dVidx)
J_KoV = -((D_K*z_K)/((lamb_extra**2) * psifac)) * (K_o[int(400/dt),:]*dVodx)
J_NaoV = -((D_Na*z_Na)/((lamb_extra**2) * psifac)) * (Na_o[int(400/dt),:]*dVodx)
J_CloV = -((D_Cl*z_Cl)/((lamb_extra**2) * psifac)) * (Cl_o[int(400/dt),:]*dVodx)

# input
Kinput = np.zeros(N)
Koutput = np.zeros(N)
comp_start = int(N*0.1)
Kinput[0:comp_start] = Imax_const*1e6
Koutput = k_dec*(K_o[int(400/dt),:] - K_o_0)*1e6


fig = plt.figure(figsize=(8,8))
plt.subplots_adjust(wspace = 0.4, hspace = 0.4)
plt.subplot(3,2,1)
plt.plot(np.linspace(0,length,N)*1e3, Kinput-Koutput, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, -Kinput+Koutput, color = '#ff7f0e')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.ylabel(r'$\mathsf{j_k^{in} - j_k^{out}}$')
plt.title('Input/Output')

plt.subplot(3,2,2)
plt.plot(np.linspace(0,length,N)*1e3, 1e6*J_K_m, label = r'$\mathsf{K^+}$', color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, 1e6*J_Na_m, label = r'$\mathsf{Na^+}$', color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, 1e6*J_Cl_m, label = r'$\mathsf{Cl_-}$', color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, 1e6*J_K_m + 1e6*J_Na_m - 1e6*J_Cl_m, label = r'$\mathsf{e^+}$', color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.ylabel(r'$\mathsf{j_kM}$')
plt.title('Transmembrane')
plt.legend(loc = 4, frameon = False)

plt.subplot(3,2,3)
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_KoV*1e6, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_NaoV*1e6, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_CloV*1e6, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_KoV*1e6 + a_o*J_NaoV*1e6 - a_o*J_CloV*1e6, color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.ylabel(r'$\mathsf{a_{E}j_{kE}^f}$')
plt.title('ECS field')


plt.subplot(3,2,4)
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_KiV*1e6, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_NaiV*1e6, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_CliV*1e6, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_KiV*1e6 + a_i*J_NaiV*1e6 - a_i*J_CliV*1e6, color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.ylabel(r'$\mathsf{a_{I}j_{kI}^f}$')
plt.title('ICS field')

plt.subplot(3,2,5)
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_KoD*1e6, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_NaoD*1e6, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_CloD*1e6, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, a_o*J_KoD*1e6 + a_o*J_NaoD*1e6 - a_o*J_CloD*1e6, color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.xlabel('x [mm]')
plt.ylabel(r'$\mathsf{a_{E}j_{kE}^d}$')
plt.title('ECS diffusion')

plt.subplot(3,2,6)
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_KiD*1e6, color = '#1f77b4')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_NaiD*1e6, color = '#ff7f0e')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_CliD*1e6, color = '#2ca02c')
plt.plot(np.linspace(0,length,N)*1e3, a_i*J_KiD*1e6 + a_i*J_NaiD*1e6 - a_i*J_CliD*1e6, color = '#d62728')
plt.xticks(np.linspace(0,length*1e3,4), [ '0' ,  '0.1',  '0.2',  '0.3'])
plt.xlim(0,0.3)
plt.xlabel('x [mm]')
plt.ylabel(r'$\mathsf{a_{I}j_{kI}^d}$')
plt.title('ICS diffusion')

plt.savefig('Halnes_2013_Fig5.png', dpi = 300)
plt.show()
