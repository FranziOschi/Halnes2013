import numpy as np
from scipy.integrate import odeint
from Parameters import *


class Astro_multi_compartment(object):

    def __init__(self, params):

        # converts dict with parameters to variables
        self.__dict__.update(params)

        # morphology
        self.X = np.linspace(0, self.length, self.N) # position along the x-axis

        # duration and time step of simulation
        self.tspan = np.arange(0,self.time, self.dt)

        # set stimulus duration and length of input zone
        self.tstart = 100.
        self.tstop = 400.
        self.comp_start = int(self.N * 0.1)

        # initialize the system
        init_Na = self.Na_0 * np.ones(self.X.shape) # mM
        init_Cl = self.Cl_0 * np.ones(self.X.shape) # mM
        init_K = self.K_0 * np.ones(self.X.shape) # mM
        init_Na_o = self.Na_o_0 * np.ones(self.X.shape) # mM
        init_Cl_o = self.Cl_o_0 * np.ones(self.X.shape) # mM
        init_K_o = self.K_o_0 * np.ones(self.X.shape) # mM
        init = np.hstack((init_Na, init_Cl, init_K, init_Na_o, init_Cl_o, init_K_o))

        # simulate spatial astro
        sol = odeint(self.spatial_astro, init, self.tspan,tcrit=[self.tstart, self.tstop])

        # transfer solution od ode system into single variables
        self.Na = sol[:,:self.N]
        self.Cl = sol[:,self.N:2*self.N]
        self.K = sol[:,2*self.N:3*self.N]
        self.Na_o = sol[:,3*self.N:4*self.N]
        self.Cl_o = sol[:,4*self.N:5*self.N]
        self.K_o = sol[:,5*self.N:6*self.N]

    def input(self, t):
        if t > self.tstart and t < self.tstop:
            max = np.ones(len(self.X))*self.Imax_const
            max[self.comp_start:] = 0
            return max
        else:
            return np.zeros(len(self.X))


    def spatial_astro(self, state, tspan):

        # initial values
        Na = state[:self.N]
        Cl = state[self.N:2*self.N]
        K = state[2*self.N:3*self.N]
        Nao = state[3*self.N:4*self.N]
        Clo = state[4*self.N:5*self.N]
        Ko = state[5*self.N:6*self.N]

        # finite difference method for calculation of pde: matrix for first derivative
        I_1a = np.zeros((self.N,self.N))
        for i in range(0,len(I_1a)-1):
            I_1a[i,i] = -1.
            I_1a[i,i+1] = 1.
        I_1a /= (self.h)
        I_1a[-1,-1] = -1./(-self.h)
        I_1a[-1,-2] = 1./(-self.h)

        I_1b = np.zeros((self.N,self.N))
        for i in range(0,len(I_1b)-1):
            I_1b[i,i] = -1.
            I_1b[i,i-1] = 1.
        I_1b /= (-self.h)
        I_1b[0,0] = -1./(self.h)
        I_1b[0,1] = 1./(self.h)
        I_1b[-1,-1] = -1./(-self.h)
        I_1b[-1,-2] = 1./(-self.h)

        # membrane potential
        X_oZ_o = -((self.O_m*self.C_m*self.V_0)/(self.a_o)) - self.F*(self.K_o_0 + self.Na_o_0 - self.Cl_o_0)
        V = -((self.a_o)/(self.C_m * self.O_m))*((self.F*(Ko + Nao - Clo)) + X_oZ_o)

        # transmemrane currents
        E_K = self.psifac*np.log(Ko/K)
        E_K_0 = self.psifac*np.log(self.K_o_0/self.K_0)
        E_Na = self.psifac*np.log(Nao/Na)
        E_Cl = self.psifac*np.log(Cl/Clo)
        delta_V_mil = (V - E_K)*1000. #in mV
        E_K_0_mil = E_K_0*1000.
        V_mil = V*1000.
        f_Kir = np.sqrt(Ko/self.K_o_0) * ((1. + np.exp(18.4/42.4))/(1. + np.exp((delta_V_mil + 18.5)/(42.5)))) * \
                ((1. + np.exp(-(118.6 + E_K_0_mil)/(44.1)))/(1. + np.exp(-(118.6 + V_mil)/(44.1))))
        P = self.P_max * ((Na**1.5)/(Na**1.5 + self.K_mN_NKA**1.5)) * ((Ko)/(Ko + self.K_mK))

        # transmembrane flux densities
        J_K_m = ((self.gK*f_Kir)/self.F) * (V - E_K) - 2.*P
        J_Na_m = (self.gNa/self.F) * (V - E_Na) + 3.*P
        J_Cl_m = (-self.gCl/self.F) * (V - E_Cl)

        # diffusive flux
        J_KiD = -(self.D_K/(self.lamb_intra**2)) * I_1a.dot(K)
        J_NaiD = -(self.D_Na/(self.lamb_intra**2)) * I_1a.dot(Na)
        J_CliD = -(self.D_Cl/(self.lamb_intra**2)) * I_1a.dot(Cl)
        J_KoD = -(self.D_K/(self.lamb_extra**2)) * I_1a.dot(Ko)
        J_NaoD = -(self.D_Na/(self.lamb_extra**2)) * I_1a.dot(Nao)
        J_CloD = -(self.D_Cl/(self.lamb_extra**2)) * I_1a.dot(Clo)

        # intra- and extracellular resistivity
        r_o = (self.psifac*(self.lamb_extra**2))/(self.F*(self.D_Na*Nao+self.D_K*Ko+self.D_Cl*Clo))
        r_i = (self.psifac*(self.lamb_intra**2))/(self.F*(self.D_Na*Na+self.D_K*K+self.D_Cl*Cl))

        # current densities due to diffusion
        i_odiff = self.F*(self.z_K*J_KoD + self.z_Na*J_NaoD + self.z_Cl*J_CloD)
        i_idiff = self.F*(self.z_K*J_KiD + self.z_Na*J_NaiD + self.z_Cl*J_CliD)

        # calculate intra- and extracellular membrane voltage
        dVidx = (I_1a.dot(V) + ((r_o*self.a_i*i_idiff)/self.a_o) + (r_o*i_odiff))*((1. + ((r_o*self.a_i)/(r_i*self.a_o)))**-1)
        dVodx = (-I_1a.dot(V) + ((r_i*self.a_o*i_odiff)/self.a_i) + (r_i*i_idiff))*((1. + ((r_i*self.a_o)/(r_o*self.a_i)))**-1)

        # field flux
        J_KiV = -((self.D_K*self.z_K)/((self.lamb_intra**2) * self.psifac)) * (K*dVidx)
        J_NaiV = -((self.D_Na*self.z_Na)/((self.lamb_intra**2) * self.psifac)) * (Na*dVidx)
        J_CliV = -((self.D_Cl*self.z_Cl)/((self.lamb_intra**2) * self.psifac)) * (Cl*dVidx)
        J_KoV = -((self.D_K*self.z_K)/((self.lamb_extra**2) * self.psifac)) * (Ko*dVodx)
        J_NaoV = -((self.D_Na*self.z_Na)/((self.lamb_extra**2) * self.psifac)) * (Nao*dVodx)
        J_CloV = -((self.D_Cl*self.z_Cl)/((self.lamb_extra**2) * self.psifac)) * (Clo*dVodx)

        # axial flux densities
        J_Ki = J_KiD + J_KiV
        J_Nai = J_NaiD + J_NaiV
        J_Cli = J_CliD + J_CliV
        J_Ko = J_KoD + J_KoV
        J_Nao = J_NaoD + J_NaoV
        J_Clo = J_CloD + J_CloV

        # set boundary conditions
        J_Ki[0] = J_Ki[-1] = 0.
        J_Nai[0] = J_Nai[-1] = 0.
        J_Cli[0] = J_Cli[-1] = 0.
        J_Ko[0] = J_Ko[-1] = 0.
        J_Nao[0] = J_Nao[-1] = 0.
        J_Clo[0] = J_Clo[-1] = 0.

        # input
        self.Imax = self.input(tspan)
        ft = self.Imax - self.k_dec*(Ko - self.K_o_0)

        #define differential equations
        dKdt = -(self.O_m/self.a_i)*(J_K_m) - I_1b.dot(J_Ki)
        dNadt = -(self.O_m/self.a_i)*(J_Na_m) - I_1b.dot(J_Nai)
        dCldt = -(self.O_m/self.a_i)*(J_Cl_m) - I_1b.dot(J_Cli)
        dKodt = (self.O_m/self.a_o)*(J_K_m + ft) - I_1b.dot(J_Ko)
        dNaodt = (self.O_m/self.a_o)*(J_Na_m - ft) - I_1b.dot(J_Nao)
        dClodt = (self.O_m/self.a_o)*(J_Cl_m) - I_1b.dot(J_Clo)

        return np.hstack((dNadt, dCldt, dKdt, dNaodt, dClodt, dKodt))
