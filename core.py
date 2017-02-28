
# coding: utf-8

# In[1]:

import numpy as np
import scipy.integrate
import constants


# In[13]:

Nc      = 100
cp_core = 860
k_core  = 100
M_core  = 2e24
P_cmb   = 135e9
gamma   = 1.5
GG=6.67408e-11

r_core = np.zeros(Nc)

def core_density(r):
    r0, r1, r2, r3=12581e0, -0.000198444e3, -8.94481e-05, -2.14483e-8
    rho = r0 + r1*r/1e3 + r2*r**2/1e6 + r3*r**3/1e9
    return rho

def core_gravity(r,rho):
    arg = rho*r**2
    ans = scipy.integrate.cumtrapz(arg, r, initial=0)
    print('divid by zero')
    g   = 4*np.pi*GG*ans/r**2
    g[0] = 0
    return g

def core_pressure(r,rho,g):
    arg = rho*g
    ans = scipy.integrate.cumtrapz(arg, r, initial=0)
    p = ans[-1]-ans + P_cmb
    return p
    
def core_adiabat(N,r,Tcmb):
    N    = N-1
    T1, T2, T3 = 3.00531e-06, -2.57791e-08, 2.8425e-13
    Tcen = Tcmb / (1 + T1*r[N]/1e3 + T2*r[N]**2/1e6 + T3*r[N]**3/1e9)
    T    = Tcen * (1 + T1*r/1e3 +   T2*r**2/1e6 +   T3*r**3/1e9)
    dTdr = Tcen * (    T1       + 2*T2*r/1e3    + 3*T3*r**2/1e6)
    dTdr = dTdr/1e3
    return T, dTdr
 
def get_hstat_profiles(N, r_l, r_t):
    r_core   = constants.radius(r_l,r_t,N)
    rho_core = core_density( r_core)
    g_core   = core_gravity( r_core,rho_core)
    p_core   = core_pressure(r_core,rho_core,g_core)
    return r_core, rho_core, g_core, p_core

def core_secular(N,r,rho,T):
    Tc = T[N-1]
    arg = cp_core * rho * T * r**2
    Qs  = -4 * np.pi * scipy.integrate.cumtrapz(arg, r, initial=0) / Tc
    arg = cp_core * rho * (T/Tc - 1) * r**2
    Es  = -4 * np.pi * scipy.integrate.cumtrapz(arg, r, initial=0) / Tc
    return Qs[N-1], Es[N-1]
    
def core_condent(N,r,T,dTdr):
    arg = k_core*(r*dTdr/T)**2
    Ek  = 4 * np.pi * scipy.integrate.cumtrapz(arg, r, initial=0)   
    return Ek[N-1]

def cool_core(Qcmb, Tc):
    r_core, rho_core, g_core, p_core = get_hstat_profiles(Nc,0,constants.r_cmb)
    T_core, dTdr_core                = core_adiabat(Nc,r_core,Tc)
    Qs, Es = core_secular(Nc, r_core, rho_core, T_core)
    Ek     = core_condent(Nc, r_core, T_core, dTdr_core)

    dTcdt = Qcmb/Qs
    EJ    = Es*dTcdt - Ek
    return EJ, dTcdt


# In[ ]:



