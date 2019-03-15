# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 12:31:22 2019

@author: chaeg
"""

"""plastic work done within body of metal ball conducting into the paper. Start
with unit and INSTANTANEOUS SOURCE."""

"""check only convolving wrt to t and not wrt x"""


from scipy.special import erf
import scipy.integrate as integrate
import scipy.special as special
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

#Material Parameters

l=25 #micrometers , the actual value is much larger but want to see the shape
K_p = 0.05
kappa_p = 8.2*10**-8
K_m = 40.1
kappa_m =1.172*10**-5
height = 1
mass = 0.5
E_before = 9.81*height*mass
e_res = (0.2)**0.5
E_heat = 0.9*E_before*(e_res**2)
x_dash = 5 # in micrometers 
t_c = 5*10**-6

#constants in C & J
k= (kappa_p/kappa_m)**0.5
sigma = k*(K_m/K_p)
alpha = (sigma - 1)/(sigma + 1)
A = sigma/((1+sigma)*(np.pi*kappa_m))
#Mesh
N = np.linspace(0,99,100,dtype = float)
x=np.linspace(-l,0,10**3)
time=np.logspace(-1,1,10**3, dtype = float)
X,T = np.meshgrid(x,time)


def v_1(t):
    counter =  np.zeros(T.shape)
    for n in N:
      B = A*(t*10**-6)**(-0.5)
      A_1 = (alpha**n)*(np.exp(-1*((((k*x_dash - X +2*n*l)*(10**-6))**(2)/(4*kappa_p*t*10**-6)))))
      A_2 = (alpha**n)*(np.exp(-1*((((k*x_dash + X +2*(n+1)*l)*(10**-6))**(2)/(4*kappa_p*t*10**-6)))))
     
      counter += B*(A_1 - A_2)
      
    return counter

#dont need the E/t as all Q delivered instantaneously
#but do account for the solid angle which means more source travels away
#DO THIS MORE CAREFULLY
temp = (E_heat/2*np.pi)*v_1(T)
temp_per_t = (E_heat/(4*np.pi)*t_c)*v_1(T)
#print(temp.shape)
temp_int = integrate.cumtrapz(temp_per_t,axis =0, initial = 0, dx = 10**-3)

"""Theta defined as 0=295K"""
temperature_re_scaled = temp_int + 295
#print(temp.shape)
print(temp_int.shape)  
"""to get an average across the length of the paper strip may need to a loop
across x values"""

CP = plt.contourf(X,T, temperature_re_scaled)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature with Plastic Heating')
plt.xlabel('distance along paper/10e-6m')
plt.ylabel('time/10e-6s')
plt.savefig('Slab_Temperature_plastic.png')

