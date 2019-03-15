# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 20:00:02 2019

@author: chaeg
"""

"""Program to calculate the temperature rise as a function of time resulting 
from frictional heating, using formula from Imado et al. 2003."""

from scipy.special import erf
from scipy.special import expi
import numpy as np
import matplotlib.pyplot as plt

#kinematic parameters
g=9.81 #SI
h=1

#paper points to mew decreasing with higher contact pressure so should try values
V = (0.5*g*h*(3)**(0.5))**(0.5)
#material properties
E=200*10**9
k_s = 40.1 #thermal conductivity W/mK
k_p = 0.05
sigma = 2*10**9 #Pa
kappa_p = 8.2*10**-8 #m^2/s
kappa_m = 1.172*10**-5
b =0.001  # contact radius in m with the PAPER, to be adjusted when calculate more precise
#b=1*10**-6 # apporx. for steel on steel 
pi = 3.14
Area = pi*(b/2)**2
tau = 400*10**3 # using value for steel
#mew = (tau*Area**2)/sigma this is for static fric. from fundamentals of sliding frict.
#very sensitive to this, need to think of other ways to think/ experimental ways
x_0 = 6.75*10**-5
alpha = (x_0/0.005)**2

Pe_p = (V*b)/(2*kappa_p)
Pe_m = (V*b)/(2*kappa_m)

eta = 1 / (1 + (kappa_m*(Pe_m)**(0.5))/(kappa_p*(Pe_p)**(0.5)))



def mew_2(n):
    X= alpha**(n/2)
    Y= alpha**(1-n/2) /4
    Z = 1/(np.pi*(1+n/2))
    W = np.arcsin(alpha) -(alpha*(1-alpha))**0.5
    return X*(Y + (1/alpha)*(Z)*W)
#^^^^^this gives negative values so it can't be correct

def T(t):
    #mew = (tau*(Area*np.sin(pi*t/10*10**-6))**2)/(sigma*(np.sin(pi*t/10*10**-6))**(3/2))
    mew = (tau/sigma)*(pi*sigma/(E*2))**2 
    #from daw kwei leu 2011, approx taking 2x first term in mew_d eqn. 30
    #mew = mew_2(0.6)
    #using eqn 30 in terms of alpha w alpha = alpha max and n=0.2 from tabors value
    A=2*eta*sigma*mew*V/kappa_p
    B=(kappa_p*t/pi)**(0.5)
    C = erf(2*b/(2*(kappa_p*t)**(0.5)))
    D = (b/(2*pi*(kappa_p*t)**(0.5)))*expi((-1*b**2)/(4*kappa_p*t))
    return A*B*(C-D)

time=np.logspace(-9,-4,10**6, dtype = float)

#want the value of time when x/2 temp. reaches auto ignition 230 ie. when function = 506 

temp = []
 
for j in time:
        sum_result = T(j)
        temp.append(sum_result)
        if 251 > sum_result > 250:
            print(j)
            
#plt.plot(time,temp,'r')    
plt.semilogx(time,temp,'r')
plt.savefig('Temp_by_friction.png')
plt.title('Slab Temperature for Fricitonal Heating')
plt.xlabel('Time/s')
plt.ylabel('Temperature/K')