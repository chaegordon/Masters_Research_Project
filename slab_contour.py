# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 15:04:08 2019

@author: chaeg
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:44:40 2019

@author: chaeg
"""

"""program to calculate the temp. of paper slab with time varying flux on one surface
and a heat sink on the other"""

from scipy.special import erf
import scipy.integrate as integrate
import scipy.special as special
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


"""Can ignore the back flux from paper as will also have forward flux from the
other paper cancelling this"""

l = 4 #micrometers
k = 0.05
kappa = 8.2*10**-8
Temp_max = 438
F_max = 0.68*(5.67*10**-8)*(Temp_max)**4 
t_c = 5*10**-6



N = np.linspace(0,4,5,dtype = float)
x=np.linspace(0*l,l,10**2)
time=np.logspace(-1,0,10**3, dtype = float)
X,T = np.meshgrid(x,time)

"""ierfc is a function defined in appendix 2 of carslaw and jaeger"""

def erfc(x):
    return 1- erf(x)

def ierfc(x):
    WE= np.exp(-x**2/np.sqrt(np.pi))
    COOL = x*erfc(x)
    return WE - COOL

def dE_1(t):
    counter = np.zeros(T.shape)
    for n in N:
      A = ((2*n+1)*l - X)*(10**-6)/(2*(kappa**(0.5)))
      a1 = ((A/(t*10**-6))**(2))
      a7 = np.exp(-((A)**2)/(t*10**-6))
      a2 = (A/(2*(t*10**-6)**(1.5)))
      a8 = erfc(A*(t*10**-6)**-0.5)
      a5 = ((A/(t*10**-6))**(2))
      a4 = (np.exp(-1*(((A)**(2))/(t*10**-6))))
      a6 = 2*(np.pi)**(-0.5)
      counter += a1*a7 + a2*a8 - a6*a4*a5
    return counter

def dE_2(t):
    counter = np.zeros(T.shape)
    for n in N:
      A = ((2*n+1)*l + X)*(10**-6)/(2*(kappa**(0.5)))
      a1 = ((A/(t*10**-6))**(2))
      a7 = np.exp(-((A)**2)/(t*10**-6))
      a2 = (A/(2*(t*10**-6)**(1.5)))
      a8 = erfc(A*((t*10**-6)**(-0.5)))
      a5 = ((A/(t*10**-6))**(2))
      a4 = (np.exp(-1*(((A)**(2))/(t*10**-6))))
      a6 = 2*(np.pi)**(-0.5)
      counter += ( a1*a7 + a2*a8 - a6*a4*a5)
    return counter

#temperature solution in time independent case per unit flux
def v(t):
    counter = np.zeros(T.shape)
    for n in N:
      A= (2/k)*(kappa*(t*10**-6)**(0.5))
      B=((-1)**(n))
      C=ierfc(((2*n+1)*l-X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      D=-1*ierfc(((2*n+1)*l+X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      counter += ( A*B*(C+D))
    return counter

def f(t):
    #room_temp = 295*np.ones(T.shape)
    flux_time = (Temp_max*(np.sin((t*10**-6)*np.pi/t_c))**(0.57))**(0.25)
   #if np.all(np.greater_equal(room_temp,flux_time)):
        #flux = np.zeros(T.shape)
    #else:
        #flux = F_max*(np.sin((t*10**-6)*np.pi/t_c))**(0.57) - 0.68*(5.67*10**-8)*(295)**4 
    return flux_time



def d_Theta(t):
    counter =  np.zeros(T.shape)
    for n in N:
      A= (2/k)*(kappa/(t*10**-6))**(0.5)
      B=((-1)**(n))
      C=ierfc(((2*n+1)*l-X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      D=-1*ierfc(((2*n+1)*l+X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      A_1 = A*B*(C+D)
      E = (2/k)*(kappa*(t*10**-6))**(0.5)
      F = dE_1(t)
      G = dE_2(t)
      A_2 = E*B*(F + G)
      counter += (A_1 + A_2)
      
    return counter

#print(f(T).shape)
#print(d_Theta(T).shape)

#accounts for the time step between points on the convolution, must be 1/len(T)
    """ potential source of error in conflating linear convolution and integral
    convolution"""
delta = 0.001  
temp = signal.fftconvolve(f(T),delta*d_Theta(T), mode = 'same')

"""Theta defined as 0=295K"""
temperature_re_scaled = temp + 295
#print(temp.shape)
  
"""to get an average across the length of the paper strip may need to a loop
across x values"""

CP = plt.contourf(X,T, temperature_re_scaled)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Void Collapse Heating')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')
plt.savefig('Slab Temperature_Paper_number_of_bubbles.png')

