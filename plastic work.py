# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 11:55:38 2019

@author: chaeg
"""

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


"""need to account for the number of hotspots recieving this Energy"""

"""need to measure final thickness accurately"""

"""may need to solve w a source inside the ball and the T at contact?"""

"""Something not right here"""

l = 25 #quoted thickness before compression 100microm
k = 0.05
kappa = 8.2*10**-8
Temp_max = 438
height = 1
mass = 0.5
E_before = 9.81*height*mass
e_res = (0.2)**0.5
E_heat = 0.9*E_before*(e_res**2)
Area = l**2
t_c = 5*10**-6
# something like this, number_bubbles = 10**-6/10**-8?
F_max = E_heat/(t_c*(Area*10**-12))
Q_dot  = E_heat/t_c

lumped_T=E_heat/((1.4*10**3)*((l*10**-6)**3)*800)
print(lumped_T)


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

temp = F_max*d_Theta(T)

#print(temp.shape)
  
"""to get an average across the length of the paper strip may need to a loop
across x values"""

CP = plt.contourf(X,T, temp)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Plastic Work')
plt.xlabel('distance along paper/10e-6m')
plt.ylabel('time/10e-6s')
plt.savefig('Slab_Work_Plastic_number_of_bubbles.png')

