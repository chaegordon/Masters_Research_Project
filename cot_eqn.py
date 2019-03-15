# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 16:47:43 2019

@author: chaeg
"""

""" program to calculate succesive roots of a transcendental eqn """

from scipy.optimize import minimize
from numpy import *
from numpy import linspace
import matplotlib.pyplot as plt

k=24.3*10**-3
h=10.45
const = k/h
L=1*10**-5
x0=linspace(1,10**4,10**6)

def f_1(x):
    return 1/(tan(x))

def f_2(x):
    return (x*k)/h

#because of how small L is have to rescale axis of the plot to x0*L
"""plt.plot(x0*L,f_1(x0*L),'r')
plt.plot(x0*L,f_2(x0),'g')

plt.show()"
"""
#this plot shows all roots of n>1 can be ignored as so large that attenuation
#is instant

def opt_fun(x, c):
    
    return (1/tan(L*x)- c*x)**2

#don't need to loop as only care about the first root, have multiplied through
    #to get c*x to help it compute

const = k/h
res = minimize(lambda x: opt_fun(x, const), x0=0.001)

# Check if the optimization was successful
print(res.success)
# >> True

# Extract the root from the minimization result
print(res.x[0])

suc=6553.048702797512

def B(x):
    return (2*sin(L*x)*cos(L**2))/(x*L + sin(x*L)*cos(x*L))

print(B(suc))
"""for i in range(len(x0)):
    res = minimize(lambda x: opt_fun(x, const), x0[i])
    # Check if the optimization was successful
      # Extract the root from the minimization result
    for i in range(20):
        print(res.x[0])
        
    if res.success == False:
        break
        print(res.success)"""
        