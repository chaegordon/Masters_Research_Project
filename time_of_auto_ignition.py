# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 12:51:24 2019

@author: chaeg
"""

"""This is a code to find the time of ignition of the paper when it collides with
the metal half plane, subject to constant gas temperature. Now that the time is
comparable to the collision time need to account for temp variation at the boundary
this should be next goal. I think also should try and do smth similar for the other
initiation mechanisms where I llike think about the temp and BC and also the adjusted
thickness of the paper"""

from scipy.special import erfc
import numpy as np
from numpy import linspace
from numpy import sum
from numpy import array
import matplotlib.pyplot as plt

#in units of nanometers

a=0.675 #for steel
#a=0.9985 #for saphire 
#a=0.902 #for glass
l=1
kappa = 0.09
x=0.5

V=7300 #for steel
#V=4610, #saphirre
#V=7400 #, if glass remains elastic


def b(x):
    return erfc(x)

def sumfunc(i,t):
  a=0.675
  kappa = 90
  #print(i)
  #kappa is 10**-6 and t is in seconds
  #trying l,x in 10 micrometers which comes from inverting mcnamara to get thickness ratio
  #and initial thickness of paper is 0.05mm
  y_1 = erfc(((1+2*i)*l+x)*4*10**-6/(2*(kappa*(t)*10**-6)**(0.5)))
  y_2 = a*erfc(((2*i+1)*l-x)*4*10**-6/(2*(kappa*(t)*10**-6)**(0.5))) 
  return sum(V*(a**(i))*(y_1 - y_2))

#writing the sum
n = linspace(0,999,10**3,dtype = float)


time=np.logspace(-6,-3,10**3, dtype = float)

#want the value of time when x/2 temp. reaches auto ignition 230 ie. when function = 506 

temp = []
 
for j in time:
        sum_result = sumfunc(n,j)
        temp.append(sum_result)
        if 300 > sum_result > 250:
            print(j)
            
    
plt.semilogx(time,temp,'r')
#plt.plot(time,temp,'r')   
"""The size of the paper makes a very big differenc to the end result, however,
it is thought that  it will assume it min thickness almost immediatedly as E of
paper is very small. I made some assumptions to allow T to be fixed at x=-l, I
think the assumption that this layer is held at the gas temp. is good
but could potenitally add the time dependence of this layer to the solution in
in time"""    
    