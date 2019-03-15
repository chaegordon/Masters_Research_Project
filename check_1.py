# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 19:15:34 2019

@author: chaeg
"""
from scipy.optimize import minimize
from numpy import *
from numpy import linspace
import matplotlib.pyplot as plt

x0=linspace(1,10,10**3)

def f_1(x):
    
    return (1/tan(x)- x)**2

def f_2(x):
    
    return (1/(x*tan(x))- 1)**2

plt.plot(x0,f_1(x0),'r')
plt.plot(x0,f_2(x0),'g')