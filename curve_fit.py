# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 20:28:14 2019

@author: chaeg
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


def func(x,a,b):
   return a*x - a*np.log(b)

"""by taking logs of the functions the data can be linearised as above however,
because of 0's in the logs the points on the plateaus of the distribution are
not defined and so this method is not suitable to my project due to the small
number of remaining data points"""

def wiebull(x,a,b):
    return 1 - np.exp(-1*(x/b)**a)

#need to renormalise for the different maximum (factor of 0.8)
def wiebull_1(x,a,b):
    return 0.8 - 0.8*np.exp(-1*(x/b)**a)
    
df = pd.read_csv('data_1.csv', index_col = 0)

x_1 = df["x'"]
y_1 = df["y'"]

x_2 = df["x''"]
y_2 = df["y''"]

x_3 = df["x'''"]
y_3 = df["y'''"]

x_plot = np.linspace(0,1.5,100)

popt_1, pcov_1 = curve_fit(wiebull_1, x_1, y_1)

popt_2, pcov_2 = curve_fit(wiebull, x_2, y_2)

popt_3, pcov_3 = curve_fit(wiebull, x_3, y_3)

print(popt_1)
#CP = plt.scatter(x,y)
"""Also renormalised here or can use wiebull_1 in the plot same difference"""
plt.plot(x_plot, wiebull_1(x_plot, *popt_1), label = 'normal impact')

plt.plot(x_plot, wiebull(x_plot, *popt_2), label = '46 degrees impact')

plt.plot(x_plot, wiebull(x_plot, *popt_3), label = '30 degrees impact')



plt.title('Drop test Data')
plt.xlabel('Height/m')
plt.ylabel('Probability of Ignition')
plt.legend()
plt.savefig('prob_v_h.png')
