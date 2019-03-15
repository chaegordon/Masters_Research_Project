#Aim to write promgram which takes data imported from a spectrometer and fits
#BB curve to it to determine the temperature of the emmision region.

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import h,c,k

#want to create  the line space of the wavelength in integers so redefine lamda 
#within the function in terms of nanometers

def B_body(lam,T):
    lam= lam*10**-9
    return 2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))


""" For future reference to import the data 
https://s3.amazonaws.com/assets.datacamp.com/blog_assets/Cheat+Sheets/Importing_Data_Python_Cheat_Sheet.pdf
"""

#T_1 and T_2 lie either side of the characteristic temperature we expect for 
#adiabatic compression

wa = np.linspace(1, 2, 10**5)   # wavelengths in nm
T1 = 3000.
T2 = 5000.
y1 = B_body(wa, T1)
y2 = B_body(wa, T2)
ytot = y1 + y2

"""http://python4esac.github.io/fitting/example_blackbody.html but be careful
as I think it fits the addition of TWO curves and I think we just want one curve
see https://www.scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
for ONE curve..."""