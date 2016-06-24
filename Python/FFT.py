# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:17:23 2016

@author: James
"""

from numpy import fft
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

site = ['US-WCr', 'US-Whs', 'US-Wkg']
year_range = ['1999-2014', '2007-2014', '2004-2014']

filename = list()
for i in xrange(len(site)):
    z = 'FLX_'+site[i]+'_FLUXNET2015_FULLSET_HH_'+year_range[i]+'_1-1.csv'
    filename.append(z)
data = pd.read_csv(filename[1],na_values=['-9999'])

n = data['NEE_VUT_MEAN'] # Number of data points
dx = 5.0 # Sampling period (in meters)
x = dx*n # x coordinates
w1 = 100.0 # wavelength (meters)
w2 = 20.0 # wavelength (meters)
fx = np.sin(2*np.pi*x/w1) + 2*np.cos(2*np.pi*x/w2) # signal
Fk = fft.fft(fx)/n # Fourier coefficients (divided by n)
nu = fft.fftfreq(n,dx) # Natural frequencies
Fk = fft.fftshift(Fk) # Shift zero freq to center
nu = fft.fftshift(nu) # Shift zero freq to center
f, ax = plt.subplots(3,1,sharex=True)
ax[0].plot(nu, np.real(Fk)) # Plot Cosine terms
ax[0].set_ylabel(r'$Re[F_k]$', size = 'x-large')
ax[1].plot(nu, np.imag(Fk)) # Plot Sine terms
ax[1].set_ylabel(r'$Im[F_k]$', size = 'x-large')
ax[2].plot(nu, np.absolute(Fk)**2) # Plot spectral power
ax[2].set_ylabel(r'$\vert F_k \vert ^2$', size = 'x-large')
ax[2].set_xlabel(r'$\widetilde{\nu}$', size = 'x-large')
plt.show()