# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 17:19:25 2018

@author: Niklas
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as con

# Very rough estimate of parameters to use in the radar equation

# Nomex honeycomb material properties from HRH-10 here:
# https://www.hexcel.com/Resources/DataSheets/Honeycomb

# EM properties
er = 1.15                       # Relative permittivity 1.15
c = con.c/np.sqrt(er)           # Speed of light in material
f = 60e9                        # EM frequency 60 GHz
k = 2*np.pi*f/c                 # EM wavenumber

# Acoustic properties
Is = 100                        # 100 W/m2 acoustic intensity
vs = 3000                       # 3000 m/s acoustic wave velocity
rho = 80                        # Nomex honeycomb HRH-10 density 80 kg/m3
s0 = np.sqrt(2*Is/rho/vs**3)    # Acoustic strain amplitude (from Saleh)

# Dimensions
Lx = 10e-2
Ly = 5e-3
Lz = 5e-3

# Assume that the Phi angular function squared is equal to 1 for the point of 
# the rx antenna. This corresponds to the rx antenna being placed at the
# location where the Bragg condition gives maximum scattering

# Radar equation properties
G = 10**(20/10)                 # Gain (both T and R) of 20 dB
lr = con.c/f                    # Wavelength at receiver (this is in air)
R = 10e-2                       # Range (both T and R)

# Photoelastic constant - wild guess based on some glass in Saleh
p = 0.25


# Calculation of cross section
sigma = er**2*k**4/16/np.pi*p**2*s0**2*Lx**2*Ly**2*Lz**2

# Calculation of received power/transmitted power
Pfrac = G**2*lr**2*sigma/(4*np.pi)**3/R**4

# Display this in dB
print(10*np.log10(Pfrac))