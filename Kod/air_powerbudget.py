# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 17:39:30 2018

@author: tfy13nwi
"""

import numpy as np
import scipy.constants as con
import matplotlib.pyplot as plt

# Material properties (air)
rho_0 = 1.293               # Density (273 K)
v = 344                     # Speed of sound (293 K)
er = 1.00059                # Relative permittivity

# Photoelastic constant from scalar model
ph = -1/3*(er-1)*(er+2)/er**2

# Acoustic transducer properties SensComp 40KPT25
fa = 40e3                   # Frequency
la = v/fa                   # Wavelength
angle_a = np.radians(23)    # Beam angle (full width, -6 dB SPL)
p_p = 20e-6*10**(110/20)    # Sound pressure amplitude at ref distance
s_p = p_p/rho_0/v**2        # Sound strain amplitude at ref distance
d_ref = 30e-2               # Reference distance for SPL measurment

# EM antenna properties
fe = 18e9                   # Frequency
le = con.c/fe               # Wavelength
k = 2*np.pi/le              # Wave number
G = 10**(20/10)             # Transmitter & Receiver gain
angle_e = 4/np.sqrt(G)      # Beam angle based on (16.2.20) in Orfanidis

# Angle between wave vectors for optimal scattering
alpha = np.arccos(le/2/la)

# Range for EM transmitter and receiver
R = d_ref/np.cos(alpha)

# Approximate beam widths using a triangle approximation
wa = 2*d_ref*np.tan(angle_a/2)
we = 2*R*np.tan(angle_e/2)
Lz = min([we, wa])

T0 = 290                        # Standard temperature, 290 K
NF = 5                          # Noise figure in dB
F = 10**(NF/10)                 # Noise factor, linear
B = fa                          # Assume receiver BW same as Ac. freq.

# RCS, assuming optimal receiver placement (aka Phi_p = 1/sin(alpha))
sigma = er**2*k**2/16/np.pi*ph**2*s_p**2*wa**2*we**2*Lz**2/np.sin(alpha)

P = (4*np.pi)**3 * R**4 * con.k*T0*B *F/(G**2 * le**2 * sigma)

print(P)
