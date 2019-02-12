# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 17:39:30 2018

@author: tfy13nwi
"""

import numpy as np
import scipy.constants as con

# Material properties (air)
rho_0 = 1.293               # Density (273 K)
v = 331.1                   # Speed of sound (273 K)
er = 1.00059                # Relative permittivity

# Bulk modulus (fluid)
K = rho_0*v**2

# Photoelastic constant from scalar model
ph = 1/3*(er-1)*(er+2)/er**2

# Acoustic transducer properties SensComp 40KPT25
fa = 40e3                   # Frequency
la = v/fa                   # Wavelength
angle_a = np.radians(23)    # Beam angle (full width, -6 dB SPL)
p0 = 20e-6*10**(110/20)     # Sound pressure amplitude at ref distance
d_ref = 30e-2               # Reference distance for SPL measurment
apa = 25.1e-3               # Ac. aperture from datasheet

# EM antenna properties
fe = 24e9                   # Frequency
le = con.c/np.sqrt(er)/fe   # Wavelength
k = 2*np.pi/le              # Wave number
G = 10**(20/10)             # Transmitter & Receiver gain
angle_e = 4/np.sqrt(G)      # Beam angle based on (16.2.20) in Orfanidis

# EM aperture size based on effective aperture
ape = np.sqrt(G*le**2/4/np.pi)

# Angle between wave vectors for optimal scattering
alpha = np.arccos(le/2/la)

# Range for EM transmitter and receiver - based on EM tx, Ac. tx, EM rx placed
# on the edge of a circle with center in the intersection point
R = d_ref

# Approximate beam widths using both transducer/antenna aperture size and
# beam spread from aperture to intersection
da = apa + 2*d_ref*np.tan(angle_a/2)
de = ape + 2*R*np.tan(angle_e/2)
Lz = min([de, da])

# Overlap length based on parallelogram model with cosine law for diagonal
d = np.sqrt((de/np.sin(np.pi-alpha))**2 + (da/np.sin(np.pi-alpha))**2 -
            2*de*da/np.sin(np.pi-alpha)**2 * np.cos(np.pi-alpha))

# Receiver things
T0 = 290                        # Standard temperature, 290 K
NF = 5                          # Noise figure in dB
F = 10**(NF/10)                 # Noise factor, linear
B = fa                          # Assume receiver BW same as Ac. freq.

# RCS, assuming optimal receiver placement
sigma = er**2*k**4*ph**2*p0**2/16/np.pi/K**2 * da**2*de**2*Lz**2/np.sin(alpha)

# Power required to obtain 0 dB SNR with 1 sample
P = (4*np.pi)**3 * R**4 * con.k*T0*B*F/(G**2 * le**2 * sigma)

# Samples required to obtain 0 dB SNR with 15 dBm EM power
N = (4*np.pi)**3 * R**4 * con.k*T0*B*F/(1e-3*10**1.5 * G**2 * le**2 * sigma)

print(60*'-')
print('{0:>60s}'.format('ULTRASOUND AND MICROWAVES IN AIR'))
print(60*'-')
print('{0:>50s} {1:.2f} kHz'.format('Ac. frequency:', fa/1e3))
print('{0:>50s} {1:.2f} GHz \n'.format('Em frequency:', fe/1e9))

print('{0:>50s} {1:.2f} deg'.format('Angle between wave vectors:',
                                    np.degrees(alpha)))
print('{0:>50s} {1:.2f} m \n'.format('Maximum overlap length ' +
      '(parallelogram model):', d))

print('{0:>50s} {1:.2f} W, {2:.2f} dBm'.format('Power required for 0 dB SNR' +
      ', using 1 sample:', P, 10*np.log10(P/1e-3)))
print('{0:>50s} {1:.0f} '.format('Samples required for 0 dB SNR, using' +
      ' 15 dBm power:', np.ceil(N)))
print(60*'-')
