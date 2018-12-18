# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 16:08:55 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# (+) or (-) scattering
pm = -1

opta = 40

# Wavelengths and wavenumbers
la = 0.03
le = 2*la*np.cos(np.radians(opta))
k = 2*np.pi/le
q = 2*np.pi/la

# Geometry
r = 30*la
angles = np.linspace(35, 45, 101)
alpha = np.radians((180*(pm + 1)/2 - pm*angles))
phi = np.linspace(0, 2*np.pi, 800)
[alpha_m, phi_m] = np.meshgrid(alpha, phi)
we = 8*le
wa = 8*la
Lx = wa/np.sin(alpha_m)
Ly = we

# Cuboid phi (with xy-dimensions), no z sinc
Phi = Lx*Ly*(np.sinc(Lx/2/np.pi*(k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
             np.sinc(Ly/2/np.pi*(-k*np.sin(phi_m) + pm*q*np.sin(alpha_m))))

# Parallelogram phi (with xy-dimensions), no z sinc
Phi_p = (wa*we/np.sin(alpha_m) *
         np.sinc(wa/2/np.pi/np.sin(alpha_m) *
                 (k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(we/2/np.pi/np.tan(alpha_m) *
                 (k - k*(np.cos(phi_m) + np.sin(phi_m)*np.tan(alpha_m)) +
                  pm*q*(np.cos(alpha_m) + np.sin(alpha_m)*np.tan(alpha_m)))))

P_ang = Phi**2
P_ang_p = Phi_p**2

plt.figure()
plt.polar(phi, P_ang[:, ::10])
plt.legend(angles[::10])

plt.figure()
plt.grid()
Ptot = np.trapz(P_ang.T, x=r*phi)
plt.plot(angles, Ptot/Ptot.max())

plt.figure()
plt.polar(phi, P_ang_p[:, ::10])
plt.legend(angles[::10])

plt.figure()
plt.grid()
Ptot_p = np.trapz(P_ang_p.T, x=r*phi)
plt.plot(angles, Ptot_p/Ptot_p.max())
