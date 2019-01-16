# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 13:48:04 2019

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# (+) or (-) scattering
pm = 1

opta = 40

# Wavelengths and wavenumbers
la = 0.0028
le = 2*la*np.cos(np.radians(opta))
k = 2*np.pi/le
q = 2*np.pi/la

# Geometry
r = 30*la
angles = np.linspace(35, 45, 101)
alpha = np.radians((180*(pm + 1)/2 - pm*angles))
phi = np.linspace(0, 360, 721)
phi = np.radians(phi)
[alpha_m, phi_m] = np.meshgrid(alpha, phi)
de = 8*le
da = 8*la

Lx = da
Ly = de

# Relative permittivity
er = 1.29

# Photoelastic constant
ph = (er-1)*(er+2)/3/er**2

# Electric field amplitude at aperture center (estimated from COMSOL)
Ei0 = 44.17

# Pressure amplitude at aperture center (estimated from COMSOL)
p0 = 4483

# Bulk modulus
kbm = 400e6

# Speed of light
c0 = 299792458

# Permeability of free space
mu0 = 4*np.pi*1e-7

# Wave impedance in material (eta_0*eta)
wi = mu0*c0/np.sqrt(er)

# Cuboid phi, 2D
Phi_c = (np.sinc(Lx/2/np.pi*(k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(Ly/2/np.pi*(-k*np.sin(phi_m) + pm*q*np.sin(alpha_m))))


# Parallelogram phi, 2D
Phi_p = (np.sinc(da/2/np.pi/np.sin(alpha_m) *
                 (k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(de/2/np.pi/np.tan(alpha_m) *
                 (k - k*(np.cos(phi_m) + np.sin(phi_m)*np.tan(alpha_m)) +
                  pm*q*(np.cos(alpha_m) + np.sin(alpha_m)*np.tan(alpha_m)))))


# Poynting vector magnitude, cuboid
Sm_c = (0.5/wi * Ei0**2 * er**2*k**3*ph**2*p0**2/8/np.pi/r/kbm**2 *
        Lx**2*Ly**2 * Phi_c**2)

# Poynting vector magnitude, parallelogram
Sm_p = (0.5/wi * Ei0**2 * er**2*k**3*ph**2*p0**2/8/np.pi/r/kbm**2 *
        da**2*de**2/np.sin(alpha_m)**2 * Phi_p**2)

# Compensation factor for transforming peak amplitude to equivalent plane
# wave power. This is multiplied to the Poynting vectors
gamma = np.pi/256
Sm_c = gamma*Sm_c
Sm_p = gamma*Sm_p

# Total scattering power and max scattering power, cuboid
Ptot_c = np.trapz(Sm_c, x=r*phi, axis=0)
Pmax_c = Sm_c.max(axis=0)

# Total scattering power and max scattering power, parallelogram
Ptot_p = np.trapz(Sm_p, x=r*phi, axis=0)
Pmax_p = Sm_p.max(axis=0)

# Plot angle sweep total power in dBm for cuboid and parallelogram
plt.figure()
plt.grid()
plt.plot(angles, 10*np.log10(Ptot_c/1e-3))
plt.plot(angles, 10*np.log10(Ptot_p/1e-3))
plt.xlabel('$\\alpha$')
plt.ylabel('Ptot [dBm]')
plt.legend(['Cuboid', 'Parallelogram'])

# Plot angle sweep max Poynting in W/m2 for cuboid and parallelogram
plt.figure()
plt.grid()
plt.plot(angles, Pmax_c)
plt.plot(angles, Pmax_p)
plt.xlabel('$\\alpha$')
plt.ylabel('Pmax [W/m2]')
plt.legend(['Cuboid', 'Parallelogram'])




