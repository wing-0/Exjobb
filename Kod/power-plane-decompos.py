# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 11:50:02 2019

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# (+) or (-) scattering
pm = -1

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


# Cuboid phi, no z sinc
def PhiC(phi, alpha):
    return (np.sinc(Lx/2/np.pi*(k - k*np.cos(phi) + pm*q*np.cos(alpha))) *
            np.sinc(Ly/2/np.pi*(-k*np.sin(phi) + pm*q*np.sin(alpha))))


# Parallelogram phi, no z sinc
def PhiP(phi, alpha):
    return (np.sinc(da/2/np.pi/np.sin(alpha) *
                    (k - k*np.cos(phi) + pm*q*np.cos(alpha))) *
            np.sinc(de/2/np.pi/np.tan(alpha) *
                    (k - k*(np.cos(phi) + np.sin(phi)*np.tan(alpha)) +
                     pm*q*(np.cos(alpha) + np.sin(alpha)*np.tan(alpha)))))

# Relative permittivity
er = 1.29

# Photoelastic constant
ph = (er-1)*(er+2)/3/er**2

# Electric field amplitude at aperture center (estimated from COMSOL)
E00 = 44.17

# Pressure amplitude at aperture center (estimated from COMSOL)
p00 = 4483

# Bulk modulus
kbm = 400e6

# Speed of light
c0 = 299792458

# Permeability of free space
mu0 = 4*np.pi*1e-7

# Integrated Poynting vector magnitude, cuboid
powerC = (E00**2*de**2*er**(5/2)*k**4/(2**15*mu0*c0*r**2) *
          ph**2*p00**2*da**2/kbm**2 *
          Lx**2*Ly**2 * Ptot)

# Integrated Poynting vector magnitude, parallelogram
powerP = (E00**2*de**2*er**(5/2)*k**4/(2**15*mu0*c0*r**2) *
          ph**2*p00**2*da**2/kbm**2 *
          da**2*de**2/np.sin(alpha)**2 * Ptot_p)

plt.figure()
plt.grid()
plt.plot(angles, powerC, angles, powerP)
plt.xlabel('$\\alpha$')
plt.ylabel('P [W]')
plt.legend(['Cuboid', 'Parallelogram'])
plt.title('Only $k_x = k$, $q_{y\'} = q$ plane wave components used')