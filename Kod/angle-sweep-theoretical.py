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

# Aperture taper adjustment (use FWHM as beam diameters)
#de = de/2*np.sqrt(np.log(2))
#da = da/2*np.sqrt(np.log(2))

Lx = da#/np.sin(alpha_m)
Ly = de

# Cuboid phi, no z sinc
Phi = (np.sinc(Lx/2/np.pi*(k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
       np.sinc(Ly/2/np.pi*(-k*np.sin(phi_m) + pm*q*np.sin(alpha_m))))

# Cuboid phi with interaction area included
PhiL = Lx*Ly*Phi

# Parallelogram phi, no z sinc
Phi_p = (np.sinc(da/2/np.pi/np.sin(alpha_m) *
                 (k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(de/2/np.pi/np.tan(alpha_m) *
                 (k - k*(np.cos(phi_m) + np.sin(phi_m)*np.tan(alpha_m)) +
                  pm*q*(np.cos(alpha_m) + np.sin(alpha_m)*np.tan(alpha_m)))))

# Parallelogram phi with interaction area included
PhiL_p = (da*de/np.sin(alpha_m))*Phi_p

P_ang = Phi**2
P_ang_p = Phi_p**2
PL_ang = PhiL**2
PL_ang_p = PhiL_p**2

###############################################################################
# Stuff for comparing which alpha give larger Phi. Assymmetric due to Phi2
###############################################################################
Phi1 = np.sinc(Lx/2/np.pi*(k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m)))
Phi2 = np.sinc(Ly/2/np.pi*(-k*np.sin(phi_m) + pm*q*np.sin(alpha_m)))
plt.figure()
plt.grid()
plt.plot(np.degrees(phi), Phi1[:, ::10])
plt.gca().set_prop_cycle(None)
plt.plot(np.degrees(phi), Phi2[:, ::10], ':')
plt.gca().set_prop_cycle(None)
plt.plot(np.degrees(phi), Phi[:, ::10], '--')
plt.legend(angles[::10])

###############################################################################
# Polar plot of Phi^2
###############################################################################
plt.figure()
plt.plot(phi, P_ang[:, ::10])
plt.plot(np.radians(260)*np.ones(2), np.arange(2), 'k:')
plt.ylim([0, 1])
plt.legend(angles[::10])

###############################################################################
# "Total power" comparisons, cuboid
###############################################################################
plt.figure()
plt.grid()

# Use this for only integrating over a small angle range around 260 deg
#a = 5
#Ptot = np.trapz(P_ang.T[:, 520-a:521+a], x=r*phi[520-a:521+a])
#PLtot = np.trapz(PL_ang.T[:, 520-a:521+a], x=r*phi[520-a:521+a])

# Without interaction area
Ptot = np.trapz(P_ang.T, x=r*phi)
Pmax = P_ang.max(axis=0)
plt.plot(angles, Ptot/Ptot.max())
plt.plot(angles, Pmax/Pmax.max())

# With interaction area
PLtot = np.trapz(PL_ang.T, x=r*phi)
PLmax = PL_ang.max(axis=0)
plt.gca().set_prop_cycle(None)
plt.plot(angles, PLtot/PLtot.max(), '--')
plt.plot(angles, PLmax/PLmax.max(), '--')

plt.legend(['$\\Phi^2$ integrated over $\\phi$', 'Max value of $\\Phi^2$',
            '$A_i\\Phi^2$ integrated over $\\phi$',
            'Max value of $A_i\\Phi^2$'])
plt.xlabel('$\\alpha$')
plt.ylabel('Normalized value')

###############################################################################
# Polar plot of Phi_p^2
###############################################################################
plt.figure()
plt.plot(phi, P_ang_p[:, ::10])
plt.plot(np.radians(260)*np.ones(2), np.arange(2), 'k:')
plt.ylim([0, 1])
plt.legend(angles[::10])

###############################################################################
# "Total power" comparisons, parallelogram
###############################################################################
plt.figure()
plt.grid()

# Without interaction area
Ptot_p = np.trapz(P_ang_p.T, x=r*phi)
Pmax_p = P_ang_p.max(axis=0)
plt.plot(angles, Ptot_p/Ptot_p.max())
plt.plot(angles, Pmax_p/Pmax_p.max())

# With interaction area
PLtot_p = np.trapz(PL_ang_p.T, x=r*phi)
PLmax_p = PL_ang_p.max(axis=0)
plt.gca().set_prop_cycle(None)
plt.plot(angles, PLtot_p/PLtot_p.max(), '--')
plt.plot(angles, PLmax_p/PLmax_p.max(), '--')

plt.legend(['$\\Phi_p^2$ integrated over $\\phi$', 'Max value of $\\Phi_p^2$',
            '$A_i\\Phi_p^2$ integrated over $\\phi$',
            'Max value of $A_i\\Phi_p^2$'])
plt.xlabel('$\\alpha$')
plt.ylabel('Normalized value')

###############################################################################
# Power estimates using equivalent plane wave power
###############################################################################

# 15 dBm input power (assumed to be equal to radiated power)
Pt = 10**(15/10)*1e-3

# Beamwidth of EM beam (estimated from COMSOL)
theta = np.radians(20)

# Antenna gain, same as directivity here
Gt = 16/theta**2

# Relative permittivity
er = 1.29

# Photoelastic constant
ph = (er-1)*(er+2)/3/er**2

# Pressure amplitude (estimated from COMSOL)
p_p = 2500

# Strain amplitude (pressure amplitude divided by bulk modulus)
s0 = p_p/400e6

# Correction factor for equivalent plane wave power
gam = np.pi**2/256

# Integrated Poynting vector magnitude, cuboid
powerC = (Pt*Gt/4/np.pi/r**2 * er**2*k**4/64/np.pi**2/r**2 *
          ph**2*s0**2*Lx**2*Ly**2 * Ptot * gam)

# Integrated Poynting vector magnitude, parallelogram
powerP = (Pt*Gt/4/np.pi/r**2 * er**2*k**4/64/np.pi**2/r**2 *
          ph**2*s0**2*da**2*de**2/np.sin(alpha)**2 * Ptot_p * gam)

plt.figure()
plt.grid()
plt.plot(angles, powerC, angles, powerP)
plt.xlabel('$\\alpha$')
plt.ylabel('P [W]')
plt.legend(['Cuboid', 'Parallelogram'])

###############################################################################
# Power estimates using decomposition of beams into plane waves
###############################################################################

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
