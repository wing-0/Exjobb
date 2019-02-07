# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:23:13 2018

@author: tfy13nwi

Takes input from COMSOL data (Ez,Hx,Hy) on a circular-ish boundary with and
without photoelastic interaction. Calculates the difference fields and from
those the power flow is calculated. This is done for a sweep of the angle
between wave vectors for either (+) or (-) scattering (this can be selected)

"""

import numpy as np
import matplotlib.pyplot as plt
import readFields_PE as read

plt.close('all')

savefigs = 0

# Read data and plot (+) or (-) scattering
pm = -1
pmchar = '+'
if(pm == -1):
    pmchar = '-'

# Optimal angle for Bragg scattering
opta = 40

# Angles between wave vectors used (in filenames for both (+) and (-)) and the
# corresponding alpha values which depend on pm
angles = np.arange(35, 46)
alpha = (180*(pm + 1)/2 - pm*angles).astype('int')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data and process it into correct format
###############################################################################

# File directory
loc = ('..\\Simulation\\All parameters correct 2018-12-28\\Angle sweep\\' +
       'Results')

files_pe = [(loc + '\\anglesweep(' + pmchar + ')_' + str(a) +
             '_pe.csv') for a in angles]
files_nope = [(loc + '\\anglesweep(' + pmchar + ')_' + str(a) +
               '_nope.csv') for a in angles]
Ez, Hx, Hy, Sx, Sy, x, y = read.readEHS_PE_various(files_pe, files_nope)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Sort by angle (since COMSOL does not)
for i in range(0, len(alpha)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot using the Poynting vector normal component
###############################################################################

# Magnitude of the Poynting vector
normS = np.sqrt(Sx**2 + Sy**2)

# Normal component of the Poynting vector along boundary
perpS = Sx*np.cos(theta) + Sy*np.sin(theta)

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi
xdiff = np.diff(np.vstack((x, x[0, :])), axis=0)
ydiff = np.diff(np.vstack((y, y[0, :])), axis=0)
dist = np.linalg.norm(np.array([xdiff, ydiff]), axis=0)
arc = np.cumsum(dist, axis=0)

# Integrate normal component of Poynting vector for all angles and plot
# against the angle between wave vectors
Stot = np.trapz(perpS, x=arc, axis=0)

plt.figure()
plt.grid()
plt.plot(alpha, Stot/Stot.max(), '.-')
plt.title('Total scattered power (normalized)')
plt.xlabel('$\\alpha$')
plt.ylabel('$P_\\mathrm{sc}/P_\\mathrm{sc, max}$')

# Adds degree sign to x ticks
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.xticks(aloc, alab)

# Save as pgf
#if(savefigs):
#    plt.rcParams['axes.unicode_minus'] = False
#    plt.savefig('../Text/Report/fig/angle-sweep-power(' + pmchar + ').pgf')

# dBm version
plt.figure()
plt.grid()
plt.plot(alpha, 10*np.log10(Stot/1e-3), '.-')
plt.title('Total scattered power, ($' + pmchar + '$) case')
plt.xlabel('$\\alpha$')
plt.ylabel('$P_\\mathrm{sc}$ [dBm]')

# Adds degree sign to x ticks
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.xticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/angle-sweep-dBm(' + pmchar + ').pgf')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Comparison with analytical results
###############################################################################

# Wavelengths and wavenumbers
la = 0.0028
le = 2*la*np.cos(np.radians(opta))
k = 2*np.pi/le
q = 2*np.pi/la

# Geometry
r = 30*la
angles_an = np.linspace(35, 45, 101)
angles_an = 180*(pm + 1)/2 - pm*angles_an
alpha_an = np.radians(angles_an)
phi_an = np.linspace(0, 360, 721)
phi_an = np.radians(phi_an)
[alpha_m, phi_m] = np.meshgrid(alpha_an, phi_an)

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

# Interaction region dimensions (most simple case)
de = 8*le
da = 8*la
Lx = da
Ly = de

# Interaction region dimensions (refined)
#de = 10*le
#da = 10*la
#Lx = da
#Ly = de

# Cuboid phi, 2D
Phi_c = (np.sinc(Lx/2/np.pi*(k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(Ly/2/np.pi*(-k*np.sin(phi_m) + pm*q*np.sin(alpha_m))))


# Parallelogram phi, 2D
Phi_p = (np.sinc(da/2/np.pi/np.sin(alpha_m) *
                 (k - k*np.cos(phi_m) + pm*q*np.cos(alpha_m))) *
         np.sinc(de/2/np.pi/np.tan(alpha_m) *
                 (k - k*(np.cos(phi_m) + np.sin(phi_m)*np.tan(alpha_m)) +
                  pm*q*(np.cos(alpha_m) + np.sin(alpha_m)*np.tan(alpha_m)))))

# Normal component of Poynting vector, cuboid
# This is the same as the magnitude for a circular boundary
Sn_c = (0.5/wi * Ei0**2 * er**2*k**3*ph**2*p0**2/8/np.pi/r/kbm**2 *
        Lx**2*Ly**2 * Phi_c**2)

# Normal component of Poynting vector, parallelogram
# This is the same as the magnitude for a circular boundary
Sn_p = (0.5/wi * Ei0**2 * er**2*k**3*ph**2*p0**2/8/np.pi/r/kbm**2 *
        da**2*de**2/np.sin(alpha_m)**2 * Phi_p**2)

# Compensation factor for transforming peak amplitude to equivalent plane
# wave power. This is multiplied to the Poynting vectors
gamma = np.pi/32
Sn_c = gamma*Sn_c
Sn_p = gamma*Sn_p

# Total scattering power, cuboid
Ptot_c = np.trapz(Sn_c, x=r*phi_an, axis=0)

# Total scattering power, parallelogram
Ptot_p = np.trapz(Sn_p, x=r*phi_an, axis=0)

# Plot total power in dBm for sim, cuboid and parallelogram
plt.figure()
plt.grid()
plt.plot(alpha, 10*np.log10(Stot/1e-3), '.-')
plt.plot(angles_an, 10*np.log10(Ptot_c/1e-3), '--')
plt.plot(angles_an, 10*np.log10(Ptot_p/1e-3), ':')
plt.title('Total scattered power, ($' + pmchar + '$) case')
plt.xlabel('$\\alpha$')
plt.ylabel('$P_\\mathrm{sc}$ [dBm]')
plt.legend(['Simulated', 'Cuboid', 'Parallelogram'])

# Adds degree sign to x ticks
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.xticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/angle-sweep-power-an(' + pmchar + ').pgf')

# Plot maximum Poynting magnitude in W/m2 for sim, cuboid and parallelogram
plt.figure()
plt.grid()
plt.plot(alpha, normS.max(axis=0), '.-')
plt.plot(angles_an, Sn_c.max(axis=0), '--')
plt.plot(angles_an, Sn_p.max(axis=0), ':')
plt.title('Maximum of the Poynting vector, ($' + pmchar + '$) case')
plt.xlabel('$\\alpha$')
plt.ylabel('$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> \\right|$' +
           '$_\\mathrm{max}$ [W/m$^2$]')
plt.legend(['Simulated', 'Cuboid', 'Parallelogram'])

# Adds degree sign to x ticks
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.xticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/angle-sweep-peak-an(' + pmchar + ').pgf')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plotting using propagation angle and not observation angle
###############################################################################

# Calculate propagation angles for all points and sort data by this angle
# instead of by the observation angle
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
for i in range(0, len(alpha)):
    ind = propang[:, i].argsort()
    propang[:, i] = propang[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]
    normS[:, i] = normS[ind, i]

# Average of propagation angle weighted by Poynting vector magnitude
wavgang = np.average(propang, axis=0, weights=normS)

# Plot the weighted average of the propagation angle for each alpha
plt.figure()
plt.grid()
plt.plot(alpha, np.degrees(wavgang), '.-')
plt.title('Weighted avg. of propagation angle (weighted by ' +
          '$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> \\right|$), ($' +
          pmchar + '$) case')
plt.xlabel('$\\alpha$')
plt.ylabel('$\\overline{\\phi}_\\mathrm{prop}$')

# Adds degree sign to x and y ticks
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.xticks(aloc, alab)
aloc = plt.yticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in aloc]
plt.yticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/angle-sweep-angles(' + pmchar + ').pgf')

# Plot Poynting vector magnitude for all PROPAGATION angles
plt.figure()
plt.polar(propang, normS/normS.max())
plt.polar(np.radians((180 - pm*2*opta))*np.ones(2),
          np.arange(2), 'k:')
plt.ylim([0, 1])
plt.title('Normalized Poynting vector magnitude ' +
          '$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> \\right|$ / ' +
          '$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> ' +
          '\\right|_\\mathrm{max}$, ($' + pmchar + '$) case\n')
plt.xlabel('$\\phi_\\mathrm{prop}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$',
           bbox_to_anchor=(1.1, 0.5), loc='center left')
plt.gca().set_rlabel_position(70)

# Fixes degree sign not displaying correctly in LaTeX
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in np.degrees(aloc)]
plt.xticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/angle-sweep-polar(' + pmchar + ').pgf')

# "Clean up" the data by removing all points with a Poynting vector magnitude
# below 5 % of the mean
sig = normS > 0.05*normS.mean(axis=0)

plt.figure()

for i in range(0, len(alpha)):
    plt.plot(np.degrees(propang[sig[:, i], i]), normS[sig[:, i], i] /
             normS.max())

plt.title('Poynting vector (time avg.), normalized magnitude\nValues below' +
          ' 5 % of mean discarded')
plt.xlabel('$\\phi_\\mathrm{prop}$')
plt.ylabel('$\\left| \\left<\\mathbf{S}\\right> \\right|$ / ' +
           '$\\left| \\left<\\mathbf{S}\\right> \\right|_\\mathrm{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')
