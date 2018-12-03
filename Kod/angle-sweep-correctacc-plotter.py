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
import scipy.interpolate as interp

plt.close('all')

# Read data and plot (+) or (-) scattering
pm = 1
pmchar = '+'
if(pm == -1):
    pmchar = '-'

# Optimal angle for Bragg scattering
opta = 40

# Angles between wave vectors used (in filenames for both (+) and (-)) and the
# corresponding alpha values which depend on pm
angles = np.array([35, 38, 40, 42, 45])
alpha = (180*(pm + 1)/2 - pm*angles).astype('int')

# %% Load data and process it into correct format

# File directory
loc = (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation'
       r'\Angle sweep verification\Correct acceleration\Results')

# Import data without and with photoelasticity (with taper)
files_pe = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
             ' (' + pmchar + ')_pe.csv') for a in angles]
files_nope = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' (' + pmchar + ')_nope.csv') for a in angles]
Ez, Hx, Hy, Sx, Sy, x, y = read.readEHS_PE(files_pe, files_nope)

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

# %% Plot using norm of the power flow

# Norm of the power flow (with taper)
normS = np.sqrt(Sx**2 + Sy**2)

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi (with taper)
xdiff = np.diff(np.vstack((x, x[0, :])), axis=0)
ydiff = np.diff(np.vstack((y, y[0, :])), axis=0)
dist = np.linalg.norm(np.array([xdiff, ydiff]), axis=0)
arc = np.cumsum(dist, axis=0)

# Integrate norm of the power flow for all angles and plot against the angle
# between wave vectors (with taper)
Stot = np.trapz(normS, arc, axis=0)

# Interpolation (2nd order spline)
f = interp.interp1d(alpha, Stot/Stot.max(), kind='quadratic')
aa = np.linspace(alpha.min(), alpha.max(), 100)

plt.figure()
plt.plot(alpha, Stot/Stot.max(), '.-', aa, f(aa), ':')
plt.title('Total scattered power (normalized)')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('P/P$_{max}$')

# %% Plotting using propagation angle and not observation angle

# TODO: don't just use one maximum power - instead, take a number of datapoints
#       near the maximum and do some averaging
# Calculate propagation angles for all points and find angle for maximum
# power flow norms (with taper)
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
maxang = np.diag(propang[normS.argmax(axis=0)])

# Re-sort data by propagation angle (with taper)
for i in range(0, len(alpha)):
    ind = propang[:, i].argsort()
    propang[:, i] = propang[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]
    normS[:, i] = normS[ind, i]


# Plot the angle where scattering is maximum for each alpha
plt.figure()
plt.plot(alpha, np.degrees(maxang), '.-')
plt.title('Angle for maximum scattering')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('Propagation angle $\\phi$ [$^\\circ$]')

# Plot norm of the power flow for all PROPAGATION angles (with taper)
plt.figure()
plt.polar(propang, normS/normS.max())
plt.title('Norm of the time avg. power flow (normalized, tapered ports)')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')
