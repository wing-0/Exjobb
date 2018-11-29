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
from readFields_PE import readEHS_PE

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
angles = np.array([35, 40, 42, 45])
alpha = (180*(pm + 1)/2 - pm*angles).astype('int')

# %% Load data and process it into correct format

# File directory
loc = (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation'
       r'\Angle sweep verification\Results')

# Import data without and with photoelasticity (with taper)
files_pe = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
             ' (' + pmchar + ')_pe.csv') for a in angles]
files_nope = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' (' + pmchar + ')_nope.csv') for a in angles]
Ez, Hx, Hy, Sx, Sy, x, y = readEHS_PE(files_pe, files_nope)

# Import data without and with photoelasticity (no taper)
files_pe_n = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' notaper (' + pmchar + ')_pe.csv') for a in angles]
files_nope_n = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
                 ' notaper (' + pmchar + ')_nope.csv') for a in angles]
Ez_n, Hx_n, Hy_n, Sx_n, Sy_n, x_n, y_n = readEHS_PE(files_pe_n, files_nope_n)

# Import data without and with photoelasticity (wider apertures)
files_pe_w = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' notaperwide (' + pmchar + ')_pe.csv') for a in angles]
files_nope_w = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
                 ' notaperwide (' + pmchar + ')_nope.csv') for a in angles]
Ez_w, Hx_w, Hy_w, Sx_w, Sy_w, x_w, y_w = readEHS_PE(files_pe_w, files_nope_w)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
# (with taper)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
# (no taper)
theta_n = np.mod(np.arctan2(y_n, x_n), 2*np.pi)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
# (wider apertures)
theta_w = np.mod(np.arctan2(y_w, x_w), 2*np.pi)

# Sort by angle (since COMSOL does not) (with taper)
for i in range(0, len(alpha)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

# Sort by angle (since COMSOL does not) (no taper)
for i in range(0, len(alpha)):
    ind = np.argsort(theta_n[:, i])
    theta_n[:, i] = theta_n[ind, i]
    x_n[:, i] = x_n[ind, i]
    y_n[:, i] = y_n[ind, i]
    Sx_n[:, i] = Sx_n[ind, i]
    Sy_n[:, i] = Sy_n[ind, i]

# Sort by angle (since COMSOL does not) (wider apertures)
for i in range(0, len(alpha)):
    ind = np.argsort(theta_w[:, i])
    theta_w[:, i] = theta_w[ind, i]
    x_w[:, i] = x_w[ind, i]
    y_w[:, i] = y_w[ind, i]
    Sx_w[:, i] = Sx_w[ind, i]
    Sy_w[:, i] = Sy_w[ind, i]

# %% Plot using norm of the power flow

# Norm of the power flow (with taper)
normS = np.sqrt(Sx**2 + Sy**2)

# Plot norm of the power flow for all angles (with taper)
plt.figure()
plt.polar(theta, normS/normS.max())

plt.title('Norm of the time avg. power flow (normalized, tapered ports)')
plt.xlabel('Observation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')

# Norm of the power flow (no taper)
normS_n = np.sqrt(Sx_n**2 + Sy_n**2)

# Plot norm of the power flow for all angles (no taper)
plt.figure()
plt.polar(theta_n, normS_n/normS_n.max())

plt.title('Norm of the time avg. power flow (normalized, non-tapered ports)')
plt.xlabel('Observation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi (with taper)
xdiff = np.diff(np.vstack((x, x[0, :])), axis=0)
ydiff = np.diff(np.vstack((y, y[0, :])), axis=0)
dist = np.linalg.norm(np.array([xdiff, ydiff]), axis=0)
arc = np.cumsum(dist, axis=0)

# Integrate norm of the power flow for all angles and plot against the angle
# between wave vectors (with taper)
Stot = np.trapz(normS, arc, axis=0)

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi (no taper)
xdiff_n = np.diff(np.vstack((x_n, x_n[0, :])), axis=0)
ydiff_n = np.diff(np.vstack((y_n, y_n[0, :])), axis=0)
dist_n = np.linalg.norm(np.array([xdiff_n, ydiff_n]), axis=0)
arc_n = np.cumsum(dist_n, axis=0)

# Integrate norm of the power flow for all angles and plot against the angle
# between wave vectors (no taper)
Stot_n = np.trapz(normS_n, arc_n, axis=0)

plt.figure()
plt.plot(alpha, Stot/Stot.max(), '.-', alpha, Stot_n/Stot_n.max(), '.--')
plt.title('Total scattered power (normalized)')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('P/P$_{max}$')
plt.legend(['Tapered ports', 'Non-tapered ports'])

# %% Other plotting

# TODO: don't just use one maximum power - instead, take a number of datapoints
#       near the maximum and do some averaging
# Calculate propagation angles for all points and find angle for maximum
# power flow norms (with taper)
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
maxang = np.diag(propang[normS.argmax(axis=0)])

# Calculate propagation angles for all points and find angle for maximum
# power flow norms (no taper)
propang_n = np.mod(np.arctan2(Sy_n, Sx_n), 2*np.pi)
maxang_n = np.diag(propang_n[normS_n.argmax(axis=0)])

# Plot the angle where scattering is maximum for each alpha (both taper cases)
plt.figure()
plt.plot(np.degrees(alpha), maxang, '.-', np.degrees(alpha), maxang_n, '.--')
plt.title('Angle for maximum scattering')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('Propagation angle $\\phi$ [$^\\circ$]')
plt.legend(['Tapered ports', 'Non-tapered ports'])

plt.figure()
plt.polar(propang, normS/normS.max())

plt.title('Norm of the time avg. power flow (normalized, tapered ports)')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')

plt.figure()
plt.polar(propang_n, normS_n/normS_n.max())

plt.title('Norm of the time avg. power flow (normalized, non-tapered ports)')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')


