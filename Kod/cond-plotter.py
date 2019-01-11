# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:23:13 2018

@author: tfy13nwi

Takes input from COMSOL data (Ez,Hx,Hy) on a circular-ish boundary with and
without photoelastic interaction. Calculates the difference fields and from
those the power flow is calculated. This is done for the case of a homogenoeus
domain and one with a circular conductive defect in the center

"""

import numpy as np
import matplotlib.pyplot as plt
import readFields_PE as read

plt.close('all')

savefigs = 0

dBm = 1

# Conductivities
sigs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# %% Load data and process it into correct format

# File directory
loc = ('..\\Simulation\\All parameters correct 2018-12-28\\Cond defect\\' +
       'Results')

# Import data without and with photoelasticity (with taper)
files_pe = [loc + '\\cond_a_40_d_le_centered_sig_' + str(a) + '_phon_1.csv'
            for a in sigs]
files_nope = [loc + '\\cond_a_40_d_le_centered_sig_' + str(a) + '_phon_0.csv'
              for a in sigs]
Ez, Hx, Hy, Sx, Sy, x, y = read.readEHS_PE(files_pe, files_nope)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Sort by angle (since COMSOL does not)
for i in range(0, len(sigs)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

# Convert conductivities into an array
sigs = np.array(sigs)

# %% Plot using Poynting vector magnitude

# Magnitude of the Poynting vector
normS = np.sqrt(Sx**2 + Sy**2)

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi
xdiff = np.diff(np.vstack((x, x[0, :])), axis=0)
ydiff = np.diff(np.vstack((y, y[0, :])), axis=0)
dist = np.linalg.norm(np.array([xdiff, ydiff]), axis=0)
arc = np.cumsum(dist, axis=0)

# Integrate Poynting vector magnitude for all angles and plot against sigma
Stot = np.trapz(normS, x=arc, axis=0)

plt.figure()
plt.grid()
plt.plot(sigs, Stot/Stot.max(), '.-')
plt.title('Total scattered power (normalized)')
plt.xlabel('$\sigma$ [S/m]')
plt.ylabel('$P_\\mathrm{sc}/P_\\mathrm{sc, max}$')

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/fig/cond-power.pgf')

# dBm version
if(dBm):
    plt.figure()
    plt.grid()
    plt.plot(sigs, 10*np.log10(Stot/1e-3), '.-')
    plt.title('Total scattered power')
    plt.xlabel('$\sigma$ [S/m]')
    plt.ylabel('$P_\\mathrm{sc}$ [dBm]')

    # Save as pgf
    if(savefigs):
        plt.rcParams['axes.unicode_minus'] = False
        plt.savefig('../Text/Report/fig/cond-dBm.pgf')

# %% Plotting using propagation angle and not observation angle

# Calculate propagation angles for all points and sort data by this angle
# instead of by the observation angle
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
for i in range(0, len(sigs)):
    ind = propang[:, i].argsort()
    propang[:, i] = propang[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]
    normS[:, i] = normS[ind, i]


# Plot the Poynting vector magnitude for all PROPAGATION angles
plt.figure()
plt.polar(propang, normS/normS.max())
plt.polar(np.radians(260)*np.ones(2), np.arange(2), 'k:')
plt.ylim([0, 1])
plt.title('Poynting vector (time avg.), normalized magnitude ' +
          '$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> \\right|$ / ' +
          '$\\left| \\left<\\mathbf{S}_\\mathrm{sc}\\right> ' +
          '\\right|_\\mathrm{max}$\n')
plt.xlabel('$\\phi_\\mathrm{prop}$')
plt.legend(sigs, title='$\sigma$ [S/m]', bbox_to_anchor=(1.1, 0.5),
           loc='center left')
plt.gca().set_rlabel_position(70)

# Fixes degree sign not displaying correctly in LaTeX
aloc = plt.xticks()[0]
alab = [str(int(a)) + '$^\\circ$' for a in np.degrees(aloc)]
plt.xticks(aloc, alab)

# Save as pgf
if(savefigs):
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('../Text/Report/cond-polar.pgf')

# "Clean up" the data by removing all points with a Poynting vector magnitude
# below 5 % of the mean
sig = normS > 0.05*normS.mean(axis=0)

plt.figure()

for i in range(0, len(sigs)):
    plt.plot(np.degrees(propang[sig[:, i], i]), normS[sig[:, i], i] /
             normS.max())

plt.title('Poynting vector (time avg.), normalized magnitude\nValues below' +
          ' 5 % of mean discarded')
plt.xlabel('$\\phi_\\mathrm{prop}$')
plt.ylabel('$\\left| \\left<\\mathbf{S}\\right> \\right|$ / ' +
           '$\\left| \\left<\\mathbf{S}\\right> \\right|_\\mathrm{max}$')
plt.legend(sigs, title='$\sigma$ [S/m]')


