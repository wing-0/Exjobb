# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:23:13 2018

@author: tfy13nwi

Takes input from COMSOL data (Ez,Hx,Hy) on a circular-ish boundary with and
without photoelastic interaction. Calculates the difference fields and from
those the power flow is calculated. This is done for the case of a homogenoeus
domain and one with a circular mechanically different defect in the center

"""

import numpy as np
import matplotlib.pyplot as plt
import readFields_PE as read

plt.close('all')

# Bulk modulus
Ks = np.arange(200, 410, 20)

# %% Load data and process it into correct format

# File directory
loc = ('..\\Simulation\\All parameters correct 2018-12-28\\Mech defect\\' +
       'Results')

# Import data without and with photoelasticity (with taper)
files_pe = [loc + '\\mech_a_40_d_le_centered_K1_' + str(a) + '_phon_1.csv'
            for a in Ks]
files_nope = [loc + '\\mech_a_40_d_le_centered_K1_' + str(a) + '_phon_0.csv'
              for a in Ks]
Ez, Hx, Hy, Sx, Sy, x, y = read.readEHS_PE(files_pe, files_nope)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Sort by angle (since COMSOL does not)
for i in range(0, len(Ks)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

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
plt.plot(Ks, Stot/Stot.max(), '.-')
plt.title('Total scattered power (normalized)')
plt.xlabel('K [MPa]')
plt.ylabel('P/P$_{max}$')

# %% Plotting using propagation angle and not observation angle

# Calculate propagation angles for all points and sort data by this angle
# instead of by the observation angle
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
for i in range(0, len(Ks)):
    ind = propang[:, i].argsort()
    propang[:, i] = propang[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]
    normS[:, i] = normS[ind, i]


# Plot the Poynting vector magnitude for all PROPAGATION angles
plt.figure()
plt.polar(propang, normS/normS.max())
plt.title('Poynting vector (time avg.), normalized magnitude')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('$\\left| \\left<\\mathbf{S}\\right> \\right|$ / ' +
           '$\\left| \\left<\\mathbf{S}\\right> \\right|_{max}$')
plt.legend(Ks, title='K [MPa]')

# "Clean up" the data by removing all points with a Poynting vector magnitude
# below 5 % of the mean
sig = normS > 0.05*normS.mean(axis=0)

plt.figure()

for i in range(0, len(Ks)):
    plt.plot(np.degrees(propang[sig[:, i], i]), normS[sig[:, i], i] /
             normS.max())

plt.title('Poynting vector (time avg.), normalized magnitude\nValues below' +
          ' 5 % of mean discarded')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('$\\left| \\left<\\mathbf{S}\\right> \\right|$ / ' +
           '$\\left| \\left<\\mathbf{S}\\right> \\right|_{max}$')
plt.legend(Ks, title='K [MPa]')


