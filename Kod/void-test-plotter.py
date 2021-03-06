# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:23:13 2018

@author: tfy13nwi

Takes input from COMSOL data (Ez,Hx,Hy) on a circular-ish boundary with and
without photoelastic interaction. Calculates the difference fields and from
those the power flow is calculated. This is done for the case of a homogenoeus
domain and one with a circular void in the center

"""

import numpy as np
import matplotlib.pyplot as plt
import readFields_PE as read

plt.close('all')

# Material used
mat = 'glass fiber'

# Optimal angle for Bragg scattering
opta = 40

# List of void sizes
voids = ['0', 'le']

# %% Load data and process it into correct format

# File directory
loc = ('..\\Simulation\\Void test\\Results\\')

# Import data without and with photoelasticity (with taper)
files_pe = [(loc + mat + ' a=' + str(opta) + ' void=' + a + '_pe.csv')
            for a in voids]
files_nope = [(loc + mat + ' a=' + str(opta) + ' void=' + a + '_nope.csv')
              for a in voids]
Ez, Hx, Hy, Sx, Sy, x, y = read.readEHS_PE(files_pe, files_nope)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Sort by angle (since COMSOL does not)
for i in range(0, len(voids)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

# %% Plot using norm of the power flow

# Norm of the power flow
normS = np.sqrt(Sx**2 + Sy**2)

# Calculate distances between all points (x,y) and construct an axis based on
# arc length from theta = 0 to theta = 2*pi
xdiff = np.diff(np.vstack((x, x[0, :])), axis=0)
ydiff = np.diff(np.vstack((y, y[0, :])), axis=0)
dist = np.linalg.norm(np.array([xdiff, ydiff]), axis=0)
arc = np.cumsum(dist, axis=0)

# Integrate norm of the power flow for all angles and plot against void size
Stot = np.trapz(normS, arc, axis=0)

plt.figure()
plt.plot(voids, Stot/Stot.max(), '.-')
plt.title('Total scattered power (normalized)')
plt.xlabel('Size of void')
plt.ylabel('P/P$_{max}$')

# %% Plotting using propagation angle and not observation angle

# Calculate propagation angles for all points and sort data by this angle
# instead of by the observation angle
propang = np.mod(np.arctan2(Sy, Sx), 2*np.pi)
for i in range(0, len(voids)):
    ind = propang[:, i].argsort()
    propang[:, i] = propang[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]
    normS[:, i] = normS[ind, i]

# Calculate propagation angle giving maximum scattering for all void. This is
# done by taking the average of the propagation angles where the norm of the
# power flow is larger than 95 % of the maximum value
maxang = np.zeros(len(voids))
for i in range(0, len(voids)):
    large = normS[:, i] > 0.95*normS[:, i].max()
    maxang[i] = np.mean(propang[large, i])

# Plot the angle where scattering is maximum for each void size
plt.figure()
plt.plot(voids, np.degrees(maxang), '.-')
plt.title('Mean angle of the 5 % with highest energy')
plt.xlabel('Size of void')
plt.ylabel('Propagation angle $\\phi$ [$^\\circ$]')

# Plot norm of the power flow for all PROPAGATION angles
plt.figure()
plt.polar(propang, normS/normS.max())
plt.title('Norm of the time avg. power flow (normalized)')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend(voids, title='void')

# "Clean up" the data by removing all points with a norm of the power flow
# below 5 % of the mean
sig = normS > 0.05*normS.mean(axis=0)

plt.figure()

for i in range(0, len(voids)):
    plt.plot(np.degrees(propang[sig[:, i], i]), normS[sig[:, i], i])

plt.title('Norm of the time avg. power flow (normalized)\nValues below 5 %' +
          ' of mean discarded')
plt.xlabel('Propagation angle $\\phi$')
plt.ylabel('S/S$_{max}$')
plt.legend(voids, title='void')


