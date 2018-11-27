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

plt.close('all')

# Read data and plot (+) or (-) scattering
pm = 1
pmchar = '+'
if(pm == -1):
    pmchar = '-'

# Optimal angle for Bragg scattering
opta = 40

# File directory
loc = (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation'
       r'\Angle sweep verification\Results')

# Angles between wave vectors used (in filenames for both (+) and (-)) and the
# corresponding alpha values which depend on pm
angles = np.array([25, 30, 35, 40, 42, 45, 50, 55])
alpha = (180*(pm + 1)/2 - pm*angles).astype('int')

# Files for data without and with photoelasticity
files_pe = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
             ' (' + pmchar + ')_pe.csv') for a in angles]
files_nope = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' (' + pmchar + ')_nope.csv') for a in angles]

# Lists for all data
xt = []
yt = []
Ezt = []
Hxt = []
Hyt = []

# Read data from both files
for f_pe, f_nope in zip(files_pe, files_nope):
    data_pe = np.genfromtxt(f_pe, delimiter=',', comments='%', dtype=str)
    data_nope = np.genfromtxt(f_nope, delimiter=',', comments='%', dtype=str)
    xt.append(data_pe[:, 0].astype(float))
    yt.append(data_pe[:, 1].astype(float))
    temp = np.array([
            np.char.replace(data_pe[:, 2], 'i', 'j').astype('complex'),
            np.char.replace(data_nope[:, 2], 'i', 'j').astype('complex')
            ]).T
    Ezt.append(temp)
    temp = np.array([
            np.char.replace(data_pe[:, 3], 'i', 'j').astype('complex'),
            np.char.replace(data_nope[:, 3], 'i', 'j').astype('complex')
            ]).T
    Hxt.append(temp)
    temp = np.array([
            np.char.replace(data_pe[:, 4], 'i', 'j').astype('complex'),
            np.char.replace(data_nope[:, 4], 'i', 'j').astype('complex')
            ]).T
    Hyt.append(temp)

# Calculate difference fields for all angles
Ez = np.array([Ezt[i][:, 0] - Ezt[i][:, 1] for i in range(0, len(Ezt))]).T
Hx = np.array([Hxt[i][:, 0] - Hxt[i][:, 1] for i in range(0, len(Hxt))]).T
Hy = np.array([Hyt[i][:, 0] - Hyt[i][:, 1] for i in range(0, len(Hyt))]).T

# Convert x,y-data to correct array format
x = np.array(xt).T
y = np.array(yt).T

# Calculate time avg. power flow
Sx = -0.5*(Ez*np.conj(Hy)).real
Sy = 0.5*(Ez*np.conj(Hx)).real

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Calculate normal vector at each angle
nx = np.cos(theta)
ny = np.sin(theta)

# Sort by angle (since COMSOL does not)
for i in range(0, len(alpha)):
    ind = np.argsort(theta[:, i])
    theta[:, i] = theta[ind, i]
    x[:, i] = x[ind, i]
    y[:, i] = y[ind, i]
    Sx[:, i] = Sx[ind, i]
    Sy[:, i] = Sy[ind, i]

# Norm of the power flow
normS = np.sqrt(Sx**2 + Sy**2)

# Plot norm of the power flow for all angles
for i in range(0, len(alpha)):
    plt.polar(theta[:, i], normS[:, i])

plt.title('Norm of the time avg. power flow')
plt.ylabel('S [W/m$^2$]')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')

# Integrate norm of the power flow for all angles and plot against the angle
# between wave vectors (boundary approximated to be perfectly circular)
r = 0.9
arc = r*theta
Stot = np.trapz(normS, arc, axis=0)

plt.figure()
plt.plot(alpha, Stot, '.-')
plt.title('Total scattered power')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('P [W/m]')

plt.figure()

# Find significant power flows (those larger than the mean value) and plot the
# observation and propagation angles for them
for i in range(0, len(alpha)):
    sig = normS[:, i] > np.mean(normS[:, i])
    plt.plot(np.degrees(theta[sig, i]), np.degrees(np.mod(
            np.arctan2(Sy[sig, i], Sx[sig, i]), 2*np.pi)))

# Reset plot color cycle
plt.gca().set_prop_cycle(None)

yrange = plt.ylim()

# Find observation angles for maximum power flow norms and plot dashed lines
# for those observation angles
maxth = np.degrees(np.diag(theta[normS.argmax(axis=0)]))
for i in range(0, len(alpha)):
    plt.plot(maxth[i]*np.ones(2), yrange, ':')

plt.title('Scattering propagation direction')
plt.xlabel('Observation angle $\\phi$ [$^\\circ$]')
plt.ylabel('Propagation angle $\\phi$ [$^\\circ$]')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')
plt.ylim(yrange)

plt.figure()
plt.plot(alpha, maxth, '.-')
plt.title('Angle for maximum scattering')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('Observation angle $\\phi$ [$^\\circ$]')





