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
import scipy.interpolate as interp

plt.close('all')

# Read data and plot (+) or (-) scattering
pm = -1
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

# Files for data without and with photoelasticity (with taper)
files_pe = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
             ' (' + pmchar + ')_pe.csv') for a in angles]
files_nope = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
               ' (' + pmchar + ')_nope.csv') for a in angles]

# Files for data without and with photoelasticity (no taper)
files_not_pe = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
                 ' notaper (' + pmchar + ')_pe.csv') for a in angles]
files_not_nope = [(loc + '\\a=' + str(a) + '_opt=' + str(opta) +
                   ' notaper (' + pmchar + ')_nope.csv') for a in angles]

# Lists for all data (with taper)
xt = []
yt = []
Ezt = []
Hxt = []
Hyt = []

# Lists for all data (no taper)
xt_n = []
yt_n = []
Ezt_n = []
Hxt_n = []
Hyt_n = []

# Read data from both files
for f_pe, f_not_pe, f_nope, f_not_nope in zip(files_pe, files_not_pe,
                                              files_nope, files_not_nope):
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
    
    data_not_pe = np.genfromtxt(f_not_pe, delimiter=',', comments='%',
                                dtype=str)
    data_not_nope = np.genfromtxt(f_not_nope, delimiter=',', comments='%',
                                  dtype=str)
    xt_n.append(data_not_pe[:, 0].astype(float))
    yt_n.append(data_not_pe[:, 1].astype(float))
    temp = np.array([
            np.char.replace(data_not_pe[:, 2], 'i', 'j').astype('complex'),
            np.char.replace(data_not_nope[:, 2], 'i', 'j').astype('complex')
            ]).T
    Ezt_n.append(temp)
    temp = np.array([
            np.char.replace(data_not_pe[:, 3], 'i', 'j').astype('complex'),
            np.char.replace(data_not_nope[:, 3], 'i', 'j').astype('complex')
            ]).T
    Hxt_n.append(temp)
    temp = np.array([
            np.char.replace(data_not_pe[:, 4], 'i', 'j').astype('complex'),
            np.char.replace(data_not_nope[:, 4], 'i', 'j').astype('complex')
            ]).T
    Hyt_n.append(temp)

# Calculate difference fields for all angles (with taper)
Ez = np.array([Ezt[i][:, 0] - Ezt[i][:, 1] for i in range(0, len(Ezt))]).T
Hx = np.array([Hxt[i][:, 0] - Hxt[i][:, 1] for i in range(0, len(Hxt))]).T
Hy = np.array([Hyt[i][:, 0] - Hyt[i][:, 1] for i in range(0, len(Hyt))]).T

# Calculate difference fields for all angles (no taper)
Ez_n = np.array([Ezt_n[i][:, 0] - Ezt_n[i][:, 1] for i in range(0,
                 len(Ezt_n))]).T
Hx_n = np.array([Hxt_n[i][:, 0] - Hxt_n[i][:, 1] for i in range(0,
                 len(Hxt_n))]).T
Hy_n = np.array([Hyt_n[i][:, 0] - Hyt_n[i][:, 1] for i in range(0,
                 len(Hyt_n))]).T

# Convert x,y-data to correct array format (with taper)
x = np.array(xt).T
y = np.array(yt).T

# Convert x,y-data to correct array format (no taper)
x_n = np.array(xt_n).T
y_n = np.array(yt_n).T

# Calculate time avg. power flow (with taper)
Sx = -0.5*(Ez*np.conj(Hy)).real
Sy = 0.5*(Ez*np.conj(Hx)).real

# Calculate time avg. power flow (with taper)
Sx_n = -0.5*(Ez_n*np.conj(Hy_n)).real
Sy_n = 0.5*(Ez_n*np.conj(Hx_n)).real

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
# (with taper)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
# (no taper)
theta_n = np.mod(np.arctan2(y_n, x_n), 2*np.pi)

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

# %% Plot using norm of the power flow

# Norm of the power flow (with taper)
normS = np.sqrt(Sx**2 + Sy**2)

# Plot norm of the power flow for all angles (with taper)
plt.figure()
for i in range(0, len(alpha)):
    plt.polar(theta[:, i], normS[:, i])

# Reset plot color cycle
plt.gca().set_prop_cycle(None)

# Norm of the power flow (no taper)
normS_n = np.sqrt(Sx_n**2 + Sy_n**2)

# Plot norm of the power flow for all angles (no taper)
for i in range(0, len(alpha)):
    plt.polar(theta[:, i], normS_n[:, i], ':')

plt.title('Norm of the time avg. power flow')
plt.ylabel('S [W/m$^2$]')
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

# %% Other plotting

# Find significant power flows (those larger than the mean value) and plot the
# observation and propagation angles for them (with taper)
plt.figure()
for i in range(0, len(alpha)):
    sig = normS[:, i] > np.mean(normS[:, i])
    plt.plot(np.degrees(theta[sig, i]), np.degrees(np.mod(
            np.arctan2(Sy[sig, i], Sx[sig, i]), 2*np.pi)))

# Reset plot color cycle
plt.gca().set_prop_cycle(None)

# Find significant power flows (those larger than the mean value) and plot the
# observation and propagation angles for them (no taper)
for i in range(0, len(alpha)):
    sig = normS_n[:, i] > np.mean(normS_n[:, i])
    plt.plot(np.degrees(theta_n[sig, i]), np.degrees(np.mod(
            np.arctan2(Sy_n[sig, i], Sx_n[sig, i]), 2*np.pi)), '--')

# Reset plot color cycle
plt.gca().set_prop_cycle(None)

yrange = plt.ylim()

# Find observation angles for maximum power flow norms and plot dotted lines
# for those observation angles (with taper)
maxth = np.degrees(np.diag(theta[normS.argmax(axis=0)]))
for i in range(0, len(alpha)):
    plt.plot(maxth[i]*np.ones(2), yrange, ':')

# Reset plot color cycle
plt.gca().set_prop_cycle(None)

# Find observation angles for maximum power flow norms and plot dotted lines
# for those observation angles (no taper)
maxth_n = np.degrees(np.diag(theta_n[normS_n.argmax(axis=0)]))
for i in range(0, len(alpha)):
    plt.plot(maxth_n[i]*np.ones(2), yrange, '-.')

plt.title('Scattering propagation direction')
plt.xlabel('Observation angle $\\phi$ [$^\\circ$]')
plt.ylabel('Propagation angle $\\phi$ [$^\\circ$]')
plt.legend([str(a) + '$^\\circ$' for a in alpha], title='$\\alpha$')
plt.ylim(yrange)

# Plot the angle where scattering is maximum for each alpha (both taper cases)
plt.figure()
plt.plot(alpha, maxth, '.-', alpha, maxth_n, '.--')
plt.title('Angle for maximum scattering')
plt.xlabel('$\\alpha$ [$^\\circ$]')
plt.ylabel('Observation angle $\\phi$ [$^\\circ$]')





