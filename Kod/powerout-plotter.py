# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:21:17 2018

@author: tfy13nwi

Takes input from COMSOL data (Ez,Hx,Hy) on a circular boundary both with and
without photoelastic interaction. Calculates the difference fields and from
those the power flow is calculated in the normal direction to the boundary.
The results are then plotted

"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# Files for data without and with photoelasticity
files = [
        (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation\Circular geometry'
         r'\data_verification_nope.csv'),
        (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation\Circular geometry'
         r'\data_verification_pe.csv')
         ]

# Angle between EM and Ac. wavevectors
alpha = 30

# Lists for all data
xt = []
yt = []
Ezt = []
Hxt = []
Hyt = []

# Read data from both files
for f in files:
    data = np.genfromtxt(f, delimiter=',', comments='%', dtype=str)
    xt.append(data[:, 0].astype(float))
    yt.append(data[:, 1].astype(float))
    Ezt.append(np.array([s.replace('i', 'j') for s in data[:, 2]])
               .astype(complex))
    Hxt.append(np.array([s.replace('i', 'j') for s in data[:, 3]])
               .astype(complex))
    Hyt.append(np.array([s.replace('i', 'j') for s in data[:, 4]])
               .astype(complex))

# Define coordinates according to both files
if((False not in (xt[0] == xt[1])) and (False not in (yt[0] == yt[1]))):
    x = xt[0]
    y = yt[0]

# Calculate difference fields
Ez = Ezt[1] - Ezt[0]
Hx = Hxt[1] - Hxt[0]
Hy = Hyt[1] - Hyt[0]

# Calculate Poynting vector (time avg.)
Sx = -0.5*(Ez*np.conj(Hy))
Sy = 0.5*(Ez*np.conj(Hx))

# Calculate angle from x and y values (angle is counterclockwise from x-axis)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

# Calculate normal vector at each angle
nx = np.cos(theta)
ny = np.sin(theta)

# Calculate power flow normal to the boundary
P = Sx.real*nx + Sy.real*ny

# Sort by angle (since COMSOL does not)
ind = np.argsort(theta)
theta = theta[ind]
P = P[ind]

# Polar plot
plt.figure()
plt.polar(theta, P)

# Arrows for incident wavevectors (red: Ac, black: EM)
yr = np.array([np.min(P), np.max(P)])
plt.polar(np.radians(180)*np.ones(2), yr, 'k:')
plt.polar([np.radians(180), np.radians(190)], 0.1*yr, 'k')
plt.polar([np.radians(180), np.radians(170)], 0.1*yr, 'k')
plt.polar(np.radians(180+alpha)*np.ones(2), yr, 'r:')
plt.polar([np.radians(180+alpha), np.radians(190+alpha)], 0.1*yr, 'r')
plt.polar([np.radians(180+alpha), np.radians(170+alpha)], 0.1*yr, 'r')

# Calculate total scattered power (line integral of power flow)
# This is not very useful on its own, but if multiple angles alpha are swept
# it can be used to verify that maximum scattering is obtained at the Bragg
# condition angle
r = 0.75
arc = r*theta
Ptot = np.trapz(P, arc)
