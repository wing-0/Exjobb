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
import csv

# Files for data without and with photoelasticity
files = [
        (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation\Circular geometry'
         ' test\powerout_noPE.csv'),
        (r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation\Circular geometry '
         'test\powerout_PE.csv')
         ]

# Lists for all data
xt = []
yt = []
Ezt = []
Hxt = []
Hyt = []

# Read data from both files
for f in files:
    r = csv.reader(open(f, 'r'))

    x = []
    y = []
    Ez = []
    Hx = []
    Hy = []

    for l in r:
        if(not l[0][0] == '%'):
            x.append(float(l[0]))
            y.append(float(l[1]))
            Ez.append(complex(l[2].replace('i', 'j')))
            Hx.append(complex(l[3].replace('i', 'j')))
            Hy.append(complex(l[4].replace('i', 'j')))
    x = np.array(x)
    y = np.array(y)
    Ez = np.array(Ez)
    Hx = np.array(Hx)
    Hy = np.array(Hy)

    xt.append(x)
    yt.append(y)
    Ezt.append(Ez)
    Hxt.append(Hx)
    Hyt.append(Hy)

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

# Plot
plt.plot(np.degrees(theta), P)
plt.xlim(0, 360)
plt.xticks(np.linspace(0, 360, 5))


# Line for theoretical scattering maximum (not a general calculation!)
alpha = 30
phi_m = np.degrees(np.arctan(-1/np.tan(np.radians(alpha))))
yr = np.array([np.min(P), np.max(P)])
plt.plot((phi_m + 360)*np.ones(2), yr, ':')

# Polar plot
plt.figure()
plt.polar(theta, P)

# Arrows for incident wavevectors (red: Ac, black: EM)
plt.polar(np.radians(270 - alpha)*np.ones(2), yr, 'k:')
plt.polar([np.radians(270 - alpha), np.radians(280 - alpha)], 0.1*yr, 'k')
plt.polar([np.radians(270 - alpha), np.radians(260 - alpha)], 0.1*yr, 'k')
plt.polar(np.radians(270)*np.ones(2), yr, 'r:')
plt.polar([np.radians(270), np.radians(280)], 0.1*yr, 'r')
plt.polar([np.radians(270), np.radians(260)], 0.1*yr, 'r')







