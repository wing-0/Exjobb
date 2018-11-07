# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 17:42:04 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

filename = r'C:\Users\tfy13nwi\Documents\Exjobb\Simulation\Circular geometry test\dataexport_test_reggrid.csv'

# 0 if data is spreadsheet, 1 if grid
isgrid = True

r = csv.reader(open(filename,'r'))

x = []
y = []
v = []
pos = 0

for l in r:
    # Grid data processing
    if(isgrid):
        if(l[0] == '% Grid'):
            pos = 1
            continue
        elif(l[0] == '% Data'):
            pos = 3
            continue
        elif(l[0][0] == '%'):
            continue
        if(pos == 1):
            for k in l:
                x.append(float(k))
            pos = 2
        elif(pos == 2):
            for k in l:
                y.append(float(k))
            pos = 0
        elif(pos == 3):
            vpart = []
            for k in l:
                vpart.append(complex(k.replace('i','j')))
            v.append(vpart)
    # Spreadsheet data processing
    else:
        if(l[0][0] == '%'):
            continue
        x.append(float(l[0]))
        y.append(float(l[1]))
        v.append(complex(l[2].replace('i','j')))

x = np.array(x)
y = np.array(y)
v = np.array(v)

# Grid plotting
if(isgrid):
    [xx,yy] = np.meshgrid(x,y)
    
    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(xx,yy,v.real,cmap='RdBu')
    
    mval = np.nanmax(np.abs(v.real))
    plt.clim(-mval,mval)
    plt.colorbar()
    
    # Boundaries for specific circular geometry, not general
    th = np.linspace(0,2*np.pi)
    plt.plot(0.1*np.cos(th),0.1*np.sin(th),'k',linewidth=0.5)
    plt.plot(0.15*np.cos(th),0.15*np.sin(th),'k',linewidth=0.5)
    
    plt.figure()
    plt.imshow(v.real,interpolation='bilinear',cmap='RdBu',origin='lower')
    
# Spreadsheet plotting
else:
    plt.figure()
    plt.tricontourf(x,y,v.real,50,cmap='RdBu')
    
    mval = np.max(np.abs(v.real))
    plt.clim(-mval,mval)
    plt.colorbar()
    
    
    
    