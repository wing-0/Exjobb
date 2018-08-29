# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:25:53 2018

@author: tfy13nwi

Calculates oblique penetration length of an optical wave into an acoustic wave 
as described in Saleh2007 (p.812). Also calculates a "more realistic" value
given a beamwidth of both waves. The ratio between them is calculated for
various angles and beamwidth ratios.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lg

# Ratio of acoustic to electromagnetic beam width B/w
Bw = np.logspace(-1,2,4)

for bw in Bw:
    w = 1
    B = bw*w
    
    # Angle of incidence (avoiding 0 and pi/2)
    theta = np.linspace(np.pi/100,np.pi/2.01)
    
    # True penetration length
    s = B/np.cos(theta) + w*np.tan(theta)
    
    # Weird length L
    L = w/np.cos(theta) + B*np.tan(theta)
    
    # Approximate penetration length
    r = L/np.sin(theta)
    
    plt.semilogy(theta,r/s)
    plt.xlabel(r'$\theta$')
    plt.ylabel('True/approx. length')
    
    nm = lg.norm(r-s)
    print('B/w = {:<7.1f} norm = {:0.2f}'.format(bw,nm))
    
plt.legend(Bw,title='Beam width ratio')

# plt.savefig('acousto-optics-length.png',dpi=400)