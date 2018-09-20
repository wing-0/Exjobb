# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:21:46 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

# Phi angular function
def pphi(pm,theta,phi,alpha):
    L = 1e-2
    k = 2*np.pi/5e-3
    q = 2*np.pi/3e-3
    A = np.sinc(L/2/np.pi*(k-k*np.sin(theta)*np.cos(phi) \
                           + np.sign(pm)*q*np.cos(np.radians(alpha))))
    B = np.sinc(L/2/np.pi*(-k*np.sin(theta)*np.sin(phi) \
                           + np.sign(pm)*q*np.sin(np.radians(alpha))))
    C = np.sinc(-L/2/np.pi*k*np.cos(theta))
    return A*B*C

theta = np.pi/2
phi = np.linspace(0,2*np.pi,500)

# %% Angle sweep
alpha = np.linspace(0,180,100)
mainlobep = np.zeros(100)
mainlobem = np.zeros(100)

aind = 0
for a in alpha:
    p = pphi(1,theta,phi,a)**2
    m = pphi(-1,theta,phi,a)**2
    maxindexp, = np.where(p==np.max(p))
    maxindexm, = np.where(m==np.max(m))
    mainlobep[aind] = phi[maxindexp[0]]
    mainlobem[aind] = phi[maxindexm[0]]
    if(len(maxindexp) > 1 or len(maxindexm) > 1):
        print("Multiple")
    aind = aind + 1

plt.plot(alpha,mainlobep)
plt.plot(alpha,mainlobem)

# %% Some angles
alpha = np.array([15,30,60,180])

plt.polar(np.pi*np.ones(2),np.arange(2),'k')

c = {alpha[0]:'r',alpha[1]:'g',alpha[2]:'b',alpha[3]:'y'}

for a in alpha:
    p = pphi(-1,theta,phi,a)**2
    plt.polar(phi,p,color=c[a])

plt.legend(['EM']+list(alpha))

for a in alpha:
    plt.polar((np.pi+np.radians(a))*np.ones(2),np.arange(2),color=c[a])