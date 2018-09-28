# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:21:46 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
lam = 5e-3
Lam = 2.5e-3
k = 2*np.pi/lam
q = 2*np.pi/Lam

# Phi angular function
def pphi(pm,theta,phi,alpha):
    L = 1e-2
    A = np.sinc(L/2/np.pi*(k-k*np.sin(theta)*np.cos(phi) \
                           + np.sign(pm)*q*np.cos(np.radians(alpha))))
    B = np.sinc(L/2/np.pi*(-k*np.sin(theta)*np.sin(phi) \
                           + np.sign(pm)*q*np.sin(np.radians(alpha))))
    C = np.sinc(-L/2/np.pi*k*np.cos(theta))
    return A*B*C

theta = np.pi/2
phi = np.linspace(0,2*np.pi,500)

# %% Various alpha, plot function value for all phi,alpha combinations
plt.close('all')
plt.figure()
alpha = np.linspace(0,180,500)

[phi_m,alpha_m] = np.meshgrid(phi,alpha)

p = pphi(1,theta,phi_m,alpha_m)**2
m = pphi(-1,theta,phi_m,alpha_m)**2

plt.pcolormesh(alpha_m,np.degrees(phi_m),p)
plt.xlabel(r'$\alpha$ [$^\circ$]')
plt.ylabel(r'$\varphi$ [$^\circ$]')
plt.colorbar(label=r'$\Phi^{+^2}$')
plt.yticks(np.linspace(0,360,10))
plt.clim(0,1)
plt.suptitle('$\\Phi$ angular function $\\varphi$ and $\\alpha$ dependence ' + \
             '($\\theta = \\pi/2$)')
plt.title('EM prop.dir. at $\\varphi = 180^\\circ$, acoustic prop.dir. at ' + \
          '$\\varphi = 180^\\circ + \\alpha$',fontsize=10)

#plt.figure()
#plt.pcolormesh(alpha_m,np.degrees(phi_m),np.abs(np.sin(phi_m/2)-lam/2/Lam) \
#               + np.abs(np.sin(np.radians(alpha_m))+np.cos(phi_m/2))\
#               ,vmin=-0.5,vmax=0.5,cmap='seismic')
#plt.yticks(np.linspace(0,360,10))
#plt.colorbar()

plt.figure()
plt.pcolormesh(alpha_m,np.degrees(phi_m),m)
plt.xlabel(r'$\alpha$ [$^\circ$]')
plt.ylabel(r'$\varphi$ [$^\circ$]')
plt.colorbar(label=r'$\Phi^{-^2}$')
plt.yticks(np.linspace(0,360,10))
plt.clim(0,1)
plt.suptitle('$\\Phi$ angular function $\\varphi$ and $\\alpha$ dependence ' + \
             '($\\theta = \\pi/2$)')
plt.title('EM prop.dir. at $\\varphi = 180^\\circ$, acoustic prop.dir. at ' + \
          '$\\varphi = 180^\\circ + \\alpha$',fontsize=10)

# %% Some angles
plt.figure()
alpha = np.array([15,30,60,180])

plt.polar(np.pi*np.ones(2),np.arange(2),'k')

c = {alpha[0]:'r',alpha[1]:'g',alpha[2]:'b',alpha[3]:'y'}

for a in alpha:
    p = pphi(-1,theta,phi,a)**2
    plt.polar(phi,p,color=c[a])

plt.legend(['EM']+list(alpha))

for a in alpha:
    plt.polar((np.pi+np.radians(a))*np.ones(2),np.arange(2),color=c[a])