# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:21:46 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
lam = 5e-3
Lam = 5e-3
k = 2*np.pi/lam
q = 2*np.pi/Lam

# Machine epsilon used to avoid "divide by zero" errors in max tracing
eps = np.finfo(float).eps

# Phi angular function
def pphi(pm,theta,phi,alpha):
    L = 1e-2
    A = np.sinc(L/2/np.pi*(k-k*np.sin(theta)*np.cos(phi) \
                           + np.sign(pm)*q*np.cos(alpha)))
    B = np.sinc(L/2/np.pi*(-k*np.sin(theta)*np.sin(phi) \
                           + np.sign(pm)*q*np.sin(alpha)))
    C = np.sinc(-L/2/np.pi*k*np.cos(theta))
    return A*B*C

theta = np.pi/2
phi = np.linspace(0,2*np.pi,500)

# %% Various alpha, plot function value for all phi,alpha combinations
plt.close('all')
namestr = 'equal'
alpha = np.radians(np.linspace(0,180,500))

[phi_m,alpha_m] = np.meshgrid(phi,alpha)

p = pphi(1,theta,phi_m,alpha_m)**2
m = pphi(-1,theta,phi_m,alpha_m)**2

plt.figure()
plt.pcolormesh(np.degrees(alpha_m),np.degrees(phi_m),p,cmap='plasma')
plt.xlabel(r'$\alpha$ [$^\circ$]')
plt.ylabel(r'$\varphi$ [$^\circ$]')
plt.colorbar(label=r'$\Phi^{+^2}$')
plt.yticks(np.linspace(0,360,10))
plt.clim(0,1)
plt.suptitle('$\\Phi$ angular function $\\varphi$ and $\\alpha$ dependence ' + \
             '($\\theta = \\pi/2$)')
plt.title('EM prop.dir. at $\\varphi = 180^\\circ$, acoustic prop.dir. at ' + \
          '$\\varphi = 180^\\circ + \\alpha$',fontsize=10)

# Plot point of maximum
phim = 2*np.arctan(np.sqrt(lam**2/(4*Lam**2-lam**2+eps)))
alpham = np.arccos(-lam/2/Lam)
plt.scatter(np.degrees(alpham),np.degrees(phim),marker='.')
plt.xlim([0,180])

#plt.savefig('Phi+'+namestr+'.png',dpi=300)

#plt.figure()
#plt.pcolormesh(np.degrees(alpha_m),np.degrees(phi_m),\
#               np.abs(np.cos(alpha_m)-q/2/k) + np.abs(np.tan(phi_m/2)+np.sqrt(q**2/(4*k**2-q**2))) ,\
#               vmin=-0.2,vmax=0.2,cmap='seismic')
#plt.yticks(np.linspace(0,360,10))
#plt.colorbar()

plt.figure()
plt.pcolormesh(np.degrees(alpha_m),np.degrees(phi_m),m,cmap='plasma')
plt.xlabel(r'$\alpha$ [$^\circ$]')
plt.ylabel(r'$\varphi$ [$^\circ$]')
plt.colorbar(label=r'$\Phi^{-^2}$')
plt.yticks(np.linspace(0,360,10))
plt.clim(0,1)
plt.suptitle('$\\Phi$ angular function $\\varphi$ and $\\alpha$ dependence ' + \
             '($\\theta = \\pi/2$)')
plt.title('EM prop.dir. at $\\varphi = 180^\\circ$, acoustic prop.dir. at ' + \
          '$\\varphi = 180^\\circ + \\alpha$',fontsize=10)

# Plot point of maximum, +2pi since negative arctan gives negative angle
phim = 2*np.arctan(-np.sqrt(lam**2/(4*Lam**2-lam**2+eps))) + 2*np.pi
alpham = np.arccos(lam/2/Lam)
plt.scatter(np.degrees(alpham),np.degrees(phim),marker='.')
plt.xlim([0,180])

#plt.savefig('Phi-'+namestr+'.png',dpi=300)

# %% Variation of max with changing wavelength fraction - polar plot (+)
plt.close('all')
plt.figure()

# Fraction EM lambda/AO lambda
lamfr = np.array([0.1,1,1.99])
alpha = np.arccos(-lamfr/2)
phi = 2*np.arctan(np.sqrt(lamfr**2/4/(1-lamfr**2/4+eps)))

# Incident EM
plt.polar(np.pi*np.ones(2),np.arange(2),'k',label='$\mathbf{\hat{k}}$')
# Arrow for incident EM
plt.polar((np.pi-0.2)*np.ones(2),np.array([0,0.1]),'k',\
          (np.pi+0.2)*np.ones(2),np.array([0,0.1]),'k')

# Colors
c = ['r','g','b']

for i in range(0,len(lamfr)):
    plt.polar((np.pi+alpha[i])*np.ones(2),np.arange(2),c[i],linestyle=':',\
              label='$\mathbf{\hat{q}}, \lambda/\Lambda = $'+str(lamfr[i]))
    plt.polar(alpha[i]*np.ones(2),np.arange(2),c[i],linestyle=':')
    # Arrow for sound
    plt.polar(np.array([alpha[i]+0.03,alpha[i]]),np.array([0.9,1]),c[i],\
          np.array([alpha[i]-0.03,alpha[i]]),np.array([0.9,1]),c[i])
    
    plt.polar(phi[i]*np.ones(2),np.arange(2),c[i],\
              label='$\mathbf{\hat{k}}_{sc}, \lambda/\Lambda = $'+str(lamfr[i]))
    # Arrow for scattered EM
    plt.polar(np.array([phi[i]+0.03,phi[i]]),np.array([0.9,1]),c[i],\
          np.array([phi[i]-0.03,phi[i]]),np.array([0.9,1]),c[i])

plt.title('Optimal scattering geometry for different $\lambda/\Lambda$ (+ case)')
plt.ylim(0,1)
plt.yticks(np.linspace(0,1,5),[])

plt.figlegend(loc=4,framealpha=1)

# %% Variation of max with changing wavelength fraction - polar plot (-)
plt.close('all')
plt.figure()

# Fraction EM lambda/AO lambda
lamfr = np.array([0.1,1,1.99])
alpha = np.arccos(lamfr/2)
phi = 2*np.arctan(-np.sqrt(lamfr**2/4/(1-lamfr**2/4+eps))) + 2*np.pi

# Incident EM
plt.polar(np.pi*np.ones(2),np.arange(2),'k',label='$\mathbf{\hat{k}}$')
# Arrow for incident EM
plt.polar((np.pi-0.2)*np.ones(2),np.array([0,0.1]),'k',\
          (np.pi+0.2)*np.ones(2),np.array([0,0.1]),'k')

# Colors
c = ['r','g','b']

for i in range(0,len(lamfr)):
    plt.polar((np.pi+alpha[i])*np.ones(2),np.arange(2),c[i],linestyle=':',\
              label='$\mathbf{\hat{q}},$ $\lambda/\Lambda = $'+str(lamfr[i]))
    plt.polar(alpha[i]*np.ones(2),np.arange(2),c[i],linestyle=':')
    # Arrow for sound
    plt.polar(np.array([alpha[i]+0.03,alpha[i]]),np.array([0.9,1]),c[i],\
          np.array([alpha[i]-0.03,alpha[i]]),np.array([0.9,1]),c[i])
    
    plt.polar(phi[i]*np.ones(2),np.arange(2),c[i],\
              label='$\mathbf{\hat{k}}_{sc},$ $\lambda/\Lambda = $'+str(lamfr[i]))
    # Arrow for scattered EM
    plt.polar(np.array([phi[i]+0.03,phi[i]]),np.array([0.9,1]),c[i],\
          np.array([phi[i]-0.03,phi[i]]),np.array([0.9,1]),c[i])

plt.title('Optimal scattering geometry for different $\lambda/\Lambda$ (- case)')
plt.ylim(0,1)
plt.yticks(np.linspace(0,1,5),[])

plt.figlegend(loc=4,framealpha=1)

# %% Polar plot of max scattering geometry - all cases
plt.close('all')
plt.figure()

# (+) or (-) case
pm = -1
pmch = '+'
if(pm is -1):
    pmch = '-'

# Fraction EM lambda/AO lambda
lamfr = 0
alpha = np.arccos(-pm*lamfr/2)
phi = 2*np.arctan(pm*np.sqrt(lamfr**2/4/(1-lamfr**2/4+eps))) + np.pi*(1-pm)

# Incident EM
plt.polar(np.pi*np.ones(2),np.arange(2),'k',label='$\mathbf{\hat{k}}$')
# Arrow for incident EM
plt.polar((np.pi-0.2)*np.ones(2),np.array([0,0.1]),'k',\
          (np.pi+0.2)*np.ones(2),np.array([0,0.1]),'k')

# Incident sound
plt.polar((np.pi+alpha)*np.ones(2),np.arange(2),':b',\
          label='$\mathbf{\hat{q}}$')
plt.polar(alpha*np.ones(2),np.arange(2),'b:')
# Arrow for sound
plt.polar(np.array([alpha+0.03,alpha]),np.array([0.9,1]),'b',\
      np.array([alpha-0.03,alpha]),np.array([0.9,1]),'b')

#Scattered EM
plt.polar(phi*np.ones(2),np.arange(2),'r',\
          label='$\mathbf{\hat{k}}_{sc}$')
# Arrow for scattered EM
plt.polar(np.array([phi+0.03,phi]),np.array([0.9,1]),'r',\
      np.array([phi-0.03,phi]),np.array([0.9,1]),'r')

plt.title('Optimal scattering geometry for $\lambda/\Lambda = $' + str(lamfr) +\
          ' (' + pmch + ' case)')
plt.ylim(0,1)
plt.yticks(np.linspace(0,1,5),[])

plt.figlegend(loc=4,framealpha=1)

# %% Some angles
plt.figure()
alpha = np.radians(np.array([15,30,60,180]))

plt.polar(np.pi*np.ones(2),np.arange(2),'k')

c = {alpha[0]:'r',alpha[1]:'g',alpha[2]:'b',alpha[3]:'y'}

for a in alpha:
    p = pphi(-1,theta,phi,a)**2
    plt.polar(phi,p,color=c[a])

plt.legend(['EM']+list(np.degrees(alpha)))

for a in alpha:
    plt.polar((np.pi+a)*np.ones(2),np.arange(2),color=c[a])