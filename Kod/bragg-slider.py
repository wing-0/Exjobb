# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 19:02:06 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

# Machine epsilon used to avoid "divide by zero" errors in max tracing
eps = np.finfo(float).eps

fig, ax = plt.subplots(subplot_kw=dict(polar=True))
plt.subplots_adjust(left=0.25, bottom=0.25)
pm = 1
lamfr_0 = 1
alpha_0 = np.arccos(-lamfr_0/2)
phi_0 = 2*np.arctan(np.sqrt(lamfr_0**2/4/(1-lamfr_0**2/4+eps)))
delta_lf = 0.1
a_s = alpha_0*np.ones(2)
a_em = phi_0*np.ones(2)
r = np.arange(2)
plt.polar(np.pi*np.ones(2),r,'k')
l_s, = plt.polar(a_s,r,'r:')
l_em, = plt.polar(a_em,r,'b')
plt.ylim(0,1)
plt.yticks(np.linspace(0,1,5),[])

axcolor = 'lightgoldenrodyellow'
axlf = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

slf = Slider(axlf, '$\lambda/\Lambda$', 0, 2, valinit=lamfr_0)

rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('+', '-'), active=0)

def update(val):
    pm_s = radio.value_selected
    pm = 1
    if(pm_s is '-'):
        pm = -1
    lamfr = slf.val
    alpha = np.arccos(-pm*lamfr/2)
    phi = 2*np.arctan(pm*np.sqrt(lamfr**2/4/(1-lamfr**2/4+eps))) + np.pi*(1-pm)
    l_s.set_xdata(alpha*np.ones(2))
    l_em.set_xdata(phi*np.ones(2))
    fig.canvas.draw_idle()
slf.on_changed(update)

def pmfunc(label):
    update(slf.val)
radio.on_clicked(pmfunc)

