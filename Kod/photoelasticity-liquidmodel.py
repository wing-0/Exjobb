# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:48:05 2018

@author: tfy13nwi
"""

import numpy as np
import matplotlib.pyplot as plt

er = np.linspace(1,10)
p = 1/3*(er-1)*(er+2)/er**2

plt.figure(figsize=(5,4))

plt.plot(er,p)
plt.xlabel('$\\varepsilon_\\mathrm{r}$')
plt.ylabel('p')

plt.rcParams['axes.unicode_minus'] = False
plt.savefig('../Text/Report/fig/photoelastic-liquid.pgf')