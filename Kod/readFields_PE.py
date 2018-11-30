# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 16:26:57 2018

@author: tfy13nwi
"""

import numpy as np


def readEHS_PE(files_pe, files_nope):

    # Lists for all data
    xt = []
    yt = []
    Ezt = []
    Hxt = []
    Hyt = []

    # Read data from both files
    for f_pe, f_nope in zip(files_pe, files_nope):
        data_pe = np.genfromtxt(f_pe, delimiter=',', comments='%', dtype=str)
        data_nope = np.genfromtxt(f_nope, delimiter=',', comments='%',
                                  dtype=str)
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

    # Calculate difference fields for all angles
    Ez = np.array([Ezt[i][:, 0] - Ezt[i][:, 1] for i in range(0, len(Ezt))]).T
    Hx = np.array([Hxt[i][:, 0] - Hxt[i][:, 1] for i in range(0, len(Hxt))]).T
    Hy = np.array([Hyt[i][:, 0] - Hyt[i][:, 1] for i in range(0, len(Hyt))]).T

    # Convert x,y-data to correct array format
    x = np.array(xt).T
    y = np.array(yt).T

    # Calculate time avg. power flow
    Sx = -0.5*(Ez*np.conj(Hy)).real
    Sy = 0.5*(Ez*np.conj(Hx)).real

    return Ez, Hx, Hy, Sx, Sy, x, y


def readEHS_PE_various(files_pe, files_nope):

    # Lists for all data
    xt = []
    yt = []
    Ezt = []
    Hxt = []
    Hyt = []

    # Read data from both files
    for f_pe, f_nope in zip(files_pe, files_nope):
        data_pe = np.genfromtxt(f_pe, delimiter=',', comments='%', dtype=str)
        data_nope = np.genfromtxt(f_nope, delimiter=',', comments='%',
                                  dtype=str)
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

    minlen = min([len(a) for a in xt])
    for i in range(0, len(xt)):
        if(len(xt[i]) > minlen):
            xt[i] = np.delete(xt[i], -1, axis=0)
            yt[i] = np.delete(yt[i], -1, axis=0)
            Ezt[i] = np.delete(Ezt[i], -1, axis=0)
            Hxt[i] = np.delete(Hxt[i], -1, axis=0)
            Hyt[i] = np.delete(Hyt[i], -1, axis=0)

    # Calculate difference fields for all angles
    Ez = np.array([Ezt[i][:, 0] - Ezt[i][:, 1] for i in range(0, len(Ezt))]).T
    Hx = np.array([Hxt[i][:, 0] - Hxt[i][:, 1] for i in range(0, len(Hxt))]).T
    Hy = np.array([Hyt[i][:, 0] - Hyt[i][:, 1] for i in range(0, len(Hyt))]).T

    # Convert x,y-data to correct array format
    x = np.array(xt).T
    y = np.array(yt).T

    # Calculate time avg. power flow
    Sx = -0.5*(Ez*np.conj(Hy)).real
    Sy = 0.5*(Ez*np.conj(Hx)).real

    return Ez, Hx, Hy, Sx, Sy, x, y
