#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 19:11:11 2020

@author: marcres
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

lx=ly=lz= 2*np.pi
V = lx**3
nx=ny=nz= 64
dx = lx/nx
ppd = 10
d= ppd*dx
Vb = 4/3.*np.pi*(d/2)**3
alpha = 0.15
nb = int(np.floor(alpha*V/Vb))
doWrite = 1


xb = np.zeros(nb)
zb = np.zeros(nb)
yb = np.zeros(nb)
domchk = np.zeros(6)
for i in xrange(nb):
    print i
    pen_chk = True
    while pen_chk:
        xb[i] = np.random.rand()*lx
        yb[i] = np.random.rand()*ly
        zb[i] = np.random.rand()*lz
        chk = np.zeros(i)
        
        domchk[0] = xb[i]-d/2 < 0
        domchk[1] = yb[i]-d/2 < 0
        domchk[2] = zb[i]-d/2 < 0
        domchk[3] = xb[i]+d/2 > lx
        domchk[4] = yb[i]+d/2 > ly
        domchk[5] = zb[i]+d/2 > lz
        
        if np.sum(domchk)==0:
            for j in xrange(i):
                chk[j] =  ((xb[j]-xb[i])**2 +
                           (yb[j]-yb[i])**2 +
                           (zb[j]-zb[i])**2)**0.5 < d*1.2
               
               
            if np.sum(chk)>=1:
               pen_chk = True
            else:
               pen_chk = False
       
       

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xb,yb,zb,s=Vb*1e3)

if doWrite:
    f = open("bub.in","write")
    for i in xrange(nb):
        f.write("%.4f %.4f %.4f %.4f \n" %(xb[i],yb[i],zb[i],d/2))
    f.close()
    
    