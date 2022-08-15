#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# to use, simply type on the terminal: python3 randomDroplet.py
#
# Import the modules
#
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#
# These inputs are just examples. 
# Tune them according to your case
#
lx      = 4.0  # domain dimension along x
ly      = 2.0  # domain dimension along y
lz      = 4.0  # domain dimension along z
nx      = 1024 # points along x
ny      = 512  # points along y 
nz      = 1024 # points along x
ppd     = 32   # grid points per diameter (of the dispersed phase)
alpha   = 0.05 # gas/liquid volume fraction
doWrite = 1    # write or not the results in bub.in
dp_dist = 1.2  # normalized intra-droplets distance
#
# Preliminary calculations
#
dx = lx/nx
dy = ly/ny
dz = lz/nz
d  = ppd*dx
V  = lx*ly*lz
Vb = 4/3.*np.pi*(d/2)**3
nb = int(np.floor(alpha*V/Vb))
#
# Generation of the random dispersed 
# phase distribution
# 
xb = np.zeros(nb)
zb = np.zeros(nb)
yb = np.zeros(nb)
domchk = np.zeros(6)
for i in range(nb):
    #print(i)
    pen_chk = True
    while pen_chk:
        #
        xb[i] = np.random.rand()*lx
        yb[i] = np.random.rand()*ly
        zb[i] = np.random.rand()*lz
        chk = np.zeros(i)
        # 
        domchk[0] = xb[i]-d/2 < 0
        domchk[1] = yb[i]-d/2 < 0
        domchk[2] = zb[i]-d/2 < 0
        domchk[3] = xb[i]+d/2 > lx
        domchk[4] = yb[i]+d/2 > ly
        domchk[5] = zb[i]+d/2 > lz
        # 
        if np.sum(domchk)==0:
            for j in range(i):
                chk[j] =  ((xb[j]-xb[i])**2 +
                           (yb[j]-yb[i])**2 +
                           (zb[j]-zb[i])**2)**0.5 < d*dp_dist
        #       
            if np.sum(chk)>=1:
               pen_chk = True
            else:
               pen_chk = False
#
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
ax.scatter(xb,yb,zb,s=Vb*1e3)
#
# print the bub.in files
#
if doWrite:
    f = open("bub.in","w")
    for i in range(nb):
        f.write("%.15f %.15f %.15f %.15f \n" %(xb[i],yb[i],zb[i],d/2))
    f.close()
    
    
