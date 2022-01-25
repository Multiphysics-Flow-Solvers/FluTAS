#
# note: use python3
#
# import needed libraries 
# 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import numpy as np
import scipy as sc
import glob
import pickle
import os
import sys as sys
#from dropFun import *
#
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
sys.setrecursionlimit(10000)
#
# definition of auxiliaries functions
# 
def readBin(filename):
    f   = open(filename,'rb')
    data = np.fromfile(f,dtype='float64')
    f.close()
    return data
#
def readBinI(filename):
    f   = open(filename,'rb')
    data = np.fromfile(f,dtype='int32')
    f.close()
    return data
#
def newTag(i,nd):
    for j in range(6):
        # print i,j
        if drop['procShared'][i,j]>-1:
           logCheck = np.logical_and((drop['myid']==drop['procShared'][i,j]), (drop['id']==drop['idShared'][i,j]))
           lid = np.squeeze(logCheck.nonzero())
           drop['procShared'][i,j]=-1
           if lid.size>0 and newId[lid] ==0:
               #print i,j, lid
               newId[lid] = nd
               newTag(lid,nd)
    return
#
def cordCorr(xi,xavg,lx,d,corrF):
    if corrF:
        if abs(xavg-xi)>0.9*d: # NOTE!!!!!!! approx factor! NOT VALID FOR VERY LONG LIGAMENTS
            sign = xavg-lx/2
            sign = sign/abs(sign)
            xf   = xi+sign*lx
            print(corrF)
        else:
            xf = xi
    else:
        xf = xi
    return xf
#
# Parameters
#
nx   = 32
ny   = 32 
nz   = 32
d0   = 0.50
maxD = 0.70 # used to check for droplet across boundaries. 
path = '../../src/data/post/tagging/' # to be modified accordingly
#
lx = 1
ly = 1
lz = 1
#
dx = lx/nx
dy = ly/ny
dz = lz/nz
#
vol0  = (np.pi/6.0)*d0**3
area0 = (np.pi/1.0)*d0**2
time  = '0000000' 
drop  = {}
#
# quantities to be plotted
#
drop_files = np.sort(glob.glob(path+'/xposfld*'))
timeSeries = {}
#
for nfile in range(len(drop_files)):
   #
   time = drop_files[nfile][-11:-4]
   print (time)
   #
   # read the binary files - volumetric
   #
   drop = {}
   drop['x_dp'] = readBin(path+'xposfld'+time+'.bin') # position-x
   drop['y_dp'] = readBin(path+'yposfld'+time+'.bin') # position-y
   drop['z_dp'] = readBin(path+'zposfld'+time+'.bin') # position-z
   drop['u_dp'] = readBin(path+'uvelfld'+time+'.bin') # velocity-x
   drop['v_dp'] = readBin(path+'vvelfld'+time+'.bin') # velocity-y
   drop['w_dp'] = readBin(path+'wvelfld'+time+'.bin') # velocity-z
   drop['vold'] = readBin(path+'voldfld'+time+'.bin') # volume of the droplet
   #
   # read the binary files - processors
   # 
   drop['id']         = readBinI(path+'dridfld'+time+'.bin')
   drop['procShared'] = readBinI(path+'procfld'+time+'.bin')
   drop['idShared']   = readBinI(path+'idshfld'+time+'.bin')
   drop['order']      = readBinI(path+'ordefld'+time+'.bin')
   drop['procShared'] = np.reshape(drop['procShared'],(len(drop['x_dp']),6),order='c')
   drop['idShared']   = np.reshape(drop['idShared'],(len(drop['x_dp']),6),order='c')
   drop['order']      = np.reshape(drop['order'],(len(drop['x_dp']),6),order='c')
   drop['myid']       = readBinI(path+'idmyfld'+time+'.bin')
   #
   nd = 0
   newId = np.zeros(len(drop['x_dp']))
   #
   # i index cycles over the disperse phase number
   #
   for i in range(len(drop['x_dp'])):
       if newId[i]==0:
           nd += 1
           newId[i] = nd
           if (drop['procShared'][i,:]>-1).any():
               newTag(i,nd)
               
   final = {}
   for i in range(1,nd+1):
       #
       corrX=corrY=corrZ= False
       #
       # Determine if the droplets is crossing a boundary, by checking the difference
       # between the min and max coordinates
       #
       spacX = abs(np.min(drop['x_dp'][newId==i])-
                   np.max(drop['x_dp'][newId==i]))/lx
       spacY = abs(np.min(drop['y_dp'][newId==i])-
                   np.max(drop['y_dp'][newId==i]))/ly
       spacZ = abs(np.min(drop['z_dp'][newId==i])-
                   np.max(drop['z_dp'][newId==i]))/lz
       #
       # Set the correction flag accordingly
       #
       corrX = spacX>maxD
       corrY = spacY>maxD
       corrZ = spacZ>maxD
       #print(corrZ)
       #
       # compute the approximate center and diameter
       #
       x_tmp = np.sum(drop['x_dp'][newId==i]*drop['vold'][newId==i])/np.sum(drop['vold'][newId==i])
       y_tmp = np.sum(drop['y_dp'][newId==i]*drop['vold'][newId==i])/np.sum(drop['vold'][newId==i])
       z_tmp = np.sum(drop['z_dp'][newId==i]*drop['vold'][newId==i])/np.sum(drop['vold'][newId==i])           
       d_tmp = (6/np.pi*np.sum(drop['vold'][newId==i]))**(1/3.)
       #
       final[i] = {'x_dp':0, 'y_dp':0, 'z_dp':0,
                   'u_dp':0, 'v_dp':0, 'w_dp':0, 'vold':0}
       #
       # some quantities are shared among processors, so we compute an 
       # weighted average based on:
       #  --> volume for volumetric quantities (pos,vel);
       # Note: j index cycles over processors 
       #
       
       for j in range(len(newId)):
           if newId[j]==i:
              # print  i,j, newId[j]
              #
              # volumetric quantities
              #
              xf = cordCorr(drop['x_dp'][j],x_tmp,lx,d_tmp,corrX)
              yf = cordCorr(drop['y_dp'][j],y_tmp,ly,d_tmp,corrY)
              zf = cordCorr(drop['z_dp'][j],z_tmp,lz,d_tmp,corrZ)
              final[i]['x_dp'] += xf*drop['vold'][j] # position-x
              final[i]['y_dp'] += yf*drop['vold'][j] # position-y
              final[i]['z_dp'] += zf*drop['vold'][j] # position-z
              final[i]['u_dp'] += drop['u_dp'][j]*drop['vold'][j] # velocity-x
              final[i]['v_dp'] += drop['v_dp'][j]*drop['vold'][j] # velocity-y
              final[i]['w_dp'] += drop['w_dp'][j]*drop['vold'][j] # velocity-z
              final[i]['vold'] += drop['vold'][j]                 # volume
              #
       #        
       # volumetric quantites 
       #
       final[i]['x_dp'] = final[i]['x_dp']/final[i]['vold']
       final[i]['y_dp'] = final[i]['y_dp']/final[i]['vold']
       final[i]['z_dp'] = final[i]['z_dp']/final[i]['vold']
       final[i]['u_dp'] = final[i]['u_dp']/final[i]['vold']
       final[i]['v_dp'] = final[i]['v_dp']/final[i]['vold']
       final[i]['w_dp'] = final[i]['w_dp']/final[i]['vold']
       final[i]['vold'] = final[i]['vold']
       #
   #
   # normalization to avoid to write many digits
   # 
   for i in range(1,nd+1):
      final[i]['vold'] = final[i]['vold']/vol0
   #
   # number of dispersed phases
   #
   f = open('num_dp_post.out','a')
   f.write("%.8f %.8f \n" %(float(time),float(nd))) 
   f.close()
   #
   # print quantities for each investigated time-step 
   # 
   f = open('stat_dp_%s.out' % time, 'w')
   for i in range(1,nd+1):
      f.write("%.1f %.6f %.6f %.6f %.6f %.4f %.4f %.4f \n" %(float(i),final[i]['vold'],final[i]['x_dp'],final[i]['y_dp'],final[i]['z_dp'],final[i]['u_dp'],final[i]['v_dp'],final[i]['w_dp']))
   f.close()
   #
   # define the time-series with all information 
   #
   timeSeries[time]=final[i]
   #
#
# example of print
#
#plt.close('all')
#x = [timeSeries[i]['z_dp'] for i in timeSeries.keys()]
#plt.plot(x,'-*')
#plt.show()

