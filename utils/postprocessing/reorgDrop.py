#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:57:46 2020

@author: marcres
"""
import matplotlib.pyplot as plt
import numpy as np
def readBin(filename):
    f   = open(filename,'rb')
    data = np.fromfile(f,dtype='float64')
    f.close()
    return data

def readBinI(filename):
    f   = open(filename,'rb')
    data = np.fromfile(f,dtype='int32')
    f.close()
    return data

def newTag(i,nd):
    for j in xrange(6):
        # print i,j
        if drop['procShared'][i,j]>-1:
           logCheck = np.logical_and((drop['myid']==drop['procShared'][i,j]), (drop['id']==drop['idShared'][i,j]))
           lid = np.squeeze(logCheck.nonzero())
           drop['procShared'][i,j]=-1
           if lid.size>0 and newId[lid] ==0:
               print i,j, lid
               newId[lid] = nd
               newTag(lid,nd)
    return
           
            
nx=ny=nz=16
dx = 1./nx
time = '0000000' 
drop = {}
drop['x'] = readBin  ('data/post/tagging/xposfld'+time+'.bin')
drop['y'] = readBin  ('data/post/tagging/yposfld'+time+'.bin')
drop['z'] = readBin  ('data/post/tagging/zposfld'+time+'.bin')
drop['u'] = readBin  ('data/post/tagging/uvelfld'+time+'.bin')
drop['v'] = readBin  ('data/post/tagging/vvelfld'+time+'.bin')
drop['w'] = readBin  ('data/post/tagging/wvelfld'+time+'.bin')
drop['vol'] = readBin('data/post/tagging/voldfld'+time+'.bin')
drop['id'] = readBinI('data/post/tagging/dridfld'+time+'.bin')
drop['procShared'] = readBinI('data/post/tagging/procfld'+time+'.bin')
drop['idShared']   = readBinI('data/post/tagging/idshfld'+time+'.bin')
drop['order']   = readBinI('data/post/tagging/ordefld'+time+'.bin')
drop['procShared'] = np.reshape(drop['procShared'],(len(drop['x']),6),order='c')
drop['idShared'] = np.reshape(drop['idShared'],(len(drop['x']),6),order='c')
drop['order'] = np.reshape(drop['order'],(len(drop['x']),6),order='c')
drop['myid']       = readBinI('data/post/tagging/idmyfld'+time+'.bin')

nd = 0
newId = np.zeros(len(drop['x']))
for i in xrange(len(drop['x'])):
    if newId[i]==0:
        nd += 1
        newId[i] = nd
        if (drop['procShared'][i,:]>-1).any():
            newTag(i,nd)
            
final = {}
for i in xrange(1,nd+1):
    final[i] = {'x':0, 'y':0, 'z':0,
                'u':0, 'v':0, 'w':0,
                'vol':0}
    for j in xrange(len(newId)):
        if newId[j]==i:
            # print  i,j, newId[j]
            final[i]['x'] += drop['x'][j]*drop['vol'][j]
            final[i]['y'] += drop['y'][j]*drop['vol'][j]
            final[i]['z'] += drop['z'][j]*drop['vol'][j]
            final[i]['u'] += drop['u'][j]*drop['vol'][j]
            final[i]['v'] += drop['v'][j]*drop['vol'][j]
            final[i]['w'] += drop['w'][j]*drop['vol'][j]
            final[i]['vol'] += drop['vol'][j]
            
    final[i]['x'] = final[i]['x']/final[i]['vol']
    final[i]['y'] /= final[i]['vol']
    final[i]['z'] /= final[i]['vol']
    final[i]['u'] /= final[i]['vol']
    final[i]['v'] /= final[i]['vol']
    final[i]['w'] /= final[i]['vol']
    print final[i]['x'], final[i]['y'], final[i]['z'], (final[i]['vol']*3/4./np.pi)**0.33333
    



