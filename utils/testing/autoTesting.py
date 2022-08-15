'''
This file performs the automatic testing for all the cases in .........
'''

import numpy as np
import os 

def readstp(fn):
  f = open(fn, 'r')
  for i in f:
    exec(i)
  f.close()
  globals().update(locals())

# define execution variables
wfold = '../../'
compiler = 'PGI'
tfold = 'utils/testing/templateTest/'
doDebugFlag = True
#conf = 'two_phase_inc_isot'
cdir = os.getcwd()

readstp(wfold+tfold+'test.stp')
os.chdir('../../src')
os.system('cp '+conf+'/* .')
os.system('make -f Makefile.mc clean; make -f Makefile.mc COMP=PGI DBG_FLAG=0')
os.system('cp flutas '+cdir+'/templateTest')
os.chdir(cdir+'/templateTest')
os.system('mpirun -np 4 flutas; pwd')
os.chdir(cdir)
