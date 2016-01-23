#!/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from future_builtins import zip

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:s:')
  except GetoptError:
    print 'randomWalk.py -t <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    if opt == '-o':
      outFile = arg
    if opt == '-n':
      nFrames = arg
    if opt == '-s':
      sat = arg

  return trajFile, outFile, int(nFrames), sat
  
def main():

  # Command line args processing
  trajFile, outFile, nFrames, saturation = getArgs(argv)

  nDims = 3

  MSD = np.zeros((nFrames-1,2))
  MSD[:,0] = np.arange(1,nFrames)
  # outFile = '%s_N%d' % (outFile, nMon)

  if saturation == 'under':
    nPoly = 175000
  else:
    nPoly = 350000

  with open(outFile, 'w') as otp:
    otp.write('# dn <dr^2>\n')

  for n in range(nPoly):
    atom = np.loadtxt(trajFile, delimiter=' ', usecols=[n])
    atom = atom.reshape((-1,nDims))
    atom = atom[:nFrames]
    for i in range(nFrames-1):
      diff = atom[i] - atom[i+1:]
      dist = diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2
      MSD[:nFrames-i-1,1] += dist
    print 'polymer %d MSD written' % (n+1) 
  MSD[:,1] = MSD[:,1] / MSD[::-1,0] / nPoly

  with open(outFile, 'a') as otp:
    np.savetxt(otp, MSD, fmt='%.5f') 
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
