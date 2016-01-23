#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *

def readFrame(file):
  line = file.readline().strip()
  nAtoms = int(line)
  atoms = np.zeros((nAtoms,4))

  line = file.readline()
  time = int(line.split()[-1])

  for i in range(nAtoms):
    line = file.readline().split()
    atoms[i,:] = np.array(line, dtype=float) 

  return time, atoms
 
def getArgs(argv, flags, types):

  flagArr = flags.split(':')
  flagArr = ''.join(flagArr)
  try:
    opts, _ = getopt(argv[1:], flags)
  except GetoptError:
    print 'randomWalk.py %s <xyz file> %s <filename> %s <nFrames>' % tuple(flagArr)
    exit(2)
  for opt, arg in opts:
    if opt == '-' + flagArr[0]:
      inFile = arg
    elif opt == '-' + flagArr[1]:
      outFile = arg
    elif opt == '-' + flagArr[2]:
      nSteps = arg
 
  # return an array or tuple of the arguments instead and perform the unpacking
  return inFile, outFile, int(nSteps)
  
def main():

  trajFile, outFile, nFrames = getArgs(argv, 'i:on:', 'ssi')

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)

      # How should the particles be masked?
      # particles = frame[mask][:,3:]
      # particles = frame[:,3:] * mask
      for n in range(63):
        mask = (frame[:,0] == (n + 1))
         
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
