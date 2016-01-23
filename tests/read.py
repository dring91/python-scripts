#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'f:n:')
  except GetoptError:
    print 'randomWalk.py -f <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-f':
      trajFile = arg
    if opt == '-n':
      nFrames = arg

  return trajFile, int(nFrames)

def readFrame(file, nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms
 
def main():
  filename, nFrames = getArgs(argv)

  pos = 0
  for m in range(nFrames-1):
    with open(filename, 'r') as file:
      file.seek(pos)
      # save atoms under a different name
      time, atoms = readFrame(file)
      pos = file.tell()
      for n in range(m+1,nFrames):
        time, atoms = readFrame(file)
        # handle function calls here
        # diffusion()
        # write result?

if __name__ == '__main__':
  main()
