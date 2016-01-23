#!/bin/python

import numpy as np
from sys import argv
from getopt import *

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <filename> -o <filename> -n <steps>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      steps = arg

  return trajFile, outFile, int(steps)
 
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
 
def main():

  inFile, outFile, nSteps = getArgs(argv)

  # read in the (final?) frame from the trajectory
  with open(inFile, 'r') as inp:
    time, particles = readFrame(inp)
  zs = particles[:,3]
  # Find the isolated particles (ie the centers)
  for z in zs:
    # print len(zs[abs(zs-z) < 0.001])
    if len(zs[abs(zs-z) < 0.001]) < 1:
      print particles[zs-z < 0.001]
      break
  # Relabel the particles based on the centers
  # write to a file

if __name__ == "__main__":
  main()
