#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
from memory_profiler import profile

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'i:o:n:')
  except GetoptError:
    print 'randomWalk.py -i <filename> -o <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      trajFile = arg
    if opt == '-o':
      outFile = arg
    if opt == '-n':
      nFrames = arg

  return trajFile, outFile, int(nFrames)

def readFrame(file, nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def COM(atoms, chainLength):
  size = atoms.shape
  atoms = atoms.reshape((size[0] / chainLength, chainLength, size[1]))
  
  return atoms.mean(1)

def diffusion(atoms_head, atoms_tail, chainLength):
  
  difference = atoms_head - atoms_tail
  difference = COM(difference, chainLength)
  dist2 = difference**2
  dist2 = dist2.mean()

  return dist2
 
@profile
def main():
  filename, outFile, nFrames = getArgs(argv)

  # Generate multiple MSDs based on the locations of the polymers
  # 1. Exclude polymers in the bulk
  # 2. Exclude polymers NOT in the bulk
  # 3. Look at polymers in the packing exclusively
  # 
  # The point is to find a reasonable time-scale for polymer relaxation
  # Polymers in the packing have z coordinate > 0 and those near the surface are less than
  # 5 sigma from the surface (most negative polymer value)
  MSD = np.zeros((nFrames-1,2))
  MSD[:,0] = np.arange(1,nFrames)
  chainLength = 1
  pos = 0
  for m in range(nFrames-1):
    with open(filename, 'r') as file:
      file.seek(pos)
      time_m, atoms_m = readFrame(file)
      # add a condition to restrict the polymers based on location (make a mask)
      mask = np.logical_and(atoms_m[:,2] < 0, atoms_m[:,2] > 5*(atoms_m[:,2].min()))
      atoms_m = atoms_m[:,1:][mask]
      pos = file.tell()
      for n in range(m+1,nFrames):
        time_n, atoms_n = readFrame(file)
        atoms_n = atoms_n[:,1:][mask]
        # handle function calls here
        MSD[n-m-1,1] += diffusion(atoms_m, atoms_n, chainLength)
  # average MSD array
  MSD[:,1] /= MSD[::-1,0]
  # write MSD array to file
  with open(outFile, 'w') as otp:
    otp.write('#  lagtime  MSD\n')
    np.savetxt(otp,MSD,'%0.5f')

if __name__ == '__main__':
  main()
