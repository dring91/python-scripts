#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt

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

def PBC(dx,x,interval):
  mask = (dx < interval[:,0])
  x = x*mask + 2*interval[:,1]*mask + x*np.logical_not(mask)
  mask = (dx >= interval[:,1])
  x = x*mask + 2*interval[:,0]*mask + x*np.logical_not(mask)
    
  return x
  
def unwrap(coords, box):
  coords[:,:2] = PBC(coords[:,:2]-coords[0,:2],coords[:,:2],box[:2,:])
      
  return coords
 
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

def writeSizes(filename, sizes, mode='w'):

  with open(filename, mode) as file:
    for row in sizes:
      file.write("%d %f\n" % tuple(row))
  
def main():

  trajFile, outFile, nFrames = getArgs(argv, 'i:o:n:', 'ssi')

  nBins = 100

  rc_max = 30
  nParticles = 63
  nPoints = 1000000
  D1, D2 = 1, 2
  # calculate porosity for a packing of spherical particles
  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)
      
      # calculate box size on a wrapped configuration!
      box = np.zeros((3,2))
      box[:,0] = frame[:,1:].min(0)
      box[:,1] = frame[:,1:].max(0)
      binSize = (box[2,1] - box[2,0]) / nBins

      # generate insertion points
      points = np.random.rand(nPoints, 3)
      points = points * (box[:,1] - box[:,0]) + box[:,0]
      test = np.zeros((nPoints, rc_max), dtype=bool)

      # find the particle centers
      COM = np.zeros((nParticles,3))
      for p in range(nParticles):
        # filter and unwrap every particle
        particle = frame[frame[:,0] == (p + 1)][:,1:]
        particle = unwrap(particle, box)

        # calculate center of mass of each particle and then calculate r
        COM[p] = particle.mean(0)
        centered = particle - COM[p]
        r2 = np.einsum('ij,ij->i',centered,centered)
        r = np.sqrt(r2.mean())

        # apply MIC on points using COM from the particles
        centerPoints = points - COM[p]
        centerPoints[:,:2] = PBC(centerPoints[:,:2], 
                                 centerPoints[:,:2], box[:2,:])

        # test if the insertion points are in the spheres
        dist2 = np.einsum('ij,ij->i',centerPoints,centerPoints)
        for rc in range(rc_max):
          test[:,rc] = np.logical_or(dist2 < (r + rc)**2, test[:,rc])

      # Calculate bulk porosity
      test = np.logical_not(test)
      sizes = np.zeros((rc_max,2))
      sizes[:,1] = test.sum(0) / float(test.sum())
      sizes[:,0] = range(rc_max)

      # plot sizes versus probability
      # plt.plot(range(rc_max), np.log(distribution), 'o')
      writeSizes(outFile, sizes)
      plt.show()
         
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
