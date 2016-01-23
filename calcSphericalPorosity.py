#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt

def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

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

def plotTest(data,dims=(0,1),plot=False,show=False):
  if plot == True:
    plt.plot(data[:,dims[0]],data[:,dims[1]], '.')
  if show == True:
    plt.gca().set_aspect('equal')
    plt.show()

def spheres(particle):
  COM = particle.mean(0)
  centered = particle - COM
  r2 = np.einsum('ij,ij->i',centered,centered)
  r = np.sqrt(r2.mean())

  return r, COM

def insert(test, r, COM, points, box, dims=(0,1), r_cut=0):
  # apply MIC on points using COM from the particles
  points = points - COM
  # plotTest(points, dims, plot)
  points[:,:2] = PBC(points[:,:2], 
                     points[:,:2], box[:2,:])
  # plotTest(points, dims, plot)

  dist2 = np.einsum('ij,ij->i',points,points)
  test = np.logical_or(dist2 < (r + r_cut)**2, test)

  return test
  
def main():

  trajFile, outFile, nFrames = getArgs(argv, 'i:o:n:', 'ssi')

  nBins = 100
  porosity = np.zeros((nBins,2))
  porosity[:,0] = (np.arange(nBins) + 0.5) / nBins

  r_cut = 0.00
  nParticles = 63
  nPoints = 1000000
  outFile = "%s_%d" % (outFile, nPoints)
  dims = (0, 1)
  plot = False
  show = False
  # calculate porosity for a packing of spherical particles
  with open(trajFile, 'r') as inp:
    for n in range(nFrames):
      t, frame = readFrame(inp)
      
      # calculate box size on a wrapped configuration!
      box = np.zeros((3,2))
      box[:,0] = frame[:,1:].min(0)
      box[:,1] = frame[:,1:].max(0)
      binSize = (box[2,1] - box[2,0]) / nBins

      # generate insertion points
      points = np.random.rand(nPoints, 3)
      points = points * (box[:,1] - box[:,0]) + box[:,0]
      test = np.zeros_like(points[:,0], dtype=bool)

      # find the particle centers
      for p in range(nParticles):
        # filter and unwrap every particle
        particle = frame[frame[:,0] == (p + 1)][:,1:]
        particle = unwrap(particle, box)

        # calculate center of mass of each particle and then calculate r
        r, COM = spheres(particle)

        # test if the insertion points are in the spheres
        test = insert(test, r, COM, points, box)

        # plotTest(particle-COM, dims, plot, show)

      # partition particles into layers and calculate porosity
      bins = ((points[:,2] - box[2,0]) / binSize).astype(int)
      porosity[:,1] = np.zeros(nBins)

      # If a value is False (accepted), add it's negation to the array
      for i in range(nBins):
        insertions = np.logical_not(test[bins == i])
        porosity[i,1] = sum(insertions) / float(len(insertions))

      mode = 'w'
      with open(outFile, mode) as otp:
        np.savetxt(otp,porosity,fmt="%.4f",header='height,h porosity,phi')
         
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
