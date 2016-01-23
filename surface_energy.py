#!/usr/python

from sys import argv
import numpy as np
import matplotlib.pyplot as plt

def energy(P, Lz):
  P = np.array(P, dtype=float)
  P = P[1:].reshape((-1,10,3))
  ave = P.mean(0)
  s = P.var(0)
  s = np.sqrt(s.sum(1))
  gamma = 0.5*Lz*(ave[:,2]-0.5*(ave[:,0]+ave[:,1]))

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

def weighted_std(x,w):

  ave = np.average(x,weights=w)
  residuals = x-ave
  SSR = sum(residuals**2)
  var = SSR/(len(x)-1)
  std = np.sqrt(var)

  return std

def main():
  
  thermoFile = argv[1]
  coordFile = argv[2]
  nFrames = int(argv[3])
  nBins = int(argv[4])

  # with open(thermoFile, 'r') as file:
  #   P = []
  #   for line in file:
  #     L = line.split()
  #     if len(L) > 0 and L[0] == '#':
  #       heading = L[1:]
  #       indices = []
  #       indices.append(heading.index('Pzz'))
  #       indices.append(heading.index('Pxx'))
  #       indices.append(heading.index('Pyy'))
  #     elif len(L) > 0:
  #       P.append([L[i] for i in indices])

  with open(coordFile, 'r') as file:
    ave = np.zeros((nBins,2))
    for frame in range(nFrames):
      time, atoms = readFrame(file)
      atoms = atoms[atoms[:,0] == 1][:,1:]
      hist = np.ones((nBins,2))
      hist[:,0] = 100*hist[:,0]
      hist[:,1] = 0*hist[:,1]

      limits = np.zeros(2)
      limits[0] = atoms[:,1].min()
      limits[1] = atoms[:,1].max()

      d = (limits[1]-limits[0])/nBins
      for atom in atoms:
        bin = int((atom[1]-limits[0])/d)
        if bin >= nBins:
          bin = 0
        if hist[bin,1] < atom[2]:
          hist[bin,1] = atom[2]
        if hist[bin,0] > atom[2]:
          hist[bin,0] = atom[2]
      hist[:,0] = hist[:,0]*(hist[:,0] != 100)

      ave[:,1] += hist[:,1]-hist[:,0]

    ave[:,0] = (np.arange(nBins)+0.5)*d + limits[0]
    ave[:,1] /= nFrames
    # stdev = ave[:,0].std(weights=ave[:,1])
    stdev = weighted_std(ave[:,0],ave[:,1])
    # plt.axis('equal')
    plt.xlim(-stdev/2**0.25,stdev/2**0.25)
    plt.plot(ave[:,0],ave[:,1])
    plt.plot(ave[:,0],18*np.exp(-(ave[:,0]/stdev)**2))
    # plt.plot(ave[:,0],ave[:,1]-18*np.exp(-(ave[:,0]/stdev)**2))
    plt.plot(stdev/2**0.25, 18*np.exp(-1/(2**0.5)),'ro')
    plt.plot(-stdev/2**0.25, 18*np.exp(-1/(2**0.5)),'ro')
    # plt.plot(ave[ave[:,0] > 0][:,0],
    plt.show()

  # create histogram along y and calculate the max in each bin
  # average the max values

if __name__ == '__main__':
  main()
