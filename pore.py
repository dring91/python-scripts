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
 
def writeHeader(filename,header):
  # Implicitly assumes a 1 line header
  with open(filename, 'w') as file:
    file.write(header)   
  
def writeFile(filename,density,time,mode):
  with open(filename, mode) as file:
    file.write('#  time: %d \n#   z   density\n' % time)
    [file.write('%f %f\n' % tuple(row)) for row in density]
    file.write('\n')

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'t:p:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <xyz file> -p <output path> -o <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      inFile = arg
    elif opt == '-p':
      path = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      nSteps = arg
 
  return inFile, path, outFile, int(nSteps)
  
def main():

  trajFile, path, outFile, nFrames = getArgs(argv)

  poreFile = path+outFile
  
  # title = '# Particle distance distribution\n'
  # writeHeader(poreFile,title+'# time  cut = 0.85  cut = 0.99\n')

  mode = 'w'

  # Pore size
  # 1. Find the distribution of lengths to nearest particles
  #    a. Calculate the distances between each polymer COM and all particle atoms
  #    b. Find the minimum distance for each polymer
  #    c. Make a histogram of these distances
  # 2. Pick a representative value from the distribution
  #    a. The mean
  #    b. The minimum

  # Insertion
  # 1. Use the porosity script, but keep the number and locations of the particles inserted (or the distances between the particles, or maybe correlation data)
  # 2. Perform insertion using a range of insertion particle sizes
  # 3. Construct a histogram of the number of accepted insertions for each particle size
  # 4. Pick a representative value (may use the locations of inserted particles as a guide)

  maxPoreSize = 30
  nBins = 50
  binSize = maxPoreSize/100.0
  hist = np.zeros((nBins,2))
  hist[:,0] = binSize*(0.5 + np.arange(nBins))
  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)

      polymers = frame[frame[:,0] == 1][:,1:]
      polymers = polymers.reshape(-1,10,3).mean(1)
      particles = frame[frame[:,0] == 3][:,1:]
      bottom = particles[:,2].min()
      polymers = polymers[polymers[:,2] > bottom]
      
      # distances = np.zeros(polymers.shape)
      for i,polymer in enumerate(polymers):
        disp = particles - polymer
        dist2 = np.einsum('ij,ij->i',disp,disp)
        minDist = np.sqrt(dist2).min()
        # distances[i] = minDist

        index = int(minDist / binSize)
        try:
          hist[index,1] += 1
        except IndexError:
          print 'Index out of bounds'
        else:
          print 'finished polymer %d' % (i+1)

  hist[:,1] /= nFrames
  writeFile(poreFile,hist,t,mode)
      
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
