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
  
def writeRDF(filename,density,time,mode):
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

  rdfFile = path+outFile
  
  title = '# Radial distribution function\n'
  writeHeader(rdfFile,title+'# time  cut = 0.85  cut = 0.99\n')

  mode = 'a'
  nBins = 200
  rdf = np.zeros(nBins)
  zc = 5.00
  binSize = zc/float(nBins)
  z = (np.arange(nBins-1)+0.5)*binSize

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)

      polymers = frame[frame[:,0] == 1][:,1:]
      box = np.array([polymers.min(0), polymers.max(0)])
      L = box[1,:] - box[0,:]
      for i in range(len(polymers) - 1):
        delta = polymers[i,:] - polymers[i+1:,:]
        delta[:,:2] = delta[:,:2] - L[:2]*np.floor(delta[:,:2]/L[:2])
        delta = np.sqrt(np.einsum('ij,ij->i',delta,delta))
        bin = np.floor(delta[delta < zc]/binSize)
        bin = bin.astype(int)
        rdf += 2*np.bincount(bin,minlength=200)
      
      # "normalize" the histogram
      rho = float(len(polymers)) / (L[0]*L[1]*L[2])
      k = np.arange(1,nBins)
      vb = ((k + 1)**3 - k**3)*binSize**3
      nid = (4.0/3)*np.pi*vb*rho
      rdf[1:] /= nid*len(polymers)

      output = np.zeros((nBins-1,2))
      output[:,0] = z
      output[:,1] = rdf[1:]
        
      plt.plot(z, rdf[1:])
      # plt.plot([0,zc], [0.2,0.2])
      plt.show()
      # writeRDF(rdfFile,output,t,mode)
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
