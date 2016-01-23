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
  
def writeChains(filename,chains,time,mode):
  with open(filename, mode) as file:
    [file.write('%d %f %f\n' % (time, com[0], com[1])) for com in chains]
    file.write('\n\n')

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

  rhoFile = path+'chains_'+outFile
  
  title = '# polymer chain trajectories\n'+'# time  COM_z   RG\n' #'# time  COM_z[0] ... COM_z[n]\n'
  writeHeader(rhoFile,title)

  mode = 'a'

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)
      # write atom number and Timestep
      # write NPs to file

      polymers = frame[frame[:,0] == 1][:,3]
      polymers = polymers.reshape((-1,10))

      COM = polymers.mean(1)
      RG = (polymers.T - COM)**2
      RG = np.sqrt(RG.mean(0))

      # Filter output using the packing base
      particles = frame[frame[:,0] == 3][:,3]
      pmin = particles.min()
      pmax = particles.max()
      COM = COM - pmin
      COM = COM[COM >= 0] / pmax
      RG = RG[COM >= 0] / pmax

      output = np.zeros((len(COM),2))
      output[:,0] = COM
      output[:,1] = RG
       
      # plt.plot(np.ones_like(polymers),polymers,'o')
      # plt.show()
      if len(output) != 0:
        # write COM, RG and surface atoms
        writeChains(rhoFile, output, t, mode)
        print('%d Chains written at timestep: %d' % (len(output),t))
      else:
        print 'No observed infiltration'
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
