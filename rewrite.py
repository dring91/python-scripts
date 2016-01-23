#!/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from pbc_tools import *

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    if opt == '-o':
      outFile = arg
    if opt == '-n':
      nFrames = arg

  return trajFile, outFile, int(nFrames)
 
def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms
  
def main():

  trajFile, outFile, nFrames = getArgs(argv)
  atomtype = 1
  nSkip = nFrames / 100

  with open(outFile+'_unwrapped.xyz', 'w') as otp:
    otp.write('')

  with open(trajFile+'.xyz','r') as inp:
    for n in range(nFrames/nSkip):
      # read and filter data
      for s in range(nSkip):
        time, frame = readFrame(inp)
      frame = frame[frame[:,0] == atomtype][:,1:]
      # build box
      box = np.zeros((2,2))
      box[:,0] = frame[:,:2].min(0)
      box[:,1] = frame[:,:2].max(0)
      # subtract prior frame
      if n == 0:
        periods = np.zeros_like(frame[:,:2])
        move = np.zeros_like(frame[:,:2])
      elif n > 0:
        move = frame[:,:2] + periods * (box[:,1] - box[:,0]) - oldFrame[:,:2]
        periods += (move < box[:,0])
        periods -= (move >= box[:,1])
        frame[:,:2] = frame[:,:2] + periods * (box[:,1] - box[:,0])
      # if differences are greater than L/2 in any frame, add 1
      oldFrame = np.copy(frame)

      with open(outFile+'_unwrapped.xyz', 'a') as otp:
        otp.write('%d\n' % len(frame))
        otp.write('Atoms. Timestep: %d\n' % time)
        np.savetxt(otp,frame,'1 %.5f %.5f %.5f')
     
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
