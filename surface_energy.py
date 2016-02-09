#!/usr/python

from sys import argv
import numpy as np
from getopt import *
from pbc_tools import *

""" Use Young's equation to get at surface energies
     - Sum interactions that correspond s-v, s-l, and l-v
     - Calculate the contact angle using Young's equation
     - Compare to the theta values found independently
"""
def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print """randomWalk.py -t <filename> -o <filename> -n <steps> """
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      steps = arg

  return trajFile, outFile, int(steps)
 
def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def interactions(potential, r_cut, box, *args):
  E = 0
  dist = lambda v: np.sqrt(v[0]**2+v[1]**2+v[2]**2)
  if len(args) == 0:
    print "No atoms"
  elif len(args) == 1:
    rows, _ = args[0].shape
    # this doesn't actually work
    for i in range(rows-1):
      for j in range(i+1,rows):
        j[:2] = PBC(i[:2]-j[:2],j[:2],box[:2,:])
        r = dist(i-j)
        if r < r_cut:
          E += potential(r)
          print E
  elif len(args) >= 2:
    # only use first two arguments from args
    # pairs = np.tile(args[0],(len(args[1]),1,1))
    # pairs = np.swapaxes(pairs,0,1)
    # pairs = pairs - args[1]
    # rs = np.sqrt(pairs[:,:,0]**2+pairs[:,:,1]**2+pairs[:,:,2]**2)
    # mask = (rs < r_cut)
    # E = potential(rs[mask]).sum()

    for i in args[0]:
      for j in args[1]:
        j[:2] = PBC(i[:2]-j[:2],j[:2],box[:2,:])
        r = dist(i-j)
        if r < r_cut:
          E += potential(r)
          print E

  return E

def main():

  inFile, outFile, nFrames = getArgs(argv)

  # define a function that calculates interactions based on LJ potential with cut-offs
  # Start out with polymers that don't have vapor
  #  thus:  cos (theta) = - gamma_sl
  # Make the interactions function compare interactions between two lists
  #  If the lists are identical, pass just one (use *args to handle this)

  NDIMS = 3
  epsilon = 1
  sigma = 1
  r_cut = 1.75
  cutOffEnergy = 4*epsilon*((sigma/r_cut)**12-(sigma/r_cut)**6)
  U_LJ = lambda r: 4*epsilon*((sigma/r)**12-(sigma/r)**12) - cutOffEnergy

  with open(outFile, 'w') as otp:
    otp.write('# Youngs law calculation of theta and surface energy\n')

  with open(inFile, 'r') as inp:
    # For each frame
    for n in range(nFrames):
      time, atoms = readFrame(inp)

      # filter out surface and polymer
      polymer = atoms[atoms[:,0] == 1][:,1:]
      surface = atoms[atoms[:,0] == 2][:,1:]
      box = makeBox(NDIMS, surface, polymer)
      
      # Calculate gamma_sl
      gamma_sl = interactions(U_LJ, r_cut, box, polymer, surface)
      
      # Calculate theta using Young's equation
      Youngs = lambda x: np.acos(-x)
      theta = Youngs(gamma_sl)

      with open(outFile, 'a') as otp:
        np.savetxt(otp, [[time, gamma_sl, theta]], fmt='%d %.5f %.5f\n')
    
if __name__ == '__main__':
  main()
