#!/usr/bin/python

from sys import argv, exit
from getopt import *
from conf_tools import *
import numpy as np

def inRange(point, interval=[0,1], left=True, right=False):
  
  if point > interval[0] and point < interval[1]:
    return True
  elif left and point == interval[0]:
    return True
  elif right and point == interval[1]:
    return True
  else:
    return False
      
def inBox(point, box=[[0,1],[0,1],[0,1]], left=True, right=False):
  
  for i in range(len(point)):
    if inRange(point[i],box[i]) == False:
      return False
  
  return True

def inSphere(point, R=[0,1], inner=True, outer=True, hemisphere=False):
  
  r = np.sqrt(np.dot(point,point))
  if hemisphere and inRange(r,R,inner,outer) and point[2] > 0:
    return True
  elif inRange(r,R,inner,outer):
    return True
  
  return False

def inCylinder(point, R, L, origin=[0,0,0], outer=True, inner=True):

  R = [0.0,R]
  L = [-0.5*L,0.5*L]
  point = np.array(point)
  origin = np.array(origin)
  r = np.sqrt(np.dot(point[1:]-origin[1:],point[1:]-origin[1:]))
  if inRange(r,R,inner,outer) and inRange(point[0]-origin[0],L):
    return True

  return False

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:l:r:z:')
  except GetoptError:
    print 'randomWalk.py -i <file> -o <file> -l <chainLength> -r <radius> -z <z-origin>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-l':
      chainLength = arg
    elif opt == '-r':
      radius = arg
    if opt == '-z':
      origin = arg
 
  return inFile, outFile, int(chainLength), float(radius), float(origin)

def cut(atoms, R, origin, chainLength):
  n = 0
  molecule = []
  drop = []
  box = [[comp-0.5*R,comp+0.5*R] for comp in origin]
  # loop over each molecule
  for atom in atoms:
    if n == chainLength:
      # unnecessary because of unwrap
      # calculate its center of mass
      molecule = np.array(molecule)
      COM = [calcCOM(molecule[:,d]) for d in range(len(molecule[0]))]
      # check whether it's in the region
      if inBox(COM, box): 
        # add it to drop if it is
        drop.extend(molecule)
      n = 0
      molecule = []
    n += 1
    molecule.append(atom)

  return drop

def main():

  inFile, outFile, chainLength, R, z = getArgs(argv)

  with open(inFile+'.conf', 'r') as inp:
    box, atoms, bonds = readConf(inp,['1'])

  atoms = np.array(atoms,dtype=float)
  conf = cut(atoms, R, [0,0,z], chainLength)

  # make new box
  box = [boundAtoms(conf,d) for d in range(3)]

  write_xyz(outFile, conf, len(conf)/chainLength, chainLength)
  # write_conf(filename,atoms,bonds,title,types,box,masses):
  write_conf(outFile, 
             conf,
             len(conf)/chainLength,
             'undersaturated infiltration simulation\n', 
             chainLength,
             [1,1], 
             box,
             [1]
            )

if __name__ == "__main__":
  main()
