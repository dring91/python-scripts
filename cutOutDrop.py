#!/usr/bin/python

from sys import argv, exit
from getopt import *
import numpy as np
from pbc_tools import *
from conf_tools import *

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

def PBC(dx,x,interval):
  if dx <= interval[0]:
    x = x + 2*interval[1]
  elif dx > interval[1]:
    x = x + 2*interval[0]

  return x  

def calcCOM(atoms):

  COM = 0
  for atom in atoms:
    COM += atom
  COM /= len(atoms)

  return COM

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

def boundAtoms(coords, dim):
  bounds = [coords[0][dim],coords[0][dim]]
  for coord in coords:
    if coord[dim] < bounds[0]:
      bounds[0] = coord[dim]
    elif coord[dim] > bounds[1]:
      bounds[1] = coord[dim]

  return bounds

def cut(atoms, R, origin, chainLength, box):
  n = 0
  molecule = []
  drop = []
  newBox = [[-0.5*R+comp,0.5*R+comp] for comp in origin]
  # loop over each molecule
  for atom in atoms:
    if n == chainLength:
      # unnecessary because of unwrap
      molecule = np.array(molecule)
      molecule = unwrap(molecule, box)
      # calculate its center of mass
      COM = [calcCOM(molecule[:,d]) for d in range(len(molecule[0]))]
      # check whether it's in the region
      if inBox(COM, newBox): # if inRange(COM[0],[-0.5*R,0.5*R]): # inCylinder(COM, R, 0.5*R, origin):
        # add it to drop if it is
        drop.extend(molecule)
      n = 0
      molecule = []
    n += 1
    molecule.append(atom)

  return drop

def main():

  inFile, outFile, chainLength, R, z = getArgs(argv)

  # try:
  #   with open(inFile+'.conf', 'r') as inp:
  #     box, atoms, bonds = readConf(inp,['1'])
  # except IOError:
  #   print 'Cannot read input file'
  #   exit(0)
  with open(inFile+'.conf', 'r') as inp:
    box, atoms, bonds = readConf(inp,['1'])

  atoms = [[float(i) for i in atom[3:]] for atom in atoms]
  box = np.array(box,dtype=float)
  drop = cut(atoms, R, [0,0,z], chainLength, box)

  # make new box
  box = [boundAtoms(drop,d) for d in range(3)]

  # write_xyz(filename, atoms, {time, mode})
  write_xyz(outFile, drop, len(drop)/chainLength, chainLength)
  # write_conf(filename, atoms, {bonds,title,types,box,masses}):
  write_conf(outFile, 
             drop,
             len(drop)/chainLength,
             'undersaturated infiltration simulation\n', 
             chainLength,
             [1,1], 
             box,
             [1]
            )

if __name__ == "__main__":
  main()
