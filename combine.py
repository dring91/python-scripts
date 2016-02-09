#!/usr/bin/python

from sys import argv, exit
from getopt import *
from writeAtoms import *
import numpy as np

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'p:s:o:')
  except GetoptError:
    print 'randomWalk.py -p <prefix file> -s <suffix file> -o <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-p':
      prefixFile = arg
    if opt == '-s':
      suffixFile = arg
    elif opt == '-o':
      outFile = arg
 
  return prefixFile, suffixFile, outFile

def formatOutput(atoms, bonds, types):
  # No modification
  polymers = [atom[:2] + [str(int(atom[2]) + 2)] + atom[3:] for atom in atoms[:types[0]]]
  # Increments atom number and molecule number
  surface  = [[str(int(atom[0])+int(polymers[-1][0])-int(atoms[types[0]][0])+1),
               str(int(atom[1])+int(polymers[-1][1])-int(atoms[types[0]][1])+1)
              ] + atom[2:] for atom in atoms[types[0]:]]
  bonds = [bond[:2] + [str(int(bond[2])+int(polymers[-1][0])-int(atoms[types[0]][0])+1),
                       str(int(bond[3])+int(polymers[-1][0])-int(atoms[types[0]][0])+1)]
                       for bond in bonds] 

  return polymers + surface, bonds

def makeBox(box1,box2):
  box = box1[:2]
  box.append([0.0,0.0])
  if box1[2][0] > box2[2][0]:
    box[2][0] = box2[2][0]
  else:
    box[2][0] = box1[2][0]
  if box1[2][1] < box2[2][1]:
    box[2][1] = box2[2][1]
  else:
    box[2][1] = box1[2][1]

  return box

def translateAtoms(atoms, offset):
  # select coordinates and convert to float
  coords = [[float(i) for i in atom[3:]] for atom in atoms]
  # add coordinates and offset
  coords = [[i + j for (i,j) in zip(line,offset)] for line in coords]
  # convert coordinates back to strings
  coords = [[str(i) for i in line] for line in coords]
  # recombine the lists
  return [atom[:3] + coord for (atom,coord) in zip(atoms,coords)]

def boundAtoms(atoms):
  coords = [[float(i) for i in atom[3:]] for atom in atoms]
  bounds = [coords[0][2],coords[0][2]]
  for coord in coords:
    if coord[2] < bounds[0]:
      bounds[0] = coord[2]
    elif coord[2] > bounds[1]:
      bounds[1] = coord[2]

  return bounds

def readConf(file, atype):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in ['Atoms','Bonds']:
      header = L[0]
    if len(L) > 0 and L[-1] in ['xhi','yhi','zhi']:
      box.append(L[:2])
    elif len(L) > 2 and L[2] in atype and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(L)
      
  return box, atoms, bonds

def main():

  inFile1, inFile2, outFile = getArgs(argv)

  nParticles = 63

  atomSet = ['3']
  with open(inFile2, 'r') as inp2:
    box2, topAtoms, bonds = readConf(inp2, atomSet)

  with open(inFile1, 'r') as inp1:
    box1, bottomAtoms, bonds = readConf(inp1,['2','1'])

  topAtoms = [[int(row[0])] + row[1:] for row in topAtoms] 
  bottomAtoms = [[int(row[0])] + row[1:] for row in bottomAtoms] 

  topAtoms.sort()
  bottomAtoms.sort()

  topAtoms = [[str(row[0])] + row[1:] for row in topAtoms] 
  bottomAtoms = [[str(row[0])] + row[1:] for row in bottomAtoms] 

  topBounds = boundAtoms(topAtoms)
  bottomBounds = boundAtoms(bottomAtoms)
  buffer = -1

  shift = topBounds[0]-bottomBounds[1]+buffer
  bottomAtoms = translateAtoms(bottomAtoms, [0,0,shift]) 

  atoms = topAtoms + bottomAtoms
  types = [len(topAtoms),len(bottomAtoms)]

  box = box1[:2] + [[str(bottomBounds[0]+shift),str(topBounds[1])]]

  write_xyz(outFile, [line[2:] for line in atoms])
  # write_conf(filename,atoms,bonds,title,types,box,masses):
  atoms, bonds = formatOutput(atoms, bonds, types)

  # # check that the bonds and atoms correspond properly
  # test_atoms = np.array(atoms)
  # test_atoms = set(test_atoms[test_atoms[:,2] == '1'][:,0])
  # test_bonds = np.array(bonds)
  # aSet = set(test_bonds[:,2])
  # bSet = set(test_bonds[:,3])
  # bSet = bSet | aSet
  # test = test_atoms - bSet
  # print test

  write_conf(outFile, 
             atoms,
             bonds,
             'saturated infiltration simulation with cylindrical capillary', 
             [3,1], 
             box,
             [1] * 3
            )

if __name__ == "__main__":
  main()
