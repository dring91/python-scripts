#!/usr/bin/python

from sys import argv, exit
from getopt import *
from writeAtoms import *
from nestedDicts import *
import numpy as np

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:')
  except GetoptError:
    print 'randomWalk.py -i <filename> -o <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
 
  return inFile, outFile
  
def readConf(file, atype):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in set(['Atoms','Bonds']):
      header = L[0]
    if len(L) > 0 and L[-1] in set(['xhi','yhi','zhi']):
      box.append(L[:2])
    elif len(L) > 2 and L[2] in set(atype) and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(L)
      
  return box, atoms, bonds
  
def PBC(dx,x,interval):
  if dx <= interval[0]:
    x = x + 2*interval[1]
  elif dx > interval[1]:
    x = x + 2*interval[0]

  return x  
  
def unwrap(info, coords, box):
  unwrapped = []
  box = [[float(i) for i in dim] for dim in box]
  # loop over each atomtype
  for atype,mols in info.iteritems():
    # loop over each molecule in an atomtype
    for mol,atoms in mols.iteritems():
      # for each atom, apply PBC(atom-atom0)
      for atom in atoms:
        line = [atom, mol, atype]
        line += [str(PBC(float(coords[atom][d])-float(coords[atoms[0]][d]),
                         float(coords[atom][d]),
                         box[d]
                        ) 
                    )
                 for d in range(2)
                ] + [coords[atom][2]]
        unwrapped.append(line)
      
  return unwrapped
 
def main():
  
  inFile, outFile = getArgs(argv)
  
  nAtoms = 63
  atomList = np.arange(63) + 1
  with open(inFile,'r') as file:
    box, atoms, bonds = readConf(file, atomList.astype(str))
  
  # substrate = [line for line in atoms if line[2] == '2']
  # atoms = [line for line in atoms if line[2] in set(['1','3'])]
  
  info,coords = makeNestedDict(atoms)

  atoms = unwrap(info,coords,box) # + substrate

  title = 'unwrapped configuration'
  types = [1,0] # [3,1]
  masses = [1] # [1,1,1]  
  write_conf(outFile,atoms,bonds,title,types,box,masses)
  write_xyz(outFile,[line[2:] for line in atoms])  
  
if __name__ == "__main__":
  main()
