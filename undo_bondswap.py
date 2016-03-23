import numpy as np
from sys import argv, exit
import os
from getopt import *
from conf_tools import write_conf, write_xyz

def PBC(x,boundaries):
  
  if x < boundaries[0]:
    x = x + 2*boundaries[1]
  elif x >= boundaries[1]:
    x = x + 2*boundaries[0]

  return x

def GetAtoms(file):

  atoms = {}
  bonds = {}
  box = []
  section = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and (L[0] == 'Atoms' or L[0] == 'Velocities' or L[0] == 'Bonds'):
      section = L[0]
    elif len(L) > 3 and (L[3] == 'xhi' or L[3] == 'yhi' or L[3] == 'zhi'):
      box.append(L[:2])
    elif len(L) > 2 and L[2] == '1' and section == 'Atoms':
      atoms.update({L[0]:L[1:6]})
    elif len(L) > 1 and L[1] == '1' and section == 'Bonds':
      bonds.update({L[2]:L[3]})
  
  return atoms, bonds, box

def getArgs(argv):
 
  try:
    opts, args = getopt(argv[1:],'i:l:')
  except GetoptError:
    print 'undo_bondswap.py -i <filename> -l <chainLength>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-l':
      chainLength = arg
 
  return inFile, int(chainLength)

def getEnds(atoms, bonds):
  # find them by counting up the number of times an atom occurs in bonds
  atom_count = dict((key,0) for (key,value) in atoms.iteritems())
  bond_atoms = [item for (key, value) in bonds.iteritems() 
                     for item in (key, value)]
  for atom in bond_atoms:
    atom_count[atom] += 1
  ends = [key for (key,value) in atom_count.iteritems() if value == 1]
  return ends

def unshuffleAtoms(atoms, bonds, ends, nMon):
  # loop over bond ends array
  # build a new list of atoms by picking a chain end and adding atoms
    # when an atom is added, remove it from the atoms list
  # if chainLength is reached and the last atom is an end, you're good
  # repeat until all bond ends are used up

  unshuffled = []
  # unshuffle all the 'forwards' chains and then flip the 'backwards' ones and unshuffle the remaining
  for tail in ends[0:5]:
    if tail in unshuffled:
      continue
    if tail in bonds.keys:
      atom = tail
      unshuffled.append(atoms[atom]) 
      for n in range(nMon):
        atom = bonds.pop(atom)
        unshuffled.append(atoms.pop(atom)) 
      print len(unshuffled)
    elif tail in bonds.values:
      print 'Value in dict reversed'

  return unshuffled

def main():
  
  inFile, nMon = getArgs(argv)

  with open(inFile+'.conf', 'r') as inp:
    atoms, bonds, box = GetAtoms(inp)

  # search for all the chain ends in bonds and store in new array
  ends = getEnds(atoms, bonds)
  print ends

  # # unshuffle bonds
  # unshuffled = unshuffleAtoms(atoms, bonds, ends, nMon)

  # # write the new atoms list and the bonds to a file
  # suffix = '_unshuffled'
  # write_conf(inFile+suffix, atoms, bonds, 'unshuffled conf', [2,1], box, [1,1])
  # write_xyz(inFile+suffix, atoms)

  print 'Finished Analyzing Trajectory'
      
if __name__ == '__main__':
  main()
