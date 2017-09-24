import numpy as np
from sys import argv, exit
import os
from getopt import *
from conf_tools import write_conf, write_traj
from copy import copy

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
    elif len(L) > 2 and L[2] in ['1','2'] and section == 'Atoms':
      try:
        atoms.update({L[0]:L[1:9]})
      except IndexError:
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
  for tail in ends:
    if tail in bonds.keys():
      atom = tail
      unshuffled.append(atoms[atom]) 
      for n in range(nMon-1):
        atom = bonds.pop(atom)
        unshuffled.append(atoms.pop(atom)) 
    elif tail in bonds.values():
      # print 'Value in dict reversed'
      continue

  return unshuffled

def main():
  
  inFile, nMon = getArgs(argv)

  with open(inFile+'.conf', 'r') as inp:
    atoms, bonds, box = GetAtoms(inp)

  # search for all the chain ends in bonds and store in new array
  ends = getEnds(atoms, bonds)

  # unshuffle bonds
  unshuffled = unshuffleAtoms(copy(atoms), copy(bonds), ends, nMon)
  atoms = np.array([[str(i+1)]+atom for i,atom in enumerate(unshuffled)] + \
          [[key]+atom for key,atom in atoms.iteritems() if atom[1] == '2']).astype(float)
  bonds = [[str(i+1),'1',key,val] for i,(key,val) in enumerate(bonds.iteritems())]

  bonded = atoms[atoms[:,2] == 1][:,3:].reshape((-1,nMon,3))
  bonded = np.diff(bonded, axis=1)
  print(np.sqrt(bonded[:,:,0]**2 + bonded[:,:,1]**2 + bonded[:,:,2]**2).mean(1))

  # write the new atoms list and the bonds to a file
  data = (atoms, bonds, box)
  options = {'masses':[1]*2,'types':{'atoms':2,'bonds':1},'title':'unshuffled configuration'}
  write_conf(inFile+'_unshuffled', *data, **options)
  write_traj(inFile+'_unshuffled', np.delete(atoms,1,1), np.array(box).astype(float))

  print 'Finished Analyzing Trajectory'
      
if __name__ == '__main__':
  main()
