import numpy as np
from sys import argv, exit
import os
from getopt import *
from conf_tools import write_conf, write_traj
from copy import copy
from argparse import ArgumentParser

def PBC(x,boundaries):
  
  if x < boundaries[0]:
    x = x + 2*boundaries[1]
  elif x >= boundaries[1]:
    x = x + 2*boundaries[0]

  return x

def GetAtoms(file, moleculetype):

  atoms = []
  molecules = {}
  bonds = {}
  box = []
  types = {}
  section = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and (L[0] == 'Atoms' or L[0] == 'Velocities' or L[0] == 'Bonds'):
      section = L[0]
    elif len(L) > 3 and (L[3] == 'xhi' or L[3] == 'yhi' or L[3] == 'zhi'):
      box.append(L[:2])
    elif len(L) > 2 and L[2] == 'types':
      types.update({L[1]:int(L[0])})
    elif len(L) > 2 and L[1] == moleculetype and section == 'Atoms':
      molecules.update({L[0]:L[1:]})
    elif len(L) > 2 and L[1] != moleculetype and section == 'Atoms':
      atoms.append(L)
    elif len(L) > 1 and L[1] == '1' and section == 'Bonds':
      bonds.update({L[2]:L[3]})
  
  return np.array(atoms), molecules, bonds, np.array(box,dtype=float), types

def getEnds(atoms, bonds):
  # find them by counting up the number of times an atom occurs in bonds
  atom_count = dict((key,0) for (key,value) in atoms.items())
  bond_atoms = [item for (key, value) in bonds.items() 
                     for item in (key, value)]
  for atom in bond_atoms:
    atom_count[atom] += 1
  ends = [key for (key,value) in atom_count.items() if value == 1]
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
      unshuffled.append([atom] + atoms[atom]) 
      for n in range(nMon-1):
        atom = bonds.pop(atom)
        unshuffled.append([atom] + atoms.pop(atom)) 
    elif tail in bonds.values():
      # print 'Value in dict reversed'
      continue

  return np.array(unshuffled)

def main():
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-l','--chainlength',type=int)
  parser.add_argument('--polymer', type=int)
  args = parser.parse_args()

  with open(args.input+'.conf', 'r') as file:
    atoms, molecules, bonds, box, types = GetAtoms(file, args.polymer)

  # search for all the chain ends in bonds and store in new array
  ends = getEnds(molecules, bonds)

  # unshuffle bonds
  unshuffled = unshuffleAtoms(copy(molecules), copy(bonds), ends, args.chainlength)
  atoms = np.concatenate((unshuffled,atoms),axis=0) 
  bonds = [[str(i+1),'1',key,val] for i,(key,val) in enumerate(bonds.items())]

  atoms[:,3:6] = atoms[:,3:6] + atoms[:,6:9] * (box[:,1] - box[:,0])
  bonded = atoms[atoms[:,2] == args.polymer][:,3:6].reshape((-1,args.chainlength,3))
  bonded = np.diff(bonded, axis=1)
  print(np.sqrt(bonded[:,:,0]**2 + bonded[:,:,1]**2 + bonded[:,:,2]**2).mean(1))

  # write the new atoms list and the bonds to a file
  atoms[:,1] += 1
  data = (atoms, bonds, box)
  options = {'masses':[1]*types['atom'],'types':types,'title':'unshuffled configuration'}
  write_conf(args.input+'_unshuffled', *data, **options)
  write_traj(args.input+'_unshuffled', np.delete(atoms[:,:6],1,1), np.array(box).astype(float))

if __name__ == '__main__':
  main()
