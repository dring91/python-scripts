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
  
  return atoms, bonds, np.array(box,dtype=float)

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
      unshuffled.append(atoms[atom]) 
      for n in range(nMon-1):
        atom = bonds.pop(atom)
        unshuffled.append(atoms.pop(atom)) 
    elif tail in bonds.values():
      # print 'Value in dict reversed'
      continue

  return unshuffled

def main():
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-l','--chainlength',type=int)
  args = parser.parse_args()

  with open(args.input+'.conf', 'r') as file:
    atoms, bonds, box = GetAtoms(file)

  # search for all the chain ends in bonds and store in new array
  ends = getEnds(atoms, bonds)

  # unshuffle bonds
  unshuffled = unshuffleAtoms(copy(atoms), copy(bonds), ends, args.chainlength)
  atoms = np.array([[str(i+1)]+atom for i,atom in enumerate(unshuffled)] + \
          [[key]+atom for key,atom in atoms.items() if atom[1] == '2']).astype(float)
  bonds = [[str(i+1),'1',key,val] for i,(key,val) in enumerate(bonds.items())]

  atoms[:,3:6] = atoms[:,3:6] + atoms[:,6:9] * (box[:,1] - box[:,0])
  bonded = atoms[atoms[:,2] == 1][:,3:6].reshape((-1,args.chainlength,3))
  bonded = np.diff(bonded, axis=1)
  print(np.sqrt(bonded[:,:,0]**2 + bonded[:,:,1]**2 + bonded[:,:,2]**2).mean(1))

  # write the new atoms list and the bonds to a file
  data = (atoms, bonds, box)
  options = {'masses':[1]*2,'types':{'atoms':2,'bonds':1},'title':'unshuffled configuration'}
  write_conf(args.input+'_unshuffled', *data, **options)
  write_traj(args.input+'_unshuffled', np.delete(atoms[:,:6],1,1), np.array(box).astype(float))

if __name__ == '__main__':
  main()
