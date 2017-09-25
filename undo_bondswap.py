import numpy as np
from sys import argv, exit
import os
from getopt import *
from conf_tools import write_conf, write_traj
from copy import copy

from conf_tools import readConf
from argparse import ArgumentParser
import matplotlib.pyplot as plt

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
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    types, box, atoms, bonds = readConf(file)

  #print(np.mean(bonds[:,2] - bonds[:,3]))
  xs, ys, zs = atoms[:,3], atoms[:,4], atoms[:,5] # - atoms[:,6:9] * (box[:,1] - box[:,0])
  diffs = np.stack((xs[bonds[:,2]] - xs[bonds[:,3]],
                    ys[bonds[:,2]] - ys[bonds[:,4]],
                    zs[bonds[:,2]] - zs[bonds[:,5]]), axis=0)
  dists = np.sqrt(diffs[:,0]**2+diffs[:,1]**2+diffs[:,2]**2)
  fx, x = np.histogram(np.abs(diffs[:,0]), bins='auto')
  fy, y = np.histogram(np.abs(diffs[:,1]), bins='auto')
  fz, z = np.histogram(np.abs(diffs[:,2]), bins='auto')
  fr, r = np.histogram(dists, bins='auto')
  
  fig, axis = plt.subplots()
  axis.plot(x[:-1]+np.diff(x)/2,fx,label='x') 
  axis.plot(y[:-1]+np.diff(y)/2,fy,label='y') 
  axis.plot(z[:-1]+np.diff(z)/2,fz,label='z') 
  axis.plot(r[:-1]+np.diff(r)/2,fr,label='r') 
  axis.legend()
  plt.show()
     
if __name__ == '__main__':
  main()
