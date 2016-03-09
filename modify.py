#!/usr/bin/python

import sys
from conf_tools import *
from pbc_tools import *
import numpy as np

MODE = 'r'

def renumber(atoms, chainLength):

  atoms[:,1] = (atoms[:,0] - 1) % chainLength
  atoms[:,1] = (atoms[:,1] - 25) 
  mask = (atoms[:,1] >= 0)
  atoms[:,1] = (atoms[:,1] + 1) * mask + abs(atoms[:,1]) * np.logical_not(mask)

  return atoms

def minImageConv(atoms, imageFlags, box):

  atoms = atoms + (box[:,1] - box[:,0]) * imageFlags

  return atoms

def main():
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  atype = ['1','2','3']
  nMon = 10

  with open(filename+'.conf', MODE) as inp:
    box, atoms, bonds = readConf(inp, atype)

  atoms = np.array(atoms)
  box = np.array(box)

  # # Use image flags written by Lammps
  # coords = atoms[:,3:6].astype(float)
  # flags = atoms[:,6:].astype(int)
  # coords = minImageConv(coords, flags, box.astype(float))
  # atoms[:,3:6] = coords.astype('|S6')

  # # unwrap atom coordinates
  # coords = atoms[:,3:6].astype(float)
  # coords = unwrap(coords, box.astype(float))
  # atoms[:,3:6] = coords.astype('|S6')

  # # renumber atoms
  # info = atoms[:,:3].astype(int)
  # info = renumber(info, nMon)
  # atoms[:,:3] = info.astype(str)

  # sort atoms
  nums = atoms[:,0].astype(int)
  atoms = atoms[nums.argsort(kind='mergesort')]

  write_conf(filename+'_sorted', atoms, bonds, 
            'atoms renumbered for bond/swap', [3,1], box, [1,1,1])
  write_xyz(filename+'_out', atoms[:,2:6])

if __name__ == '__main__':
  main()
