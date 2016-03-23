#!/usr/bin/python

import sys
from conf_tools import *
from pbc_tools import *
import numpy as np

MODE = 'r'
X = 0
Y = 1
Z = 2

def renumber(atoms, chainLength, sequential=False):

  # Renumber atoms
  atoms[:,0] = np.arange(len(atoms)) + 1
  # number of polymer units
  nPoly = len(atoms[atoms[:,2] == 1])

  if sequential:
    # sequential numbering
    atoms[:nPoly,1] = np.repeat(np.arange(nPoly/chainLength) + 1, chainLength)
    atoms[nPoly:,1] = nPoly/chainLength + 1
  else:
    # numbering for bond/swap
    atoms[:nPoly,1] = (atoms[:nPoly,0] - 1) % chainLength
    atoms[:nPoly,1] = (atoms[:nPoly,1] - chainLength/2) 
    mask = (atoms[:nPoly,1] >= 0)
    atoms[:nPoly,1] = (atoms[:nPoly,1] + 1) * mask + abs(atoms[:nPoly,1]) * np.logical_not(mask)
    atoms[nPoly:,1] = chainLength + 1

  return atoms

def minImageConv(atoms, imageFlags, box):

  atoms = atoms + (box[:,1] - box[:,0]) * imageFlags

  return atoms

def calcCOM(coords, nMon):

  coords = coords.reshape((-1,nMon,3)).mean(1)
  return np.repeat(coords, nMon, axis=0)

def main():
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  atype = ['1','2','3']
  nMon = 50

  with open(filename+'.conf', MODE) as inp:
    box, atoms, bonds = readConf(inp, atype)

  box = np.array(box, dtype=float)
  atoms = np.array(atoms)
  bonds = np.array(bonds, dtype=int)

  # sort atoms
  info = atoms[:,:3].astype(int)
  order = info[:,0].argsort(kind='mergesort')
  atoms = atoms[order]
  info = info[order]
  
  coords = atoms[:,3:6].astype(float)
  flags = atoms[:,6:].astype(int)

  # Use image flags written by Lammps
  coords = minImageConv(coords, flags, box)
  atoms[:,3:6] = coords.astype('|S6')

  # # unwrap atom coordinates in the absence of image flags
  # coords = unwrap(coords, box)
  # atoms[:,3:6] = coords.astype('|S6')

  # Adjust box configuration and discard molecules entirely outside the box
  box[X,:] = box[X,:] / 2
  com = calcCOM(coords[info[:,2] == 1], nMon)
  inside = np.copy(coords[:,X])
  inside[:len(com)] = np.logical_and(com[:,X] <= box[X,1], com[:,X] > box[X,0])
  inside[len(com):] = np.logical_and(inside[len(com):] <= box[X,1], inside[len(com):] > box[X,0])
  inside = inside.astype(bool)
  atoms = atoms[inside]
  info = info[inside]

  # Tile the current configuration
  # double the size of the atom array
  # atoms = np.concatenate((atoms,atoms),axis=0)

  # renumber atoms
  info = renumber(info, nMon)
  atoms[:,:3] = info.astype(str)

  # Generate new bonds
  bonds = np.ones((sum(inside[:len(com)]) * (nMon - 1) / nMon,4), dtype=int)
  bonds[:,0] = np.arange(len(bonds)) + 1
  mask = (info[info[:,2] == 1][:,0] % nMon).astype(bool)
  bonds[:,2] = info[mask][:,0]
  bonds[:,3] = info[mask][:,0] + 1

  # convert data to string format for output
  box = box.astype('|S8')
  bonds = bonds.astype(str)

  # write output files
  suffix = '_halved'
  write_conf(filename+suffix, atoms[:,:6], bonds, 
            'modified configuration', [3,1], box, [1,1,1])
  write_xyz(filename+suffix, atoms[:,2:6])

if __name__ == '__main__':
  main()
