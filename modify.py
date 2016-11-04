#!/usr/bin/python

import sys
from conf_tools import *
from pbc_tools import *
import numpy as np

MODE = 'r'
X = 0
Y = 1
Z = 2

def renumber(atoms, types, nMon, diffTypes=False, sequential=False):

  (nPoly, nSurf, nPart) = types

  # Renumber atoms
  atoms[:,0] = np.arange(len(atoms)) + 1
  if diffTypes is True:
    atoms[:nPoly,2] = np.arange(nPoly) + 1
    atoms[nPoly:,2] = nPoly + 1
  
  if sequential:
    # sequential numbering
    atoms[:nPart,1] = np.repeat(np.arange(nPart/4684) + 1, 4684)
    atoms[nPart:nPoly,1] = np.repeat(np.arange(nPoly/nMon) + 1, nMon)
    atoms[nPoly:,1] = nPoly/nMon + 1
  else:
    # # numbering for bond/swap
    # atoms[:nPoly,1] = (atoms[:nPoly,0] - 1) % nMon
    # atoms[:nPoly,1] = (atoms[:nPoly,1] - nMon/2) 
    # mask = (atoms[:nPoly,1] >= 0)
    # atoms[:nPoly,1] = (atoms[:nPoly,1] + 1) * mask + abs(atoms[:nPoly,1]) * np.logical_not(mask)
    # atoms[nPoly:,1] = nMon/2 + 1
    # new numbering
    if nMon % 2 == 1:
      numbers = range(1, nMon / 2 + 1) + range(nMon / 2 + 1, 0, -1)
      atoms[nPoly:,1] = nMon/2 + 2
    elif nMon % 2 == 0:
      numbers = range(1, nMon / 2 + 1) + range(nMon / 2, 0, -1)
      atoms[nPoly:,1] = nMon/2 + 1
    print numbers
    numbers = numbers * (nPoly / nMon)
    atoms[:nPoly,1] = numbers

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

  atype = '12'
  nMon = 50
  # nChains = 1000

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

  # pick a certain number of chains
  # atoms = atoms[:nMon*nChains]
  # info = info[:nMon*nChains]
  
  coords = atoms[:,3:6].astype(float)
  flags = atoms[:,6:].astype(int)

  # # Use image flags written by Lammps
  # coords = minImageConv(coords, flags, box)
  # atoms[:,3:6] = coords.astype('|S6')

  # # unwrap atom coordinates in the absence of image flags
  # coords = unwrap(coords, box)
  # atoms[:,3:6] = coords.astype('|S6')

  # # Adjust box configuration and discard molecules entirely outside the box
  # coords[:,Z] -= coords[:,Z].mean()
  # atoms[:,3:6] = coords.astype('|S6')
  # box[Z,0] = coords[:,Z].min()
  # box[Z,1] = coords[:,Z].max()
  # box /= 2
  # # box[Z,:] /= 2
  # com = calcCOM(coords, nMon)
  # inside = np.logical_and(com[:,X] <= box[X,1], com[:,X] > box[X,0])
  # inside = inside * np.logical_and(com[:,Y] <= box[Y,1], com[:,Y] > box[Y,0])
  # inside = inside * np.logical_and(com[:,Z] <= box[Z,1], com[:,Z] > box[Z,0])
  # inside = inside.astype(bool)
  # coords = coords[inside]
  # atoms = atoms[inside]
  # info = info[inside]

  # # Tile the current configuration
  # coords[:,X] -= box[X,0]
  # image = coords + (box[:,1] - box[:,0]) * np.array([1,0,0])
  # box[X,:] = box[X,:] * 2
  # coords = np.concatenate((coords,image),axis=0)
  # info = np.concatenate((info,info), axis=0)
  # atoms = np.concatenate((info.astype(str),coords.astype(str)), axis=1)

  # renumber atoms
  masks = (info[:,2] == 1, info[:,2] == 2, info[:,2] == 3)
  types = (sum(masks[0]), sum(masks[1]), sum(masks[2]))
  info = renumber(info, types, nMon) # , diffTypes=True) # , True)
  atoms[:,:3] = info.astype(str)

  # # Generate new bonds
  # poly = masks[0]
  # polyStart = info[poly][0,0]
  # nPoly = sum(poly)
  # bonds = np.ones((nPoly * (nMon - 1) / nMon,4), dtype=int)
  # bonds[:,0] = np.arange(len(bonds)) + 1
  # mask = ((info[poly][:,0] - polyStart + 1) % nMon).astype(bool)
  # bonds[:,2] = info[mask][:,0] + polyStart - 1
  # bonds[:,3] = info[mask][:,0] + polyStart

  # # create new box
  # box[:,0] = coords.min(0)
  # box[:,1] = coords.max(0)

  # convert data to string format for output
  box = box.astype('|S8')
  bonds = bonds.astype(str)

  try:
    filename = sys.argv[2]
  except IndexError:
    cont = raw_input('Warning: output filename is the same as the input file. This may overwrite data. Continue? [y]/n')
    if cont == 'no' or cont == 'n':
      sys.exit()

  # write output files
  suffix = ''
  write_conf(filename+suffix, atoms[:,:6], bonds, 
            'modified configuration', [info[-1,2],1], box, np.ones(info[-1,2],dtype=int))
  write_xyz(filename+suffix, 0, atoms[:,2:6])

if __name__ == '__main__':
  main()
