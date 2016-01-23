#!/usr/bin/python

import sys
from conf_tools import *
import numpy as np

MODE = 'r'

def renumber(atoms, chainLength):

  atoms[:,1] = (atoms[:,0] - 1) % chainLength
  atoms[:,1] = (atoms[:,1] - 25) 
  mask = (atoms[:,1] >= 0)
  atoms[:,1] = (atoms[:,1] + 1) * mask + abs(atoms[:,1]) * np.logical_not(mask)

  return atoms

def main():
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  atype = ['1','2']
  nMon = 50

  with open(filename+'.conf', MODE) as inp:
    box, atoms, bonds = readConf(inp, atype)

  atoms = np.array(atoms)
  info = atoms[:,:3].astype(int)
  info = renumber(info, nMon)
  atoms[:,:3] = info.astype(str)

  write_conf(filename+'_out', atoms, bonds, 
            'atoms renumbered for bond/swap', [2,1], box, [1,1])
  write_xyz(filename+'_out', atoms[:,2:])

if __name__ == '__main__':
  main()
