from sys import argv, exit
import numpy as np
import argparse

from conf_tools import readConf, write_conf, write_xyz

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input")
  parser.add_argument("-o","--output")
  args = parser.parse_args()
  
  # read input file
  with open(args.input,'r') as file:
    box, atoms, bonds = readConf(file)

  # convert lists to NumPy arrays
  box = np.array(box, dtype=float)
  atoms = np.array(atoms, dtype=float)
  bonds = np.array(bonds, dtype=float)

  # sort the atom array
  atoms = atoms[np.argsort(atoms[:,0])]
  
  # generate array masks
  partMask = (atoms[:,2] == 3)
  polyMask = (atoms[:,2] == 1)
  surfMask = (atoms[:,2] == 2)

  ## extract subarrays based on atomtype
  polymers = atoms[polyMask]
  surface = atoms[surfMask]

  # generate new polymer molecules
  nBeads = 350000
  nMon = 10
  nPart = 54
  nChains = nBeads / nMon
  nBonds = nChains * (nMon - 1)
  polymers[:,1] = np.repeat(np.arange(nChains)+1+nPart,nMon)
  atoms[surfMask][:,1] = nPart+nChains

  # generate new bonds
  pairs = atoms[polyMask][:,0].reshape((nChains, nMon))
  bonds = np.ones((nChains, nMon - 1, 4))
  bonds[:,:,2] = pairs[:,:-1]
  bonds[:,:,3] = pairs[:,1:]
  bonds = bonds.reshape((nBonds, 4))
  bonds[:,0] = np.arange(nBonds)+1

  ## combine the particle arrays into one
  atoms[polyMask] = polymers
  atoms[surfMask] = surface
  
  # write output to xyz and conf files
  title = 'Renumbered N=10 configuration'
  types = {"atoms":3, "bonds":1}
  masses = [1,1,1]
  write_conf(args.output, atoms, bonds, box, types, masses, title)
  write_xyz(args.output, atoms[:,2:])  
  
if __name__ == "__main__":
  main()
