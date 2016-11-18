from sys import argv, exit
import numpy as np
import argparse

from conf_tools import readConf, write_conf, write_xyz
from pbc_tools import unwrap, PBC

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
  polyMask = (atoms[:,2] == 1)

  # extract subarrays based on atomtype
  polymers = atoms[polyMask]

  # reshape the extracted array by molecule
  nBeads = 350000
  nMon = 10
  nChains = nBeads / nMon
  polymers = polymers.reshape((nChains, nMon, 6))

  # unwrap the extracted array
  polymers[:,:,3:] = unwrap(polymers[:,:,3:], box)

  # reshape the modified array so that it can be re-inserted
  polymers = polymers.reshape((nBeads, 6))

  # insert extracted array back into the original
  atoms[polyMask] = polymers
  
  # write output to xyz and conf files
  title = 'Renumbered N=10 configuration'
  types = {"atoms":3, "bonds":1}
  masses = [1,1,1]
  write_conf(args.output, atoms, bonds, box, types, masses, title)
  write_xyz(args.output, atoms[:,2:])  
  
if __name__ == "__main__":
  main()
