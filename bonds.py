from sys import argv, exit
import numpy as np
import argparse
import matplotlib.pyplot as plt

from conf_tools import readConf, write_conf, write_xyz, write_traj
from pbc_tools import unwrap, PBC

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input")
  parser.add_argument("-o","--output")
  parser.add_argument("-d","--descrip",help="Header for file describing what it contains")
  parser.add_argument("-l","--length",help="chain length of polymer molecule",type=int)
  parser.add_argument("-t","--atomtype",help="atomtype ID for layer to unwrap",type=int)
  parser.add_argument("--no-image-flags",dest="entries",action="store_const",const=6,
                      default=9,help="switch to adjust number of entries for data flags")
  args = parser.parse_args()
  
  # read input file
  with open(args.input,'r') as file:
    box, atoms, bonds = readConf(file)

  # sort the atom array
  atoms = atoms[np.argsort(atoms[:,0])]
  
  # generate array masks
  polyMask = (atoms[:,2] == args.atomtype)

  # extract subarrays based on atomtype
  polymers = atoms[polyMask]

  ## calculate the appropriate box
  box = np.zeros((3,2))
  box[:,0] = polymers[:,3:6].min(0) - 0.5
  box[:,1] = polymers[:,3:6].max(0) + 0.5

  # reshape the extracted array by molecule
  nBeads = len(polymers)
  nChains = nBeads / args.length
  polymers = polymers.reshape((nChains, args.length, args.entries))
  fractional, integral = np.modf(polymers[:,:,1])
  print sum(fractional), sum(integral)

  # calculate bond lengths
  differences = polymers[:,1:,3:6] - polymers[:,:-1,3:6]
  lengths = np.sqrt(np.einsum('ijk,ijk->ij',differences,differences))
  lengths = lengths.reshape((args.length - 1) * nChains)
  end_to_end = differences.sum(1)
  R = np.sqrt(np.einsum('ij,ij->i',end_to_end,end_to_end))

  #plt.plot(lengths)
  #plt.show()
  #plt.plot(R/np.sqrt(10))
  #plt.show()

  # check that bonds are numbered correctly
  bonds = bonds.reshape((nChains, args.length - 1, 4))
  head = polymers[:,1:,0]
  tail = polymers[:,:-1,0]
  print sum(bonds[:,:,2] - tail)
  print sum(bonds[:,:,3] - head)

  ## reshape the modified array so that it can be re-inserted
  #polymers = polymers.reshape((nBeads, args.entries))

  ## insert extracted array back into the original
  #atoms[polyMask] = polymers
  
if __name__ == "__main__":
  main()
