from sys import argv, exit
import numpy as np
import argparse

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

  if args.entries == 6:

    # reshape the extracted array by molecule
    nBeads = len(polymers)
    nChains = nBeads / args.length
    polymers = polymers.reshape((nChains, args.length, args.entries))

    # unwrap the extracted array
    polymers[:,:,3:] = unwrap(polymers[:,:,3:], box)

    # reshape the modified array so that it can be re-inserted
    polymers = polymers.reshape((nBeads, args.entries))

  elif args.entries == 9:

    polymers[:,3:6] = polymers[:,3:6] + (box[:,1] - box[:,0]) * polymers[:,6:9]

  # insert extracted array back into the original
  atoms[polyMask] = polymers
  
  # write output to xyz and conf files
  types = {"atoms":4, "bonds":1}
  masses = [1] * types['atoms']
  write_conf(args.output, atoms[:,:args.entries], bonds, box, types, masses, args.descrip)
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box)
  
if __name__ == "__main__":
  main()
