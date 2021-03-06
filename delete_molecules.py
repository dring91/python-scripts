import numpy as np
import argparse
from conf_tools import *

def main():
  'deletes molecules from a configuration'

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-l','--length', type=int)
  parser.add_argument('-s','--saturation', type=float, 
                      help='specify the saturation level; \
                            that is, tell the algorithm what fraction of molecules \
                            to remove from the configuration')
  parser.add_argument('method', choices={'random','region'},
                      help='choose the method for removing molecules')
  args = parser.parse_args()
  
  with open(args.input) as file:
    box, atoms, bonds = readConf(file)

  # degree of (under-)saturation (used to calculate the number of molecules to remove)
  #DOS = 0.8
  DOS = args.saturation
 
  # sort atoms according to atomic ID
  atoms = atoms[np.argsort(atoms[:,0])]
  # extract only polymers (type == 1)
  polymers = atoms[atoms[:,2] == 1]
  # reshape polymers so that each subarray contains exactly one molecule
  polymers = polymers.reshape((len(polymers)/args.length, args.length, atoms.shape[-1]))
  # erase polymers from atoms array
  atoms = atoms[atoms[:,2] != 1]
  print polymers.shape
  
  # calculate the number of molecules to keep in the polymers array
  nMolecules = int((1 - DOS) * polymers.shape[0])
  # slice polymers array and keep nMolecules
  if args.method == 'random':
    polymers = polymers[:nMolecules]
    print polymers.shape
  elif args.method == 'region':
    cx = box[:2,1] - box[:2,0]
    area = cx[0] * cx[1]
    height = nMolecules * args.length / area
    com = polymers[:,:,5].mean(1)
    polymers = polymers[com < height]
    print polymers.shape

  # reshape polymers
  polymers = polymers.reshape((-1, atoms.shape[-1]))

  # write modified polymers into atom matrix
  atoms = np.concatenate((polymers, atoms), axis=0)

  write_conf(args.output, atoms.astype('|S16'), bonds, box, {"atoms":2, "bonds":1}, [1]*2, '% undersaturated trajectory')
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box)
  

if __name__ == '__main__':
  main()
