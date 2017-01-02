import numpy as np
import argparse
from conf_tools import *

def main():
  'deletes molecules from a configuration'

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-n','--nMon', type=int)
  args = parser.parse_args()
  
  with open(args.input) as file:
    box, atoms, bonds = readConf(file, ['1','2','3'])
    DOS = 0.8

    # generate a few very large bins (larger than the RE of the polymer)
    # count the number of molecules in each bin by computing which bin has the most beads from each chain
    # delete a given fraction of molecules from a given bin
    
    atoms = atoms[np.argsort(atoms[:,0])]
    polymers = atoms[atoms[:,2] == 1] #[:,5]
    polymers = polymers.reshape((len(polymers)/args.nMon, args.nMon, atoms.shape[-1]))
    atoms = atoms[atoms[:,2] != 1]
    print polymers.shape
    
    # nBins = 15
    # bin_size = (box[2,1] - box[2,0]) / nBins

    # in_bin = (polymers - box[2,0])

    # try other methods first:
    # method 1: remove final molecules
    # method 2: remove every few molecules
    # method 3: remove at random

    # method 1
    nMolecules = int((1 - DOS) * polymers.shape[0])
    polymers = polymers[:nMolecules]
    print polymers.shape

    # reshape polymers
    polymers = polymers.reshape((-1, atoms.shape[-1]))

    # write modified polymers into atom matrix
    atoms = np.concatenate((polymers, atoms), axis=0)

    write_conf(args.output, atoms.astype('|S16'), bonds, box, {"atoms":3, "bonds":1}, [1]*3, '% undersaturated trajectory')
    write_traj(args.output, np.delete(atoms[:,:6],1,1), box)
    

if __name__ == '__main__':
  main()
