from sys import argv, exit
import numpy as np
import argparse

from conf_tools import readConf, write_conf, write_xyz, write_traj

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input")
  parser.add_argument("-o","--output")
  args = parser.parse_args()
  
  # read input file
  with open(args.input,'r') as file:
    box, atoms, bonds = readConf(file, atype=['1.0','2.0','3.0'])

  # sort the atom array
  atoms = atoms[np.argsort(atoms[:,0])]

  # generate new atom ids
  atoms[:,0] = np.arange(len(atoms))+1
  
  ## generate array masks
  partMask = (atoms[:,2] >= 3)
  polyMask = (atoms[:,2] == 1)
  surfMask = (atoms[:,2] == 2)

  ### extract subarrays based on atomtype
  particles = atoms[partMask]
  polymers = atoms[polyMask]
  surface = atoms[surfMask]

  ## generate new polymer molecules
  nBeads = len(polymers) # 350000
  nMon = 10
  nPart = 54
  nSites = 4684
  nChains = nBeads / nMon
  nBonds = nChains * (nMon - 1)
  particles[:,1] = np.repeat(np.arange(nPart)+1,nSites)
  polymers[:,1] = np.repeat(np.arange(nChains)+1+nPart,nMon)
  surface[:,1] = nPart+nChains+1

  ## number molecules
  #particles[:,1] += surface[-1,1]

  ### generate new bonds
  pairs = atoms[polyMask][:,0].reshape((nChains, nMon))
  bonds = np.ones((nChains, nMon - 1, 4))
  bonds[:,:,2] = pairs[:,:-1]
  bonds[:,:,3] = pairs[:,1:]
  bonds = bonds.reshape((nBonds, 4))
  #bonds = bonds[np.argsort(bonds[:,2])]
  nBonds = len(bonds)
  bonds[:,0] = np.arange(nBonds)+1

  ## reinsert subarrays back into the atom array
  atoms[partMask] = particles
  atoms[polyMask] = polymers
  atoms[surfMask] = surface

  # write output to xyz and conf files
  #title = 'Renumbered N=10 configuration'
  title = 'Renumbered cylinder and N=10 polymer configuration'
  types = {"atoms":3, "bonds":1}
  masses = [1] * types["atoms"]
  write_conf(args.output, atoms, bonds, box, types, masses, title)
  write_xyz(args.output, atoms[:,2:6])  
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')  
  
if __name__ == "__main__":
  main()
