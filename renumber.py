import numpy as np
import argparse

from conf_tools import readConf, write_conf, write_xyz, write_traj
from pbc_tools import PBC, unwrap

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input",help='input file to renumber')
  parser.add_argument("-o","--output",help='name of output file')
  parser.add_argument("-f","--formats",nargs='+',choices=['xyz','lammps','conf'],
                      help='file formats to write output file to')
  parser.add_argument("--rebuild",help='build bonds from scratch',action='store_true')
  parser.add_argument("-d","--descrip",help='descriptive header for output file')
  parser.add_argument("-l","--length",type=int,help="chain length of polymer")
  parser.add_argument("--no-image-flags",dest="entries",action="store_const",const=6,
                      default=9,help="switch to adjust number of entries for data flags")
  args = parser.parse_args()
  
  # read input file
  with open(args.input,'r') as file:
    box, atoms, bonds = readConf(file)

  # sort the atom and bond arrays
  atoms = atoms[np.argsort(atoms[:,0])]
  bonds = bonds[np.argsort(bonds[:,2])]

  # generate new atom and bond ids
  atoms[:,0] = np.arange(atoms.shape[0])+1
  bonds[:,0] = np.arange(bonds.shape[0])+1
  
  # create unique set of atomtypes
  atomtypes = set(atoms[:,2])
  bondtypes = set(bonds[:,1])

  ## generate array masks
  masks = []
  for atomtype in atomtypes:
    masks.append(atoms[:,2] == atomtype)

  ### extract subarrays based on atomtype
  subarrays = []
  for mask in masks:
    subarrays.append(atoms[mask])

  # unwrap atoms
  #try: subarrays[1][:,3:6] = subarrays[1][:,3:6] + subarrays[1][:,6:]*(box[:,1] - box[:,0])
  #except ValueError: subarrays[1][:,3:] = unwrap(subarrays[1][:,3:],box)

  molecules = set(atoms[:,1].astype(int))
  for i, mol in enumerate(molecules):
    atoms[:,1][atoms[:,1].astype(int) == mol] = i+1
  nBeads = len(subarrays[1])
  nChains = nBeads / args.length
  pairs = subarrays[1][:,0].reshape((nChains, args.length))
  if args.rebuild:
    nBonds = nChains * (args.length - 1)
    bonds = np.ones((nChains, args.length - 1, 4))
    bonds[:,:,2] = pairs[:,:-1]
    bonds[:,:,3] = pairs[:,1:]
    bonds = bonds.reshape((nBonds, 4))
    bonds[:,0] = np.arange(nBonds) + 1
  else:
    nBonds = len(bonds)
    bonds = bonds.reshape((nChains, args.length - 1, 4))
    bonds[:,:,2] = pairs[:,:-1]
    bonds[:,:,3] = pairs[:,1:]
    bonds = bonds.reshape((nBonds, 4))

  # generate correct box bounds
  box[2,0] = atoms[:,5].min()
  box[2,1] = atoms[:,5].max()

  ## reinsert subarrays back into the atom array
  #for array, mask in zip(subarrays, masks):
  #  atoms[mask] = array

  # write output to xyz and conf files
  types = {"atoms":len(atomtypes), "bonds":len(bondtypes)}
  masses = [1] * types["atoms"]

  if 'conf' in args.formats:
    write_conf(args.output, atoms[:,:args.entries], bonds, box, types, masses, args.descrip)
  if 'xyz' in args.formats:
    write_xyz(args.output, atoms[:,2:6])  
  if 'lammps' in args.formats:
    write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')  
  
if __name__ == "__main__":
  main()
